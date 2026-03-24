import argparse
import os
import time
import urllib.request
import urllib.error
import concurrent.futures

def download_fasta(uniprot_id):
    """
    Downloads the FASTA sequence for a given UniProt ID.
    Returns the sequence string or None if it fails.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    max_retries = 3
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req) as response:
                fasta_data = response.read().decode('utf-8')
            lines = fasta_data.strip().split('\n')
            if not lines:
                return None
            # Ignore the first line (header) and join the rest
            sequence = "".join(line.strip() for line in lines[1:])
            return sequence
        except urllib.error.URLError as e:
            time.sleep(1)
        except Exception as e:
            # 에러 메시지 줄바꿈 추가 (진행률과 겹치지 않도록)
            print(f"\nError downloading {uniprot_id}: {e}")
            break
            
    print(f"\nCould not download sequence for UniProt ID: {uniprot_id}")
    return None

def download_fasta_by_gene(gene_name):
    """
    Downloads the FASTA sequence for a given gene name by searching UniProt.
    Returns the sequence string for the first human protein hit or None if it fails.
    """
    # reviewed:true 추가: 검증된(Swiss-Prot) 데이터만 검색하도록 수정
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene_name}+AND+organism_id:9606+AND+reviewed:true&format=fasta"
    max_retries = 3
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url)
            # timeout 10초 추가
            with urllib.request.urlopen(req, timeout=10) as response:
                fasta_data = response.read().decode('utf-8')
            lines = fasta_data.strip().split('\n')
            if not lines:
                return None
            
            sequence_lines = []
            for line in lines[1:]:
                if line.startswith('>'):
                    break
                sequence_lines.append(line.strip())
            
            sequence = "".join(sequence_lines)
            return sequence if sequence else None
        except urllib.error.URLError as e:
            time.sleep(1)
        except Exception as e:
            print(f"\nError downloading gene {gene_name}: {e}")
            break
            
    print(f"\nCould not download sequence for gene: {gene_name}")
    return None

def load_local_fasta(fasta_path):
    print(f"Loading local FASTA from {fasta_path}...")
    local_db = {'uniprot': {}, 'gene': {}}
    try:
        with open(fasta_path, 'r', encoding='utf-8') as f:
            current_header = ""
            current_seq = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_header:
                        seq = "".join(current_seq)
                        parts = current_header.split()
                        id_part = parts[0][1:] # remove >
                        id_split = id_part.split('|')
                        if len(id_split) >= 2:
                            uniprot_id = id_split[1]
                            local_db['uniprot'][uniprot_id] = seq
                        
                        for part in parts[1:]:
                            if part.startswith('GN='):
                                gene_name = part[3:]
                                if gene_name not in local_db['gene']:
                                    local_db['gene'][gene_name] = seq
                                break

                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line)
            # last sequence
            if current_header:
                seq = "".join(current_seq)
                parts = current_header.split()
                id_part = parts[0][1:]
                id_split = id_part.split('|')
                if len(id_split) >= 2:
                    uniprot_id = id_split[1]
                    local_db['uniprot'][uniprot_id] = seq
                
                for part in parts[1:]:
                    if part.startswith('GN='):
                        gene_name = part[3:]
                        if gene_name not in local_db['gene']:
                            local_db['gene'][gene_name] = seq
                        break
        print(f"Loaded {len(local_db['uniprot'])} UniProt IDs and {len(local_db['gene'])} gene names from local FASTA.")
    except Exception as e:
        print(f"Error loading local FASTA: {e}")
    return local_db

def fetch_sequence(identifier_tuple):
    identifier, id_type = identifier_tuple
    if id_type == "uniprot":
        return download_fasta(identifier)
    else:
        return download_fasta_by_gene(identifier)

def main():
    parser = argparse.ArgumentParser(description="Create Boltz-2 input yaml scripts for protein-ligand prediction.")
    parser.add_argument("-i", "--input_list", required=True, nargs='+', help="Input files containing protein list (.tsv or .fasta) (e.g. file1.tsv file2.fasta)")
    parser.add_argument("-s", "--smiles_file", required=True, help="Small molecule SMILES text file (e.g. CK6)")
    parser.add_argument("-o", "--out_dirs", nargs='+', default=None, help="Output directories for YAML files (must match number of input_list files if provided)")
    parser.add_argument("-w", "--workers", type=int, default=10, help="Number of parallel workers for downloading FASTA files")
    parser.add_argument("--db", default=None, help="Path to a local FASTA database file to search before downloading from UniProt")
    parser.add_argument("--retry_log", default=None, help="Path to a log file containing failed downloads to retry")
    parser.add_argument("--error_log", default="failed_downloads.log", help="Path to save the failed downloads (default: failed_downloads.log)")
    parser.add_argument("--out_type", choices=['job', 'fasta', 'both'], default='job', help="Type of output to generate: 'job' for yaml only, 'fasta' for fasta only, 'both' for both (default: job)")
    parser.add_argument("--a3m", default=None, help="Directory containing .a3m files to include in the yaml job files")

    args = parser.parse_args()

    if args.out_dirs is not None and len(args.out_dirs) != len(args.input_list):
        print("Error: Number of output directories (-o/--out_dirs) must match the number of input files (-i/--input_list).")
        return

    retry_ids = None
    if args.retry_log and os.path.exists(args.retry_log):
        retry_ids = set()
        with open(args.retry_log, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split('\t')
                    retry_ids.add((parts[0], parts[1] if len(parts) > 1 else "uniprot"))
        print(f"Loaded {len(retry_ids)} IDs to retry from {args.retry_log}")

    # Determine output directories
    out_dirs = []
    fasta_dirs = []
    for i, input_file in enumerate(args.input_list):
        if args.out_dirs is not None:
            out_dir = args.out_dirs[i]
            fasta_dir = f"{out_dir}_fasta"
        else:
            input_basename = os.path.splitext(os.path.basename(input_file))[0]
            out_dir = f"{input_basename}_job"
            fasta_dir = f"{input_basename}_fasta"
        out_dirs.append(out_dir)
        fasta_dirs.append(fasta_dir)
        if args.out_type in ['job', 'both']:
            os.makedirs(out_dir, exist_ok=True)
        if args.out_type in ['fasta', 'both']:
            os.makedirs(fasta_dir, exist_ok=True)
    
    # Read the SMILES string
    with open(args.smiles_file, 'r', encoding='utf-8') as f:
        smiles_string = f.read().strip()
        
    print(f"Loaded SMILES from {args.smiles_file}: {smiles_string}")
    
    # Base name of smiles file for naming output files
    smiles_basename = os.path.splitext(os.path.basename(args.smiles_file))[0]
    
    seen_uniprots = {}

    # Parse input files to collect all unique protein identifiers
    items_to_download = set()
    input_data = [] # Store parsed data for later processing

    for input_file in args.input_list:
        file_ext = os.path.splitext(input_file)[1].lower()

        if file_ext == '.fasta':
            # Parse FASTA file - directly use sequences from file
            with open(input_file, 'r', encoding='utf-8') as f:
                fasta_lines = f.readlines()

            if not fasta_lines:
                print(f"Warning: FASTA file {input_file} is empty. Skipping.")
                input_data.append((input_file, 'fasta', [], None, None, None))
                continue

            # Parse FASTA file to extract IDs and sequences
            uniprot_ids = []
            current_id = None
            current_seq = []
            current_name = ""

            for line in fasta_lines:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_id and current_seq:
                        seq = "".join(current_seq)
                        seen_uniprots[current_id] = seq
                        uniprot_ids.append((current_id, current_name))

                    # Parse new header
                    # Format: >sp|P12345|PROTNAME or >P12345
                    current_name = ""
                    parts = line[1:].split()
                    if parts:
                        id_part = parts[0]
                        id_split = id_part.split('|')
                        if len(id_split) >= 2:
                            current_id = id_split[1]  # UniProt ID from >sp|P12345|...
                            if len(id_split) >= 3:
                                current_name = id_split[2]
                        else:
                            current_id = id_split[0]  # UniProt ID from >P12345
                            
                        if not current_name:
                            for p in parts[1:]:
                                if p.startswith("GN="):
                                    current_name = p[3:]
                                    break
                        current_seq = []
                else:
                    # Sequence line
                    current_seq.append(line)

            # Save last sequence
            if current_id and current_seq:
                seq = "".join(current_seq)
                seen_uniprots[current_id] = seq
                uniprot_ids.append((current_id, current_name))

            input_data.append((input_file, 'fasta', uniprot_ids, None, None, None))

        elif file_ext == '.tsv':
            # Parse TSV file (existing logic)
            with open(input_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()

            if not lines:
                print(f"Warning: TSV file {input_file} is empty. Skipping.")
                input_data.append((input_file, 'tsv', [], -1, None, -1))
                continue

            header = lines[0].strip('\n').split('\t')
            try:
                id_idx = header.index("source_uniprot")
                id_type = "uniprot"
            except ValueError:
                try:
                    id_idx = header.index("gene_name")
                    id_type = "gene"
                except ValueError:
                    print(f"Warning: Neither 'source_uniprot' nor 'gene_name' column found in TSV header for {input_file}. Skipping.")
                    input_data.append((input_file, 'tsv', [], -1, None, -1))
                    continue

            header_lower = [h.strip().lower() for h in header]
            source_idx = header_lower.index("source") if "source" in header_lower else -1

            input_data.append((input_file, 'tsv', lines, id_idx, id_type, source_idx))

            for i, line in enumerate(lines[1:], start=2):
                line = line.strip('\n')
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) <= id_idx:
                    continue

                id_str = parts[id_idx]
                if not id_str.strip():
                    continue

                if id_str.startswith("COMPLEX:"):
                    complex_ids = id_str.replace("COMPLEX:", "").split("_")
                    for c_id in complex_ids:
                        if retry_ids is None or (c_id, id_type) in retry_ids:
                            items_to_download.add((c_id, id_type))
                else:
                    if retry_ids is None or (id_str, id_type) in retry_ids:
                        items_to_download.add((id_str, id_type))
        else:
            print(f"Warning: Unknown file type '{file_ext}' for {input_file}. Supported types: .tsv, .fasta. Skipping.")
            input_data.append((input_file, 'unknown', [], None, None, None))
            continue

    items_to_download = list(items_to_download)

    if args.db:
        local_db = load_local_fasta(args.db)
        items_to_fetch_from_web = []
        for item in items_to_download:
            identifier, id_type = item
            if id_type == "uniprot" and identifier in local_db['uniprot']:
                seen_uniprots[identifier] = local_db['uniprot'][identifier]
            elif id_type == "gene" and identifier in local_db['gene']:
                seen_uniprots[identifier] = local_db['gene'][identifier]
            else:
                items_to_fetch_from_web.append(item)

        items_to_download = items_to_fetch_from_web
        
    total_items = len(items_to_download)
    if not items_to_download and not seen_uniprots:
        print("No valid protein IDs found in the provided input files.")
        return
        
    if total_items > 0:
        print(f"Found {total_items} unique queries to download. Downloading with {args.workers} workers...\n")
    else:
        print("All requested sequences were found in the local FASTA file. No downloads needed.\n")
    
    failed_items = []
    
    # Download in parallel with progress tracking
    completed_items = 0
    if total_items > 0:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
            future_to_item = {executor.submit(fetch_sequence, item): item for item in items_to_download}
            for future in concurrent.futures.as_completed(future_to_item):
                identifier, id_type = future_to_item[future]
                try:
                    seq = future.result()
                    if seq:
                        seen_uniprots[identifier] = seq
                    else:
                        failed_items.append((identifier, id_type))
                except Exception as exc:
                    print(f"\n[{identifier}] generated an exception: {exc}")
                    failed_items.append((identifier, id_type))
                
                # 진행률 계산 및 출력 (\r을 사용해 같은 줄에 덮어쓰기)
                completed_items += 1
                percent = (completed_items / total_items) * 100
                print(f"\r다운로드 진행 상황: {completed_items}/{total_items} 완료 ({percent:.1f}%)", end="", flush=True)

    with open(args.error_log, 'w', encoding='utf-8') as f:
        for identifier, id_type in failed_items:
            f.write(f"{identifier}\t{id_type}\n")

    if failed_items:
        print(f"\n\nFailed to download {len(failed_items)} sequences. Logged to {args.error_log}")
        print(f"Successfully downloaded {len(seen_uniprots)} FASTA sequences.")
    else:
        print(f"\n\nSuccessfully downloaded {len(seen_uniprots)} FASTA sequences.")
    
    # Process each input file and create YAMLs
    total_processed = 0
    all_out_dirs = set()

    for (input_file, file_type, data, param1, param2, param3), out_dir, fasta_dir in zip(input_data, out_dirs, fasta_dirs):
        processed_count = 0
        all_out_dirs.add(out_dir)

        if file_type == 'fasta':
            # Process FASTA input: data contains list of UniProt IDs
            uniprot_ids = data
            if not uniprot_ids:
                continue

            for u_id, prot_name in uniprot_ids:
                if not u_id:
                    continue

                seq = seen_uniprots.get(u_id)
                if not seq:
                    continue

                # Create prefix from protein name
                import re
                clean_name = re.sub(r'[^A-Za-z0-9_\-]', '_', prot_name) if prot_name else ""
                prefix = f"{clean_name}_" if clean_name else ""

                # Create YAML specific to Boltz-2 format
                msa_path = ""
                if args.a3m:
                    possible_a3m = os.path.join(args.a3m, f"{u_id}.a3m")
                    if os.path.exists(possible_a3m):
                        # Use forward slash for yaml path safety
                        abs_a3m = os.path.abspath(possible_a3m).replace('\\', '/')
                        msa_path = f"\n      msa: {abs_a3m}"

                yaml_content = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {seq}{msa_path}
  - ligand:
      id: B
      smiles: '{smiles_string}'
"""

                # Save YAML file
                if args.out_type in ['job', 'both']:
                    out_filename = f"{prefix}{u_id}_{smiles_basename}.yaml"
                    out_filepath = os.path.join(out_dir, out_filename)

                    with open(out_filepath, 'w', encoding='utf-8') as out_f:
                        out_f.write(yaml_content)

                # Save FASTA file
                if args.out_type in ['fasta', 'both']:
                    fasta_filename = f"{prefix}{u_id}.fasta"
                    fasta_filepath = os.path.join(fasta_dir, fasta_filename)

                    with open(fasta_filepath, 'w', encoding='utf-8') as fasta_f:
                        fasta_f.write(f">{u_id}\n{seq}\n")

                processed_count += 1
                total_processed += 1

        elif file_type == 'tsv':
            # Process TSV input: data contains lines, param1=id_idx, param2=id_type
            lines = data
            id_idx = param1
            id_type = param2
            source_idx = param3

            if not lines or id_idx == -1:
                continue

            for i, line in enumerate(lines[1:], start=2):
                line = line.strip('\n')
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) <= id_idx:
                    continue

                id_str = parts[id_idx]
                if not id_str.strip():
                    continue

                prot_name = ""
                if source_idx != -1 and len(parts) > source_idx:
                    prot_name = parts[source_idx].strip()
                import re
                clean_name = re.sub(r'[^A-Za-z0-9_\-]', '_', prot_name) if prot_name else ""
                prefix = f"{clean_name}_" if clean_name else ""

                ids_to_process = []
                if id_str.startswith("COMPLEX:"):
                    complex_ids = id_str.replace("COMPLEX:", "").split("_")
                    ids_to_process.extend(complex_ids)
                else:
                    ids_to_process.append(id_str)

                for u_id in ids_to_process:
                    if not u_id:
                        continue

                    seq = seen_uniprots.get(u_id)
                    if not seq:
                        continue

                    # Create YAML specific to Boltz-2 format
                    msa_path = ""
                    if args.a3m:
                        possible_a3m = os.path.join(args.a3m, f"{u_id}.a3m")
                        if os.path.exists(possible_a3m):
                            # Use forward slash for yaml path safety
                            abs_a3m = os.path.abspath(possible_a3m).replace('\\', '/')
                            msa_path = f"\n      msa: {abs_a3m}"

                    yaml_content = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {seq}{msa_path}
  - ligand:
      id: B
      smiles: '{smiles_string}'
"""

                    # Save YAML file
                    if args.out_type in ['job', 'both']:
                        out_filename = f"{prefix}{u_id}_{smiles_basename}.yaml"
                        out_filepath = os.path.join(out_dir, out_filename)

                        with open(out_filepath, 'w', encoding='utf-8') as out_f:
                            out_f.write(yaml_content)

                    # Save FASTA file
                    if args.out_type in ['fasta', 'both']:
                        fasta_filename = f"{prefix}{u_id}.fasta"
                        fasta_filepath = os.path.join(fasta_dir, fasta_filename)

                        with open(fasta_filepath, 'w', encoding='utf-8') as fasta_f:
                            fasta_f.write(f">{u_id}\n{seq}\n")

                    processed_count += 1
                    total_processed += 1

        out_msg = f"File '{input_file}': Generated {processed_count} files"
        if args.out_type in ['job', 'both']:
            out_msg += f" in '{out_dir}'"
        if args.out_type in ['fasta', 'both']:
            out_msg += f" and '{fasta_dir}'"
        print(out_msg)
        
    print(f"\nFinished processing a total of {total_processed} target sequences.")
    if args.out_type in ['job', 'both']:
        for out_dir in sorted(list(all_out_dirs)):
            print(f"You can run predictions using: boltz predict {out_dir} --use_msa_server")

if __name__ == "__main__":
    main()
