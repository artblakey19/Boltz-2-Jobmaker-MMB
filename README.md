# jobmaker.py

`jobmaker.py`는 단백질-리간드 결합 예측 모델인 **Boltz-2**의 입력 파일(YAML 형식)을 자동 생성해주는 파이썬 스크립트입니다.
표적 단백질 목록이 포함된 TSV 파일과 리간드의 SMILES 문자열을 입력받아, UniProt에서 단백질의 FASTA 서열을 다운로드하고 예측 작업에 필요한 디렉토리와 파일들을 구성합니다.

## 주요 기능

* **자동 FASTA 다운로드**: UniProt REST API를 사용하여 UniProt ID 또는 유전자 이름(Gene Name)으로 단백질 서열을 자동으로 다운로드합니다.
* **병렬 처리 지원**: 여러 단백질 서열을 빠르게 다운로드하기 위해 다중 워커(멀티스레딩)를 지원합니다.
* **로컬 DB 캐싱**: 로컬 FASTA 파일을 사용하여 이미 다운로드된 서열은 웹 요청 없이 빠르게 불러옵니다.
* **실패 복구 메커니즘**: 다운로드 실패 시 로그 파일(`failed_downloads.log`)을 생성하며, 이후 해당 로그를 바탕으로 실패한 항목만 재시도할 수 있습니다.
* **복합체(Complex) 지원**: 멀티머 형태의 단백질 복합체(예: `COMPLEX:P12345_Q98765`) 처리를 지원합니다.
* **Custom MSA 연동**: 사용자가 로컬에 가지고 있는 `.a3m` 다중 서열 정렬(MSA) 파일 디렉토리를 지정하면 생성되는 YAML 파일에 자동으로 경로를 추가합니다.
* **유연한 출력 옵션**: Boltz-2 설정 파일인 `job` 파일(YAML), `fasta` 파일, 혹은 둘 다(`both`) 생성할지 선택할 수 있습니다.

## 사용 방법

```bash
python jobmaker.py -t <단백질_TSV_파일들> -s <리간드_SMILES_파일> [옵션들]
```

### 필수 인수

* `-t` 또는 `--tsv_files`: 표적 단백질 정보가 들어있는 TSV 파일 목록 (여러 개 지정 가능). 공백으로 구분합니다.
* `-s` 또는 `--smiles_file`: 리간드의 SMILES 문자열이 적힌 텍스트 파일.

### 선택 인수

* `-o` 또는 `--out_dirs`: 결과물이 저장될 출력 디렉토리 이름 목록 (TSV 파일 개수와 같아야 함). 지정하지 않으면 입력 파일 이름을 기반으로 자동 생성됩니다 (`<파일명>_job`, `<파일명>_fasta`).
* `-w` 또는 `--workers`: FASTA 다운로드에 사용할 병렬 워커 수 (기본값: 10).
* `-f` 또는 `--fasta`: 웹에서 다운로드하기 전 먼저 검색할 로컬 FASTA 파일 경로. 이 파일을 활용하면 불필요한 API 호출을 줄일 수 있습니다.
* `--retry_log`: 다운로드 실패 로그 파일 경로. 이 파일을 지정하면 전체 목록을 다시 검색하지 않고 해당 파일에 있는 ID만 다시 다운로드를 시도합니다.
* `--error_log`: 다운로드 실패 시 기록할 로그 파일의 이름 (기본값: `failed_downloads.log`).
* `--out_type`: 생성할 출력 파일의 종류. `job` (YAML만), `fasta` (FASTA만), `both` (둘 다) 중 선택할 수 있습니다. (기본값: `job`).
* `--a3m`: 생성될 YAML 작업 파일에 포함할 `.a3m` 파일들이 있는 디렉토리 절대(또는 상대) 경로. 해당 경로 내의 파일명은 `<단백질ID>.a3m` 형식이어야 매칭됩니다.

## 입력 파일 형식

### 1. 단백질 TSV 파일 (`-t`)

TSV 파일의 첫 번째 줄(헤더)에는 반드시 `source_uniprot` 또는 `gene_name` 컬럼이 포함되어 있어야 합니다. 스크립트는 이 컬럼의 값을 바탕으로 단백질 서열을 검색하고 다운로드합니다.

**(예시)**

```tsv
source_uniprot	other_info
P04637	...
Q01196	...
COMPLEX:P04637_Q01196	...
```

### 2. SMILES 파일 (`-s`)

리간드의 SMILES 문자열이 포함된 단일 텍스트 파일이어야 합니다.

**(예시: `ligand.smi`)**

```text
CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5
```

## 출력 형식 (Boltz-2 YAML)

생성되는 YAML 파일은 Boltz-2 모델이 인식할 수 있는 구조로 자동 작성됩니다:

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTED...
      msa: /absolute/path/to/P04637.a3m  # --a3m 옵션 사용 및 파일 존재 시 추가됨
  - ligand:
      id: B
      smiles: 'CC1=C(C=C...'
```

## 사용 예시

**기본 사용 (TSV 파일과 SMILES를 이용해 YAML 작업 파일들 생성)**

```bash
python jobmaker.py -t targets1.tsv targets2.tsv -s my_ligand.txt
```

**출력 형식 변경 (FASTA와 YAML 모두 생성)**

```bash
python jobmaker.py -t targets.tsv -s my_ligand.txt --out_type both
```

**로컬 FASTA 활용 및 커스텀 A3M MSA 파일 연동**

```bash
python jobmaker.py -t targets.tsv -s my_ligand.txt -f local_db.fasta --a3m ./my_msa_files/
```

**실패한 항목(단백질 다운로드 간섭 등) 재시도**

```bash
python jobmaker.py -t targets.tsv -s my_ligand.txt --retry_log failed_downloads.log
```
## Data Source
The protein sequence data (`uniprotkb_taxonomy_id_9606_AND_reviewed_2026_03_17.fasta`) provided in this repository is sourced from the UniProt Knowledgebase (UniProtKB).

* **Organism:** Homo sapiens (Human, TaxID: 9606)
* **Dataset:** Reviewed (Swiss-Prot)
* **Date Downloaded:** 2026-03-17
* **License:** This data is distributed under the [Creative Commons Attribution 4.0 International (CC BY 4.0) License](https://creativecommons.org/licenses/by/4.0/).
