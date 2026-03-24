#!/usr/bin/env python3
"""
중복 제거 스크립트 (Deduplicate Jobs)

여러 job 폴더에서 중복된 단백질 작업(YAML 파일)을 찾아서 정리합니다.
같은 단백질-리간드 조합의 YAML 파일이 여러 폴더에 있을 경우,
가장 먼저 발견된 폴더의 파일만 남기고 나머지는 삭제합니다.

사용법:
    python deduplicate_jobs.py job_folder1 job_folder2 job_folder3 ...
    python deduplicate_jobs.py *_job --log dedup.log
    python deduplicate_jobs.py --dirs job_folder1 job_folder2 --log my_dedup.log --dry-run
"""

import argparse
import os
import sys
from datetime import datetime
from collections import defaultdict


def extract_job_key(filename):
    """
    YAML 파일명에서 job의 고유 키를 추출합니다.
    예: P12345_CK6.yaml -> (P12345, CK6)

    Returns:
        tuple: (protein_id, ligand_id) 또는 None (파싱 실패 시)
    """
    if not filename.endswith('.yaml'):
        return None

    basename = filename[:-5]  # .yaml 제거

    # '_'로 분리 (마지막 부분이 리간드 ID)
    parts = basename.rsplit('_', 1)
    if len(parts) != 2:
        return None

    protein_id, ligand_id = parts
    return (protein_id, ligand_id)


def scan_job_folders(folder_paths):
    """
    여러 job 폴더를 스캔하여 YAML 파일 정보를 수집합니다.

    Args:
        folder_paths: job 폴더 경로 리스트

    Returns:
        dict: {(protein_id, ligand_id): [(folder_path, filename), ...]}
    """
    job_map = defaultdict(list)

    for folder in folder_paths:
        if not os.path.isdir(folder):
            print(f"Warning: '{folder}' is not a valid directory. Skipping.")
            continue

        # 폴더 내의 모든 YAML 파일 스캔
        try:
            files = os.listdir(folder)
        except PermissionError:
            print(f"Warning: Permission denied for '{folder}'. Skipping.")
            continue

        for filename in files:
            if not filename.endswith('.yaml'):
                continue

            job_key = extract_job_key(filename)
            if job_key is None:
                continue

            filepath = os.path.join(folder, filename)
            job_map[job_key].append((folder, filename, filepath))

    return job_map


def deduplicate_jobs(job_map, dry_run=False):
    """
    중복된 job을 정리합니다.
    각 job_key에 대해 첫 번째 파일만 남기고 나머지는 삭제합니다.

    Args:
        job_map: scan_job_folders의 반환값
        dry_run: True이면 실제 삭제하지 않고 로그만 출력

    Returns:
        list: 삭제 기록 [(kept_folder, kept_file, deleted_folder, deleted_file), ...]
    """
    deletion_log = []
    total_duplicates = 0
    total_kept = 0

    for job_key, locations in job_map.items():
        protein_id, ligand_id = job_key

        if len(locations) <= 1:
            # 중복 없음
            total_kept += 1
            continue

        # 첫 번째 위치를 유지, 나머지 삭제
        kept_folder, kept_file, kept_path = locations[0]

        for i in range(1, len(locations)):
            del_folder, del_file, del_path = locations[i]

            deletion_log.append({
                'protein_id': protein_id,
                'ligand_id': ligand_id,
                'kept_folder': kept_folder,
                'kept_file': kept_file,
                'deleted_folder': del_folder,
                'deleted_file': del_file,
                'deleted_path': del_path
            })

            # 실제 파일 삭제
            if not dry_run:
                try:
                    os.remove(del_path)
                    total_duplicates += 1
                except Exception as e:
                    print(f"Error deleting {del_path}: {e}")
            else:
                total_duplicates += 1

        total_kept += 1

    return deletion_log, total_kept, total_duplicates


def write_log(deletion_log, log_path, total_kept, total_duplicates):
    """
    중복 제거 로그를 파일에 작성합니다.

    Args:
        deletion_log: 삭제 기록 리스트
        log_path: 로그 파일 경로
        total_kept: 유지된 파일 수
        total_duplicates: 삭제된 중복 파일 수
    """
    with open(log_path, 'w', encoding='utf-8') as f:
        # 헤더
        f.write("=" * 80 + "\n")
        f.write("Boltz-2 Job Deduplication Log\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 80 + "\n\n")

        # 요약
        f.write(f"Summary:\n")
        f.write(f"  - Unique jobs kept: {total_kept}\n")
        f.write(f"  - Duplicate jobs removed: {total_duplicates}\n\n")

        if not deletion_log:
            f.write("No duplicates found.\n")
            return

        f.write("=" * 80 + "\n")
        f.write("Detailed Deletion Log\n")
        f.write("=" * 80 + "\n\n")

        # 상세 로그
        for i, record in enumerate(deletion_log, 1):
            f.write(f"[{i}] Job: {record['protein_id']}_{record['ligand_id']}\n")
            f.write(f"    KEPT:    [{record['kept_folder']}] {record['kept_file']}\n")
            f.write(f"    DELETED: [{record['deleted_folder']}] {record['deleted_file']}\n")
            f.write(f"    Path: {record['deleted_path']}\n")
            f.write("\n")


def main():
    parser = argparse.ArgumentParser(
        description="여러 job 폴더에서 중복된 YAML 파일을 정리합니다.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
예시:
  python deduplicate_jobs.py folder1_job folder2_job folder3_job
  python deduplicate_jobs.py *_job --log dedup.log
  python deduplicate_jobs.py --dirs job1 job2 job3 --dry-run
        """
    )

    parser.add_argument(
        'folders',
        nargs='*',
        help='Job 폴더 경로들 (공백으로 구분)'
    )

    parser.add_argument(
        '--dirs',
        nargs='+',
        help='Job 폴더 경로들 (--dirs 옵션 사용 시)'
    )

    parser.add_argument(
        '--log',
        default='deduplication.log',
        help='로그 파일 경로 (기본값: deduplication.log)'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='실제로 파일을 삭제하지 않고 로그만 생성합니다'
    )

    args = parser.parse_args()

    # 폴더 경로 수집
    folder_paths = []
    if args.folders:
        folder_paths.extend(args.folders)
    if args.dirs:
        folder_paths.extend(args.dirs)

    if not folder_paths:
        print("Error: No job folders specified.")
        print("Usage: python deduplicate_jobs.py folder1 folder2 folder3 ...")
        print("       python deduplicate_jobs.py --dirs folder1 folder2 --log mylog.log")
        sys.exit(1)

    print(f"Scanning {len(folder_paths)} job folder(s)...")
    for folder in folder_paths:
        print(f"  - {folder}")
    print()

    # 1. job 폴더 스캔
    job_map = scan_job_folders(folder_paths)
    total_jobs = len(job_map)

    if total_jobs == 0:
        print("No valid YAML job files found in the specified folders.")
        sys.exit(0)

    print(f"Found {total_jobs} unique job(s) across all folders.\n")

    # 2. 중복 제거
    if args.dry_run:
        print("DRY RUN MODE: No files will be deleted.\n")

    deletion_log, total_kept, total_duplicates = deduplicate_jobs(job_map, dry_run=args.dry_run)

    # 3. 로그 작성
    write_log(deletion_log, args.log, total_kept, total_duplicates)

    # 4. 결과 출력
    print("=" * 80)
    print("Deduplication Complete")
    print("=" * 80)
    print(f"Unique jobs kept:        {total_kept}")
    print(f"Duplicate jobs removed:  {total_duplicates}")
    print(f"Log saved to:            {args.log}")

    if args.dry_run:
        print("\n*** DRY RUN: No files were actually deleted. ***")

    print()


if __name__ == "__main__":
    main()
