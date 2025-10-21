import filecmp
import difflib
import os

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
from pathlib import Path

from collections import Counter
import json

import gzip


def ordered(obj):
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    else:
        return obj


def read_json_file(json_file):
    proper_open_1 = gzip.open if json_file.suffix == ".gz" else open

    with proper_open_1(json_file.as_posix(), "rt") as f1:
        return json.load(f1)


def are_json_files_identical(file1, file2):

    data1 = read_json_file(file1)
    data2 = read_json_file(file2)
    # Load data from the second file

    return ordered(data1) == ordered(data2), []


def are_text_files_identical(expected_file, file2, outdir=Path("./")):
    diff_details = []
    common_file = expected_file.name

    if filecmp.cmp(expected_file, file2, shallow=False):
        return True, diff_details

    # Compare file content

    proper_open_e = gzip.open if expected_file.suffix == ".gz" else open

    proper_open_t = gzip.open if file2.suffix == ".gz" else open
    with proper_open_e(expected_file, "rt") as f1, proper_open_t(file2, "rt") as f2:

        f1_content, f2_content = sorted(f1.readlines()), sorted(f2.readlines())
        f1_line_count = len(f1_content)
        f2_line_count = len(f2_content)

        diff = difflib.unified_diff(
            f1_content,
            f2_content,
            fromfile=expected_file.as_posix(),
            tofile=file2.as_posix(),
        )
        diff = [line.rstrip() for line in diff]
        if len(diff) == 0:
            return True, diff_details

        diff_event_counter = {"+": 0, "-": 0}

        diff_file = outdir / f"{common_file}.diff"
        with open(diff_file, "w") as out:
            out.write("\n".join(diff) + "\n")

        # gather stat on diff
        diff_event_counter = Counter(
            (line[0] for line in diff[2:] if line[0] in ["+", "-"])
        )

        prct_inserted = 100 * diff_event_counter["+"] / f1_line_count
        prct_deleted = 100 * diff_event_counter["-"] / f2_line_count
        diff_prct = max(prct_inserted, prct_deleted)
        diff_details.append(
            f'{diff_event_counter["+"]} insertions(+), {diff_event_counter["-"]} deletions(-) - {diff_prct:.1f}% difference'
        )
        diff_details.append(f"Check out diff in {diff_file}")
        # vs code command to diff the files
        # diff_details.append(f'code --diff {expected_file} {file2}')

    return False, diff_details


def add_subdir_to_files(subdir, files):
    return [os.path.join(subdir, file) for file in files]


def compare_dir_recursively(expected_dir, tested_dir, ignored_files):
    dcmp_result = filecmp.dircmp(expected_dir, tested_dir, ignore=ignored_files)

    for subdir in dcmp_result.common_dirs:
        sub_dcmp_result = compare_dir_recursively(
            expected_dir / subdir, tested_dir / subdir, ignored_files
        )
        for files_attr in [
            "right_only",
            "left_only",
            "common_files",
            "same_files",
            "diff_files",
        ]:
            files_list = getattr(dcmp_result, files_attr)
            files_list += add_subdir_to_files(
                subdir, getattr(sub_dcmp_result, files_attr)
            )

    return dcmp_result


def get_suffix_except_gz(path: Path):
    for ext in path.suffixes[::-1]:
        if ext != ".gz":
            return ext


def compare_directories(
    expected_dir,
    tested_dir,
    ignored_files,
    diff_outdir,
    report_identical_files,
    extension_to_compare=[
        ".tsv",
        ".aln",
        ".json",
        ".gff",
        ".txt",
        ".csv",
        ".faa",
        ".fasta",
        ".yaml",
    ],
):
    # Define directory information with color
    expected_dir_info = f"- Expected directory: {expected_dir}"
    tested_dir_info = f"- Tested directory: {tested_dir}"

    # Create the panel
    print("\n===Comparison of Directories===")
    print("\n".join([expected_dir_info, tested_dir_info]))
    print("===============================")

    # Compare directories
    dcmp = compare_dir_recursively(
        expected_dir, tested_dir, ignored_files=ignored_files
    )
    ignored_files_ext = [
        common_file
        for common_file in dcmp.common_files
        if get_suffix_except_gz(Path(common_file)) not in extension_to_compare
    ]
    if ignored_files:
        print("\nFiles ignored for comparison:")
        for ignored_file in ignored_files:
            print(f"  - {ignored_file}")

    if ignored_files_ext:
        print("\nFiles ignored for comparison due to their extensions:")
        for ignored_file in ignored_files_ext:
            print(f"  - {ignored_file}")

    # Report missing files
    if dcmp.left_only:
        print("\nMissing files or directories in Tested results:")
        for file in dcmp.left_only:
            print(f"  - {file}")

    if dcmp.right_only:
        print("\nUnexpected files or directories found in Tested results:")
        for file in dcmp.right_only:
            print(f"  - {file}")

    different_files = []
    identical_files = []

    files_to_compare = list(set(dcmp.common_files) - set(ignored_files_ext))

    for common_file in files_to_compare:
        file1 = expected_dir / common_file
        file2 = tested_dir / common_file

        if get_suffix_except_gz(Path(common_file)) == ".json":
            identical_tables, details = are_json_files_identical(file1, file2)
        else:
            identical_tables, details = are_text_files_identical(
                file1, file2, outdir=diff_outdir
            )

        if identical_tables:
            identical_files.append(common_file)
        else:
            different_files.append((common_file, details))

    if report_identical_files and identical_files:
        print("\nIdentical files:")
        for file in identical_files:
            print(f"  - {file}")

    if different_files:
        print("\nDifferent files:")
        for file, details in different_files:
            print(f"  - {file}")
            for detail in details:
                print(f"{detail}")
            print()

    # Generate summary report
    print("\nSummary:")

    # Display identical files count
    print(f"{len(identical_files)} file(s) identical.")

    # Display ignored files due to extension
    print(f"{len(ignored_files_ext)} file(s) ignored due to their extension.")

    # Display unexpected files in Tested results
    print(f"{len(dcmp.right_only)} file(s) unexpectedly found in Tested results.")
    print(f"{len(dcmp.left_only)} file(s) missing in Tested results.")

    # Display different files count
    print(f"{len(different_files)} file(s) differ.")

    if different_files or dcmp.left_only or dcmp.right_only:
        print(
            "\nSome difference exist between the tested and the expected result directories"
        )
        exit(1)

    else:
        print(
            "\nNo difference have been found between the tested and the expected result directories"
        )


def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(
        description="...", formatter_class=ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-e",
        "--expected_dir",
        help="Expected result directory",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "-t", "--tested_dir", help="Tested result directory", required=True, type=Path
    )

    parser.add_argument(
        "-o",
        "--outdir",
        help="Directories where to write diff files",
        default="out_diff",
        type=Path,
    )

    parser.add_argument(
        "-i",
        "--ignored_files",
        nargs="+",
        help="File to ignore for the comparison",
        default=[
            "pangenomeGraph.json",
            "pangenomeGraph.json.gz",
            "expected_hashes.json",
        ],
    )

    args = parser.parse_args()
    return args


def main():

    args = parse_arguments()

    report_identical_files = True

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    Path.mkdir(args.outdir, exist_ok=True)

    compare_directories(
        expected_dir=args.expected_dir,
        tested_dir=args.tested_dir,
        diff_outdir=args.outdir,
        ignored_files=args.ignored_files,
        report_identical_files=report_identical_files,
    )


if __name__ == "__main__":
    main()
