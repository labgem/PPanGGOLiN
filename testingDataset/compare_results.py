import filecmp
import difflib
import os

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
from pathlib import Path

from rich.console import Console
from rich.syntax import Syntax
from rich.panel import Panel

from rich.logging import RichHandler

from rich import box
from rich.theme import Theme
from collections import Counter
import json

custom_theme = Theme({
    "identical": "bold green",
    "neutral": "bold yellow",
    "difference": "bold red"
})
console = Console(theme=custom_theme)


def are_json_files_identical(file1, file2):
    a = json.dumps(file1.as_posix(), sort_keys=True)
    b = json.dumps(file2.as_posix(), sort_keys=True)
    return a == b,[] # a normal string comparison


def are_text_files_identical(expected_file, file2, outdir=Path("./")):
    diff_details = []
    common_file = expected_file.name
    # max_lines_to_display = 10

    if filecmp.cmp(expected_file, file2, shallow=False):
        return True, diff_details


    # console.print(f"\n[bold red]{common_file} is different.[/bold red].")

    # Compare file content
    with open(expected_file, 'r') as f1, open(file2, 'r') as f2:

        f1_content, f2_content = sorted(f1.readlines()), sorted(f2.readlines())
        f1_line_count = len(f1_content)
        f2_line_count = len(f2_content)

        diff = difflib.unified_diff(f1_content, f2_content, fromfile=expected_file.as_posix(), tofile=file2.as_posix())
        diff = [line.rstrip() for line in diff]
        if len(diff) == 0:
            return True, diff_details
        # diff_lines_to_display = []
        diff_event_counter = {"+":0, "-":0}
        # truncated_diff_lines = 0
        
        diff_file = outdir / f"{common_file}.diff"
        with open(diff_file, 'w') as out:
            out.write('\n'.join(diff) + '\n')

        # gather stat on diff
        diff_event_counter = Counter((line[0] for line in diff[2:] if line[0] in ['+', '-']))
        
        prct_inserted = 100 * diff_event_counter["+"] / f1_line_count
        prct_deleted = 100 * diff_event_counter["-"] / f2_line_count
        diff_prct = max(prct_inserted, prct_deleted)
        diff_details.append(f'{diff_event_counter["+"]} insertions(+), {diff_event_counter["-"]} deletions(-) - {diff_prct:.1f}% difference')
        diff_details.append(f'code {diff_file}')
        diff_details.append(f'code --diff {expected_file} {file2}')
        # display=False
        # if display:
        #     # Display the n first line of diff with rich
        #     diff_lines_to_display = diff[:max_lines_to_display]
        #     truncated_diff_lines = len(diff_lines_to_display) - len(diff)
        
        #     if truncated_diff_lines:
        #         diff_lines_to_display.append(f'... {truncated_diff_lines} more diff lines ....')

        #     diff_content = '\n'.join( diff_lines_to_display ) 

        #     syntax = Syntax(diff_content, "diff", line_numbers=False, theme="ansi_dark")
        #     console.print(Panel(syntax, title=f"Differences in file {common_file}", expand=True, box=box.HORIZONTALS), overflow="ellipsis")
    
    return False, diff_details

def add_subdir_to_files(subdir, files):
    return [os.path.join(subdir, file) for file in  files]

def compare_dir_recursively(expected_dir, tested_dir, ignored_files):
    dcmp_result = filecmp.dircmp(expected_dir, tested_dir)

    for subdir  in dcmp_result.common_dirs:
        sub_dcmp_result = compare_dir_recursively(expected_dir/subdir, tested_dir/subdir, ignored_files)
        for files_attr in ['right_only', 'left_only', 'common_files', "same_files", "diff_files"]:
            files_list = getattr(dcmp_result, files_attr)
            files_list += add_subdir_to_files(subdir, getattr(sub_dcmp_result, files_attr))
            
        
    return dcmp_result

def compare_directories(expected_dir, tested_dir, ignored_files, diff_outdir, extension_to_compare=[".tsv"]): #, ".json", '.gff', '.txt', '.csv']):

    # Define directory information with color
    expected_dir_info = f"- Expected directory: [bold green]{expected_dir}[/bold green]"
    tested_dir_info = f"- Tested directory: [bold blue]{tested_dir}[/bold blue]"

    # Create the panel
    panel_title = "[bold cyan]Comparison of Directories:[/bold cyan]"
    panel_content = "\n".join([expected_dir_info, tested_dir_info])
    comparison_panel = Panel(panel_content, title=panel_title, title_align="center", expand=False)
    console.print(comparison_panel)

    # Compare directories
    dcmp = compare_dir_recursively(expected_dir, tested_dir, ignored_files=ignored_files)
    ignored_files_ext = [common_file for common_file in dcmp.common_files if Path(common_file).suffix not in extension_to_compare]

    if ignored_files:
        console.print("\nFiles ignored for comparison:", style="bold")
        for ignored_file in ignored_files:
            console.print(f"  - {ignored_file}")

    if ignored_files_ext:
        console.print("\nFiles ignored for comparison due to their extensions:", style="bold")
        for ignored_file in ignored_files_ext:
            console.print(f"  - {ignored_file}")

    # Report missing files
    if dcmp.left_only:
        console.print("\nMissing files or directories in Tested results:", style='neutral')
        for file in dcmp.left_only:
            console.print(f"  - {file}")

    if dcmp.right_only:
        console.print("\nUnexpected files or directories found in Tested results:", style='neutral')
        for file in dcmp.right_only:
            console.print(f"  - {file}")

    different_files = []
    identical_files = []

    
    files_to_compare = list( set(dcmp.common_files) - set(ignored_files_ext))

    for common_file in files_to_compare:
        file1 = expected_dir / common_file
        file2 = tested_dir / common_file

        # if Path(common_file).suffix == ".json":
        #     identical_tables, details = are_json_files_identical(file1, file2)
        # else:
        identical_tables, details = are_text_files_identical(file1, file2, outdir=diff_outdir)

        if identical_tables:
            identical_files.append(common_file)
        else:
            different_files.append((common_file, details))
    
    # if identical_files:
    #     console.print("\nIdentical files:", style='identical')
    #     for file in identical_files:
    #         console.print(f"  - {file}")

    if different_files:
        console.print("\nDifferent files:", style='difference')
        for file, details  in different_files:
            console.print(f"  - {file}")
            for detail in details:
                console.print(f"{detail}")
            print()

    # Generate summary report
    console.print("\n[bold]Summary:[/bold]")

    console.print(f"{len(identical_files)} file(s) identical.")

    console.print(f"{len(ignored_files_ext)} file(s) ignored because of their extension.")

    console.print(f"{len(dcmp.left_only)} file(s) unexpected in Tested results.", )
    console.print(f"{len(dcmp.right_only)} file(s) missing in Tested results.")

    console.print(f"{len(different_files)} file(s) different.")



    # Suggest updating expected directory if changes seem legitimate
    # if missing_files or different_files:
    #     code_snippet_update = f"\ncp -r {tested_dir}/* {expected_dir}\n"

    #     console.print("\nIf changes in tested directory seem legitimate, please update the expected directory with the following command:")
    #     syntax = Syntax(code_snippet_update, "bash", theme="ansi_dark")
    #     console.print(Panel(syntax, title="Update expected results", expand=False, box=box.HORIZONTALS))


    if different_files or dcmp.left_only or dcmp.right_only:
        console.print('\nSome difference exist between the tested and the expected result directories', style="difference")
        exit(1)

    else:
        console.print('\nNo difference have been found between the tested and the expected result directories', style="identical")

def parse_arguments():
    """Parse script arguments."""
    parser = ArgumentParser(description="...",
                            formatter_class=ArgumentDefaultsHelpFormatter)


    parser.add_argument('-e', '--expected_dir', help="Expected result directory", required=True, type=Path)
    
    parser.add_argument('-t', '--tested_dir', help="Tested result directory", required=True, type=Path)

    parser.add_argument('-o', '--outdir', help="Directories where to write diff files", default='out_diff', type=Path)

    parser.add_argument('-i', '--ignored_files', nargs="+", help="File to ignore for the comparison", default=['pangenomeGraph.json'])


    
    args = parser.parse_args()
    return args

def main():

    args = parse_arguments()

    logging.basicConfig(level=logging.INFO, format="%(message)s", datefmt="[%X]", handlers=[RichHandler()])
    
    Path.mkdir(args.outdir, exist_ok=True)


    compare_directories(expected_dir=args.expected_dir, tested_dir=args.tested_dir, diff_outdir=args.outdir, ignored_files=args.ignored_files)

if __name__ == '__main__':
    main()

