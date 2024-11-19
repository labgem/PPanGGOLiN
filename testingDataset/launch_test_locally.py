#!/usr/bin/env python3

"""
Description

:Example:
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
from pathlib import Path
import yaml


def parse_yaml_file(yaml_file):

    # Load the YAML file
    with open(yaml_file, "r") as stream:

        workflow = yaml.safe_load(stream)

    return workflow


def create_symbolic_links(source_dir, target_dir):
    # Convert paths to Path objects
    source_dir = Path(source_dir)
    target_dir = Path(target_dir)

    # Ensure target directory exists
    target_dir.mkdir(parents=True, exist_ok=True)

    # Loop through each item in the source directory
    for item in source_dir.iterdir():
        target_item_path = target_dir / item.name

        # Create symbolic link for each item
        try:
            target_item_path.symlink_to(item.resolve())
        except FileExistsError:
            pass


def parse_arguments(default_ci_yaml, testing_datadir):
    """Parse script arguments."""
    parser = ArgumentParser(
        description="...", formatter_class=ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--ci_yaml",
        help="increase output verbosity",
        default=default_ci_yaml,
        type=Path,
    )

    parser.add_argument(
        "--data_dir",
        help="Directory where dataset files are located",
        default=testing_datadir,
        type=Path,
    )

    parser.add_argument(
        "-o",
        "--outdir",
        help="increase output verbosity",
        default="local_CI",
        type=Path,
    )

    parser.add_argument(
        "-c",
        "--cpu",
        type=int,
        default=4,
        help="Use this amount of cpu when number of cpu is specified in the command.",
    )

    parser.add_argument(
        "-v", "--verbose", help="increase output verbosity", action="store_true"
    )

    parser.add_argument(
        "--skip_msa",
        action="store_true",
        help="Skip msa command as it takes quite some time to complete.",
    )

    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Force writing in output directory if exists.",
    )

    args = parser.parse_args()
    return args


def main():

    script_path = Path(__file__).resolve()
    ppanggolin_main_dir = script_path.parent.parent
    default_ci_yaml = ppanggolin_main_dir / ".github/workflows/main.yml"
    testing_datadir = ppanggolin_main_dir / "testingDataset"

    args = parse_arguments(default_ci_yaml, testing_datadir)

    if args.verbose:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
        logging.info("Mode verbose ON")

    else:
        logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)

    if not args.outdir.is_dir():
        logging.debug(f"Create output directory {args.outdir.absolute().as_posix()}")
        Path.mkdir(args.outdir)

    elif not args.force:
        raise FileExistsError(
            f"{args.outdir} already exists. Use -f if you want to overwrite the files in the directory"
        )

    # setup test dir
    # execution_dir = args.outdir # / args.data_dir.name
    create_symbolic_links(args.data_dir, args.outdir)

    workflow = parse_yaml_file(args.ci_yaml)

    excluded_steps = [
        "Install ppanggolin",
        "Get core number on linux",
        "Get core number on macos",
    ]

    test_script = args.outdir / "launch_test_command.sh"
    with open(test_script, "w") as fl:

        fl.write("#!/bin/bash\nset -e\n")

        # Iterate over jobs and steps
        for job in workflow["jobs"].values():

            if "steps" not in job:
                continue

            for step in job["steps"]:
                if "run" not in step:
                    continue

                if step["name"] in excluded_steps:
                    logging.info(f"Ignoring: {step}")
                    continue

                # Execute the command line
                command = step["run"].strip()

                # process the command
                command = command.replace("$NUM_CPUS", f"{args.cpu}")
                command = command.replace("cd ", "# cd ")

                if args.skip_msa:
                    command = command.replace("ppanggolin msa", "# ppanggolin msa")
                if command == "pytest":
                    command = f"pytest {ppanggolin_main_dir}"

                # log the step name
                logging.info(f'Executing: {step["name"]}')
                logging.debug(f" {command}\n")

                # write the command in the script
                fl.write(f'\n# {step["name"]}\n\n')
                fl.write(command)

    print(f"(cd {args.outdir}; bash {test_script.name})")


if __name__ == "__main__":
    main()
