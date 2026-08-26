import argparse
import importlib.util
from pathlib import Path

from ppanggolin import SUBCOMMAND_TO_SUBPARSER
from ppanggolin.main import build_parser


def test_build_parser_is_available_without_sys_argv():
    parser = build_parser()

    assert parser is not None
    assert parser.prog == "ppanggolin"
    assert "all" in parser._subparsers._group_actions[0].choices
    assert "info" in parser._subparsers._group_actions[0].choices
    assert "utils" in parser._subparsers._group_actions[0].choices


def test_build_parser_does_not_execute_analysis(monkeypatch):
    import ppanggolin.workflow.all as workflow_all

    def boom(*args, **kwargs):
        raise AssertionError(
            "analysis launch should not be called while building the parser"
        )

    monkeypatch.setattr(workflow_all, "launch", boom)

    parser = build_parser()
    assert parser is not None


def test_subparsers_own_command_metadata():
    parser = argparse.ArgumentParser(prog="ppanggolin")
    subparsers = parser.add_subparsers(dest="subcommand")

    all_parser = SUBCOMMAND_TO_SUBPARSER["all"](subparsers)
    info_parser = SUBCOMMAND_TO_SUBPARSER["info"](subparsers)

    assert all_parser.category == "Basic"
    assert all_parser.description == "Easy workflow to run all possible analysis"
    assert info_parser.category == "Output"
    assert (
        info_parser.description
        == "Prints information about a given pangenome graph file."
    )


def test_top_level_help_lists_subcommands():
    parser = build_parser()
    help_text = parser.format_help()

    assert "subcommands:" in help_text
    assert "all" in help_text
    assert "info" in help_text
    assert "utils" in help_text
    assert "All of the following subcommands have their own set of options" in help_text
    assert "Basic:" in help_text


def test_generated_reference_contains_cli_data(tmp_path):
    script_path = (
        Path(__file__).resolve().parents[2]
        / "docs"
        / "_scripts"
        / "generate_command_reference.py"
    )
    spec = importlib.util.spec_from_file_location(
        "generate_command_reference", script_path
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    output_path = tmp_path / "command_reference.md"
    rendered = module.generate_reference(output_path)

    assert "# Command-line reference" in rendered
    assert "## `ppanggolin all`" in rendered
    assert "## `ppanggolin info`" in rendered
    assert "## `ppanggolin utils`" in rendered
    assert "## Basic" in rendered
    assert "Easy workflow to run all possible analysis" in rendered
    assert "## Utility command" in rendered
    assert "Helper side commands" in rendered
    assert "--cpu" in rendered
    assert "--pangenome" in rendered
    assert "False" in rendered
    assert "Choices" in rendered
    assert "bool" in rendered
    assert "float" in rendered
    assert "Input arguments for ppanggolin all" in rendered
    assert "Optional arguments for ppanggolin all" in rendered
    assert "Common arguments for ppanggolin utils" in rendered
    assert "ppanggolin_context_<date>_<pid>" in rendered
