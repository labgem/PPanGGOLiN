import importlib.util
from pathlib import Path

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
    assert "--cpu" in rendered
    assert "--pangenome" in rendered
    assert "False" in rendered
    assert "Choices" in rendered
