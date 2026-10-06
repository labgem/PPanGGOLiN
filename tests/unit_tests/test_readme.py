"""Tests for the README rendering on both GitHub and the Sphinx documentation."""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
README = ROOT / "README.md"

_spec = importlib.util.spec_from_file_location("github_readme", ROOT / "docs" / "_ext" / "github_readme.py")
github_readme = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(github_readme)


def test_readme_has_no_myst_only_syntax():
    """README.md is rendered by GitHub, so it must only use GitHub-flavored Markdown."""
    issues = github_readme.find_myst_only_syntax(README.read_text(encoding="utf-8"))
    assert not issues, "README.md uses MyST-only syntax:\n" + "\n".join(
        f"  line {line}: {issue}" for line, issue in issues
    )


@pytest.mark.parametrize("kind", ["NOTE", "TIP", "IMPORTANT", "WARNING", "CAUTION"])
def test_github_alert_to_myst(kind):
    text = f"Before\n\n> [!{kind}]\n> First line\n>\n> Second line\n\nAfter\n"
    expected = f"Before\n\n:::{{{kind.lower()}}}\nFirst line\n\nSecond line\n:::\n\nAfter\n"
    assert github_readme.github_alerts_to_myst(text) == expected


def test_github_alert_keeps_code_block():
    text = "> [!TIP]\n> Run:\n> ```shell\n> panorama -h\n> ```\n"
    expected = ":::{tip}\nRun:\n```shell\npanorama -h\n```\n:::\n"
    assert github_readme.github_alerts_to_myst(text) == expected


def test_regular_blockquote_untouched():
    text = "> Just a quote\n> on two lines\n"
    assert github_readme.github_alerts_to_myst(text) == text


@pytest.mark.parametrize(
    "line",
    [
        "```{note}",
        ":::{warning}",
        ":::",
        "See {ref}`user/install:Installation`",
        "(my-label)=",
        "% a MyST comment",
        "+++",
    ],
)
def test_myst_only_syntax_detected(line):
    issues = github_readme.find_myst_only_syntax(f"Intro\n\n{line}\n")
    assert [number for number, _ in issues] == [3]


def test_myst_syntax_inside_code_block_ignored():
    text = "````markdown\n```{note}\nText\n```\n{ref}`a`\n````\n\nUse `{ref}` roles in the docs.\n"
    assert github_readme.find_myst_only_syntax(text) == []
