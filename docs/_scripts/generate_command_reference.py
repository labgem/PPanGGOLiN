#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
from typing import Iterable, List

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from ppanggolin.main import build_parser


def _format_type(action: argparse.Action) -> str:
    if action.nargs in (None, "?"):
        if isinstance(action.type, type):
            return action.type.__name__
        if action.type is not None:
            return getattr(action.type, "__name__", str(action.type))
        if action.const is not None and action.default is False:
            return "bool"
        if action.default is not None:
            return type(action.default).__name__
        return "str"

    if action.nargs == "*":
        return "list"
    if action.nargs == "+":
        return "list"
    if action.nargs == "?":
        return "optional"
    return str(action.nargs)


def _format_default(action: argparse.Action) -> str:
    if action.default is None:
        return "—"
    if action.default is False and action.option_strings:
        return "False"
    if action.default is True and action.option_strings:
        return "True"
    if isinstance(action.default, Path):
        return str(action.default)
    if isinstance(action.default, (list, tuple, set)):
        return str(list(action.default))
    if isinstance(action.default, str):
        return f"`{action.default}`"
    return str(action.default)


def _format_choices(action: argparse.Action) -> str:
    if not action.choices:
        return "—"
    return ", ".join(str(choice) for choice in action.choices)


def _format_param_name(action: argparse.Action) -> str:
    names = [opt for opt in action.option_strings]
    if not names:
        if action.metavar:
            return action.metavar
        return action.dest
    return ", ".join(names)


def _describe_action(action: argparse.Action) -> dict:
    return {
        "name": _format_param_name(action),
        "type": _format_type(action),
        "default": _format_default(action),
        "required": "Yes" if action.required else "No",
        "choices": _format_choices(action),
        "help": (action.help or "").replace("\n", " ").strip(),
        "is_flag": action.option_strings and action.nargs == 0,
    }


def _iter_actions(parser: argparse.ArgumentParser) -> Iterable[argparse.Action]:
    seen = set()
    for action_group in parser._action_groups:
        for action in action_group._group_actions:
            if action.dest == "help":
                continue
            if action.dest in seen:
                continue
            seen.add(action.dest)
            yield action


def _generate_command_section(
    command_name: str, parser: argparse.ArgumentParser
) -> str:
    title = f"## `ppanggolin {command_name}`"
    description = (
        parser.description
        or "Command description is not explicitly set in this parser."
    )
    lines = [
        title,
        "",
        description.strip(),
        "",
        "### Parameters",
        "",
        "| Parameter | Type | Default | Required | Choices | Description |",
        "|---|---|---|---|---|---|",
    ]

    actions = list(_iter_actions(parser))
    if not actions:
        lines.append("No parameters are defined for this command.")
        return "\n".join(lines) + "\n"

    for action in actions:
        info = _describe_action(action)
        help_text = info["help"]
        if info["is_flag"] and help_text:
            help_text = help_text
        if not help_text:
            help_text = "—"
        lines.append(
            f"| `{info['name']}` | {info['type']} | {info['default']} | {info['required']} | {info['choices']} | {help_text} |"
        )

    return "\n".join(lines) + "\n"


def _collect_subcommands(
    parser: argparse.ArgumentParser,
) -> List[tuple[str, argparse.ArgumentParser]]:
    subparsers_action = next(
        (
            group
            for group in parser._action_groups
            if isinstance(group, argparse._ArgumentGroup)
        ),
        None,
    )
    subparsers = []
    for action in getattr(parser, "_actions", []):
        if isinstance(action, argparse._SubParsersAction):
            subparsers = action.choices.items()
            break

    commands = []
    for command_name, command_parser in subparsers:
        commands.append((command_name, command_parser))
    commands.sort(key=lambda item: item[0])
    return commands


def generate_reference(output_path: Path | str | None = None) -> str:
    parser = build_parser()
    commands = _collect_subcommands(parser)

    sections = [
        "# Command-line reference",
        "",
        "Complete reference of the PPanGGOLiN command-line interface.",
        "",
        "This page is automatically generated from the PPanGGOLiN command-line parser.",
        "",
    ]

    for command_name, command_parser in commands:
        sections.append(_generate_command_section(command_name, command_parser))
        sections.append("")

    rendered = "\n".join(sections).strip() + "\n"

    if output_path is not None:
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(rendered, encoding="utf-8")
    return rendered


if __name__ == "__main__":
    script_path = Path(__file__).resolve().parents[1] / "user" / "command_reference.md"
    generate_reference(script_path)
