#!/usr/bin/env python3

import argparse
import re
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
        value = str(action.default)
        pattern = r"^(ppanggolin_[A-Za-z0-9_]+)DATE\d{4}-\d{2}-\d{2}_HOUR\d{2}\.\d{2}\.\d{2}_PID\d+$"
        match = re.match(pattern, value)
        if match:
            prefix = match.group(1)
            return f"`{prefix}<date>_<pid>`"
        return str(action.default)
    if isinstance(action.default, (list, tuple, set)):
        return str(list(action.default))
    if isinstance(action.default, str):
        return f"`{action.default}`"
    return str(action.default)


def _format_choices(action: argparse.Action) -> str:
    if not action.choices:
        return "—"
    return ", ".join(f"`{choice}`" for choice in action.choices)


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


def _iter_action_groups(
    parser: argparse.ArgumentParser,
) -> Iterable[tuple[str, list[argparse.Action]]]:
    seen = set()
    for action_group in parser._action_groups:
        group_title = action_group.title or "Arguments"
        if group_title in {"positional arguments", "options"}:
            continue

        actions = []
        for action in action_group._group_actions:
            if action.dest == "help":
                continue
            if action.dest in seen:
                continue
            seen.add(action.dest)
            actions.append(action)

        if actions:
            yield group_title, actions


def _generate_action_table(
    title: str, actions: list[argparse.Action], command_name: str | None = None
) -> list[str]:
    heading = (
        title if command_name is None else f"{title} for ppanggolin {command_name}"
    )
    lines = [
        f"#### {heading}",
        "",
        "| Parameter | Type | Default | Description |",
        "|---|---|---|---|",
    ]

    for action in actions:
        info = _describe_action(action)
        description_parts = []
        help_text = info["help"].strip() if info["help"] else "—"
        description_parts.append(help_text)

        if info["required"] == "Yes":
            description_parts.append("Required: Yes")
        if info["choices"] != "—":
            description_parts.append(f"<br>Choices: {info['choices']}")

        description = " ".join(part for part in description_parts if part)
        lines.append(
            f"| `{info['name']}` | {info['type']} | {info['default']} | {description} |"
        )

    return lines


def _generate_command_section(
    command_name: str, parser: argparse.ArgumentParser
) -> str:
    title = f"### `ppanggolin {command_name}`"
    description = (
        getattr(parser, "description", None)
        or "Command description is not explicitly set in this parser."
    )
    lines = [
        title,
        "",
        description.strip(),
        "",
    ]

    action_groups = list(_iter_action_groups(parser))
    if not action_groups:
        lines.append("No parameters are defined for this command.")
        return "\n".join(lines) + "\n"

    for group_title, actions in action_groups:
        lines.extend(_generate_action_table(group_title, actions, command_name))
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


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

    category_order = [
        "Basic",
        "Expert",
        "Output",
        "Regions of Genomic Plasticity",
        "Analysis using reference pangenomes",
        "Utility command",
    ]
    grouped_commands: dict[str, List[tuple[str, argparse.ArgumentParser]]] = {}
    for command_name, command_parser in commands:
        category = getattr(command_parser, "category", "General")
        grouped_commands.setdefault(category, []).append((command_name, command_parser))

    for category in category_order:
        if category not in grouped_commands:
            continue
        sections.append(f"## {category}")
        sections.append("")
        for command_name, command_parser in grouped_commands[category]:
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
