"""Sphinx extension to reuse the GitHub README as the documentation landing page.

README.md is written for GitHub. When it is pulled into the docs with MyST's ``{include}``
directive, this extension:

- rewrites GitHub alerts (``> [!NOTE]``) into MyST admonitions;
- warns about MyST-only syntax, which GitHub would display as raw text.

The other GitHub-flavored features (bare URL autolinks, strikethrough, task lists, mermaid
fences) are handled by MyST extensions enabled in ``conf.py``.

This module must stay importable without Sphinx so the README can be checked by the test suite.
"""

import re
from pathlib import Path

README_NAME = "README.md"

GITHUB_ALERT = re.compile(
    r"^> \[!(?P<kind>NOTE|TIP|IMPORTANT|WARNING|CAUTION)\][ \t]*\n(?P<body>(?:^>.*(?:\n|$))*)",
    flags=re.MULTILINE,
)

FENCE = re.compile(r"^ {0,3}(?P<marker>`{3,}|~{3,}|:{3,})(?P<info>.*)$")
# Scanned left to right, so a role written inside an inline code span is consumed as code
MYST_ROLE_OR_CODE = re.compile(r"(?P<role>\{[\w:.+-]+\}`[^`]+`)|`+[^`]*`+")
MYST_INLINE_SYNTAX = {
    "MyST target": re.compile(r"^\s*\([^()\s]+\)=\s*$"),
    "MyST comment": re.compile(r"^\s*%"),
    "MyST block break": re.compile(r"^\s*\+\+\+"),
}


def github_alerts_to_myst(text: str) -> str:
    """Convert GitHub alert blockquotes into MyST colon-fence admonitions.

    Args:
        text: Markdown content using GitHub alert syntax.

    Returns:
        The same content with every alert replaced by a MyST admonition.
    """

    def _replace(match: re.Match) -> str:
        body = re.sub(r"^> ?", "", match["body"], flags=re.MULTILINE).rstrip("\n")
        return f":::{{{match['kind'].lower()}}}\n{body}\n:::\n"

    return GITHUB_ALERT.sub(_replace, text)


def find_myst_only_syntax(text: str) -> list[tuple[int, str]]:
    """Find MyST syntax that GitHub cannot render.

    Content of regular code blocks is ignored, so documenting MyST syntax inside a code block is fine.

    Args:
        text: Markdown content meant to be rendered by GitHub.

    Returns:
        A list of ``(line_number, description)`` tuples, one per offending line (1-indexed).
    """
    issues = []
    open_fence = None
    for line_number, line in enumerate(text.splitlines(), start=1):
        fence = FENCE.match(line)
        if open_fence is not None:
            if fence and fence["marker"][0] == open_fence[0] and len(fence["marker"]) >= len(open_fence):
                if not fence["info"].strip():
                    open_fence = None
            continue
        if fence:
            marker, info = fence["marker"], fence["info"].strip()
            if marker.startswith(":"):
                issues.append((line_number, f"MyST colon fence '{line.strip()}'"))
            elif info.startswith("{"):
                issues.append((line_number, f"MyST directive '{line.strip()}'"))
            open_fence = marker
            continue
        if any(match["role"] for match in MYST_ROLE_OR_CODE.finditer(line)):
            issues.append((line_number, f"MyST role in '{line.strip()}'"))
        for name, pattern in MYST_INLINE_SYNTAX.items():
            if pattern.match(line):
                issues.append((line_number, f"{name} '{line.strip()}'"))
    return issues


def convert_included_readme(app, relative_path: Path, parent_docname: str, content: list[str]) -> None:
    """Sphinx ``include-read`` handler converting the included README from GitHub to MyST syntax."""
    if relative_path.name != README_NAME:
        return

    from sphinx.util import logging

    logger = logging.getLogger(__name__)
    for line_number, issue in find_myst_only_syntax(content[0]):
        logger.warning(
            "%s will not render on GitHub, use GitHub-flavored Markdown instead",
            issue,
            location=f"{(Path(app.srcdir) / relative_path).resolve()}:{line_number}",
        )
    content[0] = github_alerts_to_myst(content[0])


def setup(app):
    """Register the extension in Sphinx."""
    app.connect("include-read", convert_included_readme)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
