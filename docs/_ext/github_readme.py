"""Sphinx extension to reuse the GitHub README as the documentation landing page.

README.md is written for GitHub. When it is pulled into the docs with MyST's ``{include}``
directive, this extension replaces the directive by the README content and:

- rewrites GitHub alerts (``> [!NOTE]``) into MyST admonitions;
- warns about MyST-only syntax, which GitHub would display as raw text.

The other GitHub-flavored features (bare URL autolinks, strikethrough, task lists, mermaid
fences) are handled by MyST extensions enabled in ``conf.py``.

This module must stay importable without Sphinx so the README can be checked by the test suite.
"""

import os
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

MYST_INCLUDE = re.compile(
    r"^(?P<fence>`{3,}|:{3,})\{include\}(?P<path>[^\n]+)\n(?P<options>(?::[^\n]*\n)*)(?P=fence)[ \t]*$",
    flags=re.MULTILINE,
)
MARKDOWN_IMAGE = re.compile(
    r"(?P<alt>!\[[^\]]*\])\((?P<target>[^)\s]+)(?P<title>[^)]*)\)"
)


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
            if (
                fence
                and fence["marker"][0] == open_fence[0]
                and len(fence["marker"]) >= len(open_fence)
            ):
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


def rebase_images(text: str, readme_dir: Path, doc_dir: Path) -> str:
    """Make relative Markdown image paths of the README relative to the including document.

    Reproduces the ``:relative-images:`` option of MyST's ``{include}`` directive.

    Args:
        text: README content.
        readme_dir: Directory containing the README.
        doc_dir: Directory containing the document that includes the README.

    Returns:
        The content with every relative image path rebased on ``doc_dir``.
    """

    def _replace(match: re.Match) -> str:
        target = match["target"]
        if re.match(r"^([a-z][a-z0-9+.-]*:|/|#)", target, flags=re.IGNORECASE):
            return match[0]
        rebased = Path(os.path.relpath(readme_dir / target, doc_dir)).as_posix()
        return f"{match['alt']}({rebased}{match['title']})"

    return MARKDOWN_IMAGE.sub(_replace, text)


def expand_readme_include(app, docname: str, source: list[str]) -> None:
    """Sphinx ``source-read`` handler inlining the README converted from GitHub to MyST syntax.

    MyST's ``{include}`` directive reads the file itself, so the README cannot be converted when it
    is included. The directive is instead replaced by the converted README content before parsing.
    """
    from sphinx.util import logging

    logger = logging.getLogger(__name__)
    doc_path = Path(app.env.doc2path(docname))

    def _replace(match: re.Match) -> str:
        readme_path = (doc_path.parent / match["path"].strip()).resolve()
        if readme_path.name != README_NAME:
            return match[0]
        app.env.note_dependency(str(readme_path))
        text = readme_path.read_text(encoding="utf-8")
        for line_number, issue in find_myst_only_syntax(text):
            logger.warning(
                "%s will not render on GitHub, use GitHub-flavored Markdown instead",
                issue,
                location=f"{readme_path}:{line_number}",
            )
        text = github_alerts_to_myst(text)
        if re.search(r"^:relative-images:", match["options"], flags=re.MULTILINE):
            text = rebase_images(text, readme_path.parent, doc_path.parent)
        return text

    source[0] = MYST_INCLUDE.sub(_replace, source[0])


def setup(app):
    """Register the extension in Sphinx."""
    app.connect("source-read", expand_readme_include)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
