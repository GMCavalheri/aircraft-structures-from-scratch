"""Small helpers shared by the validation scripts: results tables and pass/fail checks."""

from __future__ import annotations

from pathlib import Path

RESULTS = Path(__file__).resolve().parent / "results"


def rel_err(computed: float, reference: float) -> float:
    return (computed - reference) / reference


def write_table(name: str, title: str, header: list[str], rows: list[list], notes: str = ""):
    """Write a Markdown table to ``validation/results/<name>.md`` and echo it to stdout."""
    lines = [f"# {title}", "", "| " + " | ".join(header) + " |"]
    lines.append("|" + "---|" * len(header))
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    if notes:
        lines += ["", notes]
    text = "\n".join(lines) + "\n"
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / f"{name}.md").write_text(text)
    print(text)
    return text
