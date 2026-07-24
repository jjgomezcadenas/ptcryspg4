#!/usr/bin/env python3
"""Compile LaTeX documents in ``docs`` and remove auxiliary files.

Usage:
    python3 docs/build_latex.py
    python3 docs/build_latex.py xsection_fit xsections_plan
    python3 docs/build_latex.py --clean
    python3 docs/build_latex.py --no-clean xsection_fit

By default, each requested document is compiled with enough passes to resolve
cross-references and its auxiliary files are removed immediately afterwards.
PDF files are retained.
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path


DOCS_DIR = Path(__file__).resolve().parent
CLUTTER_SUFFIXES = (
    ".aux",
    ".log",
    ".out",
    ".bbl",
    ".blg",
    ".synctex.gz",
    ".fdb_latexmk",
    ".fls",
    ".toc",
    ".lof",
    ".lot",
    ".nav",
    ".snm",
    ".vrb",
    ".idx",
    ".ind",
    ".ilg",
    ".glo",
    ".gls",
    ".ist",
    ".acn",
    ".acr",
    ".alg",
    ".spl",
    ".bcf",
    ".run.xml",
)


def clean(directory: Path = DOCS_DIR) -> int:
    """Remove LaTeX auxiliary files directly under ``directory``."""
    removed = 0
    for path in directory.iterdir():
        if path.is_file() and path.name.endswith(CLUTTER_SUFFIXES):
            path.unlink()
            removed += 1
    return removed


def pdflatex(stem: str) -> tuple[bool, str]:
    process = subprocess.run(
        [
            "pdflatex",
            "-interaction=nonstopmode",
            "-halt-on-error",
            f"{stem}.tex",
        ],
        cwd=DOCS_DIR,
        capture_output=True,
        text=True,
        check=False,
    )
    return process.returncode == 0, process.stdout + process.stderr


def needs_bibtex(tex_path: Path) -> bool:
    return re.search(r"\\bibliography\{", tex_path.read_text(encoding="utf-8")) is not None


def first_error(output: str) -> str:
    for line in output.splitlines():
        if line.startswith("!"):
            return line.strip()
    return "pdflatex failed"


def build(tex_path: Path) -> tuple[bool, int, str]:
    stem = tex_path.stem
    ok, output = pdflatex(stem)
    if not ok:
        return False, 0, first_error(output)

    if needs_bibtex(tex_path):
        bibliography = subprocess.run(
            ["bibtex", stem],
            cwd=DOCS_DIR,
            capture_output=True,
            text=True,
            check=False,
        )
        if bibliography.returncode != 0:
            return False, 0, first_error(
                bibliography.stdout + bibliography.stderr)
        ok, output = pdflatex(stem)
        if not ok:
            return False, 0, first_error(output)

    ok, output = pdflatex(stem)
    if not ok:
        return False, 0, first_error(output)

    undefined = len(re.findall(
        r"undefined (?:reference|citation)", output, flags=re.IGNORECASE))
    pdf_path = DOCS_DIR / f"{stem}.pdf"
    size = pdf_path.stat().st_size if pdf_path.exists() else 0
    return True, undefined, f"{size:,} bytes"


def requested_documents(arguments: list[str]) -> list[Path]:
    if not arguments:
        return sorted(DOCS_DIR.glob("*.tex"))
    documents = [DOCS_DIR / f"{Path(item).stem}.tex" for item in arguments]
    missing = [path.name for path in documents if not path.exists()]
    if missing:
        raise FileNotFoundError(", ".join(missing))
    return documents


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("documents", nargs="*", help="document names or .tex paths")
    parser.add_argument("--clean", action="store_true", help="clean without compiling")
    parser.add_argument("--no-clean", action="store_true", help="retain auxiliary files")
    arguments = parser.parse_args()

    if arguments.clean:
        print(f"cleaned {clean()} auxiliary files in {DOCS_DIR}")
        return

    if shutil.which("pdflatex") is None:
        sys.exit("error: pdflatex not found on PATH")

    try:
        documents = requested_documents(arguments.documents)
    except FileNotFoundError as error:
        sys.exit(f"error: document not found: {error}")

    failures = 0
    for tex_path in documents:
        try:
            ok, undefined, message = build(tex_path)
        finally:
            if not arguments.no_clean:
                clean()

        if not ok:
            failures += 1
            print(f"FAIL  {tex_path.name}: {message}")
        elif undefined:
            print(f"WARN  {tex_path.name}: {undefined} undefined reference(s)")
        else:
            print(f"OK    {tex_path.name}: {message}")

    if failures:
        sys.exit(f"{failures} document(s) failed")


if __name__ == "__main__":
    main()
