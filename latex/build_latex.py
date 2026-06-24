#!/usr/bin/env python3
"""Build the LaTeX documentation and clean up the build clutter.

Compiles every `*.tex` in this directory to PDF (running bibtex when the source
uses `\\bibliography`, and enough pdflatex passes to resolve cross-references),
reports undefined references, then deletes the auxiliary files LaTeX leaves
behind so the directory keeps only the essentials: `.tex`, `.pdf`, `.bib`, and
figures.

Usage:
    python latex/build_latex.py              # build all docs, then clean
    python latex/build_latex.py 02_beam_design   # build one (name or .tex path)
    python latex/build_latex.py --clean      # only remove clutter, do not build
    python latex/build_latex.py --no-clean   # build but keep the aux files
"""

import argparse
import os
import re
import shutil
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))

# Auxiliary extensions LaTeX/bibtex/tools leave behind -- safe to delete.
CLUTTER_EXT = (
    ".aux", ".log", ".out", ".bbl", ".blg", ".synctex.gz", ".fdb_latexmk",
    ".fls", ".toc", ".lof", ".lot", ".nav", ".snm", ".vrb", ".idx", ".ind",
    ".ilg", ".glo", ".gls", ".ist", ".acn", ".acr", ".alg", ".spl", ".bcf",
    ".run.xml",
)


def clean(tex_dir):
    """Delete LaTeX build clutter; return the number of files removed."""
    removed = 0
    for fn in os.listdir(tex_dir):
        if fn.endswith(CLUTTER_EXT):
            os.remove(os.path.join(tex_dir, fn))
            removed += 1
    return removed


def _pdflatex(stem, tex_dir):
    """One pdflatex pass; return (ok, combined_output)."""
    proc = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"{stem}.tex"],
        cwd=tex_dir, capture_output=True, text=True)
    return proc.returncode == 0, proc.stdout + proc.stderr


def _needs_bibtex(tex_path):
    with open(tex_path, encoding="utf-8") as f:
        return re.search(r"\\bibliography\{", f.read()) is not None


def build_one(tex_path):
    """Build a single .tex to PDF. Return (ok, n_undefined, message)."""
    tex_dir = os.path.dirname(tex_path)
    stem = os.path.splitext(os.path.basename(tex_path))[0]

    ok, out = _pdflatex(stem, tex_dir)
    if not ok:
        return False, 0, _first_error(out)

    if _needs_bibtex(tex_path):
        subprocess.run(["bibtex", stem], cwd=tex_dir, capture_output=True, text=True)
        ok, out = _pdflatex(stem, tex_dir)
        if not ok:
            return False, 0, _first_error(out)

    # Final pass to settle cross-references; keep its log for the ref check.
    ok, out = _pdflatex(stem, tex_dir)
    if not ok:
        return False, 0, _first_error(out)

    n_undef = len(re.findall(r"undefined (?:reference|citation)", out, re.I))
    pdf = os.path.join(tex_dir, f"{stem}.pdf")
    size = os.path.getsize(pdf) if os.path.exists(pdf) else 0
    return True, n_undef, f"{size:,} bytes"


def _first_error(out):
    """Pull the first TeX error line out of a failed run for the report."""
    for line in out.splitlines():
        if line.startswith("!"):
            return line.strip()
    return "pdflatex failed (see full log by re-running with --no-clean)"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("docs", nargs="*",
                    help="specific doc(s) to build (name or .tex); default: all in latex/")
    ap.add_argument("--clean", action="store_true", help="only remove clutter, do not build")
    ap.add_argument("--no-clean", action="store_true", help="build but keep the aux files")
    args = ap.parse_args()

    if args.clean:
        print(f"cleaned {clean(_HERE)} clutter files in {_HERE}")
        return

    if shutil.which("pdflatex") is None:
        sys.exit("error: pdflatex not found on PATH")

    if args.docs:
        stems = [os.path.splitext(os.path.basename(d))[0] for d in args.docs]
        texs = [os.path.join(_HERE, f"{s}.tex") for s in stems]
        missing = [t for t in texs if not os.path.exists(t)]
        if missing:
            sys.exit("error: not found: " + ", ".join(os.path.basename(m) for m in missing))
    else:
        texs = sorted(os.path.join(_HERE, f) for f in os.listdir(_HERE)
                      if f.endswith(".tex"))

    print(f"building {len(texs)} document(s) in {_HERE}\n")
    failures = warned = 0
    for tex in texs:
        name = os.path.basename(tex)
        ok, n_undef, msg = build_one(tex)
        if not ok:
            print(f"  FAIL  {name:<28} {msg}")
            failures += 1
        elif n_undef:
            print(f"  WARN  {name:<28} OK ({msg}) -- {n_undef} undefined reference(s)")
            warned += 1
        else:
            print(f"  OK    {name:<28} {msg}")

    if not args.no_clean:
        n = clean(_HERE)
        print(f"\ncleaned {n} clutter files (kept .tex .pdf .bib and figures)")
    else:
        print("\n--no-clean: aux files kept")

    if failures:
        sys.exit(f"\n{failures} document(s) failed to build")
    if warned:
        print(f"\n{warned} document(s) built with undefined references "
              f"(a rebuild usually resolves them)")
    print("\nall documents built successfully")


if __name__ == "__main__":
    main()
