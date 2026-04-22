"""Refresh the <QSP>...</QSP> block of consumer param_all.xml files.

Reads the snippet emitted by ``qsp-codegen`` (``qsp_params_xml_snippet.xml``)
and replaces the corresponding block in each target XML file.  Keeping this
in sync with the generated snippet prevents silent default-to-0 fallbacks in
``QSPParam.cpp`` when SBML param names change.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import List, Optional

QSP_BLOCK = re.compile(r"<QSP>.*?</QSP>", re.DOTALL)


def refresh(target: Path, snippet: Path) -> bool:
    """Replace the <QSP>...</QSP> block in ``target`` with ``snippet`` contents."""
    new_block = snippet.read_text().rstrip()
    text = target.read_text()
    new_text, n = QSP_BLOCK.subn(new_block.replace("\\", r"\\"), text, count=1)
    if n == 0:
        return False
    if new_text == text:
        print(f"  {target.name}: already up to date")
        return True
    target.write_text(new_text)
    print(f"  {target.name}: refreshed <QSP> block ({len(new_block)} bytes)")
    return True


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--snippet",
        type=Path,
        required=True,
        help="Path to qsp_params_xml_snippet.xml (emitted by qsp-codegen).",
    )
    parser.add_argument(
        "--xml",
        type=Path,
        action="append",
        required=True,
        help="param_all.xml target (repeatable).",
    )
    args = parser.parse_args(argv)

    if not args.snippet.exists():
        print(f"ERROR: snippet not found: {args.snippet}", file=sys.stderr)
        return 1

    n_fail = 0
    for t in args.xml:
        if not t.exists():
            print(f"  {t}: not found, skipping")
            continue
        if not refresh(t, args.snippet):
            print(f"ERROR: no <QSP>...</QSP> block found in {t}", file=sys.stderr)
            n_fail += 1
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
