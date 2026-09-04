#!/usr/bin/env python3
from __future__ import annotations

import os
import re
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]

md_link = re.compile(r"\[[^\]]*\]\(([^)]+)\)")
html_src = re.compile(r'''(?:src|href)=["']([^"']+)["']''', re.I)

IGNORE_PREFIXES = (
    "http://", "https://", "mailto:", "#", "javascript:", "data:",
)


def local_target(md: Path, raw: str):
    raw = raw.strip().strip("<>")
    if not raw or raw.startswith(IGNORE_PREFIXES):
        return None
    # Markdown can include an optional title after whitespace.
    if ' "' in raw:
        raw = raw.split(' "', 1)[0]
    raw = raw.split("#", 1)[0].split("?", 1)[0]
    if not raw:
        return None
    raw = unquote(raw)
    return (md.parent / raw).resolve()


def audit_links():
    broken = []
    checked = 0
    for md in ROOT.rglob("*.md"):
        text = md.read_text(errors="replace")
        refs = md_link.findall(text) + html_src.findall(text)
        for raw in refs:
            target = local_target(md, raw)
            if target is None:
                continue
            checked += 1
            try:
                target.relative_to(ROOT.resolve())
            except ValueError:
                broken.append((md.relative_to(ROOT), raw, "escapes repository"))
                continue
            if not target.exists():
                broken.append((md.relative_to(ROOT), raw, str(target.relative_to(ROOT))))
    return checked, broken


def scan_patterns():
    patterns = {
        "uninitialized Eigen Matrix": re.compile(r"\bMatrix\s+\w+\s*\(\s*\d+\s*,\s*\d+\s*\)\s*;"),
        "squared norm used as norm divisor": re.compile(r"\bnorm\s*=\s*[^;]*\*[^;]*;[\s\S]{0,200}?/=\s*norm\s*;"),
        "bare natom (possible missing mol.)": re.compile(r"for\s*\([^\n;]*;\s*[^;]*<\s*natom\s*;"),
        "mol.atom typo": re.compile(r"\bmol\.atom\b"),
        "an2masses typo": re.compile(r"\ban2masses\b"),
        "nonstandard M_PI": re.compile(r"\bM_PI\b"),
        "deprecated auto_ptr": re.compile(r"\bauto_ptr\b"),
    }
    hits = {k: [] for k in patterns}
    for p in ROOT.rglob("*"):
        if not p.is_file() or p.suffix.lower() not in {".md", ".cc", ".cpp", ".cxx", ".h", ".hpp", ".c"}:
            continue
        text = p.read_text(errors="replace")
        for label, rx in patterns.items():
            for m in rx.finditer(text):
                line = text.count("\n", 0, m.start()) + 1
                snippet = text[m.start():m.start()+160].replace("\n", " ")
                hits[label].append((p.relative_to(ROOT), line, snippet))
    return hits


def audit_repo_hygiene():
    ds = [p.relative_to(ROOT) for p in ROOT.rglob(".DS_Store")]
    pyc = [p.relative_to(ROOT) for p in ROOT.rglob("*.pyc")]
    return ds, pyc


def audit_project1_indexing():
    findings = []
    for p in sorted((ROOT / "Project#01" / "hints").glob("step*-solution.md")):
        text = p.read_text(errors="replace")
        if "for(int j=0; j < i; j++)" in text and "for(int k=0; k < j; k++)" in text and "Torsional angles" in text:
            # carried-forward monotonic-index torsion loop can miss valid chains
            pos = text.find("Torsional angles")
            tail = text[pos:pos+1800]
            if "for(int j=0; j < i; j++)" in tail and "for(int k=0; k < j; k++)" in tail and "for(int l=0; l < k; l++)" in tail:
                findings.append((p.relative_to(ROOT), "torsion enumeration constrained by atom-number ordering"))
        if "cross_x /= norm;" in text and "cross_y /= norm;" in text:
            findings.append((p.relative_to(ROOT), "torsion sign helper divides by squared norm instead of sqrt(norm)"))
    return findings


def audit_project3_density_convention():
    p = ROOT / "Project#03" / "README.md"
    text = p.read_text(errors="replace")
    has_explicit = "one-spin" in text.lower() or "factor of two" in text.lower() and "density" in text.lower()
    return has_explicit


def main():
    checked, broken = audit_links()
    print(f"MARKDOWN_LOCAL_LINKS_CHECKED={checked}")
    print(f"BROKEN_LOCAL_LINKS={len(broken)}")
    for x in broken:
        print("BROKEN", *x, sep=" | ")

    hits = scan_patterns()
    for label, rows in hits.items():
        print(f"PATTERN {label}: {len(rows)}")
        for row in rows[:100]:
            print("HIT", label, *row, sep=" | ")

    ds, pyc = audit_repo_hygiene()
    print(f"DS_STORE={len(ds)}")
    for p in ds:
        print("DS", p)
    print(f"PYC={len(pyc)}")

    p1 = audit_project1_indexing()
    print(f"PROJECT1_INDEX_FINDINGS={len(p1)}")
    for p, msg in p1:
        print("P1", p, msg, sep=" | ")

    print("PROJECT3_DENSITY_CONVENTION_EXPLICIT=", audit_project3_density_convention())


if __name__ == "__main__":
    main()
