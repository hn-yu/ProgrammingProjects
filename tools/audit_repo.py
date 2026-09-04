#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path
from urllib.parse import unquote

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]

md_link = re.compile(r"\[[^\]]*\]\(([^)]+)\)")
html_src = re.compile(r'''(?:src|href)=["']([^"']+)["']''', re.I)
code_fence = re.compile(r"```(?:c\+\+|cpp|cxx|c)?\s*\n(.*?)```", re.S | re.I)

IGNORE_PREFIXES = (
    "http://", "https://", "mailto:", "#", "javascript:", "data:",
)


def local_target(md: Path, raw: str):
    raw = raw.strip().strip("<>")
    if not raw or raw.startswith(IGNORE_PREFIXES):
        return None
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
        "mol.atom typo": re.compile(r"\bmol\.atom\b"),
        "an2masses typo": re.compile(r"\ban2masses\b"),
        "nonstandard M_PI": re.compile(r"\bM_PI\b"),
        "deprecated auto_ptr": re.compile(r"\bauto_ptr\b"),
        "obvious placeholder link": re.compile(r"\]\(addlink\)"),
        "Hartree-Fock misspelling": re.compile(r"\bHartee-Fock\b"),
    }
    hits = {k: [] for k in patterns}
    for p in ROOT.rglob("*"):
        if not p.is_file() or p.suffix.lower() not in {".md", ".cc", ".cpp", ".cxx", ".h", ".hpp", ".c"}:
            continue
        text = p.read_text(errors="replace")
        for label, rx in patterns.items():
            for m in rx.finditer(text):
                line = text.count("\n", 0, m.start()) + 1
                snippet = text[m.start():m.start()+180].replace("\n", " ")
                hits[label].append((p.relative_to(ROOT), line, snippet))
    return hits


def audit_code_fences():
    findings = []
    for md in ROOT.rglob("*.md"):
        text = md.read_text(errors="replace")
        for block_idx, block in enumerate(code_fence.findall(text), 1):
            declared = set(re.findall(r"\bFILE\s*\*\s*(\w+)", block))
            opened = set(re.findall(r"\b(\w+)\s*=\s*fopen\s*\(", block))
            scanned = re.findall(r"\bfscanf\s*\(\s*(\w+)", block)
            closed = re.findall(r"\bfclose\s*\(\s*(\w+)\s*\)", block)
            known = declared | opened
            for var in scanned:
                if known and var not in known:
                    findings.append((md.relative_to(ROOT), block_idx, f"fscanf uses '{var}' but FILE/fopen variables are {sorted(known)}"))
            for var in closed:
                if known and var not in known:
                    findings.append((md.relative_to(ROOT), block_idx, f"fclose uses '{var}' but FILE/fopen variables are {sorted(known)}"))

            for var in re.findall(r"\bdouble\s*\*\*\s*(\w+)\b", block):
                scalar_assign = re.search(rf"\b{re.escape(var)}\s*\[[^\]]+\]\s*=\s*[^=]", block)
                if scalar_assign:
                    rhs = block[scalar_assign.end()-1:scalar_assign.end()+80]
                    if not re.match(r"\s*new\s+double", rhs):
                        findings.append((md.relative_to(ROOT), block_idx, f"'{var}' declared double** but assigned through one subscript"))

            if re.search(r"sqrt\s*\(\s*\d+\s*/\s*\d+\s*\)", block):
                findings.append((md.relative_to(ROOT), block_idx, "sqrt contains integer division; use floating literals"))

    return findings


def audit_repo_hygiene():
    ds = [p.relative_to(ROOT) for p in ROOT.rglob(".DS_Store")]
    pyc = [p.relative_to(ROOT) for p in ROOT.rglob("*.pyc")]
    return ds, pyc


def audit_project1_indexing():
    findings = []
    for p in sorted((ROOT / "Project#01" / "hints").glob("step*-solution.md")):
        text = p.read_text(errors="replace")
        pos = text.find("Torsional angles")
        if pos >= 0:
            tail = text[pos:pos+1800]
            if "for(int j=0; j < i; j++)" in tail and "for(int k=0; k < j; k++)" in tail and "for(int l=0; l < k; l++)" in tail:
                findings.append((p.relative_to(ROOT), "torsion enumeration constrained by atom-number ordering"))
        if "cross_x /= norm;" in text and "cross_y /= norm;" in text and "sqrt(cross_x*cross_x" not in text:
            findings.append((p.relative_to(ROOT), "torsion sign helper divides by squared norm instead of sqrt(norm)"))
    return findings


def audit_project3_density_convention():
    p = ROOT / "Project#03" / "README.md"
    text = p.read_text(errors="replace").lower()
    return "one-spin" in text or "one spin" in text


def audit_pngs():
    rows = []
    for p in sorted(ROOT.rglob("*.png")):
        try:
            im = Image.open(p).convert("RGBA")
        except Exception as exc:
            rows.append((p.relative_to(ROOT), "ERROR", str(exc)))
            continue
        pixels = list(im.getdata())
        visible = [(r,g,b,a) for r,g,b,a in pixels if a > 8]
        if not visible:
            rows.append((p.relative_to(ROOT), "EMPTY", "no visible pixels"))
            continue
        alpha_transparent = sum(a < 250 for _,_,_,a in pixels) / len(pixels)
        lum = [0.2126*r + 0.7152*g + 0.0722*b for r,g,b,a in visible]
        mean_lum = sum(lum) / len(lum)
        dark_fraction = sum(x < 80 for x in lum) / len(lum)
        light_fraction = sum(x > 200 for x in lum) / len(lum)
        if alpha_transparent > 0.01 and dark_fraction > light_fraction:
            rows.append((p.relative_to(ROOT), "DARK_TRANSPARENT", f"mean={mean_lum:.1f} dark={dark_fraction:.3f} light={light_fraction:.3f}"))
    return rows


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

    code_findings = audit_code_fences()
    print(f"CODE_FENCE_FINDINGS={len(code_findings)}")
    for row in code_findings:
        print("CODE", *row, sep=" | ")

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

    png = audit_pngs()
    print(f"PNG_READABILITY_FINDINGS={len(png)}")
    for row in png:
        print("PNG", *row, sep=" | ")


if __name__ == "__main__":
    main()
