#!/usr/bin/env python3
from __future__ import annotations

import math
import re
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
P1 = ROOT / "Project#01"
CUTOFF = 4.0

OLD_ANGLE_LOOP = r'''  cout << "\nBond angles:\n";
  for(int i=0; i < mol.natom; i++) {
    for(int j=0; j < i; j++) {
      for(int k=0; k < j; k++) {
        if(mol.bond(i,j) < 4.0 && mol.bond(j,k) < 4.0)
          printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol.angle(i,j,k)*(180.0/acos(-1.0)));
      }
    }
  }
'''

NEW_ANGLE_LOOP = r'''  constexpr double distance_cutoff = 4.0;
  cout << "\nBond angles (both adjacent distances < 4.0 bohr):\n";
  // j is the central atom.  i < k removes only the i-j-k / k-j-i duplicate;
  // it does not constrain which atom may be the center.
  for(int j=0; j < mol.natom; j++) {
    for(int i=0; i < mol.natom; i++) {
      if(i == j) continue;
      for(int k=i+1; k < mol.natom; k++) {
        if(k == j) continue;
        if(mol.bond(i,j) < distance_cutoff && mol.bond(j,k) < distance_cutoff)
          printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol.angle(i,j,k)*(180.0/acos(-1.0)));
      }
    }
  }
'''

ANGLE_OUTPUT_RE = re.compile(
    r"(?m)^Bond angles:\n"
    r"(?:[ \t]*\d+-[ \t]*\d+-[ \t]*\d+[ \t]+[-+]?\d+(?:\.\d+)?[ \t]*\n)+"
)


def read_geometry(path: Path):
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    n = int(rows[0][0])
    if len(rows) < n + 1:
        raise ValueError(f"{path}: expected {n} atoms")
    coords = [tuple(map(float, rows[i + 1][1:4])) for i in range(n)]
    return coords


def distance(coords, a, b):
    return math.sqrt(sum((coords[a][q] - coords[b][q]) ** 2 for q in range(3)))


def angle_degrees(coords, i, j, k):
    rij = distance(coords, i, j)
    rkj = distance(coords, k, j)
    dot = sum(
        (coords[i][q] - coords[j][q]) * (coords[k][q] - coords[j][q])
        for q in range(3)
    )
    cosine = dot / (rij * rkj)
    cosine = max(-1.0, min(1.0, cosine))
    return math.degrees(math.acos(cosine))


def bond_angle_output(input_path: Path):
    coords = read_geometry(input_path)
    n = len(coords)
    lines = []
    for j in range(n):  # central atom
        for i in range(n):
            if i == j:
                continue
            for k in range(i + 1, n):  # i < k removes endpoint reversal only
                if k == j:
                    continue
                if distance(coords, i, j) < CUTOFF and distance(coords, j, k) < CUTOFF:
                    lines.append(f"{i:2d}-{j:2d}-{k:2d} {angle_degrees(coords, i, j, k):10.6f}")
    return "Bond angles (both adjacent distances < 4.0 bohr):\n" + "\n".join(lines) + "\n"


def replace_angle_output(text: str, output: str):
    return ANGLE_OUTPUT_RE.sub(output, text)


def fix_project1_text():
    changed = []

    readme = P1 / "README.md"
    text = readme.read_text()
    old = "Calculate all possible bond angles. For example,"
    new = (
        "Calculate all unique bond-angle candidates, allowing every atom to serve as the central atom. "
        "The 4.0 bohr distance cutoff used in the printing hints and examples is only a simple proximity "
        "filter; it is not a chemical bond definition. For example,"
    )
    if old in text:
        text = text.replace(old, new, 1)
        readme.write_text(text)
        changed.append(readme)

    hint34 = P1 / "hints" / "hint3-4.md"
    hint34.write_text(r'''Printing every central-angle combination can produce many geometrically uninteresting angles.  The examples in this project therefore use a **4.0 bohr proximity cutoff** for the two distances adjacent to the central atom.  This cutoff is only a simple printing heuristic; it is not a general chemical-bond criterion.

The important indexing rule is that the central atom must be allowed to have **any** index.  Do not enforce `i < j < k` (or `i > j > k`) while also treating `j` as the center, because that makes the result depend on atom numbering and omits valid central atoms.

A clean loop structure is:

```c++
constexpr double distance_cutoff = 4.0;

for(int j=0; j < mol.natom; j++) {          // central atom
  for(int i=0; i < mol.natom; i++) {
    if(i == j) continue;
    for(int k=i+1; k < mol.natom; k++) {    // i < k avoids i-j-k / k-j-i duplicates
      if(k == j) continue;
      if(R[i][j] < distance_cutoff && R[j][k] < distance_cutoff) {
        ...
      }
    }
  }
}
```

If you want **chemical** bond angles rather than cutoff-defined angle candidates, first construct a molecular connectivity graph using an appropriate bonding criterion, then choose pairs of neighbors around each central atom `j`.
''')
    changed.append(hint34)

    acetaldehyde_output = bond_angle_output(P1 / "input" / "acetaldehyde.dat")

    for step in range(3, 9):
        path = P1 / "hints" / f"step{step}-solution.md"
        text = path.read_text()
        count = text.count(OLD_ANGLE_LOOP)
        if count:
            text = text.replace(OLD_ANGLE_LOOP, NEW_ANGLE_LOOP)
        text = replace_angle_output(text, acetaldehyde_output)

        # Fix a few adjacent, definite compile/numerical bugs encountered in the carried-forward examples.
        text = text.replace("mol.atom", "mol.natom")
        text = text.replace("an2masses", "masses")
        text = text.replace("for(int i=0; i < natom; i++)", "for(int i=0; i < mol.natom; i++)")
        text = text.replace("Matrix I(3,3);", "Matrix I = Matrix::Zero(3,3);")
        text = text.replace("Matrix evals = solver.eigenvalues();", "Vector evals = solver.eigenvalues();")

        # Clamp the dot product before acos in the Step 3 implementation.
        old_angle_fn = r'''double Molecule::angle(int a, int b, int c)
{
  return acos(unit(0,b,a) * unit(0,b,c) + unit(1,b,a) * unit(1,b,c) + unit(2,b,a) * unit(2,b,c));
}
'''
        new_angle_fn = r'''double Molecule::angle(int a, int b, int c)
{
  double cosine = unit(0,b,a) * unit(0,b,c)
                + unit(1,b,a) * unit(1,b,c)
                + unit(2,b,a) * unit(2,b,c);
  if(cosine < -1.0) cosine = -1.0;
  if(cosine >  1.0) cosine =  1.0;
  return acos(cosine);
}
'''
        text = text.replace(old_angle_fn, new_angle_fn)

        path.write_text(text)
        changed.append(path)

    # Keep the provided reference outputs consistent with the corrected central-atom enumeration.
    for output_path in sorted((P1 / "output").glob("*_out.txt")):
        stem = output_path.stem.removesuffix("_out")
        input_path = P1 / "input" / f"{stem}.dat"
        if not input_path.exists():
            continue
        text = output_path.read_text()
        text = replace_angle_output(text, bond_angle_output(input_path))
        output_path.write_text(text)
        changed.append(output_path)

    return changed


def whiten_png(path: Path):
    """Make formula-style PNGs readable on dark backgrounds without wrecking color figures.

    * Transparent images: turn dark/gray foreground pixels white, preserve alpha and colored pixels.
    * Opaque, essentially grayscale images: invert luminance, giving white equations on a dark background.
    * Opaque color images are left unchanged.
    """
    try:
        image = Image.open(path).convert("RGBA")
    except Exception:
        return False

    pixels = list(image.getdata())
    if not pixels:
        return False

    visible = [p for p in pixels if p[3] > 0]
    if not visible:
        return False

    transparent_fraction = sum(p[3] < 250 for p in pixels) / len(pixels)
    grayscale_fraction = sum(max(p[:3]) - min(p[:3]) <= 12 for p in visible) / len(visible)

    new_pixels = []
    changed = False

    if transparent_fraction > 0.01:
        for r, g, b, a in pixels:
            chroma = max(r, g, b) - min(r, g, b)
            luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
            if a > 0 and chroma <= 24 and luminance < 235:
                new_pixels.append((255, 255, 255, a))
                changed |= (r, g, b) != (255, 255, 255)
            else:
                new_pixels.append((r, g, b, a))
    elif grayscale_fraction >= 0.97:
        for r, g, b, a in pixels:
            new_pixels.append((255 - r, 255 - g, 255 - b, a))
        changed = True
    else:
        return False

    if changed:
        image.putdata(new_pixels)
        image.save(path, format="PNG", optimize=True)
    return changed


def fix_pngs():
    changed = []
    skipped = []
    for path in sorted(ROOT.rglob("*.png")):
        if ".git" in path.parts:
            continue
        if whiten_png(path):
            changed.append(path)
        else:
            skipped.append(path)
    print(f"PNG conversion: changed {len(changed)}, left {len(skipped)} color/unchanged images alone")
    return changed


def main():
    text_changes = fix_project1_text()
    png_changes = fix_pngs()
    print(f"Project #01 text/reference files touched: {len(text_changes)}")
    print(f"PNG files touched: {len(png_changes)}")


if __name__ == "__main__":
    main()
