#!/usr/bin/env python3
from __future__ import annotations

import math
import re
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
P1 = ROOT / "Project#01"
CUTOFF = 4.0


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected exactly one match, found {count}")
    return text.replace(old, new, 1)


def fix_top_readme():
    path = ROOT / "README.md"
    text = path.read_text()
    text = text.replace("Hartee-Fock", "Hartree-Fock")
    text = text.replace("[wiki](addlink)", "[wiki](https://github.com/CrawfordGroup/ProgrammingProjects/wiki)")
    text = text.replace(
        "These input files can be found in the `inputs` directory. \nWithin `input` there are directories for several different molecule/basis-set\ncombinations where you will find integrals, molecular geometries and other files to use as input to your programs.\n",
        "Projects that require data contain their own `input` directory. Some of those directories are further organized by molecule and basis set, where you will find integrals, molecular geometries, and other files used by the exercises.\n",
    )
    text = text.replace(
        "git clone git@github.com:CrawfordGroup/ProgrammingProjects.git",
        "git clone https://github.com/hn-yu/ProgrammingProjects.git",
    )
    for n in range(1, 15):
        old = f"https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%23{n:02d}"
        new = f"./Project%23{n:02d}"
        text = text.replace(old, new)
    path.write_text(text)


def fix_project1_readme_and_hints():
    path = P1 / "README.md"
    text = path.read_text()
    oop_anchor = '<img src="./figures/oop-angle.png" height="60">\n'
    oop_note = (
        oop_anchor
        + "\nIn the notation used by the solution code, `oop(i,j,k,l)` always treats the **third argument, `k`, as the central atom**. "
          "The order in which loop variables are nested does not change that meaning; the relevant proximity checks are i-k, k-j, and k-l.\n"
    )
    text = replace_once(text, oop_anchor, oop_note, "Project#01 README OOP note")

    tors_anchor = "Can you also determine the sign of the torsional angle?\n"
    tors_note = (
        tors_anchor
        + "\nFor unique torsions, canonicalize the **central bond** (for example, require `j < k`) rather than requiring "
          "`i < j < k < l` or the reverse. Atom numbering is arbitrary, so a monotonic-index rule can silently omit valid i-j-k-l chains.\n"
    )
    text = replace_once(text, tors_anchor, tors_note, "Project#01 README torsion note")
    path.write_text(text)

    hint = P1 / "hints" / "hint5-3.md"
    hint.write_text(
        "To print each torsional angle only once, canonicalize the **central bond** rather than imposing a global ordering on all four atom indices. "
        "A rule such as `j < i`, `k < j`, `l < k` depends on the arbitrary atom numbering and can omit valid i-j-k-l chains.\n\n"
        "One robust pattern is to require `j < k` for the central bond, then choose any distinct atom `i` close to `j` and any distinct atom `l` close to `k`:\n\n"
        "```c++\n"
        "for(int j=0; j < mol.natom; j++) {\n"
        "  for(int k=j+1; k < mol.natom; k++) {\n"
        "    if(mol.bond(j,k) >= distance_cutoff) continue;\n"
        "    for(int i=0; i < mol.natom; i++) {\n"
        "      if(i == j || i == k || mol.bond(i,j) >= distance_cutoff) continue;\n"
        "      for(int l=0; l < mol.natom; l++) {\n"
        "        if(l == i || l == j || l == k || mol.bond(k,l) >= distance_cutoff) continue;\n"
        "        // i-j-k-l is now a unique candidate for this central bond.\n"
        "      }\n"
        "    }\n"
        "  }\n"
        "}\n"
        "```\n\n"
        "As in the bond-angle exercise, the 4.0 bohr cutoff is only a simple proximity filter, not a general chemical bond definition. "
        "Also skip cases where either adjacent bond angle is 0 or 180 degrees, because the corresponding plane normal and torsion are undefined.\n"
    )

    hint7 = P1 / "hints" / "hint7-1.md"
    text = hint7.read_text()
    text = text.replace(
        "typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;\n",
        "typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;\n"
        "typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;\n",
    )
    text = text.replace(
        "- Allocate your moment of inertia tensor by a line of code like:\n```c++\nMatrix I(3,3);\n```",
        "- Allocate and zero-initialize your moment of inertia tensor. Eigen's size-only constructor does **not** initialize the entries:\n```c++\nMatrix I = Matrix::Zero(3,3);\n```",
    )
    text = text.replace("  Matrix evals = solver.eigenvalues();", "  Vector evals = solver.eigenvalues();")
    hint7.write_text(text)


OLD_TORSION_LOOP = r'''  cout << "\nTorsional angles:\n\n";
  for(int i=0; i < mol.natom; i++) {
    for(int j=0; j < i; j++) {
      for(int k=0; k < j; k++) {
        for(int l=0; l < k; l++) {
          if(mol.bond(i,j) < 4.0 && mol.bond(j,k) < 4.0 && mol.bond(k,l) < 4.0)
            printf("%2d-%2d-%2d-%2d %10.6f\n", i, j, k, l, mol.torsion(i,j,k,l)*(180.0/acos(-1.0)));
        }
      }
    }
  }
'''

NEW_TORSION_LOOP = r'''  cout << "\nTorsional angles:\n\n";
  // Canonicalize only the central bond (j < k).  A global index ordering such
  // as i > j > k > l would make the result depend on arbitrary atom numbering.
  for(int j=0; j < mol.natom; j++) {
    for(int k=j+1; k < mol.natom; k++) {
      if(mol.bond(j,k) >= distance_cutoff) continue;
      for(int i=0; i < mol.natom; i++) {
        if(i == j || i == k || mol.bond(i,j) >= distance_cutoff) continue;
        for(int l=0; l < mol.natom; l++) {
          if(l == i || l == j || l == k || mol.bond(k,l) >= distance_cutoff) continue;
          // A torsion is undefined when either plane normal vanishes.
          if(fabs(sin(mol.angle(i,j,k))) < 1e-12 ||
             fabs(sin(mol.angle(j,k,l))) < 1e-12) continue;
          printf("%2d-%2d-%2d-%2d %10.6f\n", i, j, k, l, mol.torsion(i,j,k,l)*(180.0/acos(-1.0)));
        }
      }
    }
  }
'''

OLD_NORM = r'''  double norm = cross_x*cross_x + cross_y*cross_y + cross_z*cross_z;
  cross_x /= norm;
  cross_y /= norm;
  cross_z /= norm;
  double sign = 1.0;
  double dot = cross_x*unit(0,b,c)+cross_y*unit(1,b,c)+cross_z*unit(2,b,c);
  if(dot < 0.0) sign = -1.0;
'''

NEW_NORM = r'''  double norm = sqrt(cross_x*cross_x + cross_y*cross_y + cross_z*cross_z);
  double dot = 0.0;
  if(norm > 1e-14) {
    cross_x /= norm;
    cross_y /= norm;
    cross_z /= norm;
    dot = cross_x*unit(0,b,c)+cross_y*unit(1,b,c)+cross_z*unit(2,b,c);
  }
  double sign = 1.0;
  if(dot < 0.0) sign = -1.0;
'''


def read_geometry(path: Path):
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    n = int(rows[0][0])
    return [tuple(map(float, rows[i + 1][1:4])) for i in range(n)]


def sub(a, b):
    return tuple(a[q] - b[q] for q in range(3))


def dot(a, b):
    return sum(a[q] * b[q] for q in range(3))


def cross(a, b):
    return (
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    )


def norm(a):
    return math.sqrt(dot(a, a))


def unit(coords, a, b):
    v = sub(coords[b], coords[a])  # from a to b, matching Molecule::unit
    r = norm(v)
    return tuple(x / r for x in v)


def distance(coords, a, b):
    return norm(sub(coords[a], coords[b]))


def angle(coords, a, b, c):
    x = max(-1.0, min(1.0, dot(unit(coords, b, a), unit(coords, b, c))))
    return math.acos(x)


def torsion(coords, a, b, c, d):
    eabc = cross(unit(coords, b, a), unit(coords, b, c))
    ebcd = cross(unit(coords, c, b), unit(coords, c, d))
    denom = math.sin(angle(coords, a, b, c)) * math.sin(angle(coords, b, c, d))
    cosine = dot(eabc, ebcd) / denom
    cosine = max(-1.0, min(1.0, cosine))
    tau = math.acos(cosine)
    cr = cross(eabc, ebcd)
    crn = norm(cr)
    sign = 1.0
    if crn > 1e-14:
        cr = tuple(x / crn for x in cr)
        if dot(cr, unit(coords, b, c)) < 0.0:
            sign = -1.0
    return tau * sign


def torsion_lines(input_path: Path):
    coords = read_geometry(input_path)
    n = len(coords)
    lines = []
    for j in range(n):
        for k in range(j + 1, n):
            if distance(coords, j, k) >= CUTOFF:
                continue
            for i in range(n):
                if i in (j, k) or distance(coords, i, j) >= CUTOFF:
                    continue
                for l in range(n):
                    if l in (i, j, k) or distance(coords, k, l) >= CUTOFF:
                        continue
                    if abs(math.sin(angle(coords, i, j, k))) < 1e-12:
                        continue
                    if abs(math.sin(angle(coords, j, k, l))) < 1e-12:
                        continue
                    val = math.degrees(torsion(coords, i, j, k, l))
                    lines.append(f"{i:2d}-{j:2d}-{k:2d}-{l:2d} {val:10.6f}")
    return lines


def replace_torsion_output(text: str, lines: list[str]) -> str:
    marker = "Torsional angles:\n\n"
    start = text.find(marker)
    if start < 0:
        return text
    body_start = start + len(marker)
    ends = [x for x in (text.find("\nMolecular center", body_start), text.find("\n```", body_start)) if x >= 0]
    if not ends:
        return text
    end = min(ends)
    section = marker + "\n".join(lines) + "\n"
    return text[:start] + section + text[end:]


def fix_project1_solutions_and_outputs():
    acetaldehyde_lines = torsion_lines(P1 / "input" / "acetaldehyde.dat")
    for step in range(5, 9):
        path = P1 / "hints" / f"step{step}-solution.md"
        text = path.read_text()
        if OLD_TORSION_LOOP not in text:
            raise RuntimeError(f"{path}: old torsion loop not found")
        text = text.replace(OLD_TORSION_LOOP, NEW_TORSION_LOOP, 1)
        if step == 5:
            if OLD_NORM not in text:
                raise RuntimeError("step5: old torsion normalization not found")
            text = text.replace(OLD_NORM, NEW_NORM, 1)
        text = replace_torsion_output(text, acetaldehyde_lines)
        path.write_text(text)

    for output_path in sorted((P1 / "output").glob("*_out.txt")):
        stem = output_path.stem.removesuffix("_out")
        input_path = P1 / "input" / f"{stem}.dat"
        if input_path.exists():
            text = replace_torsion_output(output_path.read_text(), torsion_lines(input_path))
            output_path.write_text(text)


def fix_project2():
    path = ROOT / "Project#02" / "README.md"
    text = path.read_text()
    anchor = '<img src="./figures/hessian-file-format.png" width="200">\n'
    note = (
        anchor
        + "\nThe in-memory Hessian has `(3N) x (3N) = (3N)^2` scalar elements. The text file groups **three scalar values per data row**, "
          "so after the initial atom-count line it contains `3N^2` rows, not `(3N)^2` rows. For water (`N=3`), that is 27 rows containing 81 Hessian values.\n"
    )
    text = replace_once(text, anchor, note, "Project#02 Hessian file-format note")
    text = text.replace("squareroot", "square root")
    path.write_text(text)

    hint = ROOT / "Project#02" / "hints" / "hint1.md"
    text = hint.read_text()
    intro = "The Hessian stored in memory should be a square matrix, while the format of the input file is rectangular.  Understanding the translation between the two takes a bit of thinking."
    repl = intro + " Each input row contains three consecutive Hessian elements; therefore an N-atom system has 3N^2 data rows containing (3N)^2 scalar values in total."
    text = replace_once(text, intro, repl, "Project#02 hint1 row-count note")
    hint.write_text(text)


def fix_project3():
    path = ROOT / "Project#03" / "README.md"
    text = path.read_text()
    text = text.replace("Compute the Inital SCF Energy", "Compute the Initial SCF Energy")

    eri_anchor = "and only the permutationally unique integrals are provided in the file, with the restriction that, for each integral, the following relationships hold:\n"
    eri_repl = (
        "and only the **non-zero, permutationally unique** integrals are provided in the file. Any packed or full ERI storage must therefore be initialized to zero before the file is read. "
        "For each integral that is present, the following relationships hold:\n"
    )
    text = replace_once(text, eri_anchor, eri_repl, "Project#03 ERI sparsity clarification")

    density_anchor = "where *m* indexes the columns of the coefficient matrices, and the summation includes only the occupied spatial MOs.\n"
    density_repl = (
        density_anchor
        + "\n**Density convention used in this project.** This `P` is the one-spin (spatial-orbital) density, `P_{mu nu} = sum_m C_{mu m} C_{nu m}`, so there is deliberately **no factor of 2** in its definition. "
          "The closed-shell factor of two is instead carried explicitly by the `2J-K` Fock build and by one-electron property formulas. If you choose the spin-summed convention `D = 2P`, the Fock and energy formulas must be changed consistently.\n"
    )
    text = replace_once(text, density_anchor, density_repl, "Project#03 density convention")
    text = text.replace(
        "- The factor 2 appearing above arises because the definition of the density used in this project differs from that used in Szabo & Ostlund.",
        "- The factor 2 appearing above is the closed-shell spin factor because this project uses the one-spin spatial density `P`, rather than the spin-summed density `D = 2P`.",
    )
    path.write_text(text)

    hint = ROOT / "Project#03" / "hints" / "hint3-3.md"
    text = hint.read_text()
    text = text.replace('fscanf(hessian, "%d %d %lf", &a, &b, &S[i][j]);', 'fscanf(input, "%d %d %lf", &a, &b, &S[i][j]);')
    text = text.replace("  double **TEI;", "  double *TEI;")
    text = text.replace(
        'while(fscanf(input, "%d %d %d %d %lf", &i, &j, &k, &l, &val) != EOF) {',
        'while(fscanf(input, "%d %d %d %d %lf", &i, &j, &k, &l, &val) == 5) {'
    )
    text = text.replace(
        "However the two-electron integrals are not provided in this convenient ordering, you could go through the file and attempt to reorder them so that a loop structure works, or you could put the previously ignored index information you're reading with your fscanf to good use. This prevents the use of hard-wired loops that assume a specific ordering of the data, and it also prevents having to calculate how many elements are in the file to be read.\n",
        "However the two-electron integrals are not provided in this convenient ordering, and the file omits integrals whose value is zero. Initialize your packed TEI storage to zero, then use the explicit indices on each line to place every non-zero value. This avoids hard-wired loops that assume a specific ordering or line count.\n"
    )
    text = text.replace(
        "The `EOF` above is the [end-of-file](http://en.wikipedia.org/wiki/End-of-file) condition, so `!= EOF` means the above while loop scans the file until it reaches the end.\n",
        "`fscanf` returns the number of successfully converted fields, so `== 5` accepts only complete ERI records and stops cleanly at end-of-file or malformed input.\n"
    )
    hint.write_text(text)


def fix_hygiene_and_png():
    gitignore = ROOT / ".gitignore"
    text = gitignore.read_text()
    if ".DS_Store" not in text.splitlines():
        if text and not text.endswith("\n"):
            text += "\n"
        text += ".DS_Store\n"
        gitignore.write_text(text)

    for p in ROOT.rglob(".DS_Store"):
        p.unlink()

    png = P1 / "figures" / "determinant.png"
    image = Image.open(png).convert("RGBA")
    out = []
    for r, g, b, a in image.getdata():
        if a > 0:
            out.append((255-r, 255-g, 255-b, a))
        else:
            out.append((r, g, b, a))
    image.putdata(out)
    image.save(png, format="PNG", optimize=True)


def main():
    fix_top_readme()
    fix_project1_readme_and_hints()
    fix_project1_solutions_and_outputs()
    fix_project2()
    fix_project3()
    fix_hygiene_and_png()
    print("repository-wide fixes applied")


if __name__ == "__main__":
    main()
