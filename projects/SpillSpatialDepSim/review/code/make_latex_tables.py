"""
Script: review/code/make_latex_tables.py
Purpose: Write the numeric cells of the 3x3 and 3x4 beta-MSE tables of the
         manuscript (or of the Chapter 2 .qmd, which uses the same raw LaTeX tables) from review/tables/beta_mse_table_<grid>_<version>.csv,
         keeping the published LaTeX layout (row labels, multirow, rules).
Author: Andrew Walther (review written with Claude Code)
Created: 2026-10-25
Usage (from review/):
  python3 code/make_latex_tables.py <in.tex> <out.tex> <published|corrected>
In "published" mode the output must equal the input (checks the formatting).
"""
import csv, re, sys

# per grid: table label and the multiplier applied to each of the six columns
# (SRS MSE, SRS SD, SRS bias, BSS MSE, BSS SD, BSS bias), as in the paper
SCALE = {"3x3": ("InterventionEstimates3x3", (1, 10, 1, 1, 10, 1)),
         "3x4": ("InterventionEstimates3x4", (1, 100, 1, 10, 100, 1))}
COLS = ("SRS_MSE", "SRS_SD", "SRS_bias", "BSS_MSE", "BSS_SD", "BSS_bias")
# a data row: two label cells (which may contain an escaped "\&"), then psi, then the cells
ROW = re.compile(r"^(\s*(?:(?:\\&|[^&\n])*&){2}\s*)(0\.[5-8])(\s*&.*\\\\)\s*$")


def rows_for(grid, version):
    """Return the 16 formatted rows ('0.5 & ... \\\\') in table order."""
    _, mult = SCALE[grid]
    with open(f"tables/beta_mse_table_{grid}_{version}.csv") as f:
        recs = list(csv.DictReader(f))
    out = []
    for r in recs:
        cells = [f"{float(r[c]) * m:.3f}" for c, m in zip(COLS, mult)]
        cells = [c if c != "-0.000" else "0.000" for c in cells]
        out.append(f"{float(r['Psi']):.1f} & " + " & ".join(cells) + r" \\")
    return out


def splice(tex, grid, version):
    """Replace the 16 data rows of one table, keeping each row's label prefix."""
    label, _ = SCALE[grid]
    key = r"\label{table:" + label + "}"          # manuscript
    if key not in tex:
        key = "{#tbl-est-" + grid + "}"             # Chapter 2 .qmd
    end = tex.index(key)
    start = tex.rindex(r"\begin{table}", 0, end)
    block = tex[start:end].split("\n")
    new = iter(rows_for(grid, version))
    n = 0
    for i, line in enumerate(block):
        m = ROW.match(line)
        if m:
            block[i] = m.group(1) + next(new)
            n += 1
    assert n == 16, f"{grid}: found {n} data rows, expected 16"
    return tex[:start] + "\n".join(block) + tex[end:]


if __name__ == "__main__":
    src, dst, version = sys.argv[1:4]
    tex = open(src).read()
    for g in SCALE:
        tex = splice(tex, g, version)
    open(dst, "w").write(tex)
