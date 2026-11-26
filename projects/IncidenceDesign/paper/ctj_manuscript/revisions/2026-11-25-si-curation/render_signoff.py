#!/usr/bin/env python3
"""
Script: render_signoff.py
Purpose: Write a readable copy of a sign-off file, with each LaTeX current or
         proposed block shown as rendered text (curated-SI view) and the
         changes marked.
Author: Andrew Walther
Created: 2026-11-25
Dependencies: Python 3 standard library, pandoc

Usage (from this folder):
    python3 render_signoff.py 2g-merged-signoff.md

Writes <name>-readable.md beside the input. The LaTeX file stays the source
of truth for approval and for the apply script; this copy is for reading.

Rendering rules
- \\iffullsi A\\else B\\fi prints B and \\iffullsi A\\fi prints nothing, as in
  the curated SI.
- \\ref{label} prints the number from the last build's .aux files (SI labels
  from the curated SI, M- labels from the article, SI- labels from the SI).
- \\cite{key} prints "Surname year" taken from the citation key.
- Pandoc renders the rest to plain text, with math as Unicode.
- In each pair, words removed from the current block are ~~struck~~ and
  words added in the proposed block are **bold**.
- ```qmd blocks (Chapter 3) are rendered from Quarto markdown; chapter
  labels print as [label].
- Table rows (lines with & and \\\\) are left as LaTeX, since a row
  fragment has no plain-text form.
"""

import difflib
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]  # paper/ctj_manuscript
SI_AUX = PAPER / "supplementary_info/.build/Supplementary_Information.aux"
CTJ_AUX = PAPER / "manuscript_ctj/.build/CTJ_Manuscript.aux"


# Label lookup ----

def read_labels(aux):
    """Return {label: printed number} from a LaTeX .aux file."""
    labels = {}
    for line in aux.read_text(errors="replace").splitlines():
        m = re.match(r"\\newlabel\{([^}]*)\}\{\{([^}]*)\}", line)
        if m:
            labels[m.group(1)] = m.group(2)
    return labels


SI_LABELS = read_labels(SI_AUX)
FULL_LABELS = read_labels(SI_AUX.with_name("Supplementary_Information_Full.aux"))
CTJ_LABELS = read_labels(CTJ_AUX)


def resolve_ref(label, target, full=False):
    """Printed number for \\ref{label} in the SI ("si") or the article ("ctj").

    full=True reads unabridged-SI numbers, for text only that build prints.
    """
    if target == "si":
        if label.startswith("M-"):
            return CTJ_LABELS.get(label[2:], "??")
        return (FULL_LABELS if full else SI_LABELS).get(label, "??")
    if label.startswith("SI-"):
        return SI_LABELS.get(label[3:], "??")
    return CTJ_LABELS.get(label, "??")


# LaTeX preprocessing ----

def curated_view(tex, full=False):
    """Keep what the curated SI prints from each \\iffullsi switch.

    A fragment may start inside a switch (an orphan \\else or \\fi before any
    \\iffullsi) or end inside one (no closing \\fi); both are handled as the
    surrounding source would print them. full=True keeps the unabridged branch.
    """
    def pick(body):
        if "\\else" in body:
            a, b = body.split("\\else", 1)
            return a if full else b
        return body if full else ""

    def skip_space(t, i):
        if t[i:i + 1] == "\\" and t[i + 1:i + 2] in (" ", "\n"):
            return i + 1  # "\fi\ " leaves its control space as a plain space
        return i

    # Fragment starting inside a switch
    first_if = tex.find("\\iffullsi")
    first_fi = tex.find("\\fi")
    if first_fi >= 0 and (first_if < 0 or first_fi < first_if):
        tex = pick(tex[:first_fi]) + tex[skip_space(tex, first_fi + 3):]

    out, i = [], 0
    while True:
        j = tex.find("\\iffullsi", i)
        if j < 0:
            out.append(tex[i:])
            return "".join(out)
        out.append(tex[i:j])
        k = tex.find("\\fi", j)
        body = tex[j + len("\\iffullsi"):k if k >= 0 else len(tex)]
        out.append(pick(body))
        if k < 0:  # fragment ends inside a switch
            return "".join(out)
        i = skip_space(tex, k + len("\\fi"))


def balance(tex):
    """Drop unmatched closing braces and \\end{..}, and close open braces."""
    out, depth, i = [], 0, 0
    while i < len(tex):
        c = tex[i]
        if c == "\\":
            out.append(tex[i:i + 2])
            i += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            if depth == 0:
                i += 1
                continue
            depth -= 1
        out.append(c)
        i += 1
    tex = "".join(out) + "}" * depth
    m = re.match(r"\s*\\end\{(\w+\*?)\}", tex)
    if m and "\\begin{%s}" % m.group(1) not in tex[:m.start()]:
        tex = tex[m.end():]
    return tex


def cite_text(keys):
    """'besag_spatial_1974, lee_x_2013' -> 'Besag 1974; Lee 2013'."""
    parts = []
    for key in keys.split(","):
        bits = key.strip().split("_")
        parts.append(f"{bits[0].title()} {bits[-1]}" if len(bits) > 1 else key.strip())
    return "; ".join(parts)


def preprocess(tex, target, full=False):
    tex = curated_view(tex, full)
    tex = re.sub(r"\\caption\{", "{", tex)
    tex = re.sub(r"\\label\{[^}]*\}", "", tex)
    tex = balance(tex)
    # A comma with no space between digits is a thousands separator. Pandoc
    # spaces it inside math ("100, 000"), so hold it as THOUSANDS until render()
    tex = re.sub(r"(\d),(\d{3})(?!\d)", "\\1" + THOUSANDS + "\\2", tex)
    tex = re.sub(r"\\(?:eq)?ref\{([^}]*)\}",
                 lambda m: resolve_ref(m.group(1), target, full), tex)
    tex = re.sub(r"\\cite[pt]?\{([^}]*)\}",
                 lambda m: f"({cite_text(m.group(1))})", tex)
    return tex


UNABRIDGED = "[unabridged only] "
THOUSANDS = "\u01c2"  # placeholder character, restored to "," after pandoc


def pandoc(tex, fmt="latex"):
    res = subprocess.run(
        ["pandoc", "-f", fmt, "-t", "plain", "--wrap=none"],
        input=tex, capture_output=True, text=True)
    return res.stdout.strip() if res.returncode == 0 else ""


def render(tex, target):
    """Rendered plain text of a LaTeX fragment, or None if it is a table row.

    A fragment that prints nothing in the curated SI is rendered from the
    unabridged branch and labeled as such.
    """
    if re.search(r"&.*\\\\", tex):
        return None
    if target == "qmd":
        # Chapter (Quarto) text: chapter labels print as [label], citations as (key)
        tex = re.sub(r"\\(?:eq)?ref\{([^}]*)\}", r"[\1]", tex)
        tex = re.sub(r"@sec-([\w-]+)", r"[sec-\1]", tex)
        tex = re.sub(r"\[@([^\]]*)\]", lambda m: "(" + cite_text(m.group(1).replace("@", "").replace(";", ",")) + ")", tex)
        return pandoc(tex, "markdown") or None
    text = pandoc(preprocess(tex, target)).replace(THOUSANDS, ",")
    if text:
        return text
    if "\\iffullsi" in tex:
        full = pandoc(preprocess(tex, target, full=True)).replace(THOUSANDS, ",")
        if full:
            return UNABRIDGED + full
    return None


# Word-level marking ----

def md_escape(word):
    return re.sub(r"([*_~`\[\]<>|])", r"\\\1", word)


def tokens(text):
    """Words, with paragraph breaks kept as their own token."""
    out = []
    for k, para in enumerate(re.split(r"\n\s*\n", text)):
        if k:
            out.append("\n\n")
        out.extend(para.split())
    return out


def join(words):
    s = " ".join(words)
    return re.sub(r" ?\n\n ?", "\n\n", s)


def wrap(words, mark):
    """Mark a run of words, one mark per paragraph so markdown stays valid."""
    paras, cur = [], []
    for w in words:
        if w == "\n\n":
            paras.append(cur)
            cur = []
        else:
            cur.append(md_escape(w))
    paras.append(cur)
    return "\n\n".join(f"{mark}{' '.join(p)}{mark}" if p else "" for p in paras)


def mark_pair(cur_text, new_text):
    """Return (current with ~~removed~~ words, proposed with **added** words)."""
    a, b = tokens(cur_text), tokens(new_text)
    sm = difflib.SequenceMatcher(a=a, b=b, autojunk=False)
    ca, pb = [], []
    for op, i1, i2, j1, j2 in sm.get_opcodes():
        if op == "equal":
            ca.append(join([md_escape(w) if w != "\n\n" else w for w in a[i1:i2]]))
            pb.append(join([md_escape(w) if w != "\n\n" else w for w in b[j1:j2]]))
        else:
            if i2 > i1:
                ca.append(wrap(a[i1:i2], "~~"))
            if j2 > j1:
                pb.append(wrap(b[j1:j2], "**"))
    fix = lambda parts: re.sub(r" ?\n\n ?", "\n\n", " ".join(p for p in parts if p))
    return fix(ca), fix(pb)


def quote(text):
    return "\n".join("> " + line if line else ">" for line in text.splitlines())


# Main ----

def main(path):
    src = Path(path)
    lines = src.read_text().splitlines()

    # Pass 1: find every ```latex block, its label line and its target file
    blocks, target, i = [], "si", 0
    while i < len(lines):
        if lines[i].startswith("## Batch E"):
            target = "ctj"
        elif lines[i].startswith("## ") and not lines[i].startswith("## Batch"):
            target = "si"
        if lines[i].strip() in ("```latex", "```qmd"):
            j = i + 1
            while lines[j].strip() != "```":
                j += 1
            label = next((lines[k].strip() for k in range(i - 1, -1, -1)
                          if lines[k].strip()), "")
            blocks.append(dict(start=i, end=j, label=label,
                               target="qmd" if lines[i].strip() == "```qmd" else target,
                               tex="\n".join(lines[i + 1:j])))
            i = j
        i += 1

    for b in blocks:
        b["text"] = render(b["tex"], b["target"])

    # Pair "**Current ...**" with the next "**Proposed ...**" of the same suffix
    for k, b in enumerate(blocks):
        m = re.match(r"\*\*Current(.*)\*\*$", b["label"])
        if not m:
            continue
        for c in blocks[k + 1:]:
            if c["label"] == f"**Proposed{m.group(1)}**":
                if (b["text"] and c["text"] and UNABRIDGED not in
                        b["text"][:20] + c["text"][:20]):
                    b["marked"], c["marked"] = mark_pair(b["text"], c["text"])
                break

    # Pass 2: write the readable copy
    out, i, nb, unresolved = [], 0, 0, 0
    starts = {b["start"]: b for b in blocks}
    title_done = False
    while i < len(lines):
        if i in starts:
            b = starts[i]
            if b["text"] is None:
                out.append("*(table rows or a fragment pandoc cannot render, shown as LaTeX)*\n")
                out.extend(lines[b["start"]:b["end"] + 1])
            else:
                body = b.get("marked") or md_escape_block(b["text"])
                if b["text"].startswith(UNABRIDGED):
                    out.append("*Prints nothing in the curated SI. The unabridged "
                               "build prints (its own numbering):*\n")
                    body = md_escape_block(b["text"][len(UNABRIDGED):])
                unresolved += body.count("??")
                out.append(quote(body))
            i = b["end"] + 1
            nb += 1
            continue
        out.append(lines[i])
        if not title_done and lines[i].startswith("# "):
            out.append("")
            out.append(
                f"*Readable copy of [{src.name}]({src.name}), written by "
                "`render_signoff.py`. Each block shows the text as the curated SI "
                "(or the article) prints it, with section, table and figure numbers "
                "from the last build. In each pair, ~~struck~~ words leave the "
                "current text and **bold** words are new. Approve by ID against "
                "this copy; the apply script uses the LaTeX file.*")
            title_done = True
        i += 1

    dest = src.with_name(src.stem + "-readable.md")
    dest.write_text("\n".join(out) + "\n")
    print(f"wrote {dest.name}: {nb} blocks, "
          f"{sum(b['text'] is None for b in blocks)} left as LaTeX, "
          f"{sum('marked' in b for b in blocks)} marked, "
          f"{unresolved} unresolved references")


def md_escape_block(text):
    return "\n".join(" ".join(md_escape(w) for w in line.split())
                     for line in text.splitlines())


if __name__ == "__main__":
    main(sys.argv[1])
