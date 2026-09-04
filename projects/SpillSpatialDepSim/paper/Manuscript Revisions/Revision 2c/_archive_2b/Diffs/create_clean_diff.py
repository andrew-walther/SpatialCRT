#!/usr/bin/env python3
"""
create_clean_diff.py
Generate annotated diff between V2b_after_C6.tex (prev) and Revisions_V2b.tex (next).
Marks deletions with \st{...} (strikethrough, ulem pkg) and additions with {\color{blue}...}.
Math/table/algorithm environments use \cbstart...\cbend (changebar pkg).
"""
import difflib, re, sys, os

PREAMBLE_ADDITIONS = r"""% --- diff annotation packages (added by create_clean_diff.py) ---
\usepackage[normalem]{ulem}
\usepackage{changebar}
% ---"""

ENV_STARTS = re.compile(
    r'\\begin\{(?:equation\*?|align\*?|tabular|array|cases|algorithm|algorithmic|figure|table)\}'
)
ENV_ENDS = re.compile(
    r'\\end\{(?:equation\*?|align\*?|tabular|array|cases|algorithm|algorithmic|figure|table)\}'
)

def in_math_env(lines, idx):
    depth = 0
    for i, line in enumerate(lines[:idx+1]):
        if ENV_STARTS.search(line):
            depth += 1
        if ENV_ENDS.search(line):
            depth = max(0, depth - 1)
    return depth > 0

def annotate_diff(old_lines, new_lines):
    matcher = difflib.SequenceMatcher(None, old_lines, new_lines, autojunk=False)
    out = []
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == 'equal':
            out.extend(new_lines[j1:j2])
        elif tag == 'insert':
            block = new_lines[j1:j2]
            if any(ENV_STARTS.search(l) or ENV_ENDS.search(l) for l in block):
                out.append('\\cbstart\n')
                out.extend(block)
                out.append('\\cbend\n')
            else:
                for line in block:
                    stripped = line.rstrip('\n')
                    if stripped.strip() == '' or stripped.strip().startswith('%'):
                        out.append(line)
                    else:
                        out.append('{\\color{blue}' + stripped + '}\n')
        elif tag == 'delete':
            block = old_lines[i1:i2]
            if any(ENV_STARTS.search(l) or ENV_ENDS.search(l) for l in block):
                out.append('\\cbstart\n')
                for line in block:
                    stripped = line.rstrip('\n')
                    out.append('{\\color{red}\\st{' + stripped.replace('{', '\\{').replace('}', '\\}') + '}}\n')
                out.append('\\cbend\n')
            else:
                for line in block:
                    stripped = line.rstrip('\n')
                    if stripped.strip() == '' or stripped.strip().startswith('%'):
                        pass  # skip blank/comment deletions silently
                    else:
                        out.append('{\\color{red}\\st{' + stripped + '}}\n')
        elif tag == 'replace':
            old_block = old_lines[i1:i2]
            new_block = new_lines[j1:j2]
            is_structural = any(
                ENV_STARTS.search(l) or ENV_ENDS.search(l)
                for l in old_block + new_block
            )
            if is_structural:
                out.append('\\cbstart\n')
                for line in old_block:
                    stripped = line.rstrip('\n')
                    if stripped.strip() and not stripped.strip().startswith('%'):
                        out.append('{\\color{red}\\st{' + stripped + '}}\n')
                for line in new_block:
                    out.append(line)
                out.append('\\cbend\n')
            else:
                for line in old_block:
                    stripped = line.rstrip('\n')
                    if stripped.strip() and not stripped.strip().startswith('%'):
                        out.append('{\\color{red}\\st{' + stripped + '}}\n')
                for line in new_block:
                    stripped = line.rstrip('\n')
                    if stripped.strip() and not stripped.strip().startswith('%'):
                        out.append('{\\color{blue}' + stripped + '}\n')
                    else:
                        out.append(line)
    return out

def inject_packages(content):
    # Inject after \usepackage{xcolor} or before \begin{document}
    if '\\usepackage[normalem]{ulem}' in content:
        return content  # already present
    insert_after = '\\usepackage{xcolor}'
    idx = content.find(insert_after)
    if idx != -1:
        end = content.find('\n', idx) + 1
        return content[:end] + PREAMBLE_ADDITIONS + '\n' + content[end:]
    # fallback: insert before \begin{document}
    idx = content.find('\\begin{document}')
    return content[:idx] + PREAMBLE_ADDITIONS + '\n' + content[idx:]

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PREV = os.path.join(SCRIPT_DIR, 'snapshots', 'V2b_after_C6.tex')
NEXT = os.path.join(
    os.path.dirname(SCRIPT_DIR),
    'Walther_SpatialCRT_LaTeX_Revisions_V2b',
    'Revisions_V2b.tex'
)
OUT_TEX = os.path.join(SCRIPT_DIR, 'CLEAN_diff_annotated.tex')

with open(PREV, encoding='utf-8') as f:
    old_lines = f.readlines()
with open(NEXT, encoding='utf-8') as f:
    new_lines = f.readlines()

print(f"prev: {len(old_lines)} lines ({PREV})")
print(f"next: {len(new_lines)} lines ({NEXT})")

annotated = annotate_diff(old_lines, new_lines)
result = ''.join(annotated)
result = inject_packages(result)

with open(OUT_TEX, 'w', encoding='utf-8') as f:
    f.write(result)
print(f"Written: {OUT_TEX}")
