#!/usr/bin/env python3
"""Wrap debugMes(msg, level); statements so msg is only built when it prints."""
import re, sys

def scan_call(src, i):
    """src[i] == '(' ; return index just past the matching ')' and top-level comma positions."""
    depth, j, commas = 0, i, []
    while j < len(src):
        c = src[j]
        if c == '"':
            j += 1
            while src[j] != '"':
                j += 2 if src[j] == '\\' else 1
        elif c == "'":
            j += 1
            while src[j] != "'":
                j += 2 if src[j] == '\\' else 1
        elif c == '(':
            depth += 1
        elif c == ')':
            depth -= 1
            if depth == 0:
                return j + 1, commas
        elif c == ',' and depth == 1:
            commas.append(j)
        j += 1
    raise ValueError("unbalanced at %d" % i)

def strip_strings(s):
    return re.sub(r'"(\\.|[^"\\])*"', '""', s)

SIDE = re.compile(r'\+\+|--|(?<![=!<>])=(?!=)|\.(next|remove|add|put|poll|pop|push|set)\s*\(')

path = sys.argv[1]
lo = int(sys.argv[2]) if len(sys.argv) > 2 else 0
hi = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9
exclude = set(int(x) for x in sys.argv[4].split(",")) if len(sys.argv) > 4 and sys.argv[4] else set()
eligible = 0
src = open(path).read()
out, pos, wrapped, skipped = [], 0, 0, []
for m in re.finditer(r'\bdebugMes\s*\(', src):
    start = m.start()
    if start < pos:
        continue
    line_start = src.rfind("\n", 0, start) + 1
    prefix = src[line_start:start]
    if re.search(r'\bvoid\s*$', prefix) or prefix.strip().startswith("//") or "//" in prefix:
        continue
    open_paren = m.end() - 1
    end, commas = scan_call(src, open_paren)
    k = end
    while src[k] in " \t":
        k += 1
    if src[k] != ';' or not commas:
        skipped.append((src.count("\n", 0, start) + 1, "not a plain statement"))
        continue
    msg = src[open_paren + 1:commas[-1]]
    level = src[commas[-1] + 1:end - 1].strip()
    if SIDE.search(strip_strings(msg)):
        skipped.append((src.count("\n", 0, start) + 1, "possible side effect"))
        continue
    idx = eligible
    eligible += 1
    if not (lo <= idx < hi) or idx in exclude:
        continue
    stmt = src[start:k + 1]
    out.append(src[pos:start])
    out.append("{ if (BFLY_GLOBALS.VERBOSE_LEVEL >= (%s)) %s }" % (level, stmt))
    pos = k + 1
    wrapped += 1
out.append(src[pos:])
open(path, "w").write("".join(out))
print("%s: wrapped %d of %d eligible, skipped %d" % (path, wrapped, eligible, len(skipped)))
for ln, why in skipped:
    print("  skip line %d: %s" % (ln, why))
