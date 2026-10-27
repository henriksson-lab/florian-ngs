#!/usr/bin/env python3
"""Build the small GitHub Pages site into ``_site/``."""
from __future__ import annotations

import html
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
NGS_101 = ROOT / "NGS_101"
if not (NGS_101 / "lib").is_dir():
    raise SystemExit("missing NGS_101 dependency; run: git submodule update --init")
sys.path.insert(0, str(NGS_101 / "lib"))

from mdrender import render  # noqa: E402
from page import ASSETS, MD_STYLE  # noqa: E402

SITE = ROOT / "_site"
SCHEMATIC = ROOT / "florian-PTA-rnaseq.html"


def index_page() -> str:
    return f"""<!doctype html>
<title>florian-PTA-rnaseq</title>
{ASSETS}
<div class="wrap">
<h1>florian-PTA-rnaseq</h1>
<p>Work-in-progress single-cell RNA-seq design using an Atrandi SPC/PTA back end.</p>
<p><a href="florian-PTA-rnaseq.html"><b>Protocol schematic</b></a><br>
<a href="01_primers.html">Research and design notes</a><br>
<a href="https://henriksson-lab.github.io/NGS_101/">NGS 101 chemistry reference</a></p>
<h2>Proposed pilots</h2>
<p><a href="codex-rnaseq-i.html"><b>Codex RNAseq I: PTA-free 5-prime RNA capture
with directional adapter ligation</b></a><br>
Proposed, experimentally unvalidated design using BciVI and the original FakeD adapter.</p>
</div>
"""


def note_page() -> str:
    source = ROOT / "01_primers.md"
    body = render(source.read_text(encoding="utf-8"))
    return f"""<!doctype html>
<title>florian-PTA-rnaseq research notes</title>
{ASSETS}
{MD_STYLE}
<div class="wrap">
<nav class="crumb"><a href="index.html">florian-PTA-rnaseq</a> &nbsp;/&nbsp;
{html.escape(source.name)}</nav>
<article class="md">{body}</article>
</div>
"""


def main() -> None:
    subprocess.run([sys.executable, str(ROOT / "tools" / "build_page.py")], check=True)
    subprocess.run([sys.executable, str(ROOT / "tools" / "build_codex_rnaseq_i.py")], check=True)
    if SITE.exists():
        shutil.rmtree(SITE)
    SITE.mkdir()
    (SITE / "index.html").write_text(index_page(), encoding="utf-8")
    (SITE / "01_primers.html").write_text(note_page(), encoding="utf-8")
    shutil.copy2(SCHEMATIC, SITE / SCHEMATIC.name)
    shutil.copy2(ROOT / "codex-rnaseq-i.html", SITE / "codex-rnaseq-i.html")
    (SITE / "reference").mkdir()
    shutil.copy2(ROOT / "reference" / "codex-rnaseq-i-original.md",
                 SITE / "reference" / "codex-rnaseq-i-original.md")
    print(f"wrote {SITE} ({len(list(SITE.glob('*.html')))} pages)")


if __name__ == "__main__":
    main()
