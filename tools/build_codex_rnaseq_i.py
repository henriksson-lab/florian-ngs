#!/usr/bin/env python3
"""Render the local proposed pilot using the pinned NGS_101 model/rendering API."""
from __future__ import annotations

import html
from pathlib import Path

import codex_rnaseq_i as c
from chemdraw import oligo, panel, workflow_panel
from page import caveat, head, info, table

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'codex-rnaseq-i.html'
RAW = 'reference/codex-rnaseq-i-original.md'


def oligos():
    entries = (
        ('RT oligo-dT', c.legacy.OLIGO_DT, '', 'Existing; unmodified DNA'),
        ('Existing TSO', c.H + 'N' * 8 + 'WW' + c.legacy.TSO_G_TAIL,
         c.legacy.TSO_5_MOD, 'Existing vendor biotin, ribo-G and LNA-G notation'),
        ('H linear primer', c.H, '', 'Existing canonical TSO handle'),
        ('Q PCR primer', c.Q, '', 'Candidate; thermodynamic review required'),
        ('BciVI tailed random primer', c.RANDOM_PRIMER, '',
         'Candidate Q–B–site–spacer–N8; B is incorporated, not ordered separately'),
        ('FakeTruseq-TSO / custom Read 1', c.FAKE_TRUSEQ, '', 'Existing ordered sequence'),
        ('FakeD top', c.FAKED_TOP, '', 'Original adapter; terminal single-base 3′-T'),
        ('FakeD bottom', c.FAKED_BOTTOM, '/5Phos/', 'Original adapter; only bottom is 5′-phosphorylated'),
        ('P5 indexing primer', c.P5, '', 'Existing indexing primer'),
    )
    rows = []
    for name, sequence, mods, status in entries:
        # Raw modified oligo notation is intentionally selectable, never a DNA duplex.
        rendered = oligo(name, [c.seg(name, sequence, placeholder=True)], mods=mods)
        rows.append(rendered + '<p><small>' + html.escape(status) + '</small></p>')
    return '<seq>' + '\n'.join(rows) + '</seq>'


def build():
    workflow, graph = c.build_model()
    # Build_model checks graph provenance, terminal fates and active sequencing primers.
    raw_architecture = html.escape(c.il.P7 + ' + selected 10-base i7 + ' + c.il.TRUSEQ_READ2)
    handles = table(['Component', '5′ → 3′', 'Role'], [
        ['Q', f'<code>{c.Q}</code>', 'Candidate PCR handle'],
        ['B', f'<code>{c.B}</code>', 'Candidate extra handle; not ordered separately'],
        ['BciVI site', f'<code>{c.SITE}</code>', 'Recognition site'],
        ['Fixed spacer', f'<code>{c.SPACER}</code>', 'Sixth base T fixes the retained 3′-A'],
        ['Random octamer', '<code>NNNNNNNN</code>', 'Priming sequence; no new counting UMI'],
    ])
    text = head(c.TITLE, home='index.html', include_header=False) + f'''
<div class="wrap">
<h1>{html.escape(c.TITLE)}</h1>
<p><a href="index.html">Overview</a> · <a href="{RAW}">Original input description (verbatim)</a></p>
{caveat('''<b>Proposed pilot design, not experimentally validated.</b> Enzyme amounts,
extension times and amplification cycles remain to be established for the input and
isolation system. Q/B handles require thermodynamic review. Inferred regions and
predicted reaction states are marked throughout.''')}
{info('''A PTA-free 5′ RNA-capture design using the existing oligo-dT and UMI-bearing
TSO, H-directed linear amplification, tailed random extension, selective PCR and
BciVI-directed FakeD ligation. FakeD replaces all four Atrandi barcoding rounds in
this one-ligation prototype. Cells remain separate through indexing; no additional
cell barcode is supplied.''')}
<h2>Compartment constraints</h2>
<p>DNA fragments larger than 400 bp cannot be removed from their compartments before
barcoding. Retained long templates and unwanted amplification products must be tolerated.
Exchange of short oligos and cleaved handles is assumed; the retention constraint does
not establish clearance efficiency. Products below the effective retention threshold may
also be lost. RNA-removal chemistry remains to be specified.</p>
<h2>Candidate oligos</h2>
<p>All sequences run 5′ → 3′. Use unmodified DNA with 5′-OH except where marked.
N means mixed A/C/G/T; W means A/T. Existing oligos are reused; their role in this pilot
remains untested.</p>
{oligos()}
<h3>Full-footprint P7-i7 primer architecture</h3>
<p><code>{raw_architecture}</code></p>
<p>This is an architecture, not a literal order containing random bases. Choose a defined
10-base i7 per cell/sample and a suitable index set for the pool. The four existing index
members are insufficient for more than four compartments. The full Read 2 arm has a
24-base FakeD annealing footprint, including its terminal T at the sealed A:T junction.</p>
<h2>Random-primer and cleavage geometry</h2>
{handles}
<p>Architecture: <code>Q–B–GTATCC–GACTGT–N8</code>. BciVI cuts at offsets (6/5)
beyond GTATCC. The cuts lie in the fixed spacer, outside the random octamer.
Discard the Q/B/site handle; its rejoining would compete with FakeD ligation.</p>
<p><a href="https://www.neb.com/en/products/r0596-bcivi">NEB BciVI specification</a>:
GTATCC(6/5), single-base 3′ extension, heat inactivation at 80 °C for 20 minutes.
NEB cautions that these ends can be more difficult to ligate than blunt fragments.
See also <a href="https://www.neb.com/en/tools-and-resources/usage-guidelines/heat-inactivation">NEB heat-inactivation guidance</a>.</p>
{panel(c.cleavage_map().rows(), caption='INFERRED — Cuts in random-primer orientation; the random octamer remains intact.')}
{panel(c.cut_scene().rows(), caption='INFERRED — BciVI target fragment; sense orientation is opposite to the random-primer orientation.')}
<p>Both nicks can seal: substrate 3′-A-OH to FakeD bottom 5′-P, and FakeD top 3′-T-OH
to restriction-generated substrate 5′-P. The TSO/target-primer end remains blunt and
unphosphorylated. Do not globally end-repair or dA-tail: that would erase the intended
end asymmetry. Proceed with single-base-overhang ligation for this amplified pilot;
more substrate does not establish the fraction successfully ligated.</p>
<h2>Proposed pilot workflow</h2>
<p>The drawings follow one representative 5′-captured interval. The cDNA placeholder has
no assumed molecular length. H-directed copies may extend beyond the illustrated random
priming site. The nominal GGG tail is shown for geometry; actual tail length may vary.</p>
<ol>
<li>Reverse-transcribe using the existing oligo-dT and modified UMI-bearing TSO.</li>
<li>Exchange residual RT oligos and remove RNA while retaining first-strand cDNA;
the compatible RNA-removal method remains unspecified.</li>
<li>Use H alone for TSO-directed linear amplification by thermal extension.</li>
<li>Exchange unused H primer.</li>
<li>Extend the candidate tailed random octamer back toward H.</li>
<li>Exchange unused tailed random primer.</li>
<li>Selective PCR with FakeTruseq-TSO and Q, using a blunt-producing polymerase and
unphosphorylated primers. Long inverted handles are intended to suppress handle-handle
products; suppression remains unvalidated.</li>
<li>Digest with BciVI, heat-inactivate at 80 °C for 20 minutes, then exchange the short
released handle. Internal restriction-site losses are accepted for this transcript pilot.</li>
<li>Ligate the pre-annealed original T-overhang FakeD duplex; exchange free adapter.</li>
<li>Index each compartment with P5 and a distinct full-footprint P7-i7 primer.</li>
<li>Pool only after indexing, then perform conventional purification, size assessment
and P5/P7 library quantification.</li>
</ol>
{workflow_panel(workflow)}
<h2>Predicted final library</h2>
{panel(c.duplex(c.final_library()).rows(), caption='INFERRED — Indexed library; I10 denotes a selected compartment index, not an oligo order.')}
{c.sp.section(c.final_library(), c.PRIMERS, required_roles=c.ROLES,
    intro='Planned acquired roles: custom Read 1, i7 Index 1 and standard Read 2. '
          'Run lengths remain to be selected; Index 2 is not acquired.')}
<h2>Sequencing and counting</h2>
<p>Read 2 begins immediately with the eight random-primer-derived bases. Trim those
8 bases before alignment; there is no additional fixed AG prefix. The terminal T in
the restored Read 2 primer site pairs with the substrate A at the sealed junction.</p>
<p>Custom Read 1 starts at the original RT N8 UMI, followed by the W2 spacer and the
template-switch G tail. Extract N8 as the counting UMI; W2 is a spacer. Use a parsing
rule that accommodates tail-length variation. The nominal model places cDNA after GGG;
that cycle boundary is conditional on the actual G-tail length.</p>
<p>Random priming creates no new counting UMI. Deduplicate technical descendants
within cell/sample and gene using the RT UMI. Demultiplex compartment identity using
the chosen i7 index. Verify intended directionality and handle-handle suppression
experimentally; total library yield alone cannot establish either.</p>
<h2>Shared viral-boundary branch</h2>
<p>The same tailed random primer and handle-cleavage/FakeD interface may be used after
outward-facing viral-primer linear amplification. Viral-target primers are not specified
by this RNA oligo list, so that branch is not modelled as a complete protocol here.</p>
<p>Products bearing the engineered cassette are intentionally cleaved there. Internal
BciVI sites can also cut viral/genomic sequence and lose boundary-containing products.
Check known viral target sequences and assess boundary loss separately from accepted
transcript losses. Additional amplification does not restore lost boundaries.</p>
</div></body></html>
'''
    OUT.write_text(text, encoding='utf-8')
    print(f'wrote {OUT} ({len(text.encode())} bytes)')


if __name__ == '__main__':
    build()
