#!/usr/bin/env python3
"""Generate florian-PTA-rnaseq.html -- RNA-seq read out through Atrandi's PTA workflow."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
NGS_101 = ROOT / "NGS_101"
if not (NGS_101 / "lib").is_dir():
    raise SystemExit("missing NGS_101 dependency; run: git submodule update --init")
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(NGS_101 / "lib"))
sys.path.insert(0, str(NGS_101 / "smart-seq__10.1038+nbt.2282" / "tools"))
sys.path.insert(0, str(NGS_101 / "atrandi-wgs__10.1101+2025.06.20.660799" / "tools"))

import atrandi
import florian_pta_rnaseq as f
import seqprimers as sp
import illumina as il
import smartseq as ss
from chemdraw import (Row, annotation_rows, oligo, panel, primer_row, revcomp,
                      strand_row, tm)
from florian_pta_rnaseq import _seg as seg

# 3' nt of our i7 primer that pair with FakeD top; the rest is the graft
I7_FOOT = len(atrandi.ATRANDI_I7_ANNEAL) - len(f.I7_GRAFT)
from page import MD_STYLE, caveat, head, info, legend, table

OUT = HERE.parent / "florian-PTA-rnaseq.html"
PRE = 4


def chunks(*cs_):
    return [tuple(list(c) + [None, False])[:3] for c in cs_]


def row(*cs_, prefix="", suffix=""):
    return Row(chunks=chunks(*cs_), prefix=prefix, suffix=suffix)


def con_panel(con, cls="long"):
    rows = [strand_row(con, "top"), strand_row(con, "bottom")] + \
        annotation_rows(con, prefix_width=PRE)
    return panel(rows, cls=cls)


def preamble() -> str:
    return f"""<div class="wrap">
<h1>florian-PTA-rnaseq &mdash; RNA-seq through Atrandi's PTA workflow</h1>

{info("""An attempt to read single-cell mRNA out through the Atrandi SPC / PTA pipeline.
The design is a graft of two protocols modelled in NGS 101: a
<a href="https://henriksson-lab.github.io/NGS_101/smart-seq3xpress__10.1038+s41587-022-01311-4/smart-seq3xpress.html"><b>Smart-seq3xpress</b></a> front end that tags each
mRNA's 5' end with a UMI, and an <a href="https://henriksson-lab.github.io/NGS_101/atrandi-wgs__10.1101+2025.06.20.660799/atrandi_wgs.html"><b>Atrandi</b></a>
back end &mdash; except that <b>all four rounds of split-pool barcoding are replaced by a
single pre-annealed duplex</b>, <code>FakeD</code>, ligated on in one step.""")}

<p class="research-notes"><a href="index.html">Overview</a> &middot;
<a href="01_primers.html">Research notes</a> &middot;
<a href="https://henriksson-lab.github.io/NGS_101/">NGS 101</a></p>

{legend("""<b>Where each side of the molecule comes from.</b> The
<s5>read-1 side</s5> is installed <i>by PCR</i>, off the <tso>TSO</tso> sequence, so only
TSO-carrying molecules get it. The <t7>read-2 side</t7> is <i>ligated</i> on as
<code>FakeD</code>, which carries no barcode. That asymmetry is the whole design: one end
is sequence-selected, the other is not.""")}

{panel([row((" -> ".join(f.WORKFLOW[:3]),)), row(("   -> " + " -> ".join(f.WORKFLOW[3:]),))],
       cls="small", caption="The confirmed order of operations.")}

{caveat("""<b>Status: oligos and workflow order.</b> PTA runs on first-strand cDNA with
<b>both random hexamer and a TSO-handle spike-in</b>. The spike-in copies the tagged 3' end
that random priming alone would lose (Part 2). Primer concentration and terminator chemistry
remain to be specified. The
index PCR uses our own primers (shared with the scWGS protocol), so the final library
below is <i>predicted</i> from them. Everything marked inferred is inferred. The three
consequential design questions are at the bottom of this page and in
<a href="01_primers.html"><code>01_primers.md</code></a>.""")}
"""


def oligos_section() -> str:
    rows = [
        oligo("Smartseq3expressTSO",
              [seg("", "AGAGACAG", "me"), seg("", ss.SS3_TSO_TAG, "tso"),
               seg("", "N" * f.TSO_N_LEN, "umi", placeholder=True),
               seg("", "W" * f.TSO_W_LEN, "umi", placeholder=True),
               seg("", "rGrG+G", "r3", placeholder=True)], mods=f.TSO_5_MOD),
        oligo("Oligo dT (18T)", [seg("", f.OLIGO_DT, None)]),
        oligo("FakeTruseq-TSO",
              [seg("", f.FAKE_TRUSEQ_R1_PART, "s5"), seg("", "AGAGACAG", "me"),
               seg("", ss.SS3_TSO_TAG, "tso")]),
        oligo("FakeD top", [seg("", revcomp(f.FAKED_BOT), "t7"),
                            seg("", f.FAKED_OVERHANG, None)]),
        oligo("FakeD bottom", [seg("", il.STEM, "t7"),
                               seg("", f.FAKED_BOT[len(il.STEM):], "t7")]),
    ]
    return f"""<h2>Oligos</h2>
<seq>
{chr(10).join(rows)}
</seq>
{info(f"""Reused rather than retyped: the TSO's <me>AGAGACAG</me> is the 3' end of the Tn5
mosaic end and <tso>{ss.SS3_TSO_TAG}</tso> is Smart-seq3's own 11-nt tag &mdash; together
they are exactly <code>smartseq.SS3_TSO_HANDLE</code>. <code>FakeD bottom</code> opens with
the canonical <t7>{il.STEM}</t7> stem. <b>The one sequence that is not canonical is
<code>FakeTruseq-TSO</code>'s <s5>read-1 handle</s5></b>, which is one base short; see
below.""")}
"""


def part1_rt() -> str:
    """Steps 1-5: oligo-dT priming through the tagged first strand."""
    t = f.tso()
    prime = panel(f.scene_priming().rows(), cls="small",
                  caption="(2) Oligo-dT priming. 18 T and no handle, so nothing marks the "
                          "3' end -- and with no VN anchor the primer can sit anywhere in "
                          "the polyA tract.")
    synth = panel(f.scene_synthesis().rows(), cls="small",
                  caption="(3) First-strand synthesis. On reaching the 5' end the RT adds a "
                          "few untemplated C -- the overhang the TSO needs.")
    switch = panel(f.scene_switching().rows(), cls="small",
                   caption="(4) Template switching. The UMI is attached here -- before any "
                           "amplification, and only at a genuine transcript 5' end.")
    return f"""<h2>Part 1 &mdash; tagging the 5' end</h2>
<h3>(1) The TSO, as ordered</h3>
{con_panel(t, cls="small")}
{info(f"""{len(t)} nt. <b>The UMI is 10 nt, not Smart-seq3's 8</b> &mdash; 8 hand-mixed
<umi>N</umi> then two <umi>W</umi> (A or T only). The W positions look deliberate: they
stop the UMI ending in G, which would otherwise be indistinguishable from the
template-switch <r3>GGG</r3> and leave the UMI boundary ambiguous. Every downstream offset
depends on this being 10. The 5' biotin is Smart-seq3's &mdash; it suppresses the TSO
concatemers that otherwise dominate a template-switching reaction.""")}

<h3>(2&ndash;4) Priming, synthesis, template switching</h3>
{prime}
{synth}
{switch}

<h3>(5) The tagged first-strand cDNA</h3>
{con_panel(f.first_strand())}
{info("""Only the <b>5' end</b> of the transcript is marked, and it is marked once, with a
<umi>UMI</umi> that was never copied. The 3' end carries nothing: the oligo-dT has no
handle, so that side is not addressable by sequence and has to be supplied later by
ligation. That asymmetry is what the rest of the protocol is built around.""")}
"""


def part2_pta() -> str:
    return f"""<h2>Part 2 &mdash; amplification</h2>
<h3>(6) PTA with random hexamer + TSO-handle spike-in</h3>
<seq>
{oligo("PTA TSO-handle spike-in", [seg("constant TSO handle", f.PTA_TSO_SPIKE, "tso")])}
{oligo("PTA random hexamer", [seg("random hexamer", f.PTA_RANDOM_HEXAMER, None,
                                  placeholder=True)])}
</seq>
{panel(f.scene_pta_priming().rows(), cls="long",
       caption="Both primer classes act in PTA: random hexamers cover the transcript body; "
               "the constant-handle spike-in recovers the tagged 5' end.")}
{info("""The PTA primer pool contains <b>both</b> species. Random hexamer provides broad
transcript coverage. The 19-nt constant-handle primer copies the TSO-tagged end and carries
two protected bonds at its 3' end, like the PTA random primer. It is the constant handle,
not the full UMI-bearing TSO: a newly drawn random UMI would not match the UMI already
written into that first strand.""")}
{caveat("""<b>The spike-in is required for a 5' assay.</b>
Primary Template-directed Amplification uses a strand-displacing polymerase with
irreversible terminators, so products come off the <i>original</i> template rather than
off each other &mdash; quasi-linear instead of exponential, which is why it is so much more
even than MDA on single cells. But random priming alone does not preserve this end. If the
TSO-handle primer is omitted, this should no longer be treated as a 5'-capture design;
the alternative is to redesign selection around the transcript's 3' end.""")}
{info("""<b>Why random hexamer alone is insufficient.</b> The first strand carries
the tag only as its complement, at its <b>3' end</b>. A polymerase copies from a primer
toward the template's 5' end, so every random-primed copy <i>starts</i> at its primer and
never includes the bases 3' of it. No copy carries the <tso>handle</tso>, unless a hexamer
happens to land in the last few bases. Copies of copies don't help, because their 5' ends
are random primer sites. On a sense second strand the tag would sit at the <i>5'</i> end,
and every primer within one amplicon length would run into it. The explicit
<code>{f.PTA_TSO_SPIKE}</code> spike-in solves this by priming that terminal handle directly
and making a sense copy that random hexamers can subsequently cover. Verify preservation
with tagged:internal qPCR before and after PTA.""")}
"""


def part3_faked() -> str:
    duplex = f.scene_faked().rows()
    frag = panel([Row(chunks=[("~~~~~ enzymatic fragmentation ~~~~~", None, False)]), Row(),
                  *f.scene_fragmentation().rows()],
                 cls="small", caption="(7) NEBNext enzymatic fragmentation. One cDNA gives many "
                                      "fragments; exactly one of them spans the tagged 5' end.")
    tail = panel(f.scene_dA_tailing().rows(), cls="small",
                 caption="(8) End repair and dA-tailing. Each 3' end gets a single unpaired A "
                         "-- the partner for FakeD's 3'-T.")
    lig = panel(f.scene_ligated().rows(), cls="small",
                caption="(10) After ligation. The same adapter is now on BOTH ends of every "
                        "fragment -- there is no Y shape to make them differ.")
    return f"""<h2>Part 3 &mdash; faking the barcode cassette</h2>
<h3>(7) NEBNext enzymatic fragmentation</h3>
{frag}
<h3>(8) End repair and dA-tailing</h3>
{tail}

<h3>(9) FakeD: bottom strand phosphorylated, then annealed</h3>
{panel(duplex, cls="small",
       caption="A plain dA-tail-compatible adapter: 23 bp paired, one unpaired 3'-T.")}
{info(f"""As ordered, <b>only <code>{f.FAKED_BOT_NAME}</code> carries a 5'-phosphate</b>
&mdash; and that is the right choice. Its 5'-P seals to the insert's dA 3'-OH. The top's
3'-T seals to the insert's 5'-P, which end-prep kinased. Both nicks close. The top's
5'-OH sits at the distal blunt end, where it leaves <b>no phosphate at all</b>, so FakeD
cannot blunt-ligate to itself. Two T overhangs can't pair either, so no adapter dimers. Tm of the paired region is about
{tm(revcomp(f.FAKED_BOT)):.0f}&nbsp;&deg;C, so it stays annealed at ligation temperature.
<code>FakeD top</code> is exactly the 3' 24 nt of <code>TRUSEQ_READ2</code>;
<code>FakeD bottom</code> is exactly the first 23 nt of the i7 index-read primer. Both
checked against <code>lib/illumina.py</code>.""")}

<h3>(10) Ligation &mdash; onto both ends of every fragment</h3>
{lig}
{con_panel(f.ligated(True))}
{con_panel(f.ligated(False))}
{info("""Two products, differing only in whether the fragment happened to carry the
<tso>TSO</tso>. Both have the identical <t7>FakeD</t7> arm at each end. Nothing has
selected anything yet &mdash; that happens at the index PCR.""")}

<h3>(11) What FakeD replaces</h3>
{panel([
    row(("[TruSeq Read 2, 27 nt]", "t7"), ("[BC-D]", "cbc"), ("[L]", "r2"),
        ("[BC-C]", "cbc"), ("[L]", "r2"), ("[BC-B]", "cbc"), ("[L]", "r2"),
        ("[BC-A]", "cbc"), prefix="Atrandi: ", suffix="  -- 71 nt, 4 rounds"),
    Row(),
    row(("[TruSeq Read 2, 24 nt]", "t7"),
        prefix="FakeD:   ", suffix="  -- 23 bp, 1 ligation, no barcode"),
], cls="small", caption="The D side of an Atrandi library, and its stand-in.")}
{table(["", "Atrandi, 4 rounds", "FakeD, 1 ligation"],
       [["read-2 arm", f"{f.ATRANDI_READ2_ARM_LEN} nt",
         f"{len(f.FAKED_TOP)} nt &mdash; 10 nt short"],
        ["barcodes", f"4 &times; {atrandi.BARCODE_LEN} nt", "<b>none</b>"],
        ["linkers", f"3 &times; {atrandi.LINKER_LEN} nt", "none"],
        ["total", f"<b>{f.ATRANDI_D_SIDE_LEN} nt</b>", f"<b>{f.FAKED_DUPLEX_LEN} bp</b>"],
        ["shape", f"Y &mdash; {len(atrandi.LIGADAPT_ARM)}-nt ss 3' arm",
         "<b>blunt duplex</b>"]])}
{info(f"""Because <code>FakeD</code> stops 10 nt short of the full Read-2 arm, the i7
primer has to <b>graft <code>{f.I7_GRAFT}</code> back on</b> as a non-templated 5' tail.
Our i7 primer is truncated like Atrandi's, so only its 3' {I7_FOOT} nt pair with FakeD
(Tm &asymp; {tm(f.FAKED_TOP[:I7_FOOT]):.0f}&nbsp;&deg;C, below Atrandi's 54&nbsp;&deg;C anneal &mdash;
the first cycles will be inefficient until products carry the full 27-nt site). The same trick the
Zhang v2 adaptor uses, drawn in
<a href="https://github.com/henriksson-lab/NGS_101/blob/main/lenticrispr-v1-screening__10.1126%2Bscience.1247005/90_combined_screen_research.md">pooled-screen research notes</a>.""")}
"""


def part4_index() -> str:
    sel = panel(f.scene_selection().rows(), cls="long",
                caption="(12) The selective step. The forward primer only has a site on "
                        "fragments that kept the TSO.")
    return f"""<h2>Part 4 &mdash; index PCR</h2>
<h3>(12) FakeTruseq-TSO selects the tagged fragments</h3>
{sel}
{info("""The forward primer's 3' 19 nt anneal to the <tso>TSO handle</tso>; its 5' 32 nt
are <s5>TruSeq Read 1</s5> and are <b>not templated</b> &mdash; they are installed, exactly
the way the i7 primer installs the missing 10 nt of the Read-2 arm. An internal fragment
has no site for this primer, so it never acquires a read-1 handle and never acquires
<p5>P5</p5>.""")}

<h3>(13) The final library</h3>
{con_panel(f.final_library())}
{info(f"""Built with our index primers (<code>atrandi-wgs__10.1101+2025.06.20.660799/ref/our_index_primers.tsv</code>)
&mdash; Atrandi's universal P5 primer onto the <s5>read-1 stub</s5> and a P7 primer carrying a
{atrandi.OUR_I7_INDEX_LEN}-nt IDT UDP i7 index onto the <t7>FakeD arm</t7>. Note <p5>P5</p5> ends in <code>ACAC</code> and TruSeq Read 1
begins with <code>ACAC</code>: those four bases are <b>shared, not repeated</b>, which is
why the universal primer is 58&nbsp;nt and not 62.""")}

<h3>(14) ...and what happens to the internal fragments</h3>
{con_panel(f.ligated(False))}
{info("""They are amplifiable &mdash; <code>FakeD top</code> matches the 5' end of both
strands &mdash; but they are <b>dead ends</b>: no read-1 handle, so no <p5>P5</p5>, and a
molecule with <p7>P7</p7> at both ends cannot form a cluster. So they do not reach the
data. They do consume the PCR, which is open question (3).""")}
"""


def sequencing() -> str:
    lib = f.final_library()
    return f"""<h2>Sequencing</h2>
<h3>(15) Read 1 with the custom primer</h3>
{panel(f.scene_read1().rows(), cls="long",
       caption="The custom Read-1 primer (FakeTruseq-TSO) anneals to the bottom strand and "
               "ends on the tag, so Read 1 starts on the UMI.")}
{caveat(f"""<b>Read 1 spends {f.READ1_PREFIX_LEN} cycles before it reaches any cDNA</b>
&mdash; <me>AGAGACAG</me> 8, <tso>tag</tso> 11, <umi>UMI</umi> {f.TSO_UMI_LEN} and
<r3>GGG</r3> 3. Smart-seq3 spends only {f.SS3_READ1_PREFIX_LEN}, because its Nextera primer
hides <me>AGAGACAG</me> inside itself. At 50 cycles that is
{50 - f.READ1_PREFIX_LEN} bases of transcript instead of {50 - f.SS3_READ1_PREFIX_LEN}.
<b>Read 2 has no overhead at all</b>, so it is the one to make long if coverage matters;
read 1 only has to reach far enough past the UMI to map. With the stock TruSeq primer
instead (open question 1), add 19 cycles.""")}
{sp.section(f.final_library(), f.SEQ_PRIMERS,
            content_primers=tuple(p for p in f.SEQ_PRIMERS if not p.expect_mismatch))}
"""


def open_questions() -> str:
    _j = f.scene_read1_junction()
    return f"""<h2>Open questions</h2>
{panel([*_j[0].rows(), Row(), *_j[1].rows()], cls="small", caption="The same junction, both ways. Smart-seq3's TSO completes the "
                        "read-1 primer site exactly; the TruSeq version is one base short.")}
{caveat(f"""<b>(1) <code>FakeTruseq-TSO</code> is one base short of TruSeq Read 1 &mdash;
checked against how Smart-seq3 solves the same problem.</b> Smart-seq3 reads with the stock
Nextera primer <code>S5 + ME</code>, and the TSO's first 8 nt <b>are</b> the last 8 nt of
<me>ME</me> &mdash; so the TSO <i>completes</i> the primer site, <me>AGAGACAG</me> is never
read, and read 1 starts on the tag. That is why the Smart-seq3 TSO begins with
<me>AGAGACAG</me> at all. The TruSeq analogue needs the library to present all 33 nt of
<code>TRUSEQ_READ1</code>; it presents 32, so the stock primer's 3'-terminal <b>T faces an
A</b> &mdash; a mismatch at the one position that stops extension. If it extends anyway,
read 1 starts a base late and the <tso>{ss.SS3_TSO_TAG}</tso> tag lands at <b>offset 7, not
8</b>. The T is not arbitrary: in TruSeq it <i>is</i> the dA-ligation overhang base, and
this handle is installed by PCR, so there is no dA and no natural T &mdash; understandable,
but the stock primer does not know that. <b>Fix: one base.</b>
<code>TRUSEQ_READ1 + TSO handle</code> = {len(f.FAKE_TRUSEQ_TSO_FIXED)} nt, and the stock
primer is then an exact prefix. <b>Or, better, keep the 51-nt oligo and sequence with it</b>
&mdash; <code>FakeTruseq-TSO</code> is itself an exact custom read-1 primer, ending on the tag,
so read 1 starts on the <umi>UMI</umi>: {f.READ1_PREFIX_LEN_CUSTOM} cycles before cDNA instead
of {f.READ1_PREFIX_LEN}, and a random rather than constant start. The missing T makes it nearly
orthogonal to the stock primer, so the two can share NovaSeq&nbsp;X well CP1 when pooled with
the scWGS libraries (a custom primer there <i>replaces</i> Illumina's). Do not do both: with the
T added, the stock primer matches too and competes 19 bases out of register.""")}
{caveat("""<b>(1c) Read 1 starts 100% low-diversity with the stock primer.</b> TSO selection
means <i>every</i> cluster begins with the same 19 nt <me>AGAGACAG</me><tso>ATTGCGCAATG</tso>.
Smart-seq3 tolerates its 11-nt tag because most of its reads are internal; here none are.
Sequenced alone that fails template registration. Dilute heavily with scWGS libraries or PhiX,
or use the custom primer above.""")}
{info("""<b>Resolved by our index primers: Read 2 and Index 1.</b> Our i7 primer's 5' tail
grafts back <code>GTGACTGGAG</code>, so the stock TruSeq Read-2 and i7 index primers both find
complete sites. Read 2 starts on cDNA right after the dA. <b>Quantify by P5/P7 qPCR, not
Qubit:</b> the i7 primer's 3' 17&nbsp;nt are <code>FakeD top</code>'s 5' 17&nbsp;nt, so it
alone amplifies FakeD&ndash;FakeD fragments, which never cluster but do count by mass.""")}
{caveat(f"""<b>(1b) A larger, separate cost of the Nextera &rarr; TruSeq swap: 10 extra
read-1 cycles.</b> Smart-seq3 spends {f.SS3_READ1_PREFIX_LEN} cycles before cDNA (tag 11 +
UMI 8 + GGG 3); this design spends <b>{f.READ1_PREFIX_LEN}</b> (<me>AGAGACAG</me> 8 + tag
11 + UMI 10 + GGG 3). The difference is {f.NEXTERA_ABSORBS} cycles because TruSeq has no
<me>ME</me> to hide <me>AGAGACAG</me> inside the sequencing primer, plus
{f.TSO_UMI_LEN - ss.SS3_UMI_LEN} from the longer UMI. <b>Adding the T does not recover
them</b> &mdash; only a custom read-1 primer ending <code>...CCGATCTAGAGACAG</code>
would.""")}
{info("""<b>(2) &mdash; resolved: only FakeD bottom is phosphorylated, which is the right
design.</b> Both junctions with the insert still seal: bottom 5'-P to the insert's dA, and
top 3'-T to the insert's kinased 5'-P. The top's 5'-OH leaves FakeD's distal blunt end
unligatable, doing the job of Atrandi's <code>/5AmMC6/</code> block, so FakeD cannot
dimerise. Two T overhangs can't pair either. Residual risk, if a ~46&nbsp;bp peak appears
anyway: kinase carried over from end-prep, which the 65&nbsp;&deg;C step should have
inactivated.""")}
{caveat("""<b>(3) A blunt duplex puts the same handle on both ends.</b> With
fragmentation and dA-tailing before ligation, <code>FakeD</code> goes onto <i>both ends of
every fragment</i> &mdash; and working through the geometry, such a fragment reads
<code>FakeD top &middot; insert &middot; dA &middot; revcomp(FakeD top)</code>, so
<code>FakeD top</code> matches the 5' end of <b>both</b> strands and amplifies it
exponentially on its own. Fragmentation makes this worse: one cDNA gives many fragments and
only one carries the TSO. <b>The index PCR does filter them</b> &mdash; no read-1 handle
means no P5, and P7 at both ends does not cluster &mdash; so the <i>final library is
clean</i>. The cost is upstream: that PCR's product is <b>dominated by FakeD&ndash;FakeD</b>,
so Qubit, TapeStation or qPCR on it measures mostly material that gets discarded, and any
cycle count or input mass taken from those numbers is wrong by an unknown ratio.
<b>Testable in one reaction:</b> run it with <code>FakeD top</code> alone and compare
yield. If it is comparable, give <code>FakeD</code> a Y shape &mdash; NEBNext's own adapter
is a USER-cut hairpin for exactly this reason &mdash; or block <code>FakeD top</code> from
priming.""")}
{info("""Still unsettled: the spike-in concentration, the PTA terminator chemistry, and
whether the 10-nt UMI length is intended. The primer logic is now explicit: random hexamer
for coverage plus TSO-handle spike-in for 5' retention. Without the latter, pursue a 3'
capture design instead.""")}
</div>
"""


def main() -> None:
    parts = [head("florian-PTA-rnaseq", home="index.html", include_header=False),
             MD_STYLE, preamble(), oligos_section(),
             part1_rt(), part2_pta(), part3_faked(), part4_index(), sequencing(),
             open_questions()]
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
