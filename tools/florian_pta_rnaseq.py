"""
florian-PTA-rnaseq: single-cell RNA-seq read out through Atrandi's PTA/SPC workflow.

The idea is a graft of two protocols the repo already models:

* the **front end is Smart-seq3xpress** -- oligo-dT primes RT, the polymerase template-
  switches onto a TSO that carries a UMI, so the 5' end of each mRNA is tagged. Every
  sequence here is the same `AGAGACAGATTGCGCAATG` handle `smartseq` already defines.
* the **back end is Atrandi** -- but with the four rounds of split-pool barcoding
  replaced by a single pre-annealed duplex, `FakeD`, that is ligated on in one step.
  That is what "Fake D adaptor" means: it occupies the place of the round-D cassette
  without carrying any barcode.

So the molecule gets its read-1 side from the TSO (by PCR, with `FAKE_TRUSEQ_TSO`) and
its read-2 side from the ligated `FakeD` duplex. Only the TSO end is sequence-selected,
which is the point: `FAKE_TRUSEQ_TSO` amplifies TSO-carrying molecules and equips them
with a TruSeq handle, mimicking what Atrandi's own library prep would have produced.

Everything below is transcribed from the ordered oligos. What is *derived* -- the duplex
geometry, the finished library, the sizes -- is computed, and `selftest.py` pins it.
Open design questions are collected in `../01_primers.md`; three of them are consequential
and flagged inline here with `OPEN:`.
"""

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
sys.path.insert(0, str(NGS_101 / "atrandi-wgs__10.1101+2025.06.20.660799" / "tools"))

import atrandi
import illumina as il
import seqprimers as sp
import rt
from chemdraw import Construct, Scene, Segment, complement_segments, revcomp
from endprep import dA_tailed_scene, repair_and_dA_tail

# --------------------------------------------------------------- the ordered oligos
# Smartseq3expressTSO
#   /5BiosG/AGAGACAGATTGCGCAATG(N:25252525)(N)(N)(N)(N)(N)(N)(N)(W:5050)(W)rGrG+G
TSO_HANDLE = "AGAGACAG" + "ATTGCGCAATG"     # identical to smartseq.SS3_TSO_HANDLE
TSO_N_LEN = 8                                # (N:25252525) + 7x (N)
TSO_W_LEN = 2                                # (W:5050) + (W) -- A or T only
TSO_UMI_LEN = TSO_N_LEN + TSO_W_LEN          # 10 nt total
TSO_G_TAIL = rt.TSO_G_TAIL_LNA               # rGrG+G
TSO_5_MOD = "/5BiosG/"                       # 5' biotin, as in Smart-seq3

# Oligo dT (18T) -- note: no 5' handle at all, unlike Smart-seq3's oligo-dT
OLIGO_DT = "T" * 18

# FakeTruseq-TSO: PCR primer. Appends a TruSeq Read-1 handle to TSO-carrying cDNA.
FAKE_TRUSEQ_TSO = "ACACTCTTTCCCTACACGACGCTCTTCCGATCAGAGACAGATTGCGCAATG"

# FakeD: annealed duplex, ligated on in place of all four barcoding rounds.
FAKED_TOP = "TTCAGACGTGTGCTCTTCCGATCT"      # == il.TRUSEQ_READ2[-24:]
FAKED_BOT = "GATCGGAAGAGCACACGTCTGAA"       # == il.INDEX1_PRIMER[:23]

# -------------------------------------------------------------- derived, not transcribed
# The R1 handle FakeTruseq actually carries, and the canonical one it is one base short of.
FAKE_TRUSEQ_R1_PART = FAKE_TRUSEQ_TSO[: len(FAKE_TRUSEQ_TSO) - len(TSO_HANDLE)]
# OPEN (1): this is TRUSEQ_READ1 minus its terminal T. A 33-nt Read-1 sequencing primer
# therefore has a 3'-terminal mismatch against this library. See 01_primers.md.
R1_T_DROPPED = FAKE_TRUSEQ_R1_PART == il.TRUSEQ_READ1[:-1]

FAKED_DUPLEX_LEN = len(revcomp(FAKED_BOT))   # paired region
FAKED_OVERHANG = FAKED_TOP[FAKED_DUPLEX_LEN:]  # the 3'-T that takes the insert's dA
# As ordered (afake_adD_top_v2 / afake_adD_bottom_v2): ONLY the bottom strand is 5'-phosphorylated.
# That is exactly what TA ligation needs -- bottom 5'-P seals to the insert's dA 3'-OH, top's
# 3'-T seals to the insert's kinased 5'-P -- and the top's 5'-OH leaves FakeD's distal blunt
# end with no phosphate at all, so FakeD cannot blunt-ligate to itself.
FAKED_TOP_NAME, FAKED_BOT_NAME = "afake_adD_top_v2", "afake_adD_bottom_v2"
FAKED_TOP_MOD5, FAKED_BOT_MOD5 = "", "/5Phos/"
# OPEN (3): FakeD is a blunt duplex, not Y-shaped, so both ends of a fragment get the
# same handle and FAKED_TOP alone can amplify untagged cDNA. See 01_primers.md.

# What FakeD stands in for, on the Atrandi side: the read-2 arm plus the barcode block.
ATRANDI_READ2_ARM_LEN = 27                   # atrandi.ATRANDI_I7_ANNEAL
ATRANDI_BARCODE_BLOCK_LEN = 4 * 8 + 3 * 4    # 4 barcodes of 8, 3 linkers of 4
ATRANDI_D_SIDE_LEN = ATRANDI_READ2_ARM_LEN + ATRANDI_BARCODE_BLOCK_LEN

# The i7 indexing primer must graft back the 10 nt FakeD leaves off the Read-2 arm.
I7_GRAFT = il.TRUSEQ_READ2[: len(il.TRUSEQ_READ2) - len(FAKED_TOP)]   # "GTGACTGGAG"


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def tso() -> Construct:
    """The TSO as ordered, 5' biotin through the rGrG+G tail."""
    return Construct([
        _seg("ME 3' end", "AGAGACAG", "me"),
        _seg("Smart-seq3 tag", "ATTGCGCAATG", "tso"),
        _seg("UMI (N8)", "N" * TSO_N_LEN, "umi", placeholder=True),
        _seg("W2", "W" * TSO_W_LEN, "umi", placeholder=True,
             note="A or T only; Smart-seq3 itself stops at 8 N"),
        _seg("rGrG+G", "GGG", "r3", placeholder=True),
    ], name="Smartseq3expressTSO")


def faked_duplex() -> Construct:
    """FakeD, annealed. A 23-bp duplex with a single 3'-T on the top strand."""
    paired = revcomp(FAKED_BOT)
    return Construct([
        _seg("TruSeq Read 2 (3' 23 nt)", paired, "t7", bottom=FAKED_BOT[::-1]),
        _seg("", FAKED_OVERHANG, None,
             note="3'-T overhang; pairs with the insert's dA tail"),
    ], name="FakeD")


def final_library() -> Construct:
    """The predicted finished library, top strand, after a P5/P7 indexing PCR.

    The index PCR uses our own primers (atrandi.OUR_P5 / OUR_I7, shared with the scWGS
    protocol): a universal P5 onto the TruSeq Read-1 stub and a 10-nt UDP-indexed P7 onto
    the FakeD arm, which grafts back the I7_GRAFT that FakeD leaves off.
    """
    return Construct([
        _seg("Illumina P5", il.P5, "p5"),
        # P5 ends ACAC and TruSeq Read 1 begins ACAC: the four bases are shared, not
        # repeated -- which is why NEBNEXT_UNIVERSAL_PRIMER is 58 nt and not 62.
        _seg("TruSeq Read 1", FAKE_TRUSEQ_R1_PART[4:], "s5",
             note="one base short of canonical: see OPEN (1)"),
        _seg("ME 3' end", "AGAGACAG", "me"),
        _seg("Smart-seq3 tag", "ATTGCGCAATG", "tso"),
        _seg("UMI", "N" * TSO_N_LEN + "W" * TSO_W_LEN, "umi", placeholder=True),
        _seg("rGrG+G", "GGG", "r3", placeholder=True),
        _seg("cDNA", "X" * 20, "w1", placeholder=True),
        _seg("dA", "A", None, placeholder=True, inferred=True,
             note="pairs with FakeD's 3'-T; revcomp(FAKED_TOP) == 'A' + FAKED_BOT"),
        _seg("FakeD (Read 2)", FAKED_BOT, "t7"),
        _seg("Read 2 graft", revcomp(I7_GRAFT), "t7",
             note="installed by the i7 primer's 5' tail, not templated"),
        _seg("i7", "I" * atrandi.OUR_I7_INDEX_LEN, None, placeholder=True,
             note="10-nt IDT UDP i7, e.g. UDP0005-0008"),
        _seg("Illumina P7", il.P7_RC, "p7"),
    ], name="florian-PTA-rnaseq library")


# ------------------------------------------------- intermediates, step by step
# Each is a Construct so the page's schematics are derived from one definition rather
# than drawn by counting spaces. "X" is placeholder cDNA, "A"/"T" the dA-tail junction.

def first_strand() -> Construct:
    """After template switching: the tagged first-strand cDNA, still RNA:DNA."""
    return Construct([
        _seg("ME 3' end", "AGAGACAG", "me"),
        _seg("Smart-seq3 tag", "ATTGCGCAATG", "tso"),
        _seg("UMI", "N" * TSO_N_LEN + "W" * TSO_W_LEN, "umi", placeholder=True),
        _seg("rGrG+G", "GGG", "r3", placeholder=True),
        _seg("cDNA", "X" * 30, "w1", placeholder=True),
        _seg("polyA / oligo-dT", "A" * 10, None,
             note="top: the mRNA's polyA; bottom: the oligo-dT that primed RT "
                  "(18 T as ordered; drawn short to fit)"),
    ], name="tagged first-strand cDNA")


# ------------------------------------------------- step drawings, placed by pairing
# Every panel below is a chemdraw.Scene: strands are given 5'->3' and placed by which
# segment pairs with which. No column is typed by hand, and a Scene that would draw an
# unpaired column raises instead of rendering.

POLYA_LEN = 24          # drawn polyA tract
DT_LEN = 18             # the oligo-dT, as ordered
DT_SHIFT = 4            # where in the tract it happens to sit (any value 0..6 is valid)
_BODY = 30              # drawn transcript body


def _mrna() -> list[Segment]:
    return [_seg("body", "X" * _BODY, "w1", placeholder=True),
            _seg("polyA", "A" * POLYA_LEN)]


def _cdna() -> list[Segment]:
    """First strand 5'->3': oligo-dT, the polyA it copied, the body, untemplated CCC."""
    return [_seg("oligo-dT", "T" * DT_LEN, "r1"),
            _seg("copied polyA", "T" * DT_SHIFT, "r1"),
            *complement_segments([_mrna()[0]]),
            _seg("CCC", "CCC", "r3", note="untemplated, added by the RT off the 5' end")]


def scene_priming() -> Scene:
    sc = Scene()
    sc.strand("mRNA", _mrna())
    sc.anneal("oligo-dT", [_seg("oligo-dT", "T" * DT_LEN, "r1")], to="mRNA",
              pair=("oligo-dT", "polyA"), shift=DT_SHIFT)
    sc.mark("mRNA", "polyA", "polyA: no VN anchor, so the dT can sit anywhere in it")
    return sc


def scene_synthesis() -> Scene:
    sc = Scene()
    sc.strand("mRNA", _mrna())
    sc.anneal("cDNA", _cdna(), to="mRNA", pair=("oligo-dT", "polyA"), shift=DT_SHIFT)
    sc.mark("cDNA", "CCC", "untemplated CCC, added when the RT runs off the 5' end")
    sc.arrow("cDNA", "reverse transcription")
    return sc


def scene_switching() -> Scene:
    sc = Scene()
    sc.strand("mRNA", _mrna())
    sc.anneal("cDNA", _cdna(), to="mRNA", pair=("oligo-dT", "polyA"), shift=DT_SHIFT)
    sc.anneal("TSO", tso().segments, to="cDNA", pair=("rGrG+G", "CCC"))
    sc.mark("TSO", "rGrG+G", "rGrG+G pairs with the CCC; the RT switches template")
    sc.arrow("cDNA", "RT continues onto the TSO")
    return sc


def scene_dA_tailing() -> Scene:
    insert = Construct([_seg("cDNA", "X" * 20, "w1", placeholder=True)],
                       name="amplified cDNA")
    return dA_tailed_scene(repair_and_dA_tail(insert), label="")


def scene_faked() -> Scene:
    sc = Scene()
    sc.strand("top", [_seg("paired", FAKED_TOP[:-1], "t7"), _seg("T", FAKED_TOP[-1])],
              mod5="OH")
    sc.anneal("bottom", [_seg("bottom", FAKED_BOT, "t7")], to="top", pair=("bottom", "paired"),
              mod5="p")
    sc.mark("top", "T", "3'-T overhang: ligates to the insert's 5'-P")
    sc.note("bottom", "distal end blunt, 5'-OH only: cannot ligate, so no FakeD dimers")
    return sc


def scene_ligated() -> Scene:
    """A fragment with FakeD on both ends: the insert's dA is the first base of
    revcomp(FakeD top), so it is drawn once."""
    top = [_seg("FakeD", FAKED_TOP[:-1], "t7"), _seg("T", "T"),
           _seg("cDNA", "X" * 10, "w1", placeholder=True), _seg("dA", "A"),
           _seg("FakeD'", FAKED_BOT, "t7")]
    sc = Scene()
    sc.strand("top", top, label="")
    sc.anneal("bottom", complement_segments(top), to="top", pair=("cDNA'", "cDNA"), label="")
    sc.mark("top", "T", "T:A")
    sc.mark("top", "dA", "A:T")
    return sc


def scene_fragmentation() -> Scene:
    """The tagged cDNA (top strand) with fragmentation sites; only the first piece keeps the TSO."""
    cut = lambda i: _seg(f"cut{i}", "|", placeholder=True)
    body = lambda i: _seg(f"piece{i}", "X" * 10, "w1", placeholder=True)
    sc = Scene()
    sc.strand("cDNA", [_seg("ME 3' end", "AGAGACAG", "me"), _seg("tag", "ATTGCGCAATG", "tso"),
                       _seg("UMI", "N" * TSO_N_LEN + "W" * TSO_W_LEN, "umi", placeholder=True),
                       _seg("GGG", "GGG", "r3", placeholder=True), body(1), cut(1), body(2),
                       cut(2), body(3)], label="")
    sc.mark("cDNA", "ME 3' end", "only THIS fragment keeps the TSO", ch="-", through="piece1")
    return sc


def scene_selection() -> Scene:
    """FakeTruseq-TSO on a 5' fragment. It has the TOP strand's sequence, so it anneals to the
    BOTTOM strand; its 32-nt Read-1 part has no template there and is a 5' flap."""
    lig = ligated(True)
    sc = Scene()
    sc.strand("top", lig.segments, label="")
    sc.anneal("bottom", complement_segments(lig.segments), to="top",
              pair=("ME 3' end'", "ME 3' end"), label="")
    primer = [_seg("R1 tail", FAKE_TRUSEQ_R1_PART, "s5"), _seg("handle", TSO_HANDLE, "tso")]
    sc.anneal("FakeTruseq-TSO", primer, to="bottom", pair=("handle", "ME 3' end'"),
              unpaired=["R1 tail"], label="")
    sc.mark("FakeTruseq-TSO", "R1 tail", "5' flap: TruSeq Read 1, not templated")
    sc.mark("FakeTruseq-TSO", "handle", "anneals to the TSO handle")
    sc.arrow("FakeTruseq-TSO", "extension")
    return sc


def scene_read1() -> Scene:
    """The custom Read-1 primer on the finished library: it anneals to the bottom strand."""
    lib = final_library()
    sc = Scene()
    sc.strand("top", lib.segments, label="")
    sc.anneal("bottom", complement_segments(lib.segments), to="top",
              pair=("cDNA'", "cDNA"), label="")
    sc.anneal("Read 1 primer", [_seg("R1", FAKE_TRUSEQ_R1_PART, "s5"),
                                _seg("handle", TSO_HANDLE, "tso")],
              to="bottom", pair=("handle", "ME 3' end'"), label="")
    sc.arrow("Read 1 primer", "Read 1: UMI first")
    return sc


def scene_read1_junction() -> tuple[Scene, Scene]:
    """Read-1 primer site vs library, Smart-seq3 and this design, same-sense, for comparison."""
    ss3 = Scene()
    ss3.strand("Smart-seq3", [_seg("Nextera R1", "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", "s5"),
                              _seg("t1", "A", "tso"), _seg("tag", "TTGCGCAATG", "tso"),
                              _seg("UMI", "N" * 8, "umi", placeholder=True),
                              _seg("GGG", "GGG", "r3", placeholder=True)])
    ss3.mark("Smart-seq3", "Nextera R1", "stock Nextera R1, an EXACT prefix", ch="-")
    ss3.mark("Smart-seq3", "t1", "read 1 starts here")
    ours = Scene()
    ours.strand("this design", [_seg("R1", FAKE_TRUSEQ_R1_PART, "s5"), _seg("m1", "A", "me"),
                                _seg("ME", "GAGACAG", "me"), _seg("tag", "ATTGCGCAATG", "tso"),
                                _seg("UMI", "N" * TSO_N_LEN + "W" * TSO_W_LEN, "umi",
                                     placeholder=True),
                                _seg("GGG", "GGG", "r3", placeholder=True)])
    ours.mark("this design", "R1", "32 of the stock primer's 33 nt", ch="-")
    ours.mark("this design", "m1", "primer wants T here, library has A")
    return ss3, ours


def fragment(five_prime: bool = True) -> Construct:
    """One dA-tailed fragment after NEBNext fragmentation and end repair.

    five_prime=True is the fragment that happened to span the transcript's 5' end, so it
    still carries the TSO. Every other fragment is `five_prime=False` and carries nothing
    that distinguishes it.
    """
    head = ([_seg("ME 3' end", "AGAGACAG", "me"),
             _seg("Smart-seq3 tag", "ATTGCGCAATG", "tso"),
             _seg("UMI", "N" * TSO_N_LEN + "W" * TSO_W_LEN, "umi", placeholder=True),
             _seg("rGrG+G", "GGG", "r3", placeholder=True)]
            if five_prime else [])
    return Construct([
        _seg("dA (5' end)", "A", None, placeholder=True, bottom=" ",
             note="dA tail; unpaired, takes FakeD's 3'-T"),
        *head,
        _seg("cDNA", "X" * (18 if five_prime else 40), "w1", placeholder=True),
        _seg("dA (3' end)", "A", None, placeholder=True, bottom=" "),
    ], name=("5' fragment -- carries the TSO" if five_prime
             else "internal fragment -- carries nothing"))


def ligated(five_prime: bool = True) -> Construct:
    """After FakeD ligation. FakeD lands on BOTH ends, whatever the fragment is."""
    inner = [sg for sg in fragment(five_prime).segments
             if not sg.name.startswith("dA")]
    return Construct([
        _seg("FakeD (Read 2)", FAKED_TOP[:-1], "t7"),
        _seg("T:A junction", FAKED_TOP[-1], None,
             note="FakeD's 3'-T, ligated to the fragment's dA"),
        *inner,
        _seg("dA", "A", None, placeholder=True),
        _seg("FakeD (Read 2), other end", FAKED_BOT, "t7"),
    ], name=f"{'5-prime' if five_prime else 'internal'} fragment + FakeD both ends")


# ------------------------------------------------------------------ the workflow
# Confirmed order: RT (oligo-dT + TSO) -> PTA -> NEBNext enzymatic fragmentation ->
# dA-tailing -> ligate the annealed FakeD -> index PCR. So FakeD goes onto *fragments*,
# both ends of every one of them, and the index PCR is where TSO selection happens.
WORKFLOW = ("reverse transcription (oligo-dT + TSO)", "PTA amplification",
            "NEBNext enzymatic fragmentation", "dA-tailing",
            "ligate annealed FakeD", "index PCR")

# ------------------------------------------------- read-1 overhead, vs Smart-seq3
# Smart-seq3 reads with the stock Nextera read-1 primer, S5 + ME. The TSO's first 8 nt
# ARE the last 8 nt of ME, so the TSO completes the primer site exactly and AGAGACAG is
# never read. TruSeq has no ME, so the same 8 nt now fall inside read 1.
NEXTERA_ABSORBS = 8                              # len("AGAGACAG")
READ1_PREFIX_LEN = len(TSO_HANDLE) + TSO_UMI_LEN + 3          # 32
SS3_READ1_PREFIX_LEN = 11 + 8 + 3                             # tag + UMI8 + GGG = 22

# The one-base fix: with the T restored, the stock 33-nt TruSeq Read-1 primer becomes an
# exact prefix of the library instead of mismatching at its 3' base.
FAKE_TRUSEQ_TSO_FIXED = il.TRUSEQ_READ1 + TSO_HANDLE          # 52 nt

# The alternative: keep the 51-nt oligo and sequence with FakeTruseq-TSO itself as a custom
# read-1 primer. It is an exact match ending on the tag, so read 1 starts on the UMI; and the
# missing T makes it orthogonal to the stock primer (each mismatches the other's library).
CUSTOM_READ1_PRIMER = FAKE_TRUSEQ_TSO
READ1_PREFIX_LEN_CUSTOM = TSO_UMI_LEN + 3                      # UMI + GGG = 13


# Sequencing primers, by reference: sequences live in lib/ (or above, for the custom one).
SEQ_PRIMERS = [
    sp.custom("Read 1", "FakeTruseq-TSO, used as custom Read 1 (recommended)",
              CUSTOM_READ1_PRIMER, "this design -- see open question (1)",
              "Read 1 starts on the UMI. Load with the stock primer if pooling with scWGS."),
    sp.mismatching(sp.TRUSEQ["R1"], "its 3'-terminal T faces the A of AGAGACAG -- "
                                    "FakeTruseq-TSO lacks that T (open question 1)"),
    sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"],
]
