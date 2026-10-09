#!/usr/bin/env python3
"""
Self-test for florian-PTA-rnaseq.

Every identity asserted in ../01_primers.md is encoded here, including the three that
are *problems* -- the dropped T, the missing modifications, the non-Y adapter. Those are
checked as facts about the oligos as ordered, so that if someone re-orders a corrected
oligo the test fails and says so rather than quietly agreeing.

Run:  python3 tools/selftest.py
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
sys.path.insert(0, str(NGS_101 / "smart-seq__10.1038+nbt.2282" / "tools"))
sys.path.insert(0, str(NGS_101 / "atrandi-wgs__10.1101+2025.06.20.660799" / "tools"))

import atrandi
import florian_pta_rnaseq as f
import illumina as il
import rt
import smartseq as ss
from checks import Check, run_common
from chemdraw import bases_pair, revcomp, tm

check = Check()
run_common(check)

check.section("the TSO is Smart-seq3's, reused not retyped")
check("its 5' handle is exactly smartseq's SS3_TSO_HANDLE",
      f.TSO_HANDLE, ss.SS3_TSO_HANDLE)
check("which is the ME 3' end followed by the 11-nt tag",
      f.TSO_HANDLE, "AGAGACAG" + ss.SS3_TSO_TAG)
check("the handle is 19 nt", len(f.TSO_HANDLE), 19)
check("the ME 3' end really is the last 8 nt of the Tn5 mosaic end",
      f.TSO_HANDLE[:8], __import__("nextera").ME[-8:])
check("the 3' tail is rt's LNA G-tail", f.TSO_G_TAIL, rt.TSO_G_TAIL_LNA)
check("the tail is two riboG plus one LNA G", f.TSO_G_TAIL, "rGrG+G")
check("it is 5'-biotinylated, as Smart-seq3's TSO is", f.TSO_5_MOD, "/5BiosG/")
check("the UMI is 8 N plus 2 W", (f.TSO_N_LEN, f.TSO_W_LEN), (8, 2))
check("...so 10 nt, which is 2 MORE than Smart-seq3's own UMI",
      (f.TSO_UMI_LEN, ss.SS3_UMI_LEN), (10, 8))
check("the W positions sit between the UMI and the G-tail, so the UMI cannot end in G",
      f.tso().top().index("W") > f.tso().top().rindex("N"))
check("the TSO is 32 nt as ordered", len(f.tso()), 32)

check.section("the oligo-dT is bare -- no handle, unlike Smart-seq3's")
check("18 T and nothing else", f.OLIGO_DT, "T" * 18)
check("it carries none of Smart-seq3's oligo-dT handle",
      ss.SS3_OLIGO_DT_HANDLE not in f.OLIGO_DT)
check("and no VN anchor, so it can slip within polyA", not f.OLIGO_DT.endswith("VN"))

check.section("FakeTruseq-TSO: a TruSeq handle grafted onto the TSO sequence")
check("it is 51 nt", len(f.FAKE_TRUSEQ_TSO), 51)
check("its 3' end is the TSO handle, so it primes TSO-carrying cDNA",
      f.FAKE_TRUSEQ_TSO.endswith(f.TSO_HANDLE))
check("its 5' part is 32 nt", len(f.FAKE_TRUSEQ_R1_PART), 32)
check("OPEN (1): that is TruSeq Read 1 MINUS its terminal T", f.R1_T_DROPPED)
check("...i.e. it is NOT canonical TruSeq Read 1",
      f.FAKE_TRUSEQ_R1_PART != il.TRUSEQ_READ1)
check("...and the canonical primer is exactly one base longer",
      len(il.TRUSEQ_READ1) - len(f.FAKE_TRUSEQ_R1_PART), 1)
check("so the junction reads CCGATC|AGAGACAG, with no T",
      "CCGATCAGAGACAG" in f.FAKE_TRUSEQ_TSO)
check("a 33-nt Read-1 primer would mismatch at its 3' base: template offers T, primer T",
      il.TRUSEQ_READ1[-1] == "T" and f.FAKE_TRUSEQ_TSO[len(f.FAKE_TRUSEQ_R1_PART)] == "A")
check("it is the TruSeq analogue of smartseq's Nextera forward PCR primer, which ends "
      "in the same tag",
      ss.SS3_FWD_PCR.endswith(ss.SS3_TSO_TAG) and f.FAKE_TRUSEQ_TSO.endswith(ss.SS3_TSO_TAG))
check("both forward primers end in AGAGACAG + tag; only the 5' handle differs",
      ss.SS3_FWD_PCR[-19:], f.FAKE_TRUSEQ_TSO[-19:])

check.section("claim (1) tested against how Smart-seq3 handles the SAME junction")
import nextera as nx  # noqa: E402

check("Smart-seq3 sequences with the stock Nextera read-1 primer, S5 + ME",
      nx.READ1_PRIMER, nx.S5 + nx.ME)
check("the TSO's first 8 nt ARE the last 8 nt of ME", f.TSO_HANDLE[:8], nx.ME[-8:])
SS3_LIB = nx.S5 + nx.ME + ss.SS3_TSO_TAG + "N" * ss.SS3_UMI_LEN + "GGG"
check("so in Smart-seq3 the TSO COMPLETES the primer site exactly: the stock primer "
      "is a prefix of the library", SS3_LIB.startswith(nx.READ1_PRIMER))
check("...AGAGACAG therefore sits inside the primer and is never read",
      SS3_LIB[len(nx.READ1_PRIMER):].startswith(ss.SS3_TSO_TAG))

FL_LIB = f.FAKE_TRUSEQ_R1_PART + f.TSO_HANDLE + "N" * f.TSO_UMI_LEN + "GGG"
check("the TruSeq version does NOT line up: the stock primer is not a prefix",
      not FL_LIB.startswith(il.TRUSEQ_READ1))
check("the library offers 32 of the primer's 33 nt",
      FL_LIB.startswith(il.TRUSEQ_READ1[:-1]) and len(f.FAKE_TRUSEQ_R1_PART) == 32)
check("the primer's 3'-terminal base is T", il.TRUSEQ_READ1[-1], "T")
check("...and it faces an A, so the mismatch is at the 3' TERMINUS", FL_LIB[32], "A")
check("if it extends anyway the tag lands at offset 7, not 8",
      (FL_LIB.index(ss.SS3_TSO_TAG) - 33, FL_LIB.index(ss.SS3_TSO_TAG) - 32), (7, 8))
check("the T exists in TruSeq as the dA-ligation overhang base, which a PCR-installed "
      "handle has no counterpart for", il.NEBNEXT_ARM_READ1.endswith("T"))
check("the one-base fix is 52 nt", len(f.FAKE_TRUSEQ_TSO_FIXED), 52)
check("...and with it the stock primer IS an exact prefix",
      (f.FAKE_TRUSEQ_TSO_FIXED + "N").startswith(il.TRUSEQ_READ1))
check("the fix is exactly the ordered oligo plus one T",
      f.FAKE_TRUSEQ_TSO_FIXED,
      f.FAKE_TRUSEQ_TSO[:32] + "T" + f.FAKE_TRUSEQ_TSO[32:])

check.section("the second, larger cost: read-1 overhead vs Smart-seq3")
check("Smart-seq3 spends 22 cycles before cDNA (tag 11 + UMI 8 + GGG 3)",
      f.SS3_READ1_PREFIX_LEN, 22)
check("this design spends 32 (AGAGACAG 8 + tag 11 + UMI 10 + GGG 3)",
      f.READ1_PREFIX_LEN, 32)
check("the difference is 10 cycles", f.READ1_PREFIX_LEN - f.SS3_READ1_PREFIX_LEN, 10)
check("...which decomposes as 8 (AGAGACAG now read) + 2 (longer UMI)",
      f.NEXTERA_ABSORBS + (f.TSO_UMI_LEN - ss.SS3_UMI_LEN), 10)
check("adding the T does NOT recover those 8 cycles -- AGAGACAG is in the molecule "
      "either way", f.FAKE_TRUSEQ_TSO_FIXED.endswith(f.TSO_HANDLE))

check.section("claim (3): FakeD top primes BOTH ends of a FakeD-FakeD fragment")
check("the confirmed workflow has fragmentation and dA-tailing before ligation",
      [w for w in f.WORKFLOW if "fragmentation" in w or "dA" in w],
      ["NEBNext enzymatic fragmentation", "dA-tailing"])
check("...and the index PCR last, which is where TSO selection happens",
      f.WORKFLOW[-1], "index PCR")
# Ligation geometry: the T-bearing strand becomes contiguous with the insert, so a
# fragment with FakeD at both ends reads FakeD-top ... dA + FakeD bottom. The dA IS the
# first base of revcomp(FakeD top) -- it pairs with FakeD's 3'-T -- so it is not added twice.
check("revcomp(FakeD top) == dA + FakeD bottom", revcomp(f.FAKED_TOP), "A" + f.FAKED_BOT)
JUNK = f.FAKED_TOP + "X" * 20 + revcomp(f.FAKED_TOP)
check("such a fragment's top strand starts with FakeD top", JUNK.startswith(f.FAKED_TOP))
check("...and ends with revcomp(FakeD top)", JUNK.endswith(revcomp(f.FAKED_TOP)))
check("so FakeD top matches the 5' end of the top strand", JUNK.startswith(f.FAKED_TOP))
check("AND the 5' end of the bottom strand -- hence it primes both ends alone",
      revcomp(JUNK).startswith(f.FAKED_TOP))
check("the real library is asymmetric instead: read-1 side from the TSO, read-2 from FakeD",
      f.final_library().top().startswith(il.P5)
      and not f.final_library().top().startswith(f.FAKED_TOP))
check("a FakeD-FakeD fragment carries no TruSeq Read-1 handle, so it can never get P5",
      f.FAKE_TRUSEQ_R1_PART not in JUNK)
check("which is why the index PCR filters it and the final library stays clean",
      il.P5 not in JUNK)

check.section("FakeD: the annealed duplex")
check("top is 24 nt and bottom 23 nt", (len(f.FAKED_TOP), len(f.FAKED_BOT)), (24, 23))
check("top is the 3' 24 nt of TruSeq Read 2", f.FAKED_TOP, il.TRUSEQ_READ2[-24:])
check("bottom is the first 23 nt of the i7 index-read primer",
      f.FAKED_BOT, il.INDEX1_PRIMER[:23])
check("bottom opens with the canonical GATCGGAAGAGC stem",
      f.FAKED_BOT.startswith(il.STEM))
check("the two strands pair over 23 bp", f.FAKED_DUPLEX_LEN, 23)
check("top is revcomp(bottom) plus one base",
      f.FAKED_TOP, revcomp(f.FAKED_BOT) + f.FAKED_OVERHANG)
check("leaving a single 3'-T overhang on the top strand", f.FAKED_OVERHANG, "T")
check("the duplex melts near the ligation temperature", 55 < tm(revcomp(f.FAKED_BOT)) < 65)
check("revcomp(top) is what read 1 reads through into",
      il.TRIM_SEEN_IN_READ1.startswith(revcomp(f.FAKED_TOP)))
check("the same 3'-T geometry as NEBNext's read-1 arm",
      il.NEBNEXT_ARM_READ1.endswith("T") and f.FAKED_TOP.endswith("T"))
check("OPEN (3): it is NOT Y-shaped -- no single-stranded arm, unlike Atrandi's adapter",
      len(f.FAKED_TOP) - f.FAKED_DUPLEX_LEN == 1
      and len(atrandi.LIGADAPT_TOP) - len(atrandi.LIGADAPT_BOT) + 1 == 20)
check("Atrandi's own adapter does carry a 20-nt single-stranded 3' arm",
      len(atrandi.LIGADAPT_ARM), 20)

check.section("what FakeD stands in for on the Atrandi side")
check("the real read-2 arm is 27 nt",
      f.ATRANDI_READ2_ARM_LEN, len(atrandi.ATRANDI_I7_ANNEAL))
check("the barcode block is 4 barcodes of 8 plus 3 linkers of 4 = 44 nt",
      f.ATRANDI_BARCODE_BLOCK_LEN,
      4 * atrandi.BARCODE_LEN + 3 * atrandi.LINKER_LEN)
check("so FakeD replaces 71 nt of cassette", f.ATRANDI_D_SIDE_LEN, 71)
check("...with 23 bp carrying no barcode at all", f.FAKED_DUPLEX_LEN, 23)
check("and FakeD is shorter than what it replaces", f.FAKED_DUPLEX_LEN < f.ATRANDI_D_SIDE_LEN)

check.section("the i7 indexing primer has to graft the arm back")
check("FakeD supplies only the 3' 24 nt of the 34-nt Read-2 arm",
      len(il.TRUSEQ_READ2) - len(f.FAKED_TOP), 10)
check("the missing stretch is GTGACTGGAG", f.I7_GRAFT, "GTGACTGGAG")
check("TruSeq Read 2 is that graft followed by FakeD top",
      f.I7_GRAFT + f.FAKED_TOP, il.TRUSEQ_READ2)
check("the annealing footprint FakeD does offer is long enough to prime",
      tm(f.FAKED_TOP) > 55)
check("the truncated i7 primer (ours and Atrandi's) pairs FakeD over only 17 nt",
      len(atrandi.ATRANDI_I7_ANNEAL) - len(f.I7_GRAFT), 17)

check.section("P5 and TruSeq Read 1 share ACAC -- do not double-count it")
check("P5 ends in ACAC", il.P5.endswith("ACAC"))
check("TruSeq Read 1 begins with ACAC", il.TRUSEQ_READ1.startswith("ACAC"))
check("so the universal primer is 58 nt, not 62",
      len(il.NEBNEXT_UNIVERSAL_PRIMER), len(il.P5) + len(il.TRUSEQ_READ1) - 4)
check("and it equals P5 + Read1 with the overlap collapsed",
      il.NEBNEXT_UNIVERSAL_PRIMER, il.P5 + il.TRUSEQ_READ1[4:])
check("the modelled library collapses it the same way, so ACAC appears once",
      f.final_library().top().count("ACACTCTTTCCC"), 1)

check.section("the step-by-step intermediates the page draws")
fs = f.first_strand()
check("the tagged first strand carries the TSO handle", f.TSO_HANDLE in fs.top())
check("...the UMI", "N" * f.TSO_N_LEN in fs.top())
check("...and the polyA/oligo-dT junction at the far end", "AAAAAAAA" in fs.top())
check("on the first strand itself the tag is only the complement, at its 3' end -- the end "
      "random-primed PTA copies never include",
      revcomp(fs.top()).endswith(revcomp(f.TSO_HANDLE)))
check("the UMI sits at the 5' end, not the 3'", fs.top().index("N") < fs.top().index("AAAAAAAA"))

f5, fi = f.fragment(True), f.fragment(False)
check("only the 5' fragment keeps the Smart-seq3 tag",
      (ss.SS3_TSO_TAG in f5.top(), ss.SS3_TSO_TAG in fi.top()), (True, False))
check("both fragments are dA-tailed at each end",
      (f5.top().startswith("A") and f5.top().endswith("A"),
       fi.top().startswith("A") and fi.top().endswith("A")), (True, True))
check("the dA is drawn unpaired, as an overhang must be",
      f5.get("dA (5' end)").bottom_text(), " ")

l5, li = f.ligated(True), f.ligated(False)
check("ligation puts FakeD on the 5' side of both products",
      (l5.top().startswith(f.FAKED_TOP), li.top().startswith(f.FAKED_TOP)), (True, True))
check("...and on the 3' side of both, reverse-complemented",
      (l5.top().endswith(revcomp(f.FAKED_TOP)),
       li.top().endswith(revcomp(f.FAKED_TOP))), (True, True))
check("so the two products differ ONLY in whether the TSO survived",
      (ss.SS3_TSO_TAG in l5.top(), ss.SS3_TSO_TAG in li.top()), (True, False))
check("the T:A junction base is FakeD's own 3'-T",
      l5.get("T:A junction").top, f.FAKED_OVERHANG)
check("an internal fragment has no site for FakeTruseq-TSO, so it is a dead end",
      f.TSO_HANDLE not in li.top())
check("...and therefore can never acquire P5", il.P5 not in li.top())
check("the 5' product does have the site FakeTruseq-TSO primes on",
      f.TSO_HANDLE in l5.top())
check("every intermediate's segment names are unique, or Construct would have refused",
      all(len({sg.name for sg in c.segments}) == len(c.segments)
          for c in (fs, f5, fi, l5, li)))

check.section("phosphorylation, as ordered: bottom only")
check("top as ordered is the modelled FakeD top", f.FAKED_TOP, "TTCAGACGTGTGCTCTTCCGATCT")
check("bottom as ordered is the modelled FakeD bottom", f.FAKED_BOT, "GATCGGAAGAGCACACGTCTGAA")
check("only the bottom strand is phosphorylated", (f.FAKED_TOP_MOD5, f.FAKED_BOT_MOD5),
      ("", "/5Phos/"))
check("the bottom's phosphorylated 5' end is at the ligation end (its first 12 nt are the "
      "GATCGG stem, recessed opposite the top's 3'-T)",
      f.FAKED_BOT.startswith(il.STEM) and revcomp(f.FAKED_BOT) == f.FAKED_TOP[:-1])
check("the top's 3' end is the T overhang, which seals to the insert's kinased 5'-P",
      f.FAKED_TOP[-1], "T")
check("the distal end is blunt (top 5' base pairs bottom 3' base)",
      f.FAKED_TOP[0], revcomp(f.FAKED_BOT[-1]))
check("...and carries no 5'-phosphate on either strand, so FakeD cannot blunt-dimerise",
      f.FAKED_TOP_MOD5 == "" and f.FAKED_BOT_MOD5 == "/5Phos/")
check("two ligation ends cannot join each other: T faces T", bases_pair("T", "T"), False)
_fd = f.scene_faked()
check("the drawing labels the bottom 5'-p and the top 5'-OH",
      (_fd.strands["bottom"].mod5, _fd.strands["top"].mod5), ("p", "OH"))

check.section("the predicted library")
lib = f.final_library()
check("it starts at P5", lib.top().startswith(il.P5))
check("it ends at revcomp(P7)", lib.top().endswith(il.P7_RC))
check("it carries the Smart-seq3 tag, which is what marks a 5' UMI read",
      ss.SS3_TSO_TAG in lib.top())
check("it carries the FakeD read-2 arm on the far side",
      revcomp(f.FAKED_TOP) in lib.top())
check("read 1 crosses 19 nt of handle, 10 of UMI and 3 of G-tail before cDNA",
      f.READ1_PREFIX_LEN, 32)
check("...so a 50-cycle read 1 leaves only 18 bases of transcript", 50 - f.READ1_PREFIX_LEN, 18)
check("the full 34-nt TruSeq Read 2 arm is restored by the graft",
      revcomp(il.TRUSEQ_READ2) in lib.top())
check("i7 is the 10-nt UDP index of our primers",
      lib.segments[[x.name for x in lib].index("i7")].top, "I" * 10)

check.section("our index primers (atrandi-wgs__10.1101+2025.06.20.660799/ref/our_index_primers.tsv) fit this library")
check("our P5 3' 20 nt anneal to the start of FakeTruseq-TSO",
      f.FAKE_TRUSEQ_TSO.startswith(atrandi.ATRANDI_P5_ANNEAL))
check("our P5 3' end lies inside the 31 nt FakeTruseq-TSO does carry, "
      "so the missing terminal T does not affect it",
      len(atrandi.ATRANDI_P5_ANNEAL) < len(f.FAKE_TRUSEQ_R1_PART))
for k, p in atrandi.OUR_I7.items():
    anneal = p[len(il.P7) + atrandi.OUR_I7_INDEX_LEN:]
    check(f"{k}: 3' 17 nt pair with FakeD top, 5' 10 nt are the graft",
          (anneal[:10], anneal[10:]), (f.I7_GRAFT, f.FAKED_TOP[:17]))

check.section("sequencing primers on the finished library, with our index primers")
t = lib.top()
check("stock TruSeq Read 2 primer has its full 34-nt site (graft + FakeD)",
      revcomp(il.TRUSEQ_READ2) in t)
check("...and read 2 starts on cDNA, right after the dA",
      t[t.index(revcomp(il.TRUSEQ_READ2)) - 1], "X")
check("stock Index 1 primer has its full 33-nt site", il.INDEX1_PRIMER in t)
check("stock TruSeq Read 1 primer is NOT an exact match (3'-T missing)",
      il.TRUSEQ_READ1 not in t and il.TRUSEQ_READ1[:-1] in t)
check("FakeTruseq-TSO is an exact match, so it can be the custom read-1 primer",
      f.CUSTOM_READ1_PRIMER in t)
check("...and read 1 then starts on the UMI",
      t[t.index(f.CUSTOM_READ1_PRIMER) + len(f.CUSTOM_READ1_PRIMER)], "N")
check("...spending 13 cycles before cDNA instead of 32", f.READ1_PREFIX_LEN_CUSTOM, 13)
wgs = atrandi.final_library()
check("the custom primer has no site on the scWGS library (either strand)",
      f.TSO_HANDLE not in wgs.top() + revcomp(wgs.top()))
check("the scWGS library DOES have the stock read-1 site (its last T pairs the insert's "
      "dA, which atrandi.py leaves inside the insert placeholder)",
      revcomp(atrandi.LIGADAPT_TOP), il.TRUSEQ_READ1[:-1])
check("FakeTruseq-TSO's TSO footprint and our P5's footprint have matched Tm (within 1 C)",
      abs(tm(f.TSO_HANDLE) - tm(atrandi.ATRANDI_P5_ANNEAL)) < 1)

check("every declared sequencing primer lands where declared (lib/seqprimers.py)",
      __import__("seqprimers").verify(lib, f.SEQ_PRIMERS), [])

check("the UMI sits inside read 1, not in an index read",
      lib.top().index("N") < lib.top().index("X"))

check.report()
