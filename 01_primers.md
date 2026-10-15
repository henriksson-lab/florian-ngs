# florian-PTA-rnaseq — the oligo set

An attempt at single-cell RNA-seq read out through Atrandi's SPC/PTA workflow. The design
is a graft of two protocols already modelled here:

- **front end: Smart-seq3xpress** — oligo-dT primes RT, the RT template-switches onto a
  UMI-carrying TSO, so each mRNA's 5' end is tagged. See
  [SMART-seq3xpress schematic](https://henriksson-lab.github.io/NGS_101/smart-seq3xpress__10.1038+s41587-022-01311-4/smart-seq3xpress.html).
- **back end: Atrandi** — but the **four rounds of split-pool barcoding are replaced by a
  single pre-annealed duplex, `FakeD`**, ligated on in one step. See
  [Atrandi barcode-cassette research note](https://github.com/henriksson-lab/NGS_101/blob/main/atrandi-wgs__10.1101%2B2025.06.20.660799/05_barcode_cassette_model.md).

The molecule therefore takes its **read-1 side from the TSO** (installed by PCR with
`FakeTruseq-TSO`) and its **read-2 side from the ligated `FakeD`**. Only the TSO end is
sequence-selected, which is the point of the design.

Marking: 🟢 as ordered · 🟡 inferred · 🔴 not stated · ✅ computed here.

## 🟢 The workflow

```
reverse transcription (oligo-dT 18T + Smartseq3expressTSO)
  -> PTA amplification (random hexamer + TSO-handle spike-in)
  -> NEBNext enzymatic fragmentation
  -> dA-tailing
  -> ligate the annealed FakeD duplex      (standing in for 4 rounds of barcoding)
  -> index PCR                              (FakeTruseq-TSO selects TSO-carrying molecules)
```

🟢 PTA runs on **first-strand cDNA** with both random hexamer and a TSO-handle
spike-in. The spike-in is essential if the assay is meant to retain transcript 5' ends;
without it, the design should instead be redirected toward 3' capture. Primer concentration
and terminator chemistry remain unspecified.

## Core oligos and the PTA spike-in 🟢

| | sequence | length |
|---|---|---|
| **Smartseq3expressTSO** | `/5BiosG/AGAGACAGATTGCGCAATG` + 8×N + 2×W + `rGrG+G` | ✅ 32 nt |
| **Oligo dT (18T)** | `TTTTTTTTTTTTTTTTTT` | 18 nt |
| **FakeTruseq-TSO** | `ACACTCTTTCCCTACACGACGCTCTTCCGATCAGAGACAGATTGCGCAATG` | 51 nt |
| **FakeD top** | `TTCAGACGTGTGCTCTTCCGATCT` | 24 nt |
| **FakeD bottom** | `GATCGGAAGAGCACACGTCTGAA` | 23 nt |
| **PTA TSO-handle spike-in** | `AGAGACAGATTGCGCAATG` with two protected 3'-terminal bonds | 19 nt |

## ✅ What each one is, checked against the canonical blocks

Every identity below was computed against `lib/illumina.py` and `smart-seq__10.1038+nbt.2282/tools/smartseq.py`,
not read off by eye:

| oligo | resolves to | exact? |
|---|---|---|
| TSO 5' handle `AGAGACAGATTGCGCAATG` | `smartseq.SS3_TSO_HANDLE` = ME 3' end (8 nt) + Smart-seq3 tag (11 nt) | ✅ identical |
| TSO 3' tail `rGrG+G` | `rt.TSO_G_TAIL_LNA` — two riboG then one LNA-G | ✅ identical |
| FakeTruseq-TSO, 5' 32 nt | `illumina.TRUSEQ_READ1` **minus its terminal T** | ⚠ one base short |
| FakeTruseq-TSO, 3' 19 nt | the TSO handle | ✅ identical |
| FakeD top | `illumina.TRUSEQ_READ2[-24:]` — a 3' piece of TruSeq Read 2 | ✅ identical |
| FakeD bottom | `illumina.INDEX1_PRIMER[:23]` — begins with the `GATCGGAAGAGC` stem | ✅ identical |
| `revcomp(FakeD top)` | the start of `illumina.TRIM_SEEN_IN_READ1` | ✅ identical |

## ✅ FakeD, annealed

```
top     5'- TTCAGACGTGTGCTCTTCCGATC T -3'
            ||||||||||||||||||||||| ^ 1-nt 3'-T overhang
bottom  3'- AAGTCTGCACACGAGAAGGCTAG   -5'
```

✅ 23 bp paired, Tm ≈ 59.7 °C, with a single **3'-T on the top strand** — i.e. an ordinary
dA-tail-compatible adapter, the same geometry as NEBNext's.

### What it replaces

Atrandi's D side presents, from the P7 end: the 27-nt TruSeq Read-2 arm, then
`BC-D · linker · BC-C · linker · BC-B · linker · BC-A` — ✅ 4 × 8 + 3 × 4 = **44 nt** of
cassette, **71 nt** in total. `FakeD` collapses all of that into 23 bp carrying no barcode.

⚠ It is **not Y-shaped.** Atrandi's real ligation adapter has a 20-nt single-stranded 3'
arm (`LIGADAPT_ARM`) that presents the read-1/P5 side; `FakeD` has no arm at all. That is
deliberate here — the read-1 side arrives later, by PCR off the TSO — but it has a
consequence, below.

## ✅ Re-analysis with our index primers known (2026-10-06)

Before this, the read-2 and index sides were modelled on Atrandi's primers. With the real
oligos ([source location](https://github.com/henriksson-lab/NGS_101/tree/main/atrandi-wgs__10.1101%2B2025.06.20.660799/ref))
every sequencing primer can now be checked against the finished molecule. All ✅ items are
pinned in `tools/selftest.py`.

| Read | Primer | Site on this library | Verdict |
|---|---|---|---|
| Read 2 | stock TruSeq Read 2 (34 nt) | `GTGACTGGAG` (i7 graft) + `FakeD top` | ✅ **exact.** Starts on cDNA, after the dA |
| Index 1 | stock i7 index primer (33 nt) | `FakeD bottom` + graft | ✅ **exact.** Reports the 10-nt UDP i7 |
| Index 2 | stock | P5 / Read-1 stub | ✅ primes. Its 5' end faces the missing T, but that end doesn't matter. Reads constant `ACCGAGATCT` |
| Read 1 | stock TruSeq Read 1 (33 nt) | 32 nt offered | ⚠ **unchanged: 3'-T:T mismatch** (OPEN 1) |
| Read 1 | `FakeTruseq-TSO` as custom primer | exact, 51 nt | ✅ exact; read 1 starts **on the UMI** |

**Read 2 / Index 1 — resolved.** The doubt was whether anything restores the 10 nt `FakeD`
leaves off the Read-2 arm. Our i7 primer does (its 5' tail `GTGACTGGAG`), so both stock
P7-side primers find complete sites. The one remaining cost is the 17-nt first-cycle
footprint (Tm ≈ 51 °C), covered in the predicted-library section below.

**Read 1 — not fixed by the new primers, and the new info makes it more pressing.** Our P5
primer's footprint ends at base 20 of the Read-1 site; the missing T is base 33, which only
`FakeTruseq-TSO` supplies. Three new considerations:

1. **The two protocols share index primers, so they will likely share runs.** The scWGS
   libraries need the stock Read-1 primer (their site is complete; the T pairs the insert's
   dA). 🟢 On NovaSeq X a custom Read-1 primer in well CP1 **replaces** the Illumina primer
   for the whole run, and PhiX is then not sequenced. So pooling with a custom primer means
   mixing our own TruSeq Read-1 oligo into CP1 alongside it.
2. ✅ **Read 1 is 100 % low-diversity at its start.** TSO selection is the whole point of the
   design, so *every* cluster begins with the same `AGAGACAG ATTGCGCAATG` (19 constant
   cycles with the stock primer). Smart-seq3 gets away with its 11-nt tag because most of
   its reads are internal Tn5 fragments; here none are. 🟡 Sequenced alone, this is the
   classic template-registration failure. It needs heavy dilution in the lane: scWGS
   libraries (random genomic Read-1 starts) or PhiX.
3. ✅ **`FakeTruseq-TSO` is already the ideal custom Read-1 primer.** It matches the library
   exactly and ends on the last base of the tag, so read 1 starts on the UMI: **13 cycles
   before cDNA instead of 32**, and the first 10 cycles are random (UMI). That fixes both
   the cycle cost and the diversity problem. 🟡 The missing T then becomes an *asset*: the
   custom primer has no site on scWGS libraries, and the stock primer faces a 3'-terminal
   mismatch on RNA libraries. The two are close to orthogonal, so they can share CP1.

**So the two fixes for OPEN (1) are mutually exclusive — pick one:**

| | A: add the T (52-nt oligo) | B: keep the 51-nt oligo, sequence with it |
|---|---|---|
| Read-1 primer | stock only | stock TruSeq R1 + `FakeTruseq-TSO` mixed in CP1 |
| Cycles before cDNA | 32 | **13** |
| Read-1 start diversity | constant for 19 cycles | random (UMI) |
| Pool with scWGS | trivially | ✅ yes; primers are near-orthogonal |
| Risk | low-diversity start; needs heavy PhiX/scWGS spike | 🟡 stock primer leak-extending across the T:T mismatch on RNA clusters → mixed-phase signal, lower Q. Order the primer HPLC-purified |

🟡 **Recommendation: B**, if the facility will load custom primers. Do **not** do both:
with the T added, the stock primer also matches the RNA libraries exactly and competes with
the custom primer on the same cluster, 19 bases out of register.

**The TSO side, re-checked.**
- ✅ The `FakeTruseq-TSO` footprint on the TSO handle (19 nt, Tm ≈ 53.6 °C) matches our
  P5's footprint on its stub (20 nt, Tm ≈ 54.3 °C). If both run in one PCR at Atrandi's
  54 °C, neither limits the other.
- ✅ Our P5 primes only on `FakeTruseq-TSO` product, so TSO selection is preserved through
  indexing. TSO-concatemer and strand-invasion artefacts carry the same handle and are
  selected too; nothing in the index primers changes that.
- ⚠ **Quantify by P5/P7 qPCR, not Qubit.** Our i7 primer's 3' 17 nt are the 5' 17 nt of
  `FakeD top`, so i7 *alone* amplifies FakeD–FakeD fragments exponentially (OPEN 3). These
  are P7–P7 molecules that never cluster but do count on Qubit or TapeStation, so a mass-based
  loading concentration will underload the flow cell. A P5+P7 qPCR (KAPA-style) counts only
  clusterable molecules.
- ✅ **Model bug fixed along the way.** `final_library()` and `ligated()` drew the insert's
  dA *and* `revcomp(FakeD top)`, which already begins with that A, so the junction had one A
  too many. Now dA + `FakeD bottom`.

**Index reads, briefly.** 🟡 The 10-cycle i5 read is a constant `ACCGAGATCT`, the same for
every sample. On a 2-channel instrument its two G cycles are all-dark if the lane holds only
these libraries. Run i5 at 0 cycles unless other libraries in the lane need it. The four i7s
are colour-balanced as a set of four (every position has a non-G), not necessarily in
subsets.

## 🔴 Open questions, in order of how much they matter

### ⚠ (1) `FakeTruseq-TSO` is one base short of TruSeq Read 1

> **Update 2026-10-06:** see the re-analysis above. There are two mutually exclusive fixes;
> keeping the 51-nt oligo and using it as the custom read-1 primer (B) is now recommended.

✅ Checked the way it should be checked: **against how Smart-seq3 handles the same
junction**, since that is the design being copied.

**Smart-seq3's junction is exact.** It sequences with the stock Nextera read-1 primer,
`nextera.READ1_PRIMER = S5 + ME`, and ✅ the TSO's first 8 nt **are** the last 8 nt of ME:

```
Nextera READ1_PRIMER  TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG          33 nt
SS3 library           TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAATG NNNNNNNN GGG cDNA
                      |<------- primer, an exact prefix ------->|^ read 1 starts here
```

So the TSO **completes** the primer site, `AGAGACAG` sits *inside* the primer and is never
read, and read 1 begins on the first base of the tag. That is not luck — it is why the
Smart-seq3 TSO starts with `AGAGACAG` at all.

**The TruSeq version does not line up.** TruSeq Read 1 ends in `...CCGATCT`; the library
offers `...CCGATC`:

```
stock TruSeq Read 1   ACACTCTTTCCCTACACGACGCTCTTCCGATCT           33 nt
this library          ACACTCTTTCCCTACACGACGCTCTTCCGATC AGAGACAGATTGCGCAATG ...  32 nt offered
                                                     ^ primer's base 33 is T, library's is A
```

✅ Computed: the stock primer is **not** a prefix of this library; its 3'-terminal base
faces an A, giving a **T:T mismatch at the 3' terminus** — the one position where a
mismatch actually stops extension. If it extends regardless, read 1 begins one base late
and ✅ the `ATTGCGCAATG` tag lands at **offset 7, not 8**.

**Why the T is there in the first place, and why dropping it is understandable.** 🟢 In
TruSeq that terminal T *is* the dA-ligation overhang base — `illumina.NEBNEXT_ARM_READ1` is
commented "carries the 3'-T". Florian's read-1 handle is installed **by PCR, not by
ligation**, so there is no dA on that side and no natural T. Dropping it is a reasonable
inference, not a careless one. But the stock sequencing primer does not know that.

**The fix is one base.** ✅ `TRUSEQ_READ1 + TSO_HANDLE` = **52 nt**, and the stock primer is
then an exact prefix of the library:

```
ACACTCTTTCCCTACACGACGCTCTTCCGATCTAGAGACAGATTGCGCAATG     52 nt
```

> If the 32-nt form is kept deliberately, then **a custom 32-nt read-1 primer must be
> ordered with it** and the offset recorded, because nothing else will prime cleanly.

### ✅ A second cost of the Nextera → TruSeq swap: 10 extra read-1 cycles

This one is independent of the missing T, and is the more expensive of the two.

| | read 1 spends | on |
|---|---|---|
| **Smart-seq3** | **22 cycles** | tag 11 + UMI 8 + `GGG` 3 |
| **this design** | **32 cycles** | `AGAGACAG` 8 + tag 11 + UMI 10 + `GGG` 3 |

✅ The 10-cycle difference is `8 + 2`: **8** because TruSeq has no ME, so the `AGAGACAG`
that Smart-seq3 hides inside its sequencing primer is now read as library sequence; and
**2** from the longer UMI. At 50 cycles that is 18 bases of transcript instead of 28.

⚠ Worth knowing before committing to a read length — and it is not fixable by adding the T,
because the `AGAGACAG` has to be in the molecule either way. The only way to recover those
8 cycles is to sequence with a custom read-1 primer that ends `...CCGATCTAGAGACAG`.

### ✅ (2) Resolved — only `FakeD bottom` is phosphorylated, which is the right design

🟢 As ordered (2026-10-06):

```
afake_adD_top_v2      5'-        TTCAGACGTGTGCTCTTCCGATCT   -3'   (5'-OH)
afake_adD_bottom_v2   5'-/5Phos/GATCGGAAGAGCACACGTCTGAA     -3'
```

Same sequences as modelled. ✅ Annealed (generated by `chemdraw.Scene`, which refuses unpaired columns):

```
afake_adD_bottom_v2 5'-p GATCGGAAGAGCACACGTCTGAA -3'
afake_adD_top_v2    3'- TCTAGCCTTCTCGTGTGCAGACTT OH-5'
                        ^ ligation end: 3'-T overhang
```

✅ **Both junctions with the insert are sealed.** End-prep leaves the insert with a 5'-P and a
3'-dA on each end:

| Nick | Ends meeting | Sealed? |
|---|---|---|
| bottom strand ↔ insert | FakeD bottom **5'-P** + insert dA **3'-OH** | ✅ |
| top strand ↔ insert | FakeD top **3'-T-OH** + insert **5'-P** (kinased in end-prep) | ✅ |

The top's own 5' end sits at the *distal* end, where it never meets the insert, so leaving it
unphosphorylated costs nothing.

✅ **…and FakeD cannot dimerise.** This was the concern when both strands were assumed
phosphorylated:

- **ligation end to ligation end:** two 3'-T overhangs, T faces T → no.
- **distal end to distal end:** blunt, but the only 5' end there is the top's **5'-OH**, so
  no strand can be sealed → no. This is the job Atrandi's `/5AmMC6/` does on their adapter.

So the predicted ~46 bp dimer peak should not appear. 🟡 Two residual routes, both minor:
- **Carry-over kinase.** Atrandi's inference is that kinase activity from end-prep survives
  into the ligation, which is why they use an amino block rather than a bare OH. NEBNext's
  65 °C end-prep step should inactivate it. If a ~46 bp peak does show up, this is why, and
  the fix is a 5'-blocked top strand.
- **Blunt fragments that missed dA-tailing.** These can take FakeD's distal end by blunt
  ligation, but only one strand gets sealed (bottom 3'-OH to the insert's 5'-P). Such
  molecules carry FakeD backwards and will not amplify with the i7 primer as designed.

### ⚠ (3) A blunt duplex puts the *same* handle on both ends

🟢 Confirmed workflow order: **RT → PTA → NEBNext enzymatic fragmentation → dA-tailing →
ligate the annealed `FakeD` → index PCR.**

So `FakeD` goes onto *fragments*, and onto **both ends of every one of them**. ✅ Worked
through the ligation geometry, a FakeD–FakeD fragment's top strand is

```
TTCAGACGTGTGCTCTTCCGATCT [insert] AGATCGGAAGAGCACACGTCTGAA
|<--- FakeD top ------->|         |<- revcomp(FakeD top) ->|
                                  ^ the insert's dA, paired with FakeD's 3'-T
```

— and `FakeD top` matches the 5' end of **both** strands, so ✅ it primes both ends and
amplifies such fragments **exponentially on its own**. Fragmentation makes this worse, not
better: one cDNA molecule yields many fragments and only one of them carries the TSO.

🟢 **The index PCR is what filters them**, and that does work: a FakeD–FakeD fragment never
acquires a TruSeq Read-1 handle, so it cannot get P5, and a fragment with P7 at both ends
does not cluster. **The final library is therefore clean** — this is not a correctness bug.

⚠ What it costs is upstream of that:

- the `FakeTruseq-TSO` + `FakeD top` reaction is **dominated by FakeD–FakeD product**, so
  Qubit / TapeStation / qPCR on it measures mostly material that will be discarded;
- a cycle number or an input mass for the index PCR chosen from those numbers is therefore
  **wrong by whatever that ratio is**, and the ratio is not known;
- the junk competes for polymerase and dNTPs, so the molecules that matter amplify less.

✅ **Testable in one extra reaction:** run the PCR with `FakeD top` alone and compare yield
to the two-primer reaction. Similar yield means this dominates. 🟡 If it does, the fix is
the standard one — give `FakeD` a Y shape (NEBNext's own adapter is a hairpin cut by USER
for exactly this reason), or block `FakeD top` so it cannot prime.

### 🔴 (4) The UMI is 10 nt, not Smart-seq3's 8

✅ `(N:25252525)(N)×7` = 8 N, then `(W:5050)(W)` = 2 W, so **8 N + 2 W = 10 nt**.
Smart-seq3 itself stops at 8 N (`smartseq.SS3_UMI_LEN = 8`).

🟡 The two **W** (A or T) positions immediately before `rGrG+G` look deliberate: they stop
the UMI ending in G, which would otherwise be indistinguishable from the template-switch
`GGG` and make the UMI boundary ambiguous. Worth confirming that the intended UMI length
really is 10, because every downstream offset depends on it.

### ⚠ (5) The TSO-handle spike-in preserves the 5' end during PTA

🟢 Confirmed protocol decision: PTA runs on **first-strand** cDNA with a primer mixture
containing **random hexamer plus the TSO-handle spike-in**. The spike-in is
`AGAGACAGATTGCGCAATG`, protected at its two 3'-terminal bonds. Its concentration and the
terminator chemistry remain unspecified.

**Why the strand matters for the TSO.** ✅ The first strand carries the TSO only as its
complement, at its **3' end** (pinned in the selftest):

```
first strand   5'- TTTT…T [antisense cDNA] CCC [UMI'] [handle'] -3'
                                                    ^^^^^^^^^^^^^^^ the only copy of the tag
```

A polymerase copies **from a primer toward the template's 5' end**. A random primer landing
at distance *d* from the 3' end copies everything 5' of itself, but never the *d* bases 3'
of it. Every random-primed copy therefore *starts* at its primer and lacks the template's
3' end. ✅ To carry the full 19-nt handle, a copy would have to be primed within the last
few bases of the template, by a hexamer that happens to spell `AGAGAC`. In practice no
first-generation copy carries the tag. Later generations don't help either: they copy the
copies, whose 5' ends are random primer sites, not the TSO.

This is the general "ends are lost" property of random-primed amplification. It hits the
TSO end precisely because that end is the **template's 3' end**. It is also
strand-specific, which is why the strand question mattered:

| Template | Where the TSO sits | Random-primed copies that include it |
|---|---|---|
| **first strand** (this protocol) | 3' end, as the complement | ≈ none: every copy starts downstream of it |
| second strand (sense) | 5' end, as the handle | **all** primed within one amplicon length; they run *into* it |

🟡 **Consequences if the spike-in is omitted:**

1. The tag exists in roughly **one copy per transcript**: the original first strand. PTA
   multiplies the transcript body, not the tagged end.
2. Worse, that original 3' end stays **single-stranded**: the 32-nt TSO complement lies
   3' of the last priming site, so nothing pairs with it. NEBNext Ultra II FS end-prep
   polishes ends, removing 3' overhangs, so it will likely **trim the tag off** before
   `FakeD` ligation. Even surviving molecules need a ds end to ligate.
3. `FakeTruseq-TSO` selection in the index PCR then has almost nothing to select. Expect very
   low yield of tagged molecules, and a library dominated by whatever leaks through
   OPEN (3).

🟢 **The protocol now supplies the missing end deliberately.** The PTA reaction contains
both random hexamer and a TSO-handle primer. The latter anneals to the terminal handle and
makes a sense copy that starts with the original molecule's tag and UMI; random primers can
then cover that copy efficiently. This is what converts otherwise broad, end-losing PTA
coverage into a 5'-retaining assay.

Which TSO-derived primer is used still matters for whether the UMI stays honest:

| Primer added to PTA | 3' end lands on | Risk |
|---|---|---|
| **handle only**, `AGAGACAGATTGCGCAATG`, 3'-PS | the last base before the UMI | ✅ none: every copy reads the template's own UMI |
| the full TSO (N8W2 + `rGrG+G`) | the UMI and the G-tail | ⚠ its random UMI mismatches the template's. phi29 must trim 13 nt before extending; any copy that extends without trimming carries the **primer's** UMI → UMI inflation, more "molecules" than there were |
| `FakeTruseq-TSO` | same as handle-only | ✅ works, but the 32-nt 5' tail is dead weight in PTA |

The chosen spike-in is therefore **handle-only, with two 3'-terminal phosphorothioates**
so phi29's exo doesn't eat it (same protection as the PTA random primer), no rG/LNA, no
biotin. Its Tm (≈ 53.6 °C) is far above PTA's 30 °C, so it out-competes hexamers for
that site. 🟡 Readouts worth taking from
the experiment: tagged:internal qPCR ratio (below), and in the data the number of distinct
UMIs per gene vs a no-primer control. A jump far beyond expected molecule counts suggests
UMI-bearing primers are being extended.

**The design fork is now explicit:**

- **For 5' capture:** keep the TSO-handle spike-in in the PTA primer mixture. It copies
  exactly the terminal region that random hexamers miss.
- **Without the spike-in:** do not expect faithful 5' retention. Redesign the assay to
  capture the transcript's **3' end**, which also requires giving the oligo-dT side an
  addressable handle rather than leaving it as the current bare 18-mer.
- **Verify:** qPCR before and after PTA with `FakeTruseq-TSO` + a gene-specific
  reverse primer near a transcript 5' end, normalised to an internal amplicon of the same
  gene. If the tagged:internal ratio collapses after PTA, this is the cause.

### 🔴 (6) The oligo-dT carries no handle

🟢 Smart-seq3's oligo-dT is `ACGAGCATCAGCAGCATACGA` + T30 + `VN`. Here it is a bare 18-mer.
So the 3' end of the cDNA gets **no handle**, and nothing marks where polyA priming
happened. 🟡 Consistent with the design — the 3' side is supposed to come from the ligated
`FakeD`, not from the RT primer — but it also means no `VN` anchor, so the primer can slip
within the polyA tract.

## ✅ The predicted library

🟢 The index PCR uses **our own primers**, shared with the scWGS protocol —
[source location](https://github.com/henriksson-lab/NGS_101/tree/main/atrandi-wgs__10.1101%2B2025.06.20.660799/ref), documented
in [the Atrandi protocol research note](https://github.com/henriksson-lab/NGS_101/blob/main/atrandi-wgs__10.1101%2B2025.06.20.660799/03_our_protocol.md) §6. Atrandi's i7/P5
design with a **10-nt IDT UDP i7** (UDP0005–0008) in place of Atrandi's 6-nt `AACCTG`;
P5 identical to Atrandi's, no i5 index.

```
P5 ── TruSeq Read 1 ── AGAGACAG ── ATTGCGCAATG ── UMI(8N+2W) ── GGG ── cDNA ── dA ── FakeD ── GTGACTGGAG ── i7(10) ── P7
```

✅ P5 primer: its 3' 20 nt `ACACTCTTTCCCTACACGAC` are the first 20 nt of `FakeTruseq-TSO`, so it
primes on the TSO-selected product; the missing terminal T of OPEN (1) lies outside its footprint.

✅ Read 1 spends `19 + 10 + 3 = `**32 cycles** on handle, UMI and the G-tail before it
reaches any cDNA. Budget for that: at 50 cycles only 18 bases of transcript are read.

✅ The i7 indexing primer must **graft back** the 10 nt `GTGACTGGAG` that `FakeD` leaves off
the front of the Read-2 arm. ✅ With our primer, which is truncated like Atrandi's
(`…GTGCTCTT`, no `CCGATCT`), only its **3' 17 nt** pair with `FakeD top` and the 5' 10 nt
are the non-templated graft.

⚠ 🟡 That 17-nt footprint has Tm ≈ 51 °C (`chemdraw.tm`) — below Atrandi's 54 °C anneal —
so the first cycles prime inefficiently on this side. Once a strand carries the full 27-nt
site it is no longer limiting, so expect a few extra cycles rather than failure; a 50–52 °C
anneal for the first 2–3 cycles would avoid it. The graft itself is the same trick the v2
adaptor uses recorded in the [pooled-screen research notes](https://github.com/henriksson-lab/NGS_101/blob/main/lenticrispr-v1-screening__10.1126%2Bscience.1247005/90_combined_screen_research.md).
