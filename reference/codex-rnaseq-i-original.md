Codex RNAseq I: PTA-free 5-prime RNA capture with directional adapter ligation

STATUS

Proposed pilot design, not an experimentally validated protocol. Enzyme amounts,
extension times and amplification cycle counts remain to be established for the
input and isolation system. Candidate Q/B handles require thermodynamic review.

DESIGN AND ASSUMPTIONS

Use BciVI to expose a defined single-base 3-prime-A overhang and ligate
the original single-base 3-prime-T-overhang FakeD adapter. FakeD is the
one-ligation prototype substitute for the four Atrandi barcoding rounds.
Internal restriction-site losses are accepted for this prototype. Cells remain
separate through indexing. DNA fragments larger than 400 bp cannot be removed
from their compartments before barcoding; retained long templates and unwanted
amplification products must therefore be tolerated.

The workflow assumes short oligos and cleaved handles can be cleared by exchange
within the isolation system. Their clearance efficiency is not established by
the retention constraint alone. Products below the effective retention threshold
may also be lost.

BciVI recognizes GTATCC and cuts at offsets (6/5), producing a single-base
3-prime overhang. NEB specifies heat inactivation at 80 degrees C for 20 minutes.
Sources: https://www.neb.com/en/products/r0596-bcivi
https://www.neb.com/en/tools-and-resources/usage-guidelines/heat-inactivation

Proceed with single-base-overhang ligation for this amplified pilot. Prior
amplification increases available substrate; it does not establish the fraction
successfully ligated. NEB notes that BciVI ends can be difficult to ligate.
Clear released handles before ligation because they can rejoin the cut product.

Random-primer architecture, written 5-prime to 3-prime:

    Q -- B -- GTATCC -- GACTGT -- NNNNNNNN
               site    spacer    random

The fixed spacer's sixth base, T, determines the end; the cuts do not fall in
the random octamer. Cleavage geometry in the random-primer orientation:

    5-prime -- Q -- B -- GTATCC -- GACTGT | N8 -- target -- 3-prime
    3-prime -- Q'-- B'-- CATAGG -- CTGAC | A N8'-- target'-- 5-prime

Discard the left handle. The retained target-containing fragment has a
single-base 3-prime-A overhang complementary to FakeD top's terminal 3-prime-T.
Both nicks can seal: substrate 3-prime-A-OH to FakeD bottom 5-prime-P, and
FakeD top 3-prime-T-OH to the restriction-generated substrate 5-prime-P.
The target-primer end remains blunt and unphosphorylated.

STEP-BY-STEP PILOT WORKFLOW

1. Reverse transcription
   Retain the existing oligo-dT and UMI-bearing TSO.

2. Clear residual RT oligos and RNA
   Use the isolation system's compatible exchange/RNA-removal procedure while
   retaining first-strand cDNA. The RNA-removal method remains to be specified.

3. TSO-directed linear amplification
   Use the H primer alone to copy TSO-tagged first strands by thermal extension.

4. Clear unused H primer.

5. Random extension
   Use the revised tailed random octamer to copy back toward H.

6. Clear unused tailed random primer.

7. Selective PCR
   Amplify with FakeTruseq-TSO and Q using a blunt-producing polymerase and
   unphosphorylated primers. Long inverted handles are intended to suppress
   handle-handle products; this suppression remains unvalidated.

8. BciVI digestion
   Expose the defined single-base 3-prime-A end. Heat-inactivate at 80 degrees C
   for 20 minutes and clear the short released handle. Do not globally end-repair
   or dA-tail, as that would erase the intended end asymmetry.

9. Adapter ligation
   Ligate the original, pre-annealed T-overhang FakeD duplex and clear free adapter.

10. Cell/sample indexing
    Amplify with P5 and a distinct full-footprint P7-i7 primer for each compartment.

11. Pool after indexing
    Perform conventional purification, size assessment and P5/P7 library
    quantification after pooling becomes permissible.

CANDIDATE OLIGO LIST

All sequences below are written 5-prime to 3-prime. Except where marked, use
unmodified DNA with 5-prime-OH. N means mixed A/C/G/T; W means A/T.

RT oligo-dT
TTTTTTTTTTTTTTTTTT

Existing TSO
/5BiosG/AGAGACAGATTGCGCAATGNNNNNNNNWWrGrG+G
Retain the existing vendor-specific biotin, ribo-G and LNA-G notation.

H linear primer
AGAGACAGATTGCGCAATG

Q PCR primer
GCTACGACGTTACGCTACGA

BciVI tailed random primer
GCTACGACGTTACGCTACGACTGACTGACCGTACGATCGAGTATCCGACTGTNNNNNNNN

FakeTruseq-TSO PCR primer / custom Read 1 sequencing primer
ACACTCTTTCCCTACACGACGCTCTTCCGATCAGAGACAGATTGCGCAATG

FakeD top (terminal 3-prime-T overhang)
TTCAGACGTGTGCTCTTCCGATCT

FakeD bottom
/5Phos/GATCGGAAGAGCACACGTCTGAA

P5 indexing primer
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGAC

Full-footprint P7-i7 primer architecture
CAAGCAGAAGACGGCATACGAGAT + selected 10-base i7 + GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT

The P7-i7 entry is an architecture, not a literal order containing random bases.
Use a defined index per cell/sample and select a suitable index set for the pool.
The four existing index members are insufficient for more than four compartments.
Q and B are candidates for manual thermodynamic checking. B is incorporated
within the tailed random primer and is not ordered separately.

Handle components for checking:
Q = GCTACGACGTTACGCTACGA
B = CTGACTGACCGTACGATCGA
BciVI site = GTATCC
Fixed spacer = GACTGT

SEQUENCING AND COUNTING

The indexing primer completes the full Read-2 primer site. Its terminal T is
the FakeD overhang base, paired with the substrate A at the sealed junction.
Read 2 begins immediately with the eight random-primer-derived bases; trim
these 8 bases before alignment. There is no additional fixed AG prefix.

Custom Read 1 starts at the original RT UMI. Interpret N8 as the RT UMI and W2
as a spacer, followed by the template-switch G tail. Tail-length variation needs
an appropriate parsing rule. Random priming does not create a new counting UMI.
Deduplicate technical descendants within cell/sample and gene using the RT UMI.

Cells remain separate until indexing; this design supplies no additional cell
barcode. Intended directional ligation and suppression must be checked
experimentally rather than inferred from total library yield.

SHARED VIRAL-BOUNDARY BRANCH

The same tailed random primer and handle-cleavage/FakeD interface can be used
after outward-facing viral-primer linear amplification. This RNA oligo list
does not specify that branch's viral-target primers. Products bearing the
engineered cassette are intentionally cleaved at that cassette. Internal BciVI
sites can also cut viral/genomic sequence and lose boundary-containing products;
check known viral target sequences and assess this separately from accepted
transcript losses. Additional amplification does not restore lost boundaries.
