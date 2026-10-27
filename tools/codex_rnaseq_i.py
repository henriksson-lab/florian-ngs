"""Proposed PTA-free pilot. Molecular geometry is predicted, never validated chemistry.

Ordered legacy oligos stay in florian_pta_rnaseq; canonical sequences come from
NGS_101. Candidate handles are local. N8 random priming is not a counting UMI.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(HERE), str(ROOT / 'NGS_101' / 'lib'),
               str(ROOT / 'NGS_101' / 'smart-seq__10.1038+nbt.2282' / 'tools')]

import florian_pta_rnaseq as legacy
import illumina as il
import smartseq as ss
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      complement_segments, feature, revcomp)
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, InputRole as R, MaterialKind as M,
                            ProtocolGraph, SequencingPlan, StepKind as K, Terminal)

TITLE = 'Codex RNAseq I: PTA-free 5-prime RNA capture with directional adapter ligation'
AUTHORITY = 'User-supplied proposed pilot; reference/codex-rnaseq-i-original.md (2026-10-27)'
STATUS = 'Proposed pilot design; not experimentally validated'
Q = 'GCTACGACGTTACGCTACGA'
B = 'CTGACTGACCGTACGATCGA'
SITE = 'GTATCC'
SPACER = 'GACTGT'
CUT_TOP, CUT_BOTTOM = 6, 5
H = ss.SS3_TSO_HANDLE
FAKE_TRUSEQ = legacy.FAKE_TRUSEQ_TSO
FAKED_TOP, FAKED_BOTTOM = il.TRUSEQ_READ2[-24:], il.INDEX1_PRIMER[:23]
P5 = legacy.atrandi.OUR_P5
RANDOM_PRIMER = Q + B + SITE + SPACER + 'N' * 8
TSO_ORDER = legacy.TSO_5_MOD + H + 'N' * 8 + 'W' * 2 + legacy.TSO_G_TAIL
ROLES = ('Read 1', 'Index 1 (i7)', 'Read 2')
PRIMERS = (sp.custom('Read 1', 'FakeTruseq-TSO', FAKE_TRUSEQ, AUTHORITY,
                     'Custom Read 1 starts at the original RT UMI.'),
           sp.TRUSEQ['I1'], sp.TRUSEQ['R2'])
UMI = feature('rt_umi', 'umi', 'random', note='Original RT N8; excludes W2')
INDEX = feature('compartment_i7', 'sample_index', 'unknown',
                note='One defined 10-base i7 per compartment; index set not selected')


def seg(name, bases, tag=None, **kw):
    """Every newly predicted molecular segment is visibly inferred by default."""
    return Segment(name, bases, tag, inferred=True, **kw)


def body(*, random_derived=True):
    return [seg('RT UMI N8', 'N' * 8, 'umi', placeholder=True, length_bp=8, feature=UMI),
            seg('W2 spacer', 'WW', 'tso', placeholder=True, length_bp=2),
            seg('template-switch G tail', 'GGG', 'r3', length_bp=3,
                note='Nominal GGG; actual tail length may vary'),
            seg('cDNA', '[cDNA]', 'w1', placeholder=True),
            seg('random-primer-derived N8' if random_derived else 'octamer binding site', 'N' * 8, 'r2', placeholder=True, length_bp=8,
                note='Reverse complement on sense strand; trim eight R2 bases, not a UMI')]


def sense():
    return Construct([seg('H handle', H, 'tso'), *body(random_derived=False),
                      seg('polyA / oligo-dT', 'A' * 18)], name='H-directed sense copy')


def first_strand():
    return Construct(complement_segments(sense().segments), name='tagged first strand')


def rt_binding():
    sc = Scene()
    sc.strand('mRNA', [seg('transcript', 'X' * 30, 'w1', placeholder=True),
                       seg('polyA', 'A' * 18)])
    sc.anneal('nascent cDNA', [seg('oligo-dT', legacy.OLIGO_DT),
              seg('copied transcript', 'X' * 30, 'w1', placeholder=True),
              seg('RT-added CCC', 'CCC', 'r3')],
              to='mRNA', pair=('oligo-dT', 'polyA'), mod5='OH')
    sc.anneal('TSO', [seg('H', H, 'tso'), *body()[:3]],
              to='nascent cDNA', pair=('template-switch G tail', 'RT-added CCC'),
              mod5=legacy.TSO_5_MOD)
    sc.arrow('nascent cDNA', 'RT switches onto modified TSO')
    sc.note('TSO', 'Terminal GGG represents rGrG+G; N8 UMI followed by W2 spacer')
    return sc


def random_product():
    return Construct([seg('H handle', H, 'tso'), *body(),
                      seg('spacer complement', revcomp(SPACER), 'r2'),
                      seg('site complement', revcomp(SITE), 'r2'),
                      seg('B complement', revcomp(B), 'r2'),
                      seg('Q complement', revcomp(Q), 'r2')], name='random-extension duplex')


def h_binding():
    sc = Scene()
    sc.strand('first strand', first_strand().segments, mod5='OH')
    sc.anneal('H primer', [seg('H', H, 'tso')], to='first strand',
              pair=('H', "H handle'"), mod5='OH')
    sc.arrow('H primer', 'thermal linear extension')
    return sc


def random_binding():
    sc = Scene()
    sc.strand('sense copy', sense().segments, mod5='OH')
    tails = [seg('Q', Q, 'r2'), seg('B', B, 'r2'),
             seg('BciVI site', SITE, 'r2'), seg('spacer', SPACER, 'r2')]
    sc.anneal('tailed random primer',
              [*tails, seg('N8', 'N' * 8, 'r2', placeholder=True, length_bp=8)],
              to='sense copy', pair=('N8', 'octamer binding site'),
              unpaired=[s.name for s in tails], mod5='OH')
    sc.arrow('tailed random primer', 'extend back toward H')
    sc.note('tailed random primer', 'N8 is a representative complementary pool member; priming position varies')
    return sc


def selective_product():
    return Construct([seg('FakeTruseq-TSO', FAKE_TRUSEQ, 's5'),
                      *random_product().segments[1:]], name='blunt selective amplicon')


def cleavage():
    """Compute both strands' retained boundary from the type IIS cut offsets.

    Random-primer orientation: cut top after spacer[6], bottom after spacer[5].
    The sense-oriented target retains complement(spacer[5]) at its 3-prime end.
    """
    if len(SPACER) != CUT_TOP or CUT_TOP - CUT_BOTTOM != 1:
        raise ValueError('BciVI cuts must expose one spacer base, outside random N8')
    overhang = revcomp(SPACER[CUT_BOTTOM:CUT_TOP])
    if overhang != 'A' or revcomp(FAKED_TOP) != overhang + FAKED_BOTTOM:
        raise ValueError('Spacer-derived overhang and original FakeD are incompatible')
    return Construct([seg('FakeTruseq-TSO', FAKE_TRUSEQ, 's5'), *body(),
                      seg('substrate A', overhang)], name='BciVI target fragment')


def cleavage_map():
    """Show both cuts in random-primer orientation, without hand-placed columns."""
    cleavage()  # Enforce the spacer/offset/adapter compatibility contract first.
    con = Construct([seg('Q', Q, 'r2'), seg('B', B, 'r2'), seg('site', SITE, 'r2'),
                     seg('spacer bases 1–5', SPACER[:CUT_BOTTOM]),
                     seg('spacer base 6', SPACER[CUT_BOTTOM:CUT_TOP]),
                     seg('random N8', 'N' * 8, 'r2', placeholder=True, length_bp=8),
                     seg('target', '[target]', 'w1', placeholder=True)], name='BciVI cassette')
    sc = duplex(con)
    sc.mark('sense', 'spacer base 6', 'top cut after this base (offset 6)')
    sc.mark('antisense', 'spacer bases 1–5' + "'", 'bottom cut after these five bases (offset 5)')
    sc.note('sense', 'Discard left Q/B/site handle; retain right target with bottom-strand 3′-A')
    return sc


def duplex(con, *, junction=False):
    sc = Scene()
    sc.strand('sense', con.segments, mod5='OH')
    sc.anneal('antisense', complement_segments(con.segments), to='sense',
              pair=(con.segments[0].name + "'", con.segments[0].name), mod5='OH',
              unpaired=[s.name + "'" for s in con if s.is_role_token()])
    sc.labels('sense')
    if junction:
        sc.junction('sense', 'substrate A', 'FakeD bottom')
        # Bottom is drawn reversed; these are the original substrate 5'-P and adapter T.
        sc.junction('antisense', "random-primer-derived N8'", "substrate A'")
    return sc


def cut_scene():
    con = cleavage()
    sc = Scene()
    sc.strand('sense', con.segments, mod5='OH')
    sc.anneal('antisense', complement_segments(con.segments[:-1]), to='sense',
              pair=("FakeTruseq-TSO'", 'FakeTruseq-TSO'), mod5='p',
              unpaired=[s.name + "'" for s in con if s.is_role_token()])
    sc.overhang('sense', 'substrate A', polarity='3-prime', owner='substrate',
                text="restriction-derived 3'-A-OH; opposing recessed 5'-P")
    sc.note('sense', "TSO end stays blunt, 5'-OH; no global end repair or dA-tailing")
    return sc


def adapter_geometry():
    sc = Scene()
    sc.strand('top', [seg('FakeD paired core', FAKED_TOP[:-1], 't7'),
                      seg('FakeD T', FAKED_TOP[-1], 't7')], mod5='OH')
    sc.anneal('bottom', [seg('FakeD bottom', FAKED_BOTTOM, 't7')], to='top',
              pair=('FakeD bottom', 'FakeD paired core'), mod5='p')
    sc.overhang('top', 'FakeD T', polarity='3-prime', text="adapter 3'-T-OH")
    return AdapterLigationGeometry.known('cohesive', cut_scene(), sc,
        statement='Predicted A:T pairing from supplied oligos and NEB BciVI (6/5); '
                  'both ligating 5-prime ends phosphorylated, distal blunt ends 5-prime-OH.')


def ligated():
    return Construct([*cleavage().segments, seg('FakeD bottom', FAKED_BOTTOM, 't7')],
                     name='directionally ligated target')


def p7_primer(index):
    """A literal order requires a selected index, never N10 or another placeholder."""
    if len(index) != 10 or set(index) - set('ACGT'):
        raise ValueError('Choose a defined ten-base A/C/G/T i7 index')
    return il.P7 + index + il.TRUSEQ_READ2


def final_library():
    # P5 and the R1 stub share ACAC; the sealed A is already part of the R2 site.
    return Construct([seg('P5', il.P5, 'p5'),
                      seg('R1 stub', il.TRUSEQ_READ1[4:-1], 's5'),
                      seg('H handle', H, 'tso'), *body(),
                      seg('substrate A', cleavage().segments[-1].top),
                      seg('FakeD bottom', FAKED_BOTTOM, 't7'),
                      seg('Read 2 graft', revcomp(legacy.I7_GRAFT), 't7'),
                      seg('i7', 'I' * 10, 's7', placeholder=True, length_bp=10, feature=INDEX),
                      seg('P7', il.P7_RC, 'p7')], name='predicted indexed RNA library')


def pcr_geometries():
    selective = PCRGeometry.exact(random_product(),
        PCRPrimer('FakeTruseq-TSO', FAKE_TRUSEQ, len(H)), PCRPrimer('Q', Q, len(Q)),
        statement='Candidate selective PCR: H and Q select the captured interval; '
                  'blunt polymerase and unphosphorylated primers required.')
    # N10 is used only for modelling the unpaired index tail, never for ordering.
    indexing = PCRGeometry.exact(ligated(), PCRPrimer('P5', P5, 20),
        PCRPrimer('P7 + selected i7 + full Read 2', il.P7 + 'N' * 10 + il.TRUSEQ_READ2, 24),
        statement='Architecture only: select one defined i7 per compartment. '
                  'The complete 24-base FakeD footprint includes terminal T paired to substrate A.')
    return selective, indexing


STEPS = (
    ('rt', '1. Reverse transcription and template switching', K.REVERSE_TRANSCRIPTION),
    ('clear_rt', '2. Exchange residual RT oligos and remove RNA', K.CLEANUP),
    ('linear', '3. H-directed thermal linear amplification', K.EXTENSION),
    ('clear_h', '4. Exchange unused H primer', K.CLEANUP),
    ('random', '5. Tailed random-octamer extension toward H', K.EXTENSION),
    ('clear_random', '6. Exchange unused random primer', K.CLEANUP),
    ('selective', '7. Selective PCR with FakeTruseq-TSO and Q', K.PCR),
    ('digest', '8. BciVI digestion', K.DIGESTION),
    ('heat', '8. Heat-inactivate BciVI at 80 C for 20 minutes', K.THERMAL_TREATMENT),
    ('clear_handle', '8. Exchange released handle before ligation', K.CLEANUP),
    ('ligate', '9. Ligate pre-annealed FakeD', K.LIGATION),
    ('clear_adapter', '9. Exchange free adapter', K.CLEANUP),
    ('index', '10. Index within each compartment', K.PCR),
    ('pool', '11. Pool indexed compartments', K.POOL),
    ('purify', '11. Conventional purification after pooling', K.CLEANUP),
)


def build_model():
    """Build the ordered drawing and a typed graph with explicit material outcomes."""
    selective, indexing = pcr_geometries()
    geometry = adapter_geometry()
    initial = Scene()
    initial.strand('mRNA', [seg('transcript', 'X' * 30, 'w1', placeholder=True),
                            seg('polyA', 'A' * 18)])
    wf = Workflow(MolecularState('Compartment RNA', tuple(initial.rows())), track_graph=False)
    graph = ProtocolGraph(TITLE, authority=AUTHORITY)
    current = graph.input('rna', 'RNA in one isolated compartment', M.BIOLOGICAL,
                          state=initial, primary=True, inferred=True)
    reagents = {}
    for key, label, kind, state in (
        ('dt', 'Existing oligo-dT', M.OLIGO, legacy.OLIGO_DT),
        ('tso', 'Existing modified TSO', M.OLIGO, TSO_ORDER),
        ('h', 'H primer', M.OLIGO, H), ('random_primer', 'Candidate tailed N8 primer', M.OLIGO, RANDOM_PRIMER),
        ('fake_truseq', 'FakeTruseq-TSO', M.OLIGO, FAKE_TRUSEQ), ('q', 'Candidate Q', M.OLIGO, Q),
        ('adapter', 'FakeD pre-annealed duplex; bottom 5-prime-P', M.ADAPTER, geometry),
        ('p5', 'P5 indexing primer', M.OLIGO, P5),
        ('p7', 'Full-footprint P7 + defined compartment i7', M.OLIGO, None),
        ('rt_enzyme', 'Template-switching reverse transcriptase', M.CATALYST, None),
        ('polymerase', 'Thermal extension / blunt PCR polymerase; conditions TBD', M.CATALYST, None),
        ('bcivi', 'BciVI', M.CATALYST, None), ('ligase', 'DNA ligase; conditions TBD', M.CATALYST, None),
        ('temperature', '80 C for 20 minutes', M.CONDITION, None),
        ('qc_conditions', 'Post-pooling QC sampling; conditions TBD', M.CONDITION, None),
    ):
        reagents[key] = graph.input(key, label, kind, state=state, inferred=True)
    other = graph.input('other_indexed', 'Other compartments, already independently indexed',
                        M.LIBRARY, state=final_library(), inferred=True)
    states = {'rt': first_strand(), 'linear': sense(), 'random': random_product(),
              'selective': selective_product(), 'digest': cleavage(), 'ligate': ligated(),
              'index': final_library()}
    extras = {'rt': [('dt', R.STRUCTURAL), ('tso', R.STRUCTURAL), ('rt_enzyme', R.CATALYST)],
              'linear': [('h', R.STRUCTURAL), ('polymerase', R.CATALYST)],
              'random': [('random_primer', R.STRUCTURAL), ('polymerase', R.CATALYST)],
              'selective': [('fake_truseq', R.STRUCTURAL), ('q', R.STRUCTURAL), ('polymerase', R.CATALYST)],
              'digest': [('bcivi', R.CATALYST)], 'heat': [('temperature', R.CONDITION)],
              'ligate': [('adapter', R.STRUCTURAL), ('ligase', R.CATALYST)],
              'index': [('p5', R.STRUCTURAL), ('p7', R.STRUCTURAL), ('polymerase', R.CATALYST)]}
    current_con = None
    for key, label, kind in STEPS:
        if key in states:
            current_con = states[key]
        if key in ('rt', 'clear_rt'):
            sc = Scene(); sc.strand('first-strand cDNA', current_con.segments, mod5='OH')
            sc.labels('first-strand cDNA')
        elif key in ('linear', 'clear_h'):
            sc = Scene(); sc.strand('sense copy', current_con.segments, mod5='OH')
            sc.labels('sense copy')
        elif key in ('digest', 'heat', 'clear_handle'):
            sc = cut_scene()
        else:
            sc = duplex(current_con, junction=key in ('ligate', 'clear_adapter'))
        rows = sc.rows()
        if key == 'rt':
            wf.react('INFERRED — RT priming and template-switch binding',
                     rt_binding().rows(), name='RT and TSO binding')
        if key in ('linear', 'random'):
            binding = h_binding() if key == 'linear' else random_binding()
            wf.react('INFERRED — Anneal ' + ('H' if key == 'linear' else 'tailed random octamer'),
                     binding.rows(), name='primer binding')
        if key == 'selective':
            wf.pcr('INFERRED — ' + label, selective, rows, name=current_con.name)
        elif key == 'index':
            wf.pcr('INFERRED — ' + label, indexing, rows, name=current_con.name)
        elif key == 'ligate':
            wf.ligate_adapter('INFERRED — ' + label, geometry, rows, name=current_con.name)
        else:
            wf.react('INFERRED — ' + label, rows, name=label, note=AUTHORITY)
        ports = [(current, R.TEMPLATE)] + [(reagents[r], role) for r, role in extras.get(key, [])]
        if key == 'pool':
            ports.append((other, R.TEMPLATE))
        output = graph.product(key + '_product', label,
            M.LIBRARY if key in ('index', 'pool', 'purify') else M.MOLECULAR,
            state=current_con, inferred=True)
        outputs = [output]
        if kind is K.CLEANUP:
            waste = graph.product(key + '_waste', 'Exchanged material / possible small-product loss',
                                  M.WASTE, inferred=True)
            outputs.append(waste)
        if key == 'digest':
            handle = graph.product('released_handle', 'Q/B/site handle and internal-site loss products',
                                   M.MOLECULAR, inferred=True)
            outputs.append(handle)
        if key == 'clear_handle':
            ports.append((handle, R.TEMPLATE))
        geo = {'selective': selective, 'index': indexing, 'ligate': geometry}.get(key)
        graph.react(key, label, kind, ports, outputs, geometry=geo)
        if kind is K.CLEANUP:
            graph.discard(waste, 'Intended exchange/purification; clearance efficiency and small-product retention unestablished')
        current = output
    # Retained long templates and unwanted products are tolerated inside the compartments.
    # Their fate is unresolved rather than falsely presented as successful pre-index cleanup.
    retained = graph.product('retained_background', 'Retained long templates and unwanted amplification products',
                             M.MOLECULAR, inferred=True)
    qc = graph.product('qc', 'Pooled size assessment and P5/P7 library quantification', M.QC, inferred=True)
    library = graph.product('library', 'Purified pooled library', M.LIBRARY,
                            state=final_library(), inferred=True)
    graph.react('assess', 'Assess pooled library and account for residual background', K.SAMPLING,
                ((current, R.TEMPLATE), (reagents['qc_conditions'], R.CONDITION)),
                (library, qc, retained), changes_molecule=False,
                source=AUTHORITY + '; QC sampling conditions remain unspecified')
    graph.mark_terminal(Terminal(qc.id, Disposition.QC))
    graph.mark_terminal(Terminal(retained.id, Disposition.UNRESOLVED,
        'Long templates cannot be removed before indexing; residual background after pooled purification is unmeasured'))
    graph.sequence(library, SequencingPlan(library=final_library(), primers=PRIMERS,
                                         acquired_roles=ROLES))
    graph.declare_reviewed_scope((
        'rna', 'dt', 'tso', 'h', 'random_primer', 'fake_truseq', 'q', 'adapter',
        'p5', 'p7', 'rt_enzyme', 'polymerase', 'bcivi', 'ligase', 'temperature',
        'qc_conditions', 'other_indexed', 'rt', 'clear_rt', 'linear', 'clear_h',
        'random', 'clear_random', 'selective', 'digest', 'heat', 'clear_handle',
        'ligate', 'clear_adapter', 'index', 'pool', 'purify', 'assess', 'library',
    ), AUTHORITY)
    graph.require_publishable()  # Structural completeness only; never experimental validation.
    return wf, graph
