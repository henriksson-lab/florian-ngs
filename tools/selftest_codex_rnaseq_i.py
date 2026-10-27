"""Source transcription and cross-API regression checks for the proposed pilot."""
from pathlib import Path

import codex_rnaseq_i as c


def run(check):
    check.section('Codex RNAseq I: original proposal transcription')
    source = (Path(__file__).resolve().parent.parent /
              'reference/codex-rnaseq-i-original.md').read_text(encoding='utf-8')
    lines = source.splitlines()
    for label, actual in (
        ('RT oligo-dT', c.legacy.OLIGO_DT),
        ('Existing TSO', c.TSO_ORDER), ('H linear primer', c.H),
        ('Q PCR primer', c.Q), ('BciVI tailed random primer', c.RANDOM_PRIMER),
        ('FakeTruseq-TSO PCR primer / custom Read 1 sequencing primer', c.FAKE_TRUSEQ),
        ('FakeD top (terminal 3-prime-T overhang)', c.FAKED_TOP),
        ('FakeD bottom', '/5Phos/' + c.FAKED_BOTTOM), ('P5 indexing primer', c.P5),
    ):
        check(label + ' agrees with preserved input', actual, lines[lines.index(label) + 1])
    check('P7 architecture agrees with preserved input',
          c.il.P7 + ' + selected 10-base i7 + ' + c.il.TRUSEQ_READ2,
          lines[lines.index('Full-footprint P7-i7 primer architecture') + 1])

    check.section('Codex RNAseq I: cleavage, PCR and sequencing integration')
    lib = c.final_library()
    read2 = c.sp.locate(lib, c.sp.TRUSEQ['R2'])
    # The entire terminal T is consumed by the primer footprint at the sealed junction.
    start, end = lib.span('random-primer-derived N8')
    check('Read 2 footprint ends exactly before the random octamer', read2.start, end)
    check('Read 2 then crosses N8 and reaches cDNA with no fixed prefix',
          c.sp.read_contents(lib, c.PRIMERS, acquired_roles=c.ROLES)[-1].content.startswith(
              'random-primer-derived N8 (cycles 1–8); cDNA (from cycle 9;'), True)
    read1 = c.sp.locate(lib, c.PRIMERS[0])
    check('Custom Read 1 starts exactly at RT UMI', read1.end, lib.span('RT UMI N8')[0])
    spans = c.sp.feature_spans(lib, c.PRIMERS, {'Read 1': 50, 'Index 1 (i7)': 10, 'Read 2': 50})
    check('W2 and random priming do not become counting identifiers',
          [(s.read, s.feature.role, s.cycles) for s in spans],
          [('Read 1', 'umi', '1–8'), ('Index 1 (i7)', 'sample_index', '1–10')])
    _, indexing = c.pcr_geometries()
    check('P7 indexes via the whole original FakeD including terminal T',
          indexing.reverse_hit.end - indexing.reverse_hit.start, len(c.FAKED_TOP))
    try:
        c.p7_primer('N' * 10)
    except ValueError:
        rejects_random_index = True
    else:
        rejects_random_index = False
    check('An index architecture cannot accidentally become an N10 oligo order', rejects_random_index)

    workflow, graph = c.build_model()
    check('Typed graph accounts for inputs, outputs and sequencing', graph.publication_problems(), [])
    check('Retained background has an unresolved terminal fate',
          graph.terminals['retained_background'].disposition.value, 'unresolved')
    # Each rendered transition must consume the preceding state, including primer binding.
    check('Workflow remains connected through PCR and ligation binding states',
          all(r.before is previous for r, previous in zip(
              workflow.reactions, [workflow.initial] + [r.after for r in workflow.reactions[:-1]])))
