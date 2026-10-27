# Coding in florian-PTA-rnaseq

Read `README.md` first. This is an unpublished design with a pinned `NGS_101`
submodule, not a protocol entry in that repository's public catalogue. The parent
model policy and API are in `NGS_101/AGENTS.md` and `NGS_101/DESIGN.md`.

- Initialise the submodule with `git submodule update --init`. Do not copy `lib/` here.
- Keep ordered oligos in this repo's model; import canonical sequences from the
  submodule. Mark inferred or untested designs explicitly in the model and page.
- For new reaction steps use `chemdraw.Workflow`; use `ProtocolGraph` for process
  provenance and terminal outcomes. Existing Scene panels remain legacy code.
- Give each acquired read role one active sequencing primer. Alternative mismatching
  primers may be shown in the primer table but must not enter read-content generation.
- Use only Python 3 standard library. Run `python3 tools/selftest.py` and
  `python3 build_site.py` after model, page or submodule changes.
- Do not commit source PDFs, extracted text, raw or derived sequencing data, generated
  HTML, `_site/`, or caches.
