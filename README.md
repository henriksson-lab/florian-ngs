# florian-PTA-rnaseq

Work-in-progress design for reading single-cell RNA through the Atrandi SPC/PTA workflow.
The chemistry model, design notes, and rendered schematic live here because the protocol is
unpublished and does not belong in the public NGS 101 protocol catalogue.

Shared molecular primitives and the published SMART-seq/Atrandi models come from the pinned
[`NGS_101`](https://github.com/henriksson-lab/NGS_101) Git submodule.

```sh
git submodule update --init
python3 tools/selftest.py
python3 build_site.py
python3 -m http.server -d _site
```

## Working with the current NGS_101 model

The submodule is pinned to a tested revision; run `git submodule update --init` after
cloning or pulling this repo. Do not rely on whichever `NGS_101` checkout happens to be
installed elsewhere. When advancing the pin, run both `python3 tools/selftest.py` and
`python3 build_site.py` before recording the new revision. The build matters: the shared
sequencing table now requires one active primer per read role. The page can still show
the stock Read 1 primer as a mismatching alternative; it passes only the active primers
to read-content generation.

This is an unpublished design and stays outside the public NGS_101 catalogue. Its
existing `Scene` panels and construct model are legacy code. For new reaction steps, use
the shared `Workflow` and typed `ProtocolGraph` APIs. Mark oligo inputs as
source-backed or explicitly inferred, and account for material outcomes. Declare
barcode and UMI meaning as `chemdraw.feature(...)` and derive read spans with
`seqprimers.feature_spans()` when
adding run lengths. Keep canonical Illumina, SMART-seq and Atrandi sequences imported
from the submodule. Research notes and the source text should distinguish an ordered
oligo from an inferred design; local source material belongs outside Git.

The upstream [`NGS_101/AGENTS.md`](NGS_101/AGENTS.md) and
[`NGS_101/DESIGN.md`](NGS_101/DESIGN.md) describe the shared model's current contracts.

GitHub Pages is built by `.github/workflows/pages.yml`. In repository settings, select
**Pages → Build and deployment → GitHub Actions** once the repository has been pushed.
