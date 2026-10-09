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

GitHub Pages is built by `.github/workflows/pages.yml`. In repository settings, select
**Pages → Build and deployment → GitHub Actions** once the repository has been pushed.
