# Changelog

## 0.1.1 (2026-10-25)

- First release archived on Zenodo.
- The app footer links to the source code and issue tracker.

## 0.1.0 (2026-10-24)

First public release.

- NEXUS and TNT matrix import with data checks (MorphoBank exports supported).
- Equal-weights parsimony with three gap treatments: inapplicable (Brazeau, Guillerme & Smith 2019), missing, extra state.
- Search: random addition sequence, SPR and TBR, parsimony ratchet, tree fusing, collection of equally parsimonious trees; parallel Web Workers.
- Exact indirect length calculation for Fitch characters and incremental rescoring of inapplicable characters during branch swapping.
- Strict and majority-rule consensus, bootstrap and jackknife support, per-character steps, CI and RI, ancestral state mapping.
- Exports: Newick, NEXUS trees, SVG figure, TNT script, methods paragraph.
- Validation suite and full-size benchmark on 31 published matrices (see VALIDATION.md).
