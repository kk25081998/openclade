# Changelog

## Unreleased

- **Compare gap treatments.** One button runs the same search with gaps as inapplicable, as missing data and as an extra state, and shows how the results differ:
  - a headline count of the clades found under all three treatments;
  - tree lengths, MPT counts and consensus sizes for each treatment;
  - the strict consensus tree, with branches that depend on the treatment drawn in red and marked with which treatments recover them;
  - a table of every treatment-dependent clade;
  - pairwise Robinson–Foulds distances;
  - Newick and tab-separated exports;
  - a methods-and-results paragraph for the paper.

  Each comparison keeps at least 1,000 trees per treatment, because a strict consensus of a truncated tree set is over-resolved and would exaggerate the differences. Searches that hit a limit are flagged.
- Loading a new matrix cancels a running analysis.

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
