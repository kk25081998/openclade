# OpenClade

[![tests](https://github.com/kk25081998/openclade/actions/workflows/test.yml/badge.svg)](https://github.com/kk25081998/openclade/actions/workflows/test.yml) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23634128.svg)](https://doi.org/10.5281/zenodo.23634128)

**Use it now:** [kartikkhosa.com/openclade](https://kartikkhosa.com/openclade/). No installation, no account; your data stays on your computer.

**New to it?** The [user guide](https://kartikkhosa.com/openclade/docs/) walks through preparing a matrix, choosing a gap treatment (with a worked example), weighting and character settings, search settings, branch support, synapomorphies, comparing analyses, reading and exporting the results, and troubleshooting.

Parsimony analysis of morphological character matrices in the web browser. Load a NEXUS or TNT file (MorphoBank exports work as they are), run a tree search, and get most parsimonious trees, consensus trees, bootstrap, jackknife and Bremer support, per-character statistics, character-state mapping, synapomorphy lists, and a methods paragraph, all in one page. Nothing is uploaded: the analysis runs locally in a Web Worker, so unpublished data stay on the user's computer.

**Status: v0.2 in preparation.** The scoring and search are validated (see below). TNT remains much faster on large matrices.

## Why

Morphological phylogenetics usually involves several programs: a matrix editor (Mesquite, MorphoBank), a search program (TNT, PAUP\*), and another tool for figures. The best treatment of inapplicable characters (Brazeau, Guillerme & Smith 2019) is available only through the TreeSearch R package. OpenClade puts the whole path from matrix to tree to figure in a single page that needs no installation.

## Features

- **Input.** NEXUS (incl. interleaved, CHARSTATELABELS, quoted names, `(01)`/`{01}` polymorphisms) and TNT `xread`, with data checks. Ordered characters, weights and exclusions are read from the file (`ASSUMPTIONS`/`SETS`/`PAUP` blocks, TNT `ccode`) and can be edited.
- **Three gap treatments:** inapplicable (Brazeau et al. 2019 algorithm, as implemented in MorphyLib/TreeSearch; cells coded as "a state or a gap", such as `{0-}`, are read as unknown, which differs from TreeSearch), missing, or extra state.
- **Equal or implied weights** (Goloboff 1993).
- **Constraints:** groups that must, or must not, form a clade. Several outgroup taxa, optionally with a monophyletic ingroup.
- **Search:** independent replicates of random addition, SPR and TBR swapping, the parsimony ratchet and tree drifting, then tree fusing (Goloboff 1999) and collection of all equally parsimonious trees reachable by TBR. Replicates run in parallel across CPU cores (Web Workers).
- **Convergence you can report.** A search stops when the best score has been found in a set number of replicates, and says in how many it was found.
- **Repeatable.** The same seed and settings give the same trees on any computer and any number of cores. An analysis file (matrix, settings, result) can be opened by anyone and run again, and the page says whether the result was reproduced.
- **Robustness report.** A common reviewer question is "what difference does the treatment of inapplicable characters make?" One click runs the same search under all three gap treatments; the Robustness tab extends this to any combination of gap treatment, weighting (equal, or implied with several k) and ordering. It reports which clades of the strict consensus are found in every analysis, marks the others on the tree, and writes the sentences for the paper.
- **Support:** bootstrap, jackknife and Bremer support.
- **Reading the trees:** strict and majority-rule consensus; optional collapsing of branches with no support; per-character steps, CI and RI; ancestral state mapping; synapomorphy lists (unambiguous, ACCTRAN, DELTRAN) for single trees and for the strict consensus.
- **Your own trees.** Trees from TNT, PAUP\* or R can be loaded, scored and mapped.
- **Data check.** Flags cells whose gap coding contradicts the character they depend on.
- **Exports:** Newick, NEXUS trees, SVG/PNG/PDF figures, tab-separated tables, an R script that re-scores the trees in TreeSearch, a MrBayes file, a TNT script, and a methods paragraph.
- **Fast rescoring:** exact indirect length calculation for ordinary characters; for characters scored with the inapplicable algorithm, a lower-bound screen followed by incremental re-evaluation.

## Use

Open `dist/index.html` in a browser, or build it:

```
DATASETS=path/to/nexus/files node build.js   # inlines engine + example matrices; also builds dist/docs (the user guide)
```

The guide's screenshots are generated from the built app: `npm i --no-save playwright sharp && npm run guide-images`.

At a release, `npm run build && npm run snapshot` copies the single-file build to `releases/v<version>/`, which the site serves at `/v/<version>/`. Analysis files record the version that made them, so an old file can always be opened in the version it came from.

## Validation

Run with `npm test`, `npm run test:data` (needs `DATASETS` pointing at the 31 published matrices in the TreeSearch package, `inst/datasets`), and `npm run test:ui`. GitHub Actions runs all of these on every push except the two slowest data tests (`test/indirect_check.js` and `test/exhaustive.js`), which are run before a release.

| Test | What it checks | Result |
|---|---|---|
| Hand-worked inapplicable cases | 41 single-character cases with expected lengths from the MorphyLib test suite, on two rootings | 41/41 correct |
| Agreement with TreeSearch | Best + 25 random trees per matrix, all 31 matrices, re-scored with TreeSearch: equal weights, implied weights (k = 10), ordered characters | 806/806 identical in each |
| Independent scoring | Engine vs a separate Sankoff dynamic-programming scorer: unordered characters (930 trees), and ordered characters, weights, exclusions and implied-weights fit (372 trees) | 930/930; 1,519/1,519 |
| Fast rearrangement scoring | Every candidate's shortcut score vs full rescoring, and the screening bound vs the exact score, across gap treatments, weights, ordered characters and implied weighting | about 6.8 million candidates, all identical |
| Exhaustive search check | 9-taxon subsets of all 31 matrices, both gap modes: all 135,135 trees enumerated | optimum found 124/124; all MPTs recovered 124/124 |
| Constraints and Bremer support | 8-taxon subsets, all 10,395 trees enumerated, equal and implied weights | 257/257, including 76 Bremer values |
| Synapomorphy lists and collapsing | All most parsimonious reconstructions enumerated on 400 small trees | 14,800/14,800 |
| Same result on any number of cores | 1, 2 and 3 worker threads and the single-thread driver: trees, bootstrap and Bremer values | 69/69 |
| Two-seed benchmark | Standard search on all 31 matrices (20–88 taxa) from two seeds | same best length on 31/31; 28 within 3 minutes on one core |
| TBR move generator | Engine neighbourhood vs brute-force TBR on 42 random trees | 42/42 identical |
| Published result | Vinther et al. (2008) machaeridian matrix | length 79, matching the TreeSearch documentation |
| Page tests | End-to-end tests of the built page in jsdom | 68/68 and 74/74 |
| Fixes from the independent review | One check per defect found, against hand-worked answers or all trees | 25/25 |

Not yet done: a comparison of search results with TNT (`validation/tnt`), which has to be run by hand because TNT cannot be redistributed. See `VALIDATION.md` for the full benchmark, what was tried for speed, and where OpenClade differs from TreeSearch on purpose.

## Known limitations

- Speed: TNT is far faster. Matrices of 60 or more taxa with many inapplicable characters can still need minutes.
- No step matrices, irreversible or Dollo characters, or continuous characters.
- Bootstrap replicates use one random-addition + SPR search each (a "quick" resampling strategy). Bremer values come from heuristic searches and are upper bounds.
- ACCTRAN and DELTRAN are offered for characters scored with ordinary parsimony only.
- The comparison of search results with TNT (`validation/tnt`) has to be run by hand, and has not been run yet.

## How to cite

Please cite OpenClade, and cite Brazeau, Guillerme & Smith (2019) for the inapplicable-character algorithm:

> Khosa, K. (2026). *OpenClade: parsimony analysis of morphological character matrices in the browser* (v0.2.0). Zenodo. https://doi.org/10.5281/zenodo.23634128

> Brazeau, M. D., Guillerme, T. & Smith, M. R. (2019). An algorithm for morphological phylogenetic analysis with inapplicable data. *Systematic Biology* 68: 619–631. https://doi.org/10.1093/sysbio/syy083

The same metadata is in `CITATION.cff` (GitHub's "Cite this repository" button).

## AI usage disclosure

OpenClade was developed with substantial assistance from Claude (Anthropic; Opus 5.5 in October 2026, Fable 5.2, Sonnet 5.5 and Opus 5.5 for version 0.2 in November 2026), which wrote much of the code, test suite and documentation under the author's direction. Scoring correctness is established by the independent tests described in `VALIDATION.md`, not by the tool that wrote the code.

## Credits

Inapplicable-character algorithm: Brazeau, M. D., Guillerme, T. & Smith, M. R. (2019) An algorithm for morphological phylogenetic analysis with inapplicable data. *Systematic Biology* 68: 619–631. The scoring code ports the logic of MorphyLib (M. D. Brazeau) as distributed in TreeSearch (M. R. Smith); both GPL-3. Example matrices are published datasets distributed with TreeSearch; cite the original authors when using them.

## License

GPL-3.0-or-later.
