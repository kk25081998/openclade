# OpenClade

[![tests](https://github.com/kk25081998/openclade/actions/workflows/test.yml/badge.svg)](https://github.com/kk25081998/openclade/actions/workflows/test.yml) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23634128.svg)](https://doi.org/10.5281/zenodo.23634128)

**Use it now:** [kartikkhosa.com/projects/openclade](https://kartikkhosa.com/projects/openclade/). No installation, no account; your data stays on your computer.

Parsimony analysis of morphological character matrices in the web browser. Load a NEXUS or TNT file (MorphoBank exports work as they are), run a tree search, and get most parsimonious trees, consensus trees, bootstrap or jackknife support, per-character statistics, character-state mapping, and a methods paragraph, all in one page. Nothing is uploaded: the analysis runs locally in a Web Worker, so unpublished data stay on the user's computer.

**Status: prototype (v0.1).** The scoring and search are validated (see below); the interface and speed on large matrices are still being developed.

## Why

Morphological phylogenetics usually involves several programs: a matrix editor (Mesquite, MorphoBank), a search program (TNT, PAUP\*), and another tool for figures. The best treatment of inapplicable characters (Brazeau, Guillerme & Smith 2019) is available only through the TreeSearch R package. OpenClade puts the whole path from matrix to tree to figure in a single page that needs no installation.

## Features

- NEXUS (incl. interleaved, CHARSTATELABELS, quoted names, `(01)`/`{01}` polymorphisms) and TNT `xread` input, with data checks.
- Three gap treatments: inapplicable (Brazeau et al. 2019 algorithm, as implemented in MorphyLib/TreeSearch), missing, or extra state.
- Search: random addition sequence, SPR and TBR branch swapping, parsimony ratchet, tree fusing (Goloboff 1999), and collection of all equally parsimonious trees reachable by TBR. Replicates run in parallel across CPU cores (Web Workers); time budget is shared across replicates. Seeded and reproducible for a given seed and core count.
- Fast rescoring: exact indirect length calculation for ordinary characters and incremental (per node, per character) re-evaluation of the inapplicable algorithm during branch swapping.
- Strict and majority-rule consensus; bootstrap and jackknife support.
- Per-character steps, CI and RI; ensemble CI/RI; ancestral state mapping on any MPT.
- Exports: Newick, NEXUS trees block, SVG figure, a TNT script, and an auto-generated methods paragraph.

## Use

Open `dist/index.html` in a browser, or build it:

```
DATASETS=path/to/nexus/files node build.js   # inlines engine + example matrices
```

## Validation

Run with `npm test`, `npm run test:data` (needs `DATASETS` pointing at the 31 published matrices in the TreeSearch package, `inst/datasets`), and `npm run test:ui`.

| Test | What it checks | Result |
|---|---|---|
| Hand-worked inapplicable cases | 41 single-character cases with expected lengths from the MorphyLib test suite, on two rootings | 41/41 correct |
| Root invariance | Same length on every rooting (41 test cases; random trees on 30 real matrices) | 41/41; 30/30 |
| Independent Fitch check | Engine vs a separate Sankoff dynamic-programming implementation, 930 random trees across 31 matrices | 930/930 identical |
| Fast rearrangement scoring | Every candidate's shortcut score (indirect Fitch + incremental inapplicable passes) vs full rescoring, unit and random weights, 3 gap modes | 2,028,648/2,028,648 identical |
| Parallel pipeline | Real worker threads: result, MPT count, merged bootstrap, early stop | 9/9 checks |
| TBR move generator | Engine neighbourhood vs brute-force TBR on 42 random trees | 42/42 identical |
| Exhaustive search check | 9-taxon subsets of all 31 matrices, both gap modes: all 135,135 trees enumerated | optimum found 124/124; all MPTs recovered 124/124 |
| Published result | Vinther et al. (2008) machaeridian matrix | length 79, matching the TreeSearch documentation |
| Unit and UI tests | Parsers, Newick, consensus, resampling, mapping; end-to-end page test in jsdom | all pass |

See `VALIDATION.md` for the full-size benchmark on all 31 matrices.

## Known limitations

- Speed: large matrices with many inapplicable characters (60+ taxa) still need minutes to converge; TNT is far faster. Planned: WebAssembly core, sectorial searches.
- Ordered characters, implied weighting and Bremer support are not yet implemented.
- Bootstrap replicates use one random-addition + SPR search each (a "quick" resampling strategy).

## How to cite

Please cite OpenClade, and cite Brazeau, Guillerme & Smith (2019) for the inapplicable-character algorithm:

> Khosa, K. (2026). *OpenClade: parsimony analysis of morphological character matrices in the browser* (v0.1.1). Zenodo. https://doi.org/10.5281/zenodo.23634128

> Brazeau, M. D., Guillerme, T. & Smith, M. R. (2019). An algorithm for morphological phylogenetic analysis with inapplicable data. *Systematic Biology* 68: 619–631. https://doi.org/10.1093/sysbio/syy083

The same metadata is in `CITATION.cff` (GitHub's "Cite this repository" button).

## AI usage disclosure

OpenClade was developed with substantial assistance from Claude Opus 5.5 (Anthropic), which wrote much of the initial code, test suite and documentation under the author's direction, in October 2026. Scoring correctness is established by the independent tests described in `VALIDATION.md`, not by the tool that wrote the code.

## Credits

Inapplicable-character algorithm: Brazeau, M. D., Guillerme, T. & Smith, M. R. (2019) An algorithm for morphological phylogenetic analysis with inapplicable data. *Systematic Biology* 68: 619–631. The scoring code ports the logic of MorphyLib (M. D. Brazeau) as distributed in TreeSearch (M. R. Smith); both GPL-3. Example matrices are published datasets distributed with TreeSearch; cite the original authors when using them.

## License

GPL-3.0-or-later.
