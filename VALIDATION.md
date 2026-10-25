# OpenClade validation

All tests use the 31 published morphological matrices distributed with the TreeSearch R package (`inst/datasets`; the benchmark set of Brazeau, Guillerme & Smith 2019), plus synthetic cases. Commands are in `package.json`.

## Correctness of scoring

| Test | Method | Result |
|---|---|---|
| Inapplicable algorithm | 41 hand-worked single-character cases with expected lengths from the MorphyLib test suite, on two rootings | 41/41 |
| Root invariance | 41 cases × 12 rootings; random trees on 30 matrices × every rooting | 41/41; 30/30 |
| Fitch scoring | Engine vs an independent Sankoff dynamic-programming implementation, 930 random trees, all 31 matrices, random rootings | 930/930 identical |
| Fast rescoring during search | Each candidate's shortcut score (indirect Fitch calculation + incremental inapplicable passes) vs full rescoring; unit and random character weights; three gap treatments; 7 matrices (plus large gap-heavy matrices in separate runs) | 2,028,648/2,028,648 identical |
| Returned trees | Every tree returned by every benchmark run re-scored independently | all equal to the reported length |

| **Agreement with TreeSearch** | OpenClade's best tree plus 25 random trees for each of the 31 matrices, scored independently with `TreeSearch::TreeLength()` (TreeSearch 1.7.0 from CRAN, which uses the MorphyLib C library), run in GitHub Actions by `validation/treesearch/check.R` | **806/806 identical** |

The TreeSearch comparison reruns automatically whenever the cross-check inputs change; its log and full table are in `validation/treesearch/`.

## Correctness of search

| Test | Method | Result |
|---|---|---|
| TBR move generator | Neighbourhood visited by the engine vs brute-force TBR (bisect every edge, reconnect every edge pair), 42 random trees of 5–12 taxa | 42/42 identical |
| Exhaustive optimum | 9-taxon subsets (2 per matrix) under inapplicable and missing treatments; all 135,135 trees enumerated and scored | optimum found 124/124; every MPT recovered 124/124 |
| Published result | Vinther et al. (2008) machaeridian matrix | length 79 (52 MPTs), matching the TreeSearch documentation |
| Parallel pipeline | Real worker threads: result, MPT count, merged bootstrap, early stop | 9/9 |
| Web page | End-to-end test of the built page (jsdom): load, run, map, consensus, exports, re-score exported tree, error handling, stop | 23/23 |

## Full-size benchmark

Each matrix was analysed twice (seeds 1 and 2) with a deliberately short **25-second** budget, under the inapplicable treatment. "Baseline" is the first single-threaded version; "current" adds incremental rescoring, time-aware replicate scheduling, tree fusing, and 2 parallel workers (the test machine has 2 cores; typical laptops have 4–10).

| Summary | Baseline | Current |
|---|---|---|
| Matrices where both runs agree on best length | 17/31 | 20/31 |
| Best length found (min of 2 runs): better / same / worse than baseline | — | 7 / 24 / 0 |
| Total disagreement between runs (sum of steps) | 47 | 27 |

| Matrix | Taxa | Baseline (seed 1 / 2) | Current (seed 1 / 2) |
|---|---|---|---|
| Agnarsson2004 | 62 | 778 / 778 | 778 / 778 |
| Aguado2009 | 76 | 583 / 587 | 579 / 579 |
| Aria2015 | 35 | 143 / 143 | 143 / 143 |
| Asher2005 | 23 | 345 / 345 | 345 / 345 |
| Capa2011 | 67 | 385 / 385 | 385 / 385 |
| Conrad2008 | 64 | 1764 / 1767 | 1763 / 1762 |
| DeAssis2011 | 33 | 64 / 64 | 64 / 64 |
| Dikow2009 | 88 | 1619 / 1615 | 1615 / 1618 |
| Eklund2004 | 54 | 440 / 441 | 440 / 441 |
| Geisler2001 | 68 | 1299 / 1304 | 1299 / 1295 |
| Giles2015 | 78 | 718 / 718 | 712 / 713 |
| Griswold1999 | 43 | 407 / 407 | 407 / 407 |
| Liljeblad2008 | 68 | 2884 / 2877 | 2877 / 2871 |
| Loconte1991 | 56 | 539 / 539 | 539 / 539 |
| Longrich2010 | 20 | 131 / 131 | 131 / 131 |
| OLeary1999 | 40 | 508 / 508 | 508 / 508 |
| OMeara2014 | 63 | 1215 / 1215 | 1215 / 1215 |
| Rougier2012 | 58 | 1149 / 1149 | 1149 / 1149 |
| Rousset2004 | 34 | 259 / 259 | 259 / 259 |
| Sano2011 | 36 | 223 / 224 | 223 / 223 |
| Sansom2010 | 23 | 189 / 189 | 189 / 189 |
| Schulze2007 | 52 | 165 / 164 | 164 / 164 |
| Shultz2007 | 59 | 434 / 434 | 434 / 434 |
| Sun2018 | 54 | 631 / 636 | 633 / 631 |
| Vinther2008 | 23 | 79 / 79 | 79 / 79 |
| Wetterer2000 | 63 | 559 / 560 | 559 / 560 |
| Wills2012 | 55 | 273 / 273 | 273 / 273 |
| Wilson2003 | 61 | 882 / 881 | 882 / 881 |
| Wortley2006 | 37 | 483 / 482 | 483 / 482 |
| Zanol2014 | 74 | 1322 / 1331 | 1318 / 1312 |
| Zhu2013 | 75 | 642 / 646 | 639 / 639 |

**Reading this:** matrices up to about 40 taxa converge reliably within seconds. For 50–90 taxa, 25 seconds is not enough for full convergence; runs differ by a few steps. The app's default limit is 120 s, and users should run several seeds and longer limits for publication analyses, as with any heuristic search.

## Speed work

Profiling showed 97–98% of search time was spent rescoring characters handled by the inapplicable algorithm. Incremental re-evaluation (only cells whose inputs changed, journalled and rolled back) cut a single TBR hill-climb on gap-heavy matrices roughly in half (Conrad2008: 30.5 s → 14.7 s; Zhu2013: 35.3 s → 18.1 s). A candidate-screening bound for inapplicable characters was evaluated and rejected: NA(tree) ≥ NA(part A) + NA(part B) fails in 2.8% of cases, by up to 43 steps.

## Bugs found and fixed by these tests

1. After trying TBR reconnections, the pruned subtree was restored through a different branch, silently changing the tree (found by re-scoring returned trees).
2. Moves of the outgroup taxon itself were never tried, making part of the tree space unreachable (found by exhaustive enumeration).
3. Re-rooting candidates were identified by the current rather than the original rooting, skipping one valid move and repeating another (found by the brute-force TBR comparison).
4. Early stop could return no trees on large matrices (found by the parallel test).
