# OpenClade validation

All tests use the 31 published morphological matrices distributed with the TreeSearch R package (`inst/datasets`; the benchmark set of Brazeau, Guillerme & Smith 2019), plus synthetic cases. Commands are in `package.json`. Figures below are from version 0.2.0, run on 14 November 2026.

## Correctness of scoring

| Test | Method | Result |
|---|---|---|
| Inapplicable algorithm | 41 hand-worked single-character cases with expected lengths from the MorphyLib test suite, on two rootings | 41/41 |
| Root invariance | 41 cases × 12 rootings | 41/41 |
| Fitch scoring | Engine vs an independent Sankoff dynamic-programming implementation, 930 random trees, all 31 matrices, random rootings | 930/930 identical |
| Ordered characters, weights, exclusions, implied weights | 372 random trees on all 31 matrices with 761 characters ordered at random, random weights (including 0.5) and random exclusions. Steps per character, weighted length and implied-weights fit vs an independent Sankoff scorer with \|i−j\| costs; minimum steps per character vs a brute-force search over sets of states (`test/char_model.js`) | 1,519/1,519 |
| Fast rescoring during search | Each candidate's shortcut score (indirect Fitch calculation + incremental inapplicable passes) vs full rescoring; unit and random character weights; three gap treatments; 7 matrices | 2,028,648/2,028,648 identical |
| The same, with the new character settings | 5 matrices with random ordered characters, weights and exclusions, three gap treatments, equal weights and implied weights (k = 3, 12), both ways of counting regions; SPR, TBR, ratchet weights and drifting | about 4.8 million candidates, all identical |
| Lower bound used to skip rearrangements | For every inapplicable-type column, Fitch length with gaps read as missing vs the length from the inapplicable algorithm: 217,272 column–tree pairs on all 31 matrices (`test/screen_bound.js`), and every candidate in the two rows above | bound never exceeded the exact length |
| Returned trees | Every tree returned by every benchmark run re-scored independently | all equal to the reported score |
| Changes on branches | 400 random trees × 6 characters, unordered and ordered, with polymorphic and missing cells. Every assignment of states to internal nodes is enumerated, giving all most parsimonious reconstructions (MPRs). Unambiguous changes = changes in every MPR; ACCTRAN and DELTRAN are MPRs with changes nearest the root / the tips; a collapsed branch = one with a zero-change MPR for every character (`test/reconstructions.js`) | 14,800/14,800 |
| **Agreement with TreeSearch** | OpenClade's best tree plus 25 random trees for each of the 31 matrices, scored independently with `TreeSearch::TreeLength()`, by `validation/treesearch/check.R` | see below |

### Agreement with TreeSearch

| Analysis | Result |
|---|---|
| Equal weights, inapplicable algorithm | 806/806 identical |
| Implied weights, k = 10 (fit to 6 decimals) | 806/806 identical |
| Ordered characters (odd-numbered multistate characters not scored with the inapplicable algorithm; `TreeTools::Decompose()`) | 806/806 identical |
| Ordered characters including those scored with the inapplicable algorithm, in OpenClade's TreeSearch-compatible mode | 806/806 identical |

The first row has run in GitHub Actions against TreeSearch 1.7.0 from CRAN since v0.1.1 (`crosscheck.log`). All four rows were run locally for this version against a source build of TreeSearch (development version 1.8.0.9001; see `validation/treesearch/local_crosscheck.log` for how it was built). The GitHub Actions run repeats all four against the CRAN release when these files reach `main`, and overwrites `crosscheck.log` and `openclade_vs_treesearch.csv`.

**Where OpenClade deliberately differs from TreeSearch.** For an ordered character scored with the inapplicable algorithm, TreeSearch analyses the binary characters produced by `Decompose()`, each of which keeps the gaps, so an extra applicable region costs one step *per binary character*. The cost of a region then depends on how many states the character has and on how they are numbered: states 1, 2, 3 give one more binary character than states 0, 1, 2. OpenClade's default counts an extra region once per character. The engine option `orderedRegions: 'decompose'` reproduces TreeSearch exactly (fourth row above); it is not offered in the app. The same applies under implied weights, where OpenClade (like TNT) weights an ordered character as one character and TreeSearch weights each binary character separately. The R script exported by the app states these differences when they apply.

Cells coded as "a state or a gap", such as `{0-}`, are read as unknown (`?`) unless gaps are an extra state. TreeSearch keeps the gap as one of the possibilities, which gave different lengths on 1,314 characters with such cells in 60 random matrices; neither program scores these cells exactly (OpenClade's reading equalled the exact minimum over both readings in 94% of those characters, TreeSearch's in 94% as well, and TreeSearch's was never below it). None of the 31 benchmark matrices has such cells. The page warns when a file has them.

**Exports run in the programs they are for.** The R script exported for the Vinther et al. (2008) example was run in R (50/50 trees agree with TreeSearch), as was one with three ordered characters and a constraint (18/18). The exported MrBayes files (plain, and with ordered characters and a constraint) were run in MrBayes 3.2 to the end of a short analysis; MrBayes reported the ordered characters and the constraint as set. The TNT script has not been run in TNT.

## Correctness of search

| Test | Method | Result |
|---|---|---|
| TBR move generator | Neighbourhood visited by the engine vs brute-force TBR (bisect every edge, reconnect every edge pair), 42 random trees of 5–12 taxa | 42/42 identical |
| Exhaustive optimum | 9-taxon subsets (2 per matrix) under inapplicable and missing treatments; all 135,135 trees enumerated and scored | optimum found 124/124; every MPT recovered 124/124 |
| Constraints, Bremer support, implied-weights search | 8-taxon subsets of 10 matrices, equal and implied weights; all 10,395 trees enumerated. Best tree with a clade, best tree without it, all such trees found, and the Bremer support of every clade of the strict consensus (76 values) vs the exhaustive answer (`test/constraints.js`) | 257/257 |
| Sectorial search (off by default) | 4 matrices × 2 gap treatments × equal and implied weights, in verify mode: the returned tree re-scores to the reported score and is never worse than the starting tree; where all characters are standard and weights equal, the change of length inside the sector equals the change of length of the whole tree (`test/sectorial.js`) | 16/16 |
| Fixes from the independent review | One check for each defect the review found, each against a hand-worked answer or all trees enumerated: "state or gap" cells, MATCHCHAR, lower-case MISSING/GAP symbols, TNT `piwe` before the matrix, ordered distance with unused states, gap changes as synapomorphies, implied-weights ties, NEXUS TRANSLATE with quoted names, Bremer values and bootstrap replicates under a "must not" constraint (7 taxa, all 945 trees), searches cut by the time limit, support values after re-rooting (`test/review_fixes.js`) | 25/25 |
| Published result | Vinther et al. (2008) machaeridian matrix | length 79 (52 MPTs), matching the TreeSearch documentation |
| Same result on any number of workers | Three analyses (equal weights; implied weights; constraint with ordered characters), each with bootstrap and Bremer support, on 1, 2 and 3 real worker threads and in the single-thread driver: best score, replicates used, tree set, bootstrap frequencies and Bremer values (`test/determinism.js`) | 69/69 |
| Parallel pipeline | Real worker threads: result, MPT count, merged bootstrap, early stop | 9/9 |
| Tree files, outgroups, collapsing, coding check | Newick, NEXUS (with TRANSLATE) and TNT tree files; re-rooting on several outgroup taxa; ordered characters with gaps by hand; collapsing on real trees; the coding check on a matrix with planted errors (`test/features.js`) | 34/34 |
| Web page | End-to-end tests of the built page in jsdom: the original 68 checks (load, run, map, consensus, exports, re-score exported tree, error handling, stop, gap comparison) and 74 for the new features, including a taxon or matrix name with a line break (the R script must stay free of code from the file), malformed analysis files, re-scoring exported trees under implied weights and constraints, the PDF's cross-reference table, and opening a saved analysis and reproducing it on a different number of cores | 68/68; 74/74 |
| Real browser (by hand, not part of the automated tests) | The built page driven in headless Chrome: example run, Bremer support, a six-analysis robustness comparison, SVG/PNG/PDF downloads; the PDF rasterised and inspected | no console errors |
| User guide | Every link from the app into the guide, and within it, reaches a section that exists; every image was built; the worked gap example's printed lengths (3/4, 3/3, 4/5) match the engine | 16/16 |

### Comparison with TNT

Run on 16 November 2026 with TNT 1.6 (64-bit Linux), on all 31 matrices with gaps read as missing data, equal weights and every character unordered (TNT's default reading, so a tree has the same length in both programs). Method and scripts: `validation/tnt/README.md`; results: `validation/tnt/openclade_vs_tnt.csv`.

- **TNT:** for each matrix, 1000 random addition sequences with TBR, and separately TNT's driven search (sectorial search, ratchet, drifting and fusing until the best length is hit 10 times), each followed by TBR on all trees held to collect the equally short ones. OpenClade scored every tree TNT returned, and its scores equalled the lengths TNT reported.
- **OpenClade:** replicates until the best length is found in 10 (at most 200), 20 ratchet and drift rounds each, one tree kept per collapsed tree (as TNT does), up to 5,000 trees, one thread; run with a 300-second limit (0.2.0 and 0.2.1) and a 900-second limit (the app's default from 0.2.2).

| Result | 0.2.0, 300 s | 0.2.1, 300 s | 0.2.1, 900 s |
|---|---|---|---|
| Same best length as TNT | 31 of 31 | 31 of 31 | **31 of 31** (OpenClade shorter on 0, TNT shorter on 0) |
| Strict consensus identical to TNT's | 24 of 31 | 25 of 31 | **25 of 31** |
| Searches that reached the time limit | 4 | 4 | **0** (longest 530 s) |
| OpenClade search time, median | | 43 s | 43 s |
| TNT driven search time | 1 to 54 s per matrix | | |

TNT's two searches agreed with each other on the best length for every matrix. OpenClade's times were measured with two searches sharing a 2-core machine and TNT's with one, so times are only indicative; TNT is clearly faster on the largest matrices.

**Where the strict consensus differs, OpenClade's has extra clades, because it found only some of the equally short trees.** Every collapsed tree OpenClade found was also among TNT's (checked on Capa2011, Schulze2007 and Wortley2006). There were three causes:

1. **The tree limit filled with resolutions of the same tree (fixed in 0.2.1).** OpenClade 0.2.0 kept fully resolved trees, many of which collapse to the same tree; TNT keeps one tree per collapsed tree. On Capa2011, 20,000 trees kept by 0.2.0 collapsed to 226 distinct trees, against TNT's 756. OpenClade now keeps one tree per collapsed tree (on by default in the app): on Capa2011 it finds exactly TNT's 756, and the strict consensus is identical. A test checks this against all trees enumerated for 124 small subsets of the benchmark matrices (`test/collapsed_set.js`).
2. **The time limit (not the cause at 900 s).** At 300 s the four largest matrices were cut short. At 900 s every search finished (the longest took 530 s), and each found more of the equally short trees (Zhu2013: 85 at 300 s, 308 at 900 s), but the same six strict consensus trees still differed, so the remaining gap is due to the next cause and, on Aguado2009, to the tree limit: it has about 9,600 collapsed shortest trees, and OpenClade stopped at the 5,000-tree maximum.
3. **Islands not visited (Dikow2009, Schulze2007, Wortley2006, Zanol2014, Zhu2013).** Some equally short trees are not reachable by TBR from the trees OpenClade found, and its 10 to 41 replicates landed in fewer of them than TNT's 1000. The search stops once the best length has been found 10 times, even while new equally short trees are still turning up. On Schulze2007, 40 hits instead of 10 found all 2,592 of TNT's collapsed trees; on Wortley2006, OpenClade found 255 of TNT's 356, and with up to 1000 replicates 355.

So OpenClade reaches TNT's best length on every matrix tested. Its strict consensus should be treated as possibly too resolved when the time or tree limit was reached, or when the best length was found in few of the replicates; running more replicates (Thorough effort, or a higher number of hits) is the remedy for now. A stopping rule that also waits until replicates stop finding new trees is planned.

TNT was made available with the sponsorship of the Willi Hennig Society. Goloboff, P. A. & Morales, M. E. (2023). TNT version 1.6, with a graphical interface for MacOS and Linux, including new routines in parallel. *Cladistics*. https://doi.org/10.1111/cla.12524

## Full-size benchmark

Each matrix was analysed twice (seeds 1 and 2) with the app's **Standard** effort (10 to 50 replicates, stopping when the best length has been found in 5; 10 ratchet and drift rounds each) under the inapplicable treatment, single-threaded, with a 600-second limit. The machine has 2 shared cores and ran both seeds at once, so times are indicative. The v0.1.1 column is that version's result with a 25-second budget on 2 workers.

| Summary | v0.1.1 | v0.2.0 |
|---|---|---|
| Matrices where both seeds agree on best length | 20/31 | 31/31 |
| Best length found (min of 2 seeds): better / same / worse than v0.1.1 | — | 8 / 23 / 0 |
| Matrices finished within 30 s / 180 s (slower seed) | — | 21 / 28 |

| Matrix | Taxa | v0.1.1 (seed 1 / 2) | v0.2.0 seed 1 (found in) | v0.2.0 seed 2 (found in) | Seconds (seed 1 / 2) |
|---|---|---|---|---|---|
| Agnarsson2004 | 62 | 778 / 778 | 778 (10 of 10) | 778 (10 of 10) | 24 / 24 |
| Aguado2009 | 76 | 579 / 579 | 579 (5 of 12) | 579 (5 of 12) | 19 / 20 |
| Aria2015 | 35 | 143 / 143 | 143 (5 of 10) | 143 (6 of 10) | 2 / 1 |
| Asher2005 | 23 | 345 / 345 | 345 (10 of 10) | 345 (10 of 10) | 1 / 1 |
| Capa2011 | 67 | 385 / 385 | 385 (10 of 10) | 385 (10 of 10) | 10 / 11 |
| Conrad2008 | 64 | 1763 / 1762 | 1761 (5 of 10) | 1761 (5 of 11) | 59 / 82 |
| DeAssis2011 | 33 | 64 / 64 | 64 (10 of 10) | 64 (10 of 10) | 1 / 1 |
| Dikow2009 | 88 | 1615 / 1618 | 1611 (5 of 18) | 1611 (5 of 21) | 111 / 129 |
| Eklund2004 | 54 | 440 / 441 | 440 (7 of 10) | 440 (5 of 10) | 6 / 6 |
| Geisler2001 | 68 | 1299 / 1295 | 1295 (2 of 50) | 1295 (3 of 50) | 136 / 138 |
| Giles2015 | 78 | 712 / 713 | 710 (5 of 10) | 710 (6 of 10) | 482 / 482 |
| Griswold1999 | 43 | 407 / 407 | 407 (7 of 10) | 407 (6 of 10) | 11 / 11 |
| Liljeblad2008 | 68 | 2877 / 2871 | 2868 (5 of 21) | 2868 (5 of 20) | 147 / 144 |
| Loconte1991 | 56 | 539 / 539 | 539 (10 of 10) | 539 (10 of 10) | 6 / 6 |
| Longrich2010 | 20 | 131 / 131 | 131 (10 of 10) | 131 (10 of 10) | 1 / 1 |
| OLeary1999 | 40 | 508 / 508 | 508 (10 of 10) | 508 (10 of 10) | 4 / 4 |
| OMeara2014 | 63 | 1215 / 1215 | 1215 (10 of 10) | 1215 (10 of 10) | 31 / 31 |
| Rougier2012 | 58 | 1149 / 1149 | 1149 (8 of 10) | 1149 (8 of 10) | 25 / 25 |
| Rousset2004 | 34 | 259 / 259 | 259 (10 of 10) | 259 (10 of 10) | 2 / 2 |
| Sano2011 | 36 | 223 / 223 | 223 (5 of 10) | 223 (5 of 10) | 2 / 2 |
| Sansom2010 | 23 | 189 / 189 | 189 (10 of 10) | 189 (9 of 10) | 1 / 1 |
| Schulze2007 | 52 | 164 / 164 | 164 (7 of 10) | 164 (8 of 10) | 9 / 8 |
| Shultz2007 | 59 | 434 / 434 | 434 (10 of 10) | 434 (10 of 10) | 10 / 10 |
| Sun2018 | 54 | 633 / 631 | 628 (2 of 31, time limit) | 628 (3 of 31, time limit) | 600 / 600 |
| Vinther2008 | 23 | 79 / 79 | 79 (10 of 10) | 79 (10 of 10) | 1 / 1 |
| Wetterer2000 | 63 | 559 / 560 | 559 (5 of 16) | 559 (5 of 15) | 28 / 22 |
| Wills2012 | 55 | 273 / 273 | 273 (6 of 10) | 273 (6 of 10) | 5 / 5 |
| Wilson2003 | 61 | 882 / 881 | 879 (5 of 12) | 879 (5 of 10) | 74 / 64 |
| Wortley2006 | 37 | 483 / 482 | 482 (0 of 16) | 482 (0 of 10) | 7 / 4 |
| Zanol2014 | 74 | 1318 / 1312 | 1311 (1 of 11, time limit) | 1311 (1 of 11, time limit) | 600 / 600 |
| Zhu2013 | 75 | 639 / 639 | 638 (5 of 19) | 638 (5 of 20) | 144 / 155 |

**Reading this.** "579 (5 of 12)" means the best length, 579, was found independently in 5 of the 12 replicates used. Four matrices did not meet the stopping rule, though both seeds still reached the same length:

- **Sun2018 and Zanol2014** reached the 600-second limit. They need a longer limit or more cores.
- **Geisler2001** used all 50 replicates with the best length found in only 2 or 3.
- **Wortley2006**: no single replicate found 482; it came from fusing trees of different replicates. The page reports this as "0 of N".

For these the page shows a warning and suggests Thorough. Per-replicate scores are in `test/convergence_v0.2.json`; `npm run benchmark` repeats the run.

**Why there is a minimum number of replicates.** In a first version of the rule, without a minimum, Aguado2009 with seed 1 stopped after five replicates that had all found length 580, while seed 2 found 579. Five replicates in a row had landed on the same local optimum. Standard now always runs at least 10 replicates (Thorough 20), and both seeds find 579. This is the reason to check a result with a second seed.

## Speed work

Profiling of v0.1.1 showed 97–98% of search time was spent rescoring characters handled by the inapplicable algorithm, almost all of it on rearrangements that were then rejected.

| Change | Outcome |
|---|---|
| **Lower-bound screening** (kept). The Fitch length of a column with gaps read as missing never exceeds its length under the inapplicable algorithm (tested above), and it can be computed for a rearrangement by the same indirect calculation as for ordinary characters. A rearrangement whose bound is already too long is dropped without running the inapplicable algorithm. | Starting tree, then SPR and TBR to a local optimum, same seed, same final length, v0.1.1 engine vs this one: Conrad2008 8.1 s → 1.1 s; Zhu2013 21.3 s → 1.5 s; Dikow2009 6.9 s → 0.65 s. On matrices where most characters have gaps (Zanol2014) the bound rules out fewer rearrangements and the gain is smaller. |
| The same shortcut in random addition sequences (kept) | Starting trees for 60–90 taxa in tens of milliseconds instead of seconds. |
| Tree drifting alternating with the ratchet (kept) | Similar or slightly better hit rate than the ratchet alone at the same cost. |
| Continuing a sweep after an improvement instead of restarting; climbing back with SPR before TBR (rejected) | 15–30% slower per replicate with the same hit rate (Dikow2009, Conrad2008; 12 replicates each). |
| Random sectorial search (implemented, exact, off by default) | Improved SPR-optimal trees in 9 of 16 trials, but did not raise the share of replicates reaching the best length on Dikow2009 or Conrad2008, and cost about 15% more time. May pay off above 100 taxa; `sectors` option of `replicate()`. |
| WebAssembly kernel (measured, not shipped) | A C version of the full inapplicable pass, compiled to WebAssembly, was 2.5× faster than the same pass in JavaScript, but slower than the incremental JavaScript scorer the search uses (106 vs 72 µs per tree on Zanol2014). See `experiments/wasm-kernel`. |

Earlier: incremental re-evaluation (only cells whose inputs changed, journalled and rolled back) cut a hill-climb roughly in half in v0.1.0. A different screening bound, NA(tree) ≥ NA(part A) + NA(part B), was evaluated then and rejected: it fails in 2.8% of cases.

## Bugs found and fixed by these tests

1. After trying TBR reconnections, the pruned subtree was restored through a different branch, silently changing the tree (found by re-scoring returned trees).
2. Moves of the outgroup taxon itself were never tried, making part of the tree space unreachable (found by exhaustive enumeration).
3. Re-rooting candidates were identified by the current rather than the original rooting, skipping one valid move and repeating another (found by the brute-force TBR comparison).
4. Early stop could return no trees on large matrices (found by the parallel test).
5. With gaps treated as missing data, state 0 was left out of the minimum and maximum steps of each character, so CI and RI were wrong in that treatment in v0.1.x (found by the independent minimum-steps check).
6. The MrBayes export put a non-ASCII character in a comment, which stops MrBayes (found by running the export in MrBayes).
7. The first convergence rule could stop after a few replicates that had all found the same local optimum (found by the two-seed benchmark; see above).
8. An independent review of this version, by a reviewer that had not seen the work, found the defects listed under "Found by an independent review" in `CHANGELOG.md`: among them, Bremer support and bootstrap ignored "must not" constraints, an exported R script could run text from a hostile file, and searches cut by the time limit were counted as complete. Writing the regression test for the first of these found one more: with two or more "must not" constraints, a starting tree could keep a forbidden clade, because breaking one clade could re-create another; the search could then report a tree that broke the constraint.
