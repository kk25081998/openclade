# Changelog

## 0.2.0 (2026-11-15)

### Analyses

- **Ordered characters, weights and exclusions.** Read from NEXUS (`ASSUMPTIONS`, `SETS` and `PAUP` blocks: `TYPESET`, `WTSET`, `EXSET`, `CHARSET`, `CTYPE`, `WEIGHTS`, `EXCLUDE`, `OPTIONS DEFTYPE`) and TNT (`ccode`), applied by default, and editable in step 2. Previously they were ignored with a warning. Step matrices and irreversible or Dollo characters are still unsupported and are reported.
- **Implied weighting** (Goloboff 1993), with the fit reported to five decimals. TNT `piwe` in a file turns it on.
- **Topological constraints.** Groups that must, or must not, form a clade.
- **Several outgroup taxa**, with an option to require a monophyletic ingroup.
- **Bremer support**, by a converse-constraint search for each clade of the strict consensus.
- **Collapsing of branches with no support** (minimum length zero, TNT's rule 1), as an option; the number of distinct collapsed trees is always reported.
- **Synapomorphy lists.** Changes on each branch of a tree, or common to all shortest trees on the strict consensus: unambiguous changes, ACCTRAN or DELTRAN.
- **Your own trees.** Trees from TNT, PAUP* or R (Newick, NEXUS, TNT notation) are scored against the matrix and can be mapped and listed like OpenClade's own.
- **Robustness report.** The gap-treatment comparison now takes any combination of gap treatment, weighting (equal, or implied with several k) and ordering, and exports a consensus tree labelled with the number of analyses that recover each clade.
- **Data check.** Finds, for each character with gaps, the character that controls it, and lists cells coded inconsistently with it.

### Search

- **Stops on convergence, not on the clock.** A search runs numbered replicates until the best score has been found in a set number of them (2, 5 or 10 for Quick, Standard, Thorough). The result reports "found in N of M replicates".
- **Same result on any computer.** A replicate depends only on the seed and its number, so the result no longer depends on the number of CPU cores or on the speed of the machine. The time limit (now 300 s by default) only cuts a search short, and says so.
- **Several times faster on matrices with gaps** (7 to 14 times for a hill-climb on three benchmark matrices; less where most characters have gaps). Rearrangements are first screened with a lower bound (gaps read as missing data), and the inapplicable algorithm runs only on those that could still be good enough. Random addition sequences use the same shortcut.
- Tree drifting (Goloboff 1999) alternates with the ratchet. Sectorial search is implemented but off: it did not help on the benchmark matrices.

### Output

- **Saved analyses.** One file holds the matrix, settings and result; opening it restores everything, and "Run it again and compare" reports whether the result is reproduced.
- **Figures as PNG and PDF**, as well as SVG; node numbers; Bremer values below branches.
- **Exports:** an R script that re-scores the trees in TreeSearch, a MrBayes file, collapsed trees, Bremer values, and a fuller TNT script.
- The methods paragraph covers all of the above.

### Fixes

- **CI and RI were wrong with gaps treated as missing data.** State 0 was left out when the minimum and maximum number of steps of a character were worked out, in that gap treatment only. Tree lengths and trees were not affected.
- Symbols are now ordered digits first, then letters, whatever order the file declares them in.
- The MrBayes export is plain ASCII; MrBayes stops at other characters inside comments.
- Found by an independent review before release:
  - The exported R script could run text from a hostile file: a taxon or matrix name containing a line break was written into a comment line. Names are now written on one line.
  - Bremer support and bootstrap/jackknife replicates ignored groups that the user's constraint said must not form a clade, so Bremer values could be too low (even negative) and support could go to forbidden groups.
  - A Bremer clade or support replicate whose search was cut by the time limit was counted as complete, so its value depended on the clock and on the number of cores. Such searches are now dropped and reported as not done ("?", or fewer replicates).
  - Cells coded "state or gap", such as `{0-}`, were read as the state alone. They are now read as unknown (or as state-or-extra-state with gaps as an extra state), with a warning.
  - NEXUS `MATCHCHAR` is supported, and lower-case `MISSING=`/`GAP=` symbols are recognised; both were read as extra character states.
  - TNT `piwe=` before `xread` (where TNT expects it) was ignored. `cc+1 2;` without a space is now read too.
  - Ordered characters measured distance by position among the states present: with no `SYMBOLS` list and no taxon in state 2, 0 → 3 counted 2 steps instead of 3. Without a `SYMBOLS` list, states are now 0–9 then A–Z.
  - With gaps as an extra state, changes to and from the gap were missing from the unambiguous synapomorphy list and the tree drawing.
  - Under implied weights, trees with exactly the same fit could compare as different after rounding, so one could be dropped from the set of best trees.
  - NEXUS tree files written by OpenClade could not be read back when a taxon name contained an apostrophe or comma.
  - Majority-rule percentages and state dots were missing on some branches when the tree was rooted on several outgroup taxa.
  - A malformed analysis file could leave the previous file's saved-analysis notice on screen, or show "undefined trees".
  - The TNT script now says when the analysis used constraints that it does not include.
  - Found while testing these fixes: with two or more "must not" constraints, breaking one forbidden clade in a starting tree could re-create another, so a search could return a tree that broke the constraint.

### Also in 0.2.0 (made after the 0.1.1 release)

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
- **User guide** at kartikkhosa.com/openclade/docs/ covers:
  - preparing a matrix, and what each data warning means;
  - the three gap treatments, with a worked "tail colour" example whose lengths are checked by the tests;
  - rooting the tree, search settings and branch support;
  - reading and exporting the results, and comparing gap treatments;
  - citing, troubleshooting, limitations and a glossary.
- **Clearer app layout and wording:**
  - a top bar with links to the guide;
  - steps that say what to do (Load a matrix, Choose settings, Run), each setting linked to its section of the guide;
  - Run always in view on wide screens;
  - plain-language summary numbers and a line under each tab saying what it shows;
  - a notice explaining the example analysis;
  - Hide settings, to give the results the full width;
  - wider tree drawings, so support values and change labels have room.
- **Download buttons** for every export.
- The saved SVG figure now uses real colours on a white background. It previously used the page's colour variables, so the tree was invisible when the file was opened elsewhere.
- **The data check warns when a file orders, weights or excludes characters** (NEXUS ASSUMPTIONS blocks, TNT `ccode`). (Superseded in 0.2.0, which applies these settings.)

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
