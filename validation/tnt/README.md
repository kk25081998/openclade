# Comparing OpenClade's search with TNT

The TreeSearch cross-check (`../treesearch`) shows that OpenClade *scores* trees correctly. This comparison is about the *search*: does OpenClade find trees as short as TNT does, and the same strict consensus?

TNT cannot be part of the automatic tests, because its licence does not allow it to be redistributed. You run this yourself with your own copy of TNT (free from <https://www.lillo.org.ar/phylogeny/tnt/>).

**Status: run on 16 November 2026 with TNT 1.6 (64-bit Linux, October 2026 build) and OpenClade 0.2.0.** Results are in `openclade_vs_tnt.csv` and summarised in `../../VALIDATION.md`. TNT's tree files are not committed (`out/` is ignored); rerun the scripts to regenerate them.

TNT was made available with the sponsorship of the Willi Hennig Society. Goloboff, P. A. & Morales, M. E. (2023). TNT version 1.6, with a graphical interface for MacOS and Linux, including new routines in parallel. *Cladistics*. https://doi.org/10.1111/cla.12524

## Steps

```sh
git clone --depth 1 https://github.com/ms609/TreeSearch /tmp/TreeSearch
export DATASETS=/tmp/TreeSearch/inst/datasets

node make_runs.js             # writes runs/<matrix>.tnt and runs/<matrix>.xmult.tnt for the 31 matrices
TNT=/path/to/tnt sh run_all.sh  # TNT writes out/<matrix>.tre
node compare.js 900           # OpenClade searches each matrix (up to 900 s, the app default) and compares
# SHARD=0/2 and SHARD=1/2 split the matrices between two processes
```

On Windows, open each `runs/<matrix>.tnt` in TNT with `proc runs/<matrix>.tnt;` from this folder.

## What is compared

Both programs analyse the matrix with equal weights, all characters unordered and gaps read as missing data (TNT's reading), so a tree has the same length in both.

- **TNT's searches.** Two per matrix: 1000 random addition sequences with TBR (`mult`), and TNT's driven search (`xmult = hits 10 level 10`: sectorial search, ratchet, drifting and tree fusing until the best length has been hit 10 times). Each is followed by TBR on all trees held (up to 20,000) to collect the equally short ones.
- **Best length.** TNT's trees are scored by OpenClade, so no TNT output has to be parsed (on every matrix OpenClade's score of TNT's trees equalled the length TNT reported). OpenClade then runs its own search: replicates until the best length has been found in 10 (at most 200), 20 rounds of ratchet and drifting each, up to 5,000 trees, 900 s, one thread.
- **Strict consensus.** The shortest trees from both TNT runs, and OpenClade's trees, are collapsed by the same rule (minimum length zero) before the two strict consensus trees are compared.

`compare.js` prints one line per matrix and writes `openclade_vs_tnt.csv`. A matrix where TNT is shorter means OpenClade's search missed the optimum within its limits; a matrix where OpenClade is shorter would be a surprise worth reporting.
