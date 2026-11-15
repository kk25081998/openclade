# Comparing OpenClade's search with TNT

The TreeSearch cross-check (`../treesearch`) shows that OpenClade *scores* trees correctly. This comparison is about the *search*: does OpenClade find trees as short as TNT does, and the same strict consensus?

TNT cannot be part of the automatic tests, because its licence does not allow it to be redistributed. You run this yourself with your own copy of TNT (free from <https://www.lillo.org.ar/phylogeny/tnt/>).

**Status: the TNT side has not been run yet.** The scripts in `runs/` were written without access to TNT. `compare.js` has been tested on tree files in the same format written by OpenClade itself (`node selftest.js`). If TNT rejects a command in a script, correct `make_runs.js` and regenerate.

## Steps

```sh
git clone --depth 1 https://github.com/ms609/TreeSearch /tmp/TreeSearch
export DATASETS=/tmp/TreeSearch/inst/datasets

node make_runs.js             # writes runs/<matrix>.tnt for the 31 matrices
TNT=/path/to/tnt sh run_all.sh  # TNT writes out/<matrix>.tre
node compare.js 300           # OpenClade searches each matrix (up to 300 s) and compares
```

On Windows, open each `runs/<matrix>.tnt` in TNT with `proc runs/<matrix>.tnt;` from this folder.

## What is compared

Both programs analyse the matrix with equal weights, all characters unordered and gaps read as missing data (TNT's reading), so a tree has the same length in both.

- **Best length.** TNT's trees are scored by OpenClade, so no TNT output has to be parsed. OpenClade then runs its own Thorough search.
- **Strict consensus.** TNT collapses branches with no support; OpenClade's trees are collapsed by the same rule (minimum length zero) before the two strict consensus trees are compared.

`compare.js` prints one line per matrix and writes `openclade_vs_tnt.csv`. A matrix where TNT is shorter means OpenClade's search missed the optimum within its limits; a matrix where OpenClade is shorter would be a surprise worth reporting.
