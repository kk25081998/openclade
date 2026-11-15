# WebAssembly kernel experiment (not used by the app)

Question: would moving the inapplicable-algorithm scoring to WebAssembly make the search faster?

`na.c` is the four-pass scoring of all inapplicable-type columns of a tree, written in C. `bench.js` checks it against the JavaScript scorer and times three things on one matrix: the JavaScript full pass, the WebAssembly full pass, and the JavaScript incremental scorer that the search actually uses (per rearrangement, during a real TBR climb).

```
clang --target=wasm32 -O3 -msimd128 -nostdlib -Wl,--no-entry -Wl,--export=__heap_base -Wl,--initial-memory=16777216 -o na.wasm na.c
DATASETS=path/to/TreeSearch/inst/datasets node bench.js Zanol2014
```

Result (Node 22, November 2026, microseconds per tree):

| Matrix | JS full pass | WASM full pass | JS incremental (in use) |
|---|---|---|---|
| Zanol2014 (74 taxa, 115 columns) | 267 | 106 | 72 |
| Giles2015 (78 taxa, 86 columns) | 254 | 102 | 59 |
| Conrad2008 (64 taxa, 81 columns) | 152 | 62 | 36 |

The WebAssembly full pass is 2.5 times faster than the same pass in JavaScript, and gave identical lengths on 300 random trees per matrix. It is still slower than the incremental scorer already in the engine, which only recomputes the nodes a rearrangement touches. A gain would need the incremental scorer itself ported (with its change lists and undo journal), for an expected factor of about two. That was not done for this release: the lower-bound screening in the search gave 10 to 20 times on the same matrices with no change to the build.
