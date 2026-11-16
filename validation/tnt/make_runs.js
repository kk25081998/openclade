// Writes one TNT script per matrix (runs/<name>.tnt) for comparing OpenClade's
// search with TNT's. Gaps are read as missing data, which is what TNT does, with
// equal weights and all characters unordered, so the two programs score trees
// the same way and their best lengths are directly comparable.
//   DATASETS=path/to/TreeSearch/inst/datasets node make_runs.js
// Then run TNT on each script (see README.md) and node compare.js.
//
// Two scripts per matrix: runs/<name>.tnt (1000 random addition sequences + TBR)
// and runs/<name>.xmult.tnt (TNT's driven search, xmult). Run with TNT 1.6.
const fs = require('fs');
const path = require('path');
const P = require('../../src/parser.js');
const D = process.env.DATASETS;
if (!D) { console.error('Set DATASETS to TreeSearch/inst/datasets'); process.exit(1); }
const OUT = path.join(__dirname, 'runs');
fs.mkdirSync(OUT, { recursive: true });
fs.mkdirSync(path.join(__dirname, 'out'), { recursive: true });
const names = [];
for (const f of fs.readdirSync(D).filter((x) => x.endsWith('.nex')).sort()) {
  const name = f.replace('.nex', '');
  const mx = P.parseMatrix(fs.readFileSync(path.join(D, f), 'utf8'));
  // taxa are written as t1..tN (matrix order) so that names cannot be misread
  const rows = mx.taxa.map((_, i) => 't' + (i + 1) + ' ' + mx.cells[i].map((c) => (c === '?' || c === '-' ? c : c.length === 1 ? c[0] : '[' + c.join('') + ']')).join(''));
  const head = [
    'mxram 1024;', 'nstates 32;', 'xread', `'${name}: ${mx.taxa.length} taxa, ${mx.nchar} characters. Written by OpenClade for validation.'`,
    `${mx.nchar} ${mx.taxa.length}`, ...rows, ';',
    'cc - .;',          // every character non-additive
    'hold 20000;', 'rseed 1;',
  ];
  const tail = (file) => [
    'bbreak = tbr;',    // swap on all trees held, to collect the equally short ones
    'taxname=;',        // write taxon names (t1 ... tN), not numbers
    `export - out/${file};`, // trees only, NEXUS format
    'zzz;', '',
  ];
  // 1. 1000 random addition sequences + TBR
  fs.writeFileSync(path.join(OUT, name + '.tnt'), head.concat(['mult = replic 1000 tbr hold 10;'], tail(name + '.tre')).join('\n'));
  // 2. TNT's driven search: sectorial search, ratchet, drifting and tree fusing until
  //    the best length has been found 10 times (the slack gives sectorial search memory)
  fs.writeFileSync(path.join(OUT, name + '.xmult.tnt'), head.concat(['sect: slack 60;', 'xmult = hits 10 level 10;'], tail(name + '.xmult.tre')).join('\n'));
  names.push(name);
}
const runAll = [
  '#!/bin/sh',
  '# Runs TNT on every script. Usage: TNT=/path/to/tnt sh run_all.sh',
  'cd "$(dirname "$0")" || exit 1',
  'mkdir -p out',
  'for f in runs/*.tnt; do',
  '  echo "$f"',
  '  "${TNT:-tnt}" proc "$f"\\; || echo "  TNT failed on $f"',
  'done',
  '',
].join('\n');
fs.writeFileSync(path.join(__dirname, 'run_all.sh'), runAll);
console.log('wrote', 2 * names.length, 'TNT scripts to runs/ and run_all.sh');
