// Writes one TNT script per matrix (runs/<name>.tnt) for comparing OpenClade's
// search with TNT's. Gaps are read as missing data, which is what TNT does, with
// equal weights and all characters unordered, so the two programs score trees
// the same way and their best lengths are directly comparable.
//   DATASETS=path/to/TreeSearch/inst/datasets node make_runs.js
// Then run TNT on each script (see README.md) and node compare.js.
//
// NOT YET RUN AGAINST TNT: these scripts were written without access to TNT
// (its licence does not allow redistribution). If TNT rejects a command, fix
// it here and please report it.
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
  const script = [
    'mxram 1024;', 'nstates 32;', 'xread', `'${name}: ${mx.taxa.length} taxa, ${mx.nchar} characters. Written by OpenClade for validation.'`,
    `${mx.nchar} ${mx.taxa.length}`, ...rows, ';',
    'cc - .;',          // every character non-additive
    'hold 20000;', 'rseed 1;',
    'mult = replic 1000 tbr hold 10;',   // 1000 random addition sequences + TBR
    'bbreak = tbr;',    // swap on all trees held, to collect the equally short ones
    `export - out/${name}.tre;`,         // trees only, NEXUS format
    'zzz;', '',
  ].join('\n');
  fs.writeFileSync(path.join(OUT, name + '.tnt'), script);
  names.push(name);
}
fs.writeFileSync(path.join(__dirname, 'run_all.sh'), '#!/bin/sh\n# Runs TNT on every script. Usage: TNT=/path/to/tnt sh run_all.sh\ncd "$(dirname "$0")"\nfor f in runs/*.tnt; do echo "$f"; "${TNT:-tnt}" proc "$f"\; ; done\n');
console.log('wrote', names.length, 'TNT scripts to runs/ and run_all.sh');
