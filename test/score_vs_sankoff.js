// Independent check of Fitch scoring: compares OpenClade lengths (gap = missing)
// with a separate Sankoff dynamic-programming implementation (unit costs),
// which computes the same optimum by a completely different algorithm.
const fs = require('fs');
const P = require('../src/parser.js');
const E = require('../src/engine.js');

const D = process.env.DATASETS;
const NTREES = parseInt(process.env.NTREES || '30', 10);

function sankoffLength(mx, tree) {
  const n = mx.taxa.length, nc = mx.nchar, S = mx.symbols.length;
  const sym = new Map(mx.symbols.map((s, i) => [s, i]));
  const post = E.postorder(tree, []);
  let total = 0;
  const cost = new Array(2 * n - 1);
  for (let c = 0; c < nc; c++) {
    for (let t = 0; t < n; t++) {
      const cell = mx.cells[t][c];
      const v = new Float64Array(S);
      if (cell !== '?' && cell !== '-') { v.fill(Infinity); for (const s of cell) v[sym.get(s)] = 0; }
      cost[t] = v;
    }
    for (const v of post) {
      const a = cost[tree.left[v]], b = cost[tree.right[v]];
      const minA = Math.min(...a), minB = Math.min(...b);
      const out = new Float64Array(S);
      for (let s = 0; s < S; s++) out[s] = Math.min(a[s], minA + 1) + Math.min(b[s], minB + 1);
      cost[v] = out;
    }
    total += Math.min(...cost[tree.root]);
  }
  return total;
}

let total = 0, bad = 0;
const rng = new E.RNG(2024);
for (const f of fs.readdirSync(D).filter((x) => x.endsWith('.nex')).sort()) {
  const mx = P.parseMatrix(fs.readFileSync(D + '/' + f, 'utf8'));
  const data = E.encode(mx, 'missing');
  const sc = new E.Scorer(data);
  let ok = 0;
  for (let i = 0; i < NTREES; i++) {
    const t = E.rootOn(E.randomTree(mx.taxa.length, rng), rng.int(mx.taxa.length));
    const a = sc.length(t), b = sankoffLength(mx, t);
    if (a === b) ok++; else { bad++; if (bad < 5) console.log('  MISMATCH', f, a, b); }
    total++;
  }
  console.log(f.padEnd(20), `${ok}/${NTREES} identical`);
}
console.log(`\nTOTAL: ${total - bad}/${total} random trees: Fitch engine == independent Sankoff DP`);
process.exit(bad ? 1 : 0);
