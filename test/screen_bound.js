// The search skips a rearrangement without running the inapplicable algorithm when a
// cheap lower bound already rules it out. The bound for an inapplicable-type column is
// its Fitch length with gaps read as missing data. This test checks, column by column
// on random and near-optimal trees of all 31 matrices, that the bound never exceeds the
// length the inapplicable algorithm gives. (test/indirect_check.js checks the same
// inequality for every rearrangement scored in verify mode.)
const fs = require('fs');
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const D = process.env.DATASETS;
const NTREES = parseInt(process.env.NTREES || '60', 10);
const rng = new E.RNG(5);
let pairs = 0, viol = 0, sumB = 0, sumL = 0, files = 0;

// Fitch steps of every column on the screening tips
function screenSteps(data, tree) {
  const n = data.ntax, nc = data.nchar, m = 2 * n - 1;
  const D1 = new Int32Array(m * nc), cs = new Int32Array(nc);
  D1.set(data.screenTips);
  for (const v of E.postorder(tree, [])) {
    const o = v * nc, l = tree.left[v] * nc, r = tree.right[v] * nc;
    for (let j = 0; j < nc; j++) { let x = D1[l + j] & D1[r + j]; if (x === 0) { x = D1[l + j] | D1[r + j]; cs[j]++; } D1[o + j] = x; }
  }
  return cs;
}

for (const f of fs.readdirSync(D).filter((x) => x.endsWith('.nex')).sort()) {
  const mx = P.parseMatrix(fs.readFileSync(D + '/' + f, 'utf8'));
  files++;
  // unordered, and with every third multistate character ordered (both ways of counting regions)
  const nstates = (c) => { const s = new Set(); for (const r of mx.cells) if (Array.isArray(r[c])) r[c].forEach((q) => s.add(q)); return s.size; };
  const ord = Array.from({ length: mx.nchar }, (_, c) => c % 3 === 0 && nstates(c) > 2);
  const variants = [E.encode(mx, 'inapplicable'), E.encode(Object.assign({}, mx, { ordered: ord }), 'inapplicable'), E.encode(Object.assign({}, mx, { ordered: ord }), 'inapplicable', { orderedRegions: 'decompose' })];
  const trees = [];
  for (let i = 0; i < NTREES; i++) trees.push(E.randomTree(mx.taxa.length, rng));
  const S = new E.Searcher(variants[0], { seed: 3 });
  for (let i = 0; i < 3; i++) { const t = S.wagner(); trees.push(t.clone()); trees.push(S.climb(t, null, 'spr', Date.now() + 1500).tree); }
  for (const data of variants) {
    const sc = new E.Scorer(data);
    for (const t of trees) {
      const exact = sc.colSteps(t), bound = screenSteps(data, t);
      for (const j of data.naChars) {
        pairs++; sumB += bound[j]; sumL += exact[j];
        if (bound[j] > exact[j]) { viol++; if (viol < 6) console.log('  VIOLATION', f, 'column', j, 'bound', bound[j], 'exact', exact[j]); }
      }
      for (const j of data.stdChars) if (bound[j] !== exact[j]) { viol++; if (viol < 6) console.log('  standard column differs', f, j); }
    }
  }
}
console.log(`${files} matrices, ${pairs} column-tree pairs (unordered, ordered, and TreeSearch-compatible ordered): bound exceeded the exact length ${viol} times; bound / exact = ${(sumB / sumL).toFixed(3)} on average`);
process.exit(viol ? 1 : 0);
