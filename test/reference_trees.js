// Scores the published most-parsimonious trees (TreeSearch inapplicable.trees,
// from Brazeau, Guillerme & Smith 2019) under the three gap treatments.
// Consistency check: the tree optimised under a treatment should be no longer
// than the trees optimised under the other treatments, when scored that way.
// Also checks root-invariance on real data and reports scoring speed.
const fs = require('fs');
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const D = process.env.DATASETS;
const refs = JSON.parse(fs.readFileSync(process.env.REFEDGES, 'utf8'));

function edgesToNested(edges, n) {
  const ch = new Map();
  for (const [p, c] of edges) { if (!ch.has(p)) ch.set(p, []); ch.get(p).push(c); }
  const root = n + 1;
  const rec = (x) => (x <= n ? x - 1 : ch.get(x).map(rec));
  return rec(root);
}

const out = {};
let consistent = 0, total = 0, rootOK = 0;
const rng = new E.RNG(99);
for (const f of fs.readdirSync(D).filter((x) => x.endsWith('.nex')).sort()) {
  const name = f.replace('.nex', '');
  const mx = P.parseMatrix(fs.readFileSync(D + '/' + f, 'utf8'));
  const n = mx.taxa.length;
  const r = refs[name];
  if (!r || !r.inapplicable || !r.ambiguous || !r.extraState) { console.log(name.padEnd(15), '(no reference trees for the three gap treatments)'); continue; }
  const trees = {};
  for (const k of ['inapplicable', 'ambiguous', 'extraState']) trees[k] = E.rootOn(E.fromNested(n, edgesToNested(r[k], n)), 0);
  const modes = { inapplicable: 'inapplicable', ambiguous: 'missing', extraState: 'extra' };
  const row = {};
  let ok = true;
  for (const [treeKind, mode] of Object.entries(modes)) {
    const sc = new E.Scorer(E.encode(mx, mode));
    const own = sc.length(trees[treeKind]);
    const others = Object.keys(trees).filter((k) => k !== treeKind).map((k) => sc.length(trees[k]));
    row[mode] = own;
    if (others.some((o) => o < own)) { ok = false; console.log('  inconsistent', name, mode, own, others); }
  }
  // random-tree baseline + root invariance + speed (inapplicable mode)
  const sc = new E.Scorer(E.encode(mx, 'inapplicable'));
  const rt = E.randomTree(n, rng);
  const lens = new Set();
  for (let k = 0; k < n; k++) lens.add(sc.length(E.rootOn(rt, k)));
  if (lens.size === 1) rootOK++; else console.log('  root-dependent length', name, [...lens]);
  const t0 = process.hrtime.bigint();
  const reps = 300;
  for (let i = 0; i < reps; i++) sc.length(rt);
  const us = Number(process.hrtime.bigint() - t0) / 1000 / reps;
  row.random = [...lens][0];
  row.usPerTree = Math.round(us);
  out[name] = row;
  total++; if (ok) consistent++;
  console.log(name.padEnd(15), 'ntax', String(n).padStart(3), 'nchar', String(mx.nchar).padStart(4),
    '| ref MPT length  inapp', String(row.inapplicable).padStart(5), ' missing', String(row.missing).padStart(5), ' extra', String(row.extra).padStart(5),
    '| random', String(row.random).padStart(5), '|', row.usPerTree, 'µs/tree');
}
fs.writeFileSync(process.env.OUT || 'reference_scores.json', JSON.stringify(out, null, 1));
console.log(`\n${consistent}/${total} datasets: each published MPT is optimal-or-tied under its own gap treatment`);
console.log(`${rootOK}/${total} datasets: identical length on every rooting of a random tree`);
