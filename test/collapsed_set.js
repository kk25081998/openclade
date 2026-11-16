// Keeping one tree per collapsed tree (the "collapse" search option, as in TNT).
// For 8-taxon subsets of real matrices, every tree is enumerated (10,395) and the
// distinct collapsed shortest trees are worked out. The search with collapse on must
// find exactly that set, report one tree per collapsed tree, and give the same
// strict consensus as the full set. With a small tree limit, collapse on must keep
// more distinct trees than collapse off.
const fs = require('fs');
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const D = process.env.DATASETS;
if (!D) { console.error('Set DATASETS to the folder of the 31 TreeSearch matrices.'); process.exit(1); }
const K = 8;
let pass = 0, fail = 0, cases = 0, withMany = 0;
const check = (n, c, x) => { if (c) pass++; else { fail++; console.log('  FAIL', n, x === undefined ? '' : JSON.stringify(x)); } };

function allTrees(n, f) {
  const t = E.fromNested(n, [0, [1, 2]]);
  (function rec(k) {
    if (k === n) { f(t); return; }
    const w = n + k - 1, edges = [];
    for (let v = 1; v < 2 * n - 1; v++) if (t.parent[v] >= 0) edges.push(v);
    for (const e of edges) {
      const p = t.parent[e];
      t.parent[w] = p; if (t.left[p] === e) t.left[p] = w; else t.right[p] = w;
      t.left[w] = e; t.right[w] = k; t.parent[e] = w; t.parent[k] = w;
      rec(k + 1);
      if (t.left[p] === w) t.left[p] = e; else t.right[p] = e;
      t.parent[e] = p; t.parent[w] = -1; t.parent[k] = -1; t.left[w] = -1; t.right[w] = -1;
    }
  })(3);
}
function rng(seed) { let s = seed >>> 0; return () => { s = (s + 0x6D2B79F5) >>> 0; let t = s; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
const ckey = (keys) => keys.slice().sort().join('|');
const strictOf = (sets) => sets[0].filter((c) => sets.every((s) => s.includes(c))).sort().join('|');

const files = fs.readdirSync(D).filter((f) => f.endsWith('.nex')).sort();
const rnd = rng(42);
for (const f of files) for (const mode of ['missing', 'inapplicable']) {
  const mx = P.parseMatrix(fs.readFileSync(D + '/' + f, 'utf8'));
  for (let rep = 0; rep < 2; rep++) {
    const idx = [...mx.taxa.keys()].sort(() => rnd() - 0.5).slice(0, K).sort((a, b) => a - b);
    const sub = Object.assign({}, mx, { taxa: idx.map((i) => mx.taxa[i]), cells: idx.map((i) => mx.cells[i]) });
    const data = E.encode(sub, mode), sc = new E.Scorer(data);
    let best = Infinity; const opt = [];
    allTrees(K, (t) => { const s = sc.length(t); if (s < best) { best = s; opt.length = 0; } if (s === best) opt.push(E.rootOn(t.clone ? t.clone() : E.fromNested(K, E.toNested(t)), 0)); });
    const truth = new Map(); for (const c of E.collapseAll(data, opt)) truth.set(ckey(c.keys), c.keys);
    const res = E.analyse(data, { hits: 5, minReps: 10, maxReps: 50, ratchetIter: 10, maxTrees: 5000, seed: 1, timeLimit: 60000, collapse: true }); // the Standard preset
    const got = E.collapseAll(data, res.trees);
    const gotKeys = new Set(got.map((c) => ckey(c.keys)));
    cases++; if (truth.size > 1) withMany++;
    if (process.env.DEBUG && gotKeys.size !== truth.size) {
      const off = E.analyse(data, { hits: 3, maxReps: 20, ratchetIter: 3, maxTrees: 5000, seed: 1, timeLimit: 60000 });
      const offKeys = new Set(E.collapseAll(data, off.trees).map((c) => ckey(c.keys)));
      console.log('DEBUG', f, mode, rep, 'binary MPTs', opt.length, 'off finds binary', off.trees.length, 'collapsed', offKeys.size, 'on collapsed', gotKeys.size, 'truth', truth.size);
    }
    check(`${f} ${mode} #${rep}: best length`, res.bestScore === best, [res.bestScore, best]);
    check(`${f} ${mode} #${rep}: one tree per collapsed tree`, got.length === res.trees.length, [got.length, res.trees.length]);
    check(`${f} ${mode} #${rep}: all ${truth.size} collapsed shortest trees found`, gotKeys.size === truth.size && [...truth.keys()].every((k) => gotKeys.has(k)), [gotKeys.size, truth.size]);
    check(`${f} ${mode} #${rep}: strict consensus`, strictOf(got.map((c) => c.keys)) === strictOf([...truth.values()]));
  }
}
// a small tree limit: collapse keeps more distinct trees
{
  const mx = P.parseMatrix(fs.readFileSync(D + '/Capa2011.nex', 'utf8'));
  const data = E.encode(mx, 'missing');
  const o = { hits: 3, maxReps: 10, ratchetIter: 3, maxTrees: 100, seed: 1, timeLimit: 120000 };
  const off = E.analyse(data, o), on = E.analyse(data, Object.assign({}, o, { collapse: true }));
  const a = E.collapseAll(data, off.trees).length, b = E.collapseAll(data, on.trees).length;
  check('Capa2011, 100-tree limit: collapse keeps more distinct trees', b > a && b === on.trees.length, [a, b]);
}
console.log(`${pass}/${pass + fail} checks of keeping one tree per collapsed tree passed (${cases} subsets, ${withMany} with more than one collapsed shortest tree)`);
process.exit(fail ? 1 : 0);
