// Search correctness: for random 9-taxon subsets of real matrices, enumerate
// ALL 135,135 unrooted binary trees, find the true minimum length and the full
// set of most-parsimonious trees, then check that the heuristic search finds
// the same optimum and how many of the MPTs it recovers.
const fs = require('fs');
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const D = process.env.DATASETS;
const NSUB = parseInt(process.env.NSUB || '2', 10);
const K = parseInt(process.env.K || '9', 10);

function subset(mx, idx) {
  return Object.assign({}, mx, { taxa: idx.map((i) => mx.taxa[i]), cells: idx.map((i) => mx.cells[i]) });
}

function enumerate(n, sc) {
  // stepwise insertion of tips 3..n-1 into every edge of a tree rooted on tip 0
  const t = E.fromNested(n, [0, [1, 2]]);
  let best = Infinity, count = 0, total = 0;
  const optKeys = new Set();
  function rec(k) {
    if (k === n) {
      total++;
      const s = sc.length(t);
      if (s < best) { best = s; optKeys.clear(); }
      if (s === best) optKeys.add(E.treeKey(E.rootOn(t, 0)));
      return;
    }
    const w = n + k - 1;
    const edges = [];
    for (let v = 1; v < 2 * n - 1; v++) if (t.parent[v] >= 0) edges.push(v);
    for (const e of edges) {
      // insert tip k above e
      const p = t.parent[e];
      t.parent[w] = p; if (t.left[p] === e) t.left[p] = w; else t.right[p] = w;
      t.left[w] = e; t.right[w] = k; t.parent[e] = w; t.parent[k] = w;
      rec(k + 1);
      // undo
      if (t.left[p] === w) t.left[p] = e; else t.right[p] = e;
      t.parent[e] = p; t.parent[w] = -1; t.parent[k] = -1; t.left[w] = -1; t.right[w] = -1;
    }
  }
  rec(3);
  return { best, nMPT: optKeys.size, total, optKeys };
}

const rng = new E.RNG(31337);
let runs = 0, found = 0, allMPT = 0;
const files = fs.readdirSync(D).filter((x) => x.endsWith('.nex')).sort();
for (const f of files) {
  const mx = P.parseMatrix(fs.readFileSync(D + '/' + f, 'utf8'));
  for (let s = 0; s < NSUB; s++) {
    const idx = rng.shuffle([...Array(mx.taxa.length).keys()]).slice(0, K).sort((a, b) => a - b);
    const sub = subset(mx, idx);
    for (const mode of ['inapplicable', 'missing']) {
      const data = E.encode(sub, mode);
      const sc = new E.Scorer(data);
      const ex = enumerate(K, sc);
      const res = E.analyse(data, { replicates: parseInt(process.env.REPS||'3',10), ratchetIter: parseInt(process.env.RAT||'3',10), maxTrees: 200, seed: 1 + s, timeLimit: 20000 });
      const keys = new Set(res.trees.map((t) => E.treeKey(t)));
      const recovered = [...ex.optKeys].filter((k) => keys.has(k)).length;
      const ok = res.bestScore === ex.best;
      runs++; if (ok) found++; if (recovered === ex.nMPT) allMPT++;
      if (!ok || recovered !== ex.nMPT) console.log('  NOTE', f, mode, 'exhaustive', ex.best, `(${ex.nMPT} MPTs)`, 'search', res.bestScore, `(${recovered} recovered)`);
    }
  }
  process.stdout.write('.');
}
console.log(`\nTrees enumerated per subset: 135135 (K=${K})`);
console.log(`${found}/${runs} searches found the exact optimum length`);
console.log(`${allMPT}/${runs} searches recovered every most-parsimonious tree`);
process.exit(found === runs ? 0 : 1);
