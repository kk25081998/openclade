// Topological constraints, Bremer support and implied-weights search, checked
// by brute force: on 8-taxon subsets of real matrices all 10,395 trees are
// enumerated, so the best tree with a clade, the best tree without it and the
// true Bremer support of every clade are known exactly.
const fs = require('fs');
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const D = process.env.DATASETS;
const K = 8;
let pass = 0, fail = 0;
const check = (n, c, x) => { if (c) pass++; else { fail++; if (fail < 15) console.log('  FAIL', n, x === undefined ? '' : JSON.stringify(x)); } };
const subset = (mx, idx) => Object.assign({}, mx, { taxa: idx.map((i) => mx.taxa[i]), cells: idx.map((i) => mx.cells[i]) });

// every tree of n taxa (rooted on tip 0) with its score and clade keys
function allTrees(n, sc) {
  const t = E.fromNested(n, [0, [1, 2]]), out = [];
  (function rec(k) {
    if (k === n) { out.push({ score: sc.length(t), clades: new Set(E.clades(E.rootOn(t, 0))) }); return; }
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
  return out;
}

const rng = new E.RNG(4242);
const files = (process.env.FILES || 'Vinther2008,Longrich2010,Sansom2010,Aria2015,DeAssis2011,Asher2005,Griswold1999,Sano2011,Rousset2004,OLeary1999').split(',');
let subsets = 0, bremerClades = 0;
for (const f of files) {
  const mx = P.parseMatrix(fs.readFileSync(D + '/' + f + '.nex', 'utf8'));
  for (const variant of ['equal', 'implied']) {
    const idx = rng.shuffle([...Array(mx.taxa.length).keys()]).slice(0, K).sort((a, b) => a - b);
    const sub = subset(mx, idx);
    const data = E.encode(sub, 'inapplicable', variant === 'implied' ? { concavity: 6 } : {});
    const sc = new E.Scorer(data);
    const all = allTrees(K, sc);
    subsets++;
    const best = Math.min(...all.map((x) => x.score));
    const opt = all.filter((x) => x.score === best);
    // unconstrained search (also checks the implied-weights search against the exhaustive optimum)
    const res = E.analyse(data, { hits: 3, maxReps: 10, ratchetIter: 3, maxTrees: 500, seed: 3, timeLimit: 30000 });
    check(`${f} ${variant}: search finds the exhaustive optimum`, res.bestScore === best, [res.bestScore, best]);
    check(`${f} ${variant}: search finds every optimal tree`, res.trees.length === opt.length, [res.trees.length, opt.length]);
    // a clade that is in no optimal tree, and one that is in all of them
    const allClades = new Set(); for (const x of all) for (const c of x.clades) allClades.add(c);
    const absent = [...allClades].filter((c) => !opt.some((x) => x.clades.has(c)));
    const strict = [...opt[0].clades].filter((c) => opt.every((x) => x.clades.has(c)));
    const want = absent[rng.int(absent.length)];
    const tipsOf = (key) => E.cladeTips(key, K);
    // positive constraint
    const bestWith = Math.min(...all.filter((x) => x.clades.has(want)).map((x) => x.score));
    const nWith = all.filter((x) => x.clades.has(want) && x.score === bestWith).length;
    const pos = E.analyse(data, { hits: 3, maxReps: 10, ratchetIter: 3, maxTrees: 500, seed: 5, timeLimit: 30000, constraint: { pos: [tipsOf(want)] } });
    check(`${f} ${variant}: constrained search finds the best tree with the clade`, pos.bestScore === bestWith, [pos.bestScore, bestWith]);
    check(`${f} ${variant}: every constrained tree has the clade`, pos.trees.every((t) => E.clades(t).includes(want)));
    check(`${f} ${variant}: all best trees with the clade are found`, pos.trees.length === nWith, [pos.trees.length, nWith]);
    // the same constraint given by its complement (a set that includes the root tip)
    const comp = []; { const s = new Set(tipsOf(want)); for (let i = 0; i < K; i++) if (!s.has(i)) comp.push(i); }
    const pos2 = E.analyse(data, { hits: 2, maxReps: 6, ratchetIter: 2, maxTrees: 500, seed: 5, timeLimit: 30000, constraint: { pos: [comp] } });
    check(`${f} ${variant}: constraint given as its complement`, pos2.bestScore === bestWith && pos2.trees.every((t) => E.clades(t).includes(want)), pos2.bestScore);
    // negative constraint and Bremer support for every clade of the strict consensus
    if (strict.length) {
      const exact = new Map(strict.map((c) => [c, Math.min(...all.filter((x) => !x.clades.has(c)).map((x) => x.score))]));
      const c0 = strict[0];
      const neg = E.analyse(data, { hits: 3, maxReps: 10, ratchetIter: 3, maxTrees: 500, seed: 9, timeLimit: 30000, constraint: { neg: [tipsOf(c0)] } });
      check(`${f} ${variant}: converse constraint finds the best tree without the clade`, neg.bestScore === exact.get(c0), [neg.bestScore, exact.get(c0)]);
      check(`${f} ${variant}: no tree of the converse search has the clade`, neg.trees.every((t) => !E.clades(t).includes(c0)));
      const br = E.bremer(data, { clades: strict.slice().sort(), startTrees: res.trees.slice(0, 8).map((t) => E.toNested(t)), reps: 2, ratchetIter: 2, seed: 1, timeLimit: 60000 });
      for (const [key, bound] of br.bounds) { bremerClades++; check(`${f} ${variant}: Bremer support of a clade`, bound === exact.get(key), [bound - best, exact.get(key) - best]); }
      // a clade enforced by a constraint cannot be searched without
      const en = E.bremer(data, { clades: [c0], reps: 1, ratchetIter: 0, seed: 1, timeLimit: 10000, constraint: { pos: [tipsOf(c0)] } });
      check(`${f} ${variant}: enforced clade is skipped by the Bremer search`, en.skipped.length === 1 && en.done.length === 0, en);
    }
  }
}
// incompatible constraints are reported, not silently dropped
{
  const mx = P.parseMatrix(fs.readFileSync(D + '/Vinther2008.nex', 'utf8'));
  const data = E.encode(mx, 'inapplicable');
  let err = null;
  try { E.analyse(data, { replicates: 1, ratchetIter: 0, timeLimit: 5000, constraint: { pos: [[1, 2, 3], [3, 4, 5]] } }); } catch (e) { err = e.message; }
  check('overlapping constraints raise an error', /constraint/i.test(err || ''), err);
}
console.log(`${subsets} subsets of ${K} taxa (10,395 trees each), equal and implied weights; ${bremerClades} Bremer values checked against exhaustive enumeration`);
console.log(`${pass}/${pass + fail} constraint, Bremer and implied-weights search checks passed`);
process.exit(fail ? 1 : 0);
