// Generates the trees used to cross-check OpenClade against TreeSearch (R).
// For each of the 31 matrices shipped with TreeSearch (inst/datasets): OpenClade's
// best tree from a short search plus 25 random trees, labelled by row number
// (t1..tN, file order), with OpenClade's scores for three analyses:
//   openclade_length   equal weights, inapplicable algorithm
//   openclade_iw10     implied weights (k = 10), inapplicable algorithm
//   openclade_ordered  equal weights with the characters on line 1 of
//                      trees/<name>.ordered.txt treated as ordered (multistate
//                      characters scored by Fitch parsimony)
//   openclade_ordered_na  equal weights with the characters on line 2 ordered (all
//                      odd-numbered multistate characters, including those scored with
//                      the inapplicable algorithm), in OpenClade's TreeSearch-compatible
//                      mode (encode option orderedRegions: 'decompose')
// check.R scores the same trees with TreeSearch and compares.
// Existing tree files are kept (so earlier results stay comparable); set REGEN=1
// to search again.
// usage: DATASETS=path/to/TreeSearch/inst/datasets node make_trees.js
const fs = require('fs');
const path = require('path');
const P = require('../../src/parser.js');
const E = require('../../src/engine.js');
const D = process.env.DATASETS;
const OUT = __dirname;
const LIMIT = parseInt(process.env.LIMIT || '15000', 10);
fs.mkdirSync(path.join(OUT, 'trees'), { recursive: true });
const rows = ['dataset,tree,kind,openclade_length,openclade_iw10,openclade_ordered,openclade_ordered_na'];
const rng = new E.RNG(20261025);
for (const f of fs.readdirSync(D).filter((x) => x.endsWith('.nex')).sort()) {
  const name = f.replace('.nex', '');
  const mx = P.parseMatrix(fs.readFileSync(path.join(D, f), 'utf8'));
  const n = mx.taxa.length;
  const labels = mx.taxa.map((_, i) => 't' + (i + 1));
  const data = E.encode(mx, 'inapplicable');
  const sc = new E.Scorer(data);
  const treeFile = path.join(OUT, 'trees', name + '.tre');
  let trees;
  if (fs.existsSync(treeFile) && !process.env.REGEN) {
    trees = fs.readFileSync(treeFile, 'utf8').trim().split('\n').map((nw, i) => ({ kind: i ? 'random' : 'best', t: E.fromNested(n, E.parseNewick(nw.replace(/;\s*$/, ''), labels)) }));
  } else {
    const res = E.analyse(data, { replicates: 4, ratchetIter: 5, maxTrees: 1, seed: 1, timeLimit: LIMIT });
    trees = [{ kind: 'best', t: res.trees[0] }];
    for (let i = 0; i < 25; i++) trees.push({ kind: 'random', t: E.randomTree(n, rng) });
    fs.writeFileSync(treeFile, trees.map((x) => E.toNewick(x.t, labels)).join('\n') + '\n');
    fs.writeFileSync(path.join(OUT, 'trees', name + '.taxa.txt'), mx.taxa.join('\n') + '\n');
  }
  // implied weights
  const dIW = E.encode(mx, 'inapplicable', { concavity: 10 }), scIW = new E.Scorer(dIW);
  // ordered: every odd-numbered character with three or more states (set B);
  // set A = those of set B that are not scored with the inapplicable algorithm
  const setA = [], setB = [];
  for (let c = 0; c < mx.nchar; c += 2) {
    const s = new Set();
    let gaps = 0;
    for (let t = 0; t < n; t++) { const x = mx.cells[t][c]; if (Array.isArray(x)) x.forEach((q) => s.add(q)); else if (x === '-') gaps++; }
    if (s.size >= 3) { setB.push(c + 1); if (gaps <= 2) setA.push(c + 1); }
  }
  fs.writeFileSync(path.join(OUT, 'trees', name + '.ordered.txt'), setA.join(' ') + '\n' + setB.join(' ') + '\n');
  const withOrd = (set) => Object.assign({}, mx, { ordered: Array.from({ length: mx.nchar }, (_, c) => set.includes(c + 1)) });
  const dA = E.encode(withOrd(setA), 'inapplicable'), scA = new E.Scorer(dA);
  const dB = E.encode(withOrd(setB), 'inapplicable', { orderedRegions: 'decompose' }), scB = new E.Scorer(dB);
  trees.forEach((x, i) => rows.push(`${name},${i + 1},${x.kind},${sc.length(x.t)},${(scIW.length(x.t) / dIW.scoreDiv).toFixed(6)},${scA.length(x.t)},${scB.length(x.t)}`));
  console.log(name.padEnd(14), 'best', sc.length(trees[0].t), ' ordered characters:', setA.length, '/', setB.length);
}
fs.writeFileSync(path.join(OUT, 'expected.csv'), rows.join('\n') + '\n');
