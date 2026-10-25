// Generates the trees used to cross-check OpenClade against TreeSearch (R).
// For each of the 31 matrices shipped with TreeSearch (inst/datasets), writes
// OpenClade's best tree from a short search plus 25 random trees, labelled by
// row number (t1..tN, file order), together with OpenClade's tree lengths
// under the inapplicable algorithm. check.R scores the same trees with
// TreeSearch::TreeLength and compares.
// usage: DATASETS=path/to/TreeSearch/inst/datasets node make_trees.js
const fs = require('fs');
const path = require('path');
const P = require('../../src/parser.js');
const E = require('../../src/engine.js');
const D = process.env.DATASETS;
const OUT = __dirname;
const LIMIT = parseInt(process.env.LIMIT || '15000', 10);
fs.mkdirSync(path.join(OUT, 'trees'), { recursive: true });
const rows = ['dataset,tree,kind,openclade_length'];
const rng = new E.RNG(20261025);
for (const f of fs.readdirSync(D).filter((x) => x.endsWith('.nex')).sort()) {
  const name = f.replace('.nex', '');
  const mx = P.parseMatrix(fs.readFileSync(path.join(D, f), 'utf8'));
  const n = mx.taxa.length;
  const labels = mx.taxa.map((_, i) => 't' + (i + 1));
  const data = E.encode(mx, 'inapplicable');
  const sc = new E.Scorer(data);
  const res = E.analyse(data, { replicates: 4, ratchetIter: 5, maxTrees: 1, seed: 1, timeLimit: LIMIT });
  const trees = [{ kind: 'best', t: res.trees[0] }];
  for (let i = 0; i < 25; i++) trees.push({ kind: 'random', t: E.randomTree(n, rng) });
  const lines = [];
  trees.forEach((x, i) => {
    lines.push(E.toNewick(x.t, labels));
    rows.push(`${name},${i + 1},${x.kind},${sc.length(x.t)}`);
  });
  fs.writeFileSync(path.join(OUT, 'trees', name + '.tre'), lines.join('\n') + '\n');
  fs.writeFileSync(path.join(OUT, 'trees', name + '.taxa.txt'), mx.taxa.join('\n') + '\n');
  console.log(name.padEnd(14), 'best', res.bestScore);
}
fs.writeFileSync(path.join(OUT, 'expected.csv'), rows.join('\n') + '\n');
