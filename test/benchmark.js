// Full-size benchmark: runs the default search on every published matrix
// (inapplicable gap treatment) with two independent random seeds and a fixed
// time budget, recording best length, number of MPTs and time. Agreement
// between seeds indicates the search is converging on the optimum.
const fs = require('fs');
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const D = process.env.DATASETS;
const LIMIT = parseInt(process.env.LIMIT || '25000', 10);
const out = [];
for (const f of fs.readdirSync(D).filter((x) => x.endsWith('.nex')).sort()) {
  const mx = P.parseMatrix(fs.readFileSync(D + '/' + f, 'utf8'));
  const data = E.encode(mx, 'inapplicable');
  const row = { dataset: f.replace('.nex', ''), ntax: mx.taxa.length, nchar: mx.nchar, naChars: data.naChars.length, runs: [] };
  for (const seed of [1, 2]) {
    const t0 = Date.now();
    const res = E.analyse(data, { replicates: 5, ratchetIter: 10, maxTrees: 100, seed, timeLimit: LIMIT });
    // independent re-score of every returned tree
    const sc = new E.Scorer(data);
    const allOk = res.trees.every((t) => sc.length(t) === res.bestScore);
    row.runs.push({ seed, best: res.bestScore, trees: res.trees.length, ms: Date.now() - t0, timedOut: res.timedOut, rescoreOk: allOk });
  }
  out.push(row);
  const [a, b] = row.runs;
  console.log(row.dataset.padEnd(14), String(row.ntax).padStart(3), 'taxa', String(row.nchar).padStart(4), 'chars |',
    'seed1', a.best, `(${a.trees} MPTs, ${(a.ms / 1000).toFixed(1)}s${a.timedOut ? ', time limit' : ''})`, '| seed2', b.best, `(${b.trees}, ${(b.ms / 1000).toFixed(1)}s${b.timedOut ? ', time limit' : ''})`,
    a.best === b.best ? 'AGREE' : 'DIFFER', a.rescoreOk && b.rescoreOk ? '' : 'RESCORE-FAIL');
}
fs.writeFileSync(process.env.OUT || 'benchmark.json', JSON.stringify(out, null, 1));
