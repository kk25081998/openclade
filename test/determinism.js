// The same analysis must give the same result on 1, 2 or 3 workers and in the
// single-thread driver: best score, replicates used, the set of trees,
// bootstrap frequencies and Bremer values. Uses real worker_threads.
const fs = require('fs');
const path = require('path');
const { Worker } = require('worker_threads');
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const { runParallel } = require('../src/parallel.js');
const D = process.env.DATASETS;
const workerCode = `
const { parentPort } = require('worker_threads');
globalThis.OCEngine = require(${JSON.stringify(path.join(__dirname, '../src/engine.js'))});
require(${JSON.stringify(path.join(__dirname, '../src/worker.js'))});
parentPort.on('message', (m) => globalThis.OCWorkerHandle(m, (x) => parentPort.postMessage(x)));
`;
function spawn() {
  const w = new Worker(workerCode, { eval: true });
  const wrap = { onmessage: null, postMessage: (m) => w.postMessage(m), terminate: () => w.terminate() };
  w.on('message', (m) => wrap.onmessage && wrap.onmessage(m));
  w.on('error', (e) => console.log('worker error', e));
  return wrap;
}
function run(cfg, nWorkers) {
  return new Promise((resolve) => {
    const out = {};
    runParallel(Object.assign({ engine: E, spawn, nWorkers, onTrees: (r) => { out.trees = r; }, onSupport: (s) => { out.support = s; }, onBremer: (b) => { out.bremer = b; },
      onError: (e) => { out.error = e; resolve(out); }, onDone: () => resolve(out) }, cfg));
  });
}
let pass = 0, fail = 0;
const check = (n, c, x) => { if (c) pass++; else { fail++; console.log('  FAIL', n, x === undefined ? '' : JSON.stringify(x)); } };
const sig = (mx, r) => ({
  best: r.trees.best, reps: r.trees.reps, hits: r.trees.hits, timedOut: r.trees.timedOut,
  trees: r.trees.trees.map((nt) => E.treeKey(E.fromNested(mx.taxa.length, nt))).sort().join('#'),
  support: r.support ? r.support.freq.map(([k, f]) => k + '=' + f.toFixed(6)).sort().join(' ') : '',
  bremer: r.bremer ? r.bremer.values.map((v) => v.key + '=' + v.bound).sort().join(' ') : '',
});

(async () => {
  const cases = [
    { file: 'Vinther2008', gapMode: 'inapplicable', charOpts: null, constraint: null, label: 'equal weights' },
    { file: 'Longrich2010', gapMode: 'inapplicable', charOpts: { concavity: 12 }, constraint: null, label: 'implied weights k=12' },
    { file: 'Aria2015', gapMode: 'missing', charOpts: null, constraint: { pos: [[3, 4, 5]] }, label: 'constraint, some characters ordered', order: true },
    { file: 'Capa2011', gapMode: 'missing', charOpts: null, constraint: null, label: 'one tree per collapsed tree', collapse: true },
  ];
  for (const c of cases) {
    const mx = P.parseMatrix(fs.readFileSync(D + '/' + c.file + '.nex', 'utf8'));
    if (c.order) mx.ordered = Array.from({ length: mx.nchar }, (_, i) => i % 4 === 0);
    const cfg = { mx: { taxa: mx.taxa, cells: mx.cells, nchar: mx.nchar, symbols: mx.symbols, ordered: mx.ordered, weights: mx.weights, excluded: mx.excluded },
      gapMode: c.gapMode, charOpts: c.charOpts, constraint: c.constraint,
      opts: { hits: 3, maxReps: 12, ratchetIter: 4, maxTrees: 200, seed: 11, timeLimit: 120000, mode: 'tbr', collapse: !!c.collapse },
      support: { method: 'bootstrap', replicates: 30 }, bremer: { enabled: true, reps: 2, ratchetIter: 1 } };
    const sigs = [];
    for (const W of [1, 2, 3]) { const r = await run(cfg, W); check(`${c.file}: no error on ${W} worker(s)`, !r.error, r.error); sigs.push(sig(mx, r)); }
    const data = E.encode(mx, c.gapMode, c.charOpts);
    const single = E.analyse(data, Object.assign({}, cfg.opts, { constraint: c.constraint }));
    const sc = new E.Scorer(data);
    console.log(`  ${c.file} (${c.label}): score ${(sigs[0].best / data.scoreDiv).toFixed(data.iw ? 5 : 0)}, found by ${sigs[0].hits} of ${sigs[0].reps} replicates, ${sigs[0].trees.split('#').length} trees, ${sigs[0].bremer.split(' ').length} Bremer values`);
    check(`${c.file}: not cut short by the time limit`, !sigs[0].timedOut);
    for (const key of ['best', 'reps', 'hits', 'trees', 'support', 'bremer']) {
      check(`${c.file}: same ${key} on 1 and 2 workers`, sigs[0][key] === sigs[1][key], key === 'best' || key === 'reps' || key === 'hits' ? [sigs[0][key], sigs[1][key]] : undefined);
      check(`${c.file}: same ${key} on 1 and 3 workers`, sigs[0][key] === sigs[2][key], key === 'best' || key === 'reps' || key === 'hits' ? [sigs[0][key], sigs[2][key]] : undefined);
    }
    check(`${c.file}: single-thread driver gives the same score`, single.bestScore === sigs[0].best, [single.bestScore, sigs[0].best]);
    check(`${c.file}: single-thread driver uses the same replicates`, single.reps === sigs[0].reps && single.hits === sigs[0].hits, [single.reps, single.hits]);
    check(`${c.file}: single-thread driver gives the same trees`, single.trees.map((t) => E.treeKey(t)).sort().join('#') === sigs[0].trees);
    check(`${c.file}: every tree re-scores to the reported score`, single.trees.every((t) => sc.length(t) === single.bestScore));
    if (c.constraint) {
      const w = new Uint32Array(Math.ceil(mx.taxa.length / 32)); for (const i of c.constraint.pos[0]) w[i >> 5] |= 1 << (i & 31);
      const key = Array.from(w).join('.');
      check(`${c.file}: every tree satisfies the constraint`, single.trees.every((t) => E.clades(t).includes(key)));
    }
    // a different seed is a different search, but must still report a consistent result
    const other = E.analyse(data, Object.assign({}, cfg.opts, { seed: 12, constraint: c.constraint }));
    check(`${c.file}: second seed agrees on the best score`, other.bestScore === single.bestScore, [other.bestScore, single.bestScore]);
  }
  // the stopping rule
  check('repsNeeded: best found 3 times', E.repsNeeded([5, 4, 4, 6, 4, 4], 3, 20) === 5);
  check('repsNeeded: a later, better score restarts the count', E.repsNeeded([5, 5, 5], 3, 20) === 3 && E.repsNeeded([5, 5, 4, 5, 4, 4], 3, 20) === 6);
  check('repsNeeded: waits for missing replicates', E.repsNeeded([5, undefined, 5, 5], 3, 20) === 0);
  check('repsNeeded: a minimum number of replicates', E.repsNeeded([5, 5, 5, 5, 5], 5, 50, 10) === 0 && E.repsNeeded([5, 5, 5, 5, 5, 5, 6, 6, 6, 6], 5, 50, 10) === 10 && E.repsNeeded([5, 5, 5, 5, 5, 5, 6, 4, 6, 6], 5, 50, 10) === 0);
  check('repsNeeded: fixed number of replicates', E.repsNeeded([5, 4, 4], 0, 3) === 3 && E.repsNeeded([5, 4], 0, 3) === 0);
  console.log(`${pass}/${pass + fail} reproducibility checks passed`);
  process.exit(fail ? 1 : 0);
})();
