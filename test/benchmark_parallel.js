// Full-size benchmark of the app's actual pipeline (parallel workers + fusing
// + collection) on every published matrix, two seeds, fixed time budget.
// Compare with test/benchmark_baseline.log (single-threaded v0.1 engine).
const fs = require('fs');
const path = require('path');
const { Worker } = require('worker_threads');
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const { runParallel } = require('../src/parallel.js');
const D = process.env.DATASETS;
const LIMIT = parseInt(process.env.LIMIT || '25000', 10);
const CORES = parseInt(process.env.CORES || '2', 10);
const workerCode = `
const { parentPort } = require('worker_threads');
globalThis.OCEngine = require(${JSON.stringify(path.join(__dirname, '../src/engine.js'))});
require(${JSON.stringify(path.join(__dirname, '../src/worker.js'))});
parentPort.on('message', (m) => globalThis.OCWorkerHandle(m, (x) => parentPort.postMessage(x)));`;
function spawn() {
  const w = new Worker(workerCode, { eval: true });
  const wrap = { onmessage: null, postMessage: (m) => w.postMessage(m), terminate: () => w.terminate() };
  w.on('message', (m) => wrap.onmessage && wrap.onmessage(m));
  return wrap;
}
function run(mx, seed) {
  return new Promise((resolve) => {
    const t0 = Date.now(); let res = null;
    runParallel({ engine: E, spawn, nWorkers: CORES, mx: { taxa: mx.taxa, cells: mx.cells, nchar: mx.nchar, symbols: mx.symbols }, gapMode: 'inapplicable',
      opts: { replicates: 6, ratchetIter: 10, maxTrees: 100, seed, timeLimit: LIMIT, mode: 'tbr' }, support: { method: 'none' },
      onTrees: (r) => { res = r; }, onDone: () => resolve(Object.assign(res, { ms: Date.now() - t0 })), onError: (e) => console.log('error', e) });
  });
}
(async () => {
  const out = [];
  for (const f of fs.readdirSync(D).filter((x) => x.endsWith('.nex')).sort()) {
    const mx = P.parseMatrix(fs.readFileSync(D + '/' + f, 'utf8'));
    const sc = new E.Scorer(E.encode(mx, 'inapplicable'));
    const runs = [];
    for (const seed of [1, 2]) {
      const r = await run(mx, seed);
      const ok = r.trees.every((nt) => sc.length(E.fromNested(mx.taxa.length, nt)) === r.best);
      runs.push({ seed, best: r.best, trees: r.trees.length, ms: r.ms, rescoreOk: ok });
    }
    out.push({ dataset: f.replace('.nex', ''), runs });
    const [a, b] = runs;
    console.log(f.replace('.nex', '').padEnd(14), String(mx.taxa.length).padStart(3), 'taxa |', 'seed1', a.best, `(${a.trees})`, '| seed2', b.best, `(${b.trees})`, a.best === b.best ? 'AGREE' : 'DIFFER', a.rescoreOk && b.rescoreOk ? '' : 'RESCORE-FAIL');
  }
  fs.writeFileSync(process.env.OUT || 'benchmark_parallel.json', JSON.stringify(out, null, 1));
  process.exit(0);
})();
