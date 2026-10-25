// Runs the parallel orchestrator with real worker_threads (as the browser does
// with Web Workers): checks results, support merging and early stop.
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
function run(mx, gapMode, opts, support, nWorkers, stopAfter) {
  return new Promise((resolve) => {
    const out = { progress: 0 };
    const h = runParallel({ engine: E, spawn, nWorkers, mx: { taxa: mx.taxa, cells: mx.cells, nchar: mx.nchar, symbols: mx.symbols }, gapMode, opts, support,
      onProgress: () => out.progress++, onTrees: (r) => { out.trees = r; }, onSupport: (s) => { out.support = s; }, onError: (e) => { out.error = e; },
      onDone: () => resolve(out) });
    if (stopAfter) setTimeout(() => { out.stopped = h.stop(); resolve(out); }, stopAfter);
  });
}

let pass = 0, fail = 0;
const check = (n, c, x) => { if (c) pass++; else { fail++; console.log('  FAIL', n, x === undefined ? '' : JSON.stringify(x)); } };
(async () => {
  const mx = P.parseMatrix(fs.readFileSync(D + '/Vinther2008.nex', 'utf8'));
  const data = E.encode(mx, 'inapplicable');
  const sc = new E.Scorer(data);
  const r = await run(mx, 'inapplicable', { replicates: 6, ratchetIter: 5, maxTrees: 100, seed: 1, timeLimit: 30000, mode: 'tbr' }, { method: 'bootstrap', replicates: 40 }, 2);
  check('no worker error', !r.error, r.error);
  check('parallel search finds 79', r.trees && r.trees.best === 79, r.trees && r.trees.best);
  check('52 MPTs (as single-thread)', r.trees && r.trees.trees.length === 52, r.trees && r.trees.trees.length);
  check('all returned trees re-score to best', r.trees.trees.every((nt) => sc.length(E.fromNested(23, nt)) === 79));
  check('support merged from 40 replicates', r.support && r.support.replicates === 40, r.support && r.support.replicates);
  check('support frequencies in [0,1]', r.support.freq.every(([, f]) => f >= 0 && f <= 1));
  check('progress reported', r.progress > 0);
  // early stop on a big matrix
  const big = P.parseMatrix(fs.readFileSync(D + '/Conrad2008.nex', 'utf8'));
  const s = await run(big, 'inapplicable', { replicates: 4, ratchetIter: 10, maxTrees: 100, seed: 1, timeLimit: 120000, mode: 'tbr' }, { method: 'none' }, 2, 12000);
  check('stop returns best-so-far trees', s.stopped && s.stopped.trees.length > 0 && isFinite(s.stopped.best), s.stopped && s.stopped.best);
  if (s.stopped) {
    const bsc = new E.Scorer(E.encode(big, 'inapplicable'));
    check('stopped trees re-score to reported best', s.stopped.trees.every((nt) => bsc.length(E.fromNested(big.taxa.length, nt)) === s.stopped.best));
  }
  console.log(`${pass}/${pass + fail} parallel checks passed`);
  process.exit(fail ? 1 : 0);
})();
