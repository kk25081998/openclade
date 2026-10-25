/* OpenClade — analysis worker. Runs inside a Web Worker (browser) or a
 * worker_threads Worker (Node tests) after parser.js and engine.js.
 * Messages in:  {cmd:'search'|'finish'|'resample', ...}
 * Messages out: progress, searchResult, trees, sprogress, support, error.
 * License: GPL-3.0-or-later */
(function (root) {
  'use strict';
  const E = root.OCEngine;
  function handle(msg, post) {
    try {
      const data = E.encode(msg.mx, msg.gapMode);
      if (msg.cmd === 'search') {
        let last = 0;
        const res = E.analyse(data, Object.assign({}, msg.opts, {
          searchOnly: true,
          onProgress: function (p) {
            const now = Date.now();
            if (now - last < 400) return;
            last = now;
            const snap = p.snapshot();
            post({ type: 'progress', p: { phase: p.phase, replicate: p.replicate, iteration: p.iteration, best: p.best, trees: p.trees, fraction: p.fraction, elapsed: p.elapsed },
              snap: { best: snap.best, trees: snap.trees.slice(0, 100).map((t) => E.toNested(t)), repBests: snap.repBests.map((t) => E.toNested(t)) } });
          },
        }));
        post({ type: 'searchResult', best: res.bestScore, trees: res.trees.map((t) => E.toNested(t)), repBests: res.repBests.map((t) => E.toNested(t)), evaluations: res.evaluations, timedOut: res.timedOut });
      } else if (msg.cmd === 'finish') {
        let last = 0;
        const res = E.mergeAndFinish(data, msg.results, Object.assign({}, msg.opts, {
          onProgress: function (p) { const now = Date.now(); if (now - last > 400) { last = now; post({ type: 'progress', p: { phase: p.phase, best: p.best, trees: p.trees, fraction: 1 } }); } },
        }));
        post({ type: 'trees', best: res.bestScore, trees: res.trees.map((t) => E.toNested(t)), evaluations: res.evaluations, timedOut: res.timedOut });
      } else if (msg.cmd === 'resample') {
        let last = 0;
        const r = E.resample(data, { method: msg.method, replicates: msg.replicates, seed: msg.seed, mode: 'spr', timeLimit: msg.timeLimit,
          onProgress: function (p) { const now = Date.now(); if (now - last > 400) { last = now; post({ type: 'sprogress', p: p }); } } });
        const counts = [];
        for (const [k, f] of r.freq) counts.push([k, Math.round(f * r.replicates)]);
        post({ type: 'support', counts: counts, replicates: r.replicates });
      }
    } catch (e) {
      post({ type: 'error', message: String(e && e.message || e) });
    }
  }
  root.OCWorkerHandle = handle;
  if (typeof root.postMessage === 'function' && typeof root.addEventListener === 'function' && typeof window === 'undefined') {
    root.onmessage = function (ev) { handle(ev.data, function (m) { root.postMessage(m); }); };
  }
})(typeof self !== 'undefined' ? self : globalThis);
