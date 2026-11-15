/* OpenClade — analysis worker. Runs inside a Web Worker (browser) or a
 * worker_threads Worker (Node tests) after parser.js and engine.js.
 * Messages in:  {cmd:'init', mx, gapMode, charOpts, constraint}, then
 *               {cmd:'rep', index, opts, timeLeft} | {cmd:'finish', results, opts, timeLeft}
 *               | {cmd:'resample', indices, method, seed, timeLeft}
 *               | {cmd:'bremer', clades, indices, startTrees, reps, ratchetIter, seed, timeLeft}
 * Messages out: progress, rep, trees, sprogress, support, bprogress, bremer, error.
 * License: GPL-3.0-or-later */
(function (root) {
  'use strict';
  const E = root.OCEngine;
  let st = null; // { data, S, constraint }
  function handle(msg, post) {
    try {
      if (msg.cmd === 'init') {
        const data = E.encode(msg.mx, msg.gapMode, msg.charOpts);
        const S = new E.Searcher(data, { seed: 1 });
        if (msg.constraint) S.setConstraint(msg.constraint);
        st = { data, S, constraint: msg.constraint || null };
        return;
      }
      if (!st) throw new Error('worker used before init');
      const { data, S } = st;
      const deadline = Date.now() + (msg.timeLeft || 60000);
      if (msg.cmd === 'rep') {
        let last = 0;
        S.onImprove = function (t, sc) { // lets an early stop return the best tree so far
          const now = Date.now();
          if (now - last < 500) return;
          last = now;
          post({ type: 'progress', index: msg.index, p: { phase: 'climb', best: sc }, snap: { score: sc, tree: E.toNested(E.rootOn(t, 0)) } });
        };
        const r = E.replicate(S, msg.index, msg.opts, deadline, function (p) {
          post({ type: 'progress', index: msg.index, p: { phase: p.phase, iteration: p.iteration, best: p.best } });
        });
        S.onImprove = null;
        post({ type: 'rep', index: r.index, score: r.score, trees: r.trees.map((t) => E.toNested(t)), timedOut: r.timedOut, evaluations: S.evaluations });
      } else if (msg.cmd === 'finish') {
        let last = 0;
        const fin = E.finishTrees(S, msg.results, msg.opts, deadline, function (p) {
          const now = Date.now();
          if (now - last > 400) { last = now; post({ type: 'progress', p: { phase: p.phase, best: p.best, trees: p.trees } }); }
        });
        post({ type: 'trees', best: fin.bestScore, trees: fin.trees.map((t) => E.toNested(t)), hits: fin.hits, timedOut: fin.timedOut });
      } else if (msg.cmd === 'resample') {
        const counts = new Map();
        let done = 0, last = 0;
        const ropts = { method: msg.method, seed: msg.seed, mode: 'spr' };
        if (st.constraint && st.constraint.pos) S.setConstraint({ pos: st.constraint.pos });
        for (const i of msg.indices) {
          if (Date.now() > deadline) break;
          for (const c of E.resampleRep(S, i, ropts, deadline)) counts.set(c, (counts.get(c) || 0) + 1);
          done++;
          const now = Date.now();
          if (now - last > 400) { last = now; post({ type: 'sprogress', p: { replicate: done, of: msg.indices.length } }); }
        }
        if (st.constraint) S.setConstraint(st.constraint);
        post({ type: 'support', counts: [...counts], replicates: done, asked: msg.indices.length });
      } else if (msg.cmd === 'bremer') {
        let last = 0;
        const r = E.bremer(data, {
          clades: msg.clades, indices: msg.indices, startTrees: msg.startTrees, reps: msg.reps, ratchetIter: msg.ratchetIter,
          seed: msg.seed, timeLimit: msg.timeLeft, constraint: st.constraint,
          onProgress: function (p) { const now = Date.now(); if (now - last > 300) { last = now; post({ type: 'bprogress', p: { replicate: p.replicate, of: p.of } }); } },
        });
        post({ type: 'bremer', bounds: r.bounds.map(([k, v]) => [k, isFinite(v) ? v : null]), done: r.done, skipped: r.skipped, timedOut: r.timedOut });
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
