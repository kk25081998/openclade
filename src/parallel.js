/* OpenClade — parallel analysis orchestrator (environment-agnostic).
 * spawn() must return {postMessage(msg), onmessage: fn(msg), terminate()}.
 * Phase 1: replicates split across workers (independent seeds).
 * Phase 2: results merged; one worker fuses and collects equal trees.
 * Phase 3: optional bootstrap/jackknife, replicates split across workers.
 * License: GPL-3.0-or-later */
(function (root) {
  'use strict';
  function runParallel(cfg) {
    const E = cfg.engine;
    const W = Math.max(1, Math.min(cfg.nWorkers || 1, cfg.opts.replicates));
    const T = cfg.opts.timeLimit;
    const workers = [];
    let stopped = false;
    const snaps = new Array(W).fill(null);
    const fractions = new Array(W).fill(0);
    const t0 = Date.now();
    const emit = (f, o) => { if (f) f(o); };
    function spawn(i) { const w = cfg.spawn(); workers.push(w); return w; }
    function killAll() { for (const w of workers) { try { w.terminate(); } catch (e) { /* ignore */ } } workers.length = 0; }
    function shares(total) { const out = []; for (let i = 0; i < W; i++) out.push(Math.floor(total / W) + (i < total % W ? 1 : 0)); return out; }

    // Merge worker snapshots in the main thread (used when stopped early)
    function mergeSnapshots() {
      let best = Infinity;
      for (const s of snaps) if (s && s.best < best) best = s.best;
      const seen = new Set(), trees = [];
      for (const s of snaps) if (s && s.best === best) for (const nt of s.trees) {
        const t = E.fromNested(cfg.mx.taxa.length, nt);
        const k = E.treeKey(t);
        if (!seen.has(k)) { seen.add(k); trees.push(nt); }
      }
      return { best, trees };
    }

    const reps = shares(cfg.opts.replicates);
    const results = new Array(W).fill(null);
    let remaining = W;
    for (let i = 0; i < W; i++) {
      const w = spawn(i);
      w.onmessage = (m) => {
        if (stopped) return;
        if (m.type === 'progress') {
          snaps[i] = m.snap || snaps[i];
          fractions[i] = m.p.fraction || 0;
          let best = Infinity; for (const s of snaps) if (s && s.best < best) best = s.best;
          emit(cfg.onProgress, { phase: 'search', best, workers: W, fraction: 0.85 * fractions.reduce((a, b) => a + b, 0) / W, elapsed: Date.now() - t0, detail: m.p });
        } else if (m.type === 'searchResult') {
          results[i] = { trees: m.trees, repBests: m.repBests, best: m.best };
          snaps[i] = { best: m.best, trees: m.trees, repBests: m.repBests };
          w.terminate();
          if (--remaining === 0) finishPhase();
        } else if (m.type === 'error') { emit(cfg.onError, m.message); }
      };
      w.postMessage({ cmd: 'search', mx: cfg.mx, gapMode: cfg.gapMode,
        opts: Object.assign({}, cfg.opts, { replicates: reps[i], seed: cfg.opts.seed * 7919 + i, timeLimit: T * 0.85 }) });
    }

    function finishPhase() {
      workers.length = 0;
      const w = spawn(0);
      const left = Math.max(5000, T - (Date.now() - t0));
      emit(cfg.onProgress, { phase: 'fuse', best: Math.min(...results.map((r) => r.best)), workers: W, fraction: 0.85, elapsed: Date.now() - t0 });
      w.onmessage = (m) => {
        if (stopped) return;
        if (m.type === 'progress') emit(cfg.onProgress, { phase: m.p.phase, best: m.p.best, workers: W, fraction: 0.9, elapsed: Date.now() - t0 });
        else if (m.type === 'trees') {
          w.terminate();
          emit(cfg.onTrees, { best: m.best, trees: m.trees, timedOut: m.timedOut, elapsed: Date.now() - t0, workers: W });
          if (cfg.support && cfg.support.method !== 'none') supportPhase(); else { done(); }
        } else if (m.type === 'error') emit(cfg.onError, m.message);
      };
      w.postMessage({ cmd: 'finish', mx: cfg.mx, gapMode: cfg.gapMode, results, opts: { maxTrees: cfg.opts.maxTrees, mode: cfg.opts.mode, seed: cfg.opts.seed, timeLimit: left } });
    }

    function supportPhase() {
      workers.length = 0;
      const S = cfg.support;
      const Ws = Math.max(1, Math.min(cfg.nWorkers || 1, S.replicates));
      let n = 0; const sh = []; for (let i = 0; i < Ws; i++) sh.push(Math.floor(S.replicates / Ws) + (i < S.replicates % Ws ? 1 : 0));
      const counts = new Map(); let total = 0, left = Ws;
      const prog = new Array(Ws).fill(0);
      for (let i = 0; i < Ws; i++) {
        const w = spawn(i);
        w.onmessage = (m) => {
          if (stopped) return;
          if (m.type === 'sprogress') { prog[i] = m.p.replicate; emit(cfg.onProgress, { phase: S.method, replicate: prog.reduce((a, b) => a + b, 0), of: S.replicates, fraction: prog.reduce((a, b) => a + b, 0) / S.replicates, workers: Ws }); }
          else if (m.type === 'support') {
            for (const [k, c] of m.counts) counts.set(k, (counts.get(k) || 0) + c);
            total += m.replicates; w.terminate();
            if (--left === 0) {
              const freq = []; for (const [k, c] of counts) freq.push([k, c / Math.max(1, total)]);
              emit(cfg.onSupport, { freq, replicates: total, method: S.method });
              done();
            }
          } else if (m.type === 'error') emit(cfg.onError, m.message);
        };
        w.postMessage({ cmd: 'resample', mx: cfg.mx, gapMode: cfg.gapMode, method: S.method, replicates: sh[i], seed: cfg.opts.seed * 104729 + i + 1, timeLimit: T });
        n++;
      }
      void n;
    }
    function done() { killAll(); emit(cfg.onDone, {}); }

    return {
      stop: function () {
        stopped = true;
        killAll();
        return mergeSnapshots();
      },
      nWorkers: W,
    };
  }
  const api = { runParallel };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.OCParallel = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);
