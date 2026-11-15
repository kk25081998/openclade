/* OpenClade — parallel analysis orchestrator (environment-agnostic).
 * spawn() must return {postMessage(msg), onmessage: fn(msg), terminate()}.
 *
 * Phase 1: search replicates 0, 1, 2, ... are handed to workers one at a time.
 *          Replicate i depends only on (seed, i). The search uses replicates
 *          0..m-1, where m is the first count at which the best score has been
 *          found `hits` times (engine.repsNeeded); later replicates that happen
 *          to be running are discarded. The result is therefore the same on
 *          one core or sixteen, and on a slow or a fast computer.
 * Phase 2: one worker fuses the replicate results and collects equal trees.
 * Phase 3: optional bootstrap/jackknife; replicate i goes to worker i mod W.
 * Phase 4: optional Bremer support; clade i goes to worker i mod W.
 * The time limit only cuts a phase short; results are then flagged timedOut.
 * License: GPL-3.0-or-later */
(function (root) {
  'use strict';
  function runParallel(cfg) {
    const E = cfg.engine;
    const opts = E.searchOpts(cfg.opts);
    const W = Math.max(1, Math.min(cfg.nWorkers || 1, opts.maxReps));
    const T = opts.timeLimit;
    const t0 = Date.now();
    const searchDeadline = t0 + 0.8 * T;
    const workers = [];
    let stopped = false, phase = 'search';
    const emit = (f, o) => { if (f) f(o); };
    const init = { cmd: 'init', mx: cfg.mx, gapMode: cfg.gapMode, charOpts: cfg.charOpts || null, constraint: cfg.constraint || null };
    function spawn() { const w = cfg.spawn(); workers.push(w); w.postMessage(init); return w; }
    function kill(w) { try { w.terminate(); } catch (e) { /* ignore */ } const i = workers.indexOf(w); if (i >= 0) workers.splice(i, 1); }
    function killAll() { for (const w of workers.slice()) kill(w); }
    const fail = (m) => { if (!stopped) { stopped = true; killAll(); emit(cfg.onError, m); } };

    /* ---- phase 1 ---- */
    // implied weights: replicate scores this close are the same fit (see engine setImplied)
    const tieTol = cfg.charOpts && cfg.charOpts.concavity ? E.encode(cfg.mx, cfg.gapMode, cfg.charOpts).tieTol : 0;
    const results = [];          // by replicate index
    const scores = [];           // scores[i] undefined until replicate i is finished
    const snaps = new Map();     // running replicates: best tree so far
    let next = 0, running = 0, cut = false, finalTrees = null, evaluations = 0;
    const repOpts = { seed: opts.seed, ratchetIter: opts.ratchetIter, mode: opts.mode, drift: opts.drift, keepPerRep: opts.keepPerRep };
    function report(extra) {
      let best = Infinity, done = 0, atBest = 0;
      for (const s of scores) if (s !== undefined) { done++; if (s < best) best = s; }
      for (const s of scores) if (s !== undefined && s <= best + tieTol) atBest++;
      for (const s of snaps.values()) if (s.score < best) best = s.score;
      const frac = opts.hits ? Math.min(1, Math.max(Math.min(atBest / opts.hits, opts.minReps ? done / opts.minReps : 1), done / opts.maxReps)) : done / opts.maxReps;
      emit(cfg.onProgress, Object.assign({ phase: 'search', best, hits: atBest, reps: done, workers: W, fraction: 0.85 * frac, elapsed: Date.now() - t0 }, extra || {}));
    }
    function dispatch(w) {
      if (stopped) return;
      const used = E.repsNeeded(scores, opts.hits, opts.maxReps, opts.minReps, tieTol);
      if (used) { finishPhase(used, false); return; }
      const timeLeft = searchDeadline - Date.now();
      if (next >= opts.maxReps || timeLeft <= 0) {
        kill(w);
        if (running === 0) finishPhase(0, true); // out of time before the stopping rule was met
        return;
      }
      const index = next++;
      running++;
      w.onmessage = (m) => {
        if (stopped || phase !== 'search') return;
        if (m.type === 'progress') { if (m.snap) snaps.set(index, m.snap); report(); }
        else if (m.type === 'rep') {
          running--; snaps.delete(index);
          results[index] = { index, score: m.score, trees: m.trees, timedOut: m.timedOut };
          scores[index] = m.score;
          if (m.timedOut) cut = true;
          evaluations = Math.max(evaluations, m.evaluations || 0);
          report();
          dispatch(w);
        } else if (m.type === 'error') fail(m.message);
      };
      w.postMessage({ cmd: 'rep', index, opts: repOpts, timeLeft });
    }
    function finishPhase(used, outOfTime) {
      if (phase !== 'search') return;
      phase = 'finish';
      killAll();
      const use = used ? results.slice(0, used) : results.filter(Boolean);
      const timedOut = outOfTime || cut && use.some((r) => r.timedOut);
      if (!use.length) { fail('The time limit ran out before the first tree was built. Raise the time limit.'); return; }
      const w = spawn();
      const left = Math.max(5000, T - (Date.now() - t0));
      emit(cfg.onProgress, { phase: 'fuse', best: Math.min(...use.map((r) => r.score)), workers: W, fraction: 0.85, elapsed: Date.now() - t0 });
      w.onmessage = (m) => {
        if (stopped) return;
        if (m.type === 'progress') emit(cfg.onProgress, { phase: m.p.phase, best: m.p.best, workers: W, fraction: 0.9, elapsed: Date.now() - t0 });
        else if (m.type === 'trees') {
          kill(w);
          finalTrees = { best: m.best, trees: m.trees };
          emit(cfg.onTrees, { best: m.best, trees: m.trees, hits: m.hits, reps: use.length, repScores: use.map((r) => r.score),
            timedOut: timedOut || m.timedOut, elapsed: Date.now() - t0, workers: W });
          afterTrees();
        } else if (m.type === 'error') fail(m.message);
      };
      w.postMessage({ cmd: 'finish', results: use, opts: { maxTrees: opts.maxTrees, mode: opts.mode, seed: opts.seed, fuse: opts.fuse }, timeLeft: left });
    }

    /* ---- phases 3 and 4 ---- */
    function afterTrees() {
      if (cfg.support && cfg.support.method && cfg.support.method !== 'none') supportPhase();
      else if (cfg.bremer && cfg.bremer.enabled) bremerPhase();
      else done();
    }
    function supportPhase() {
      phase = 'support';
      const S = cfg.support;
      const Ws = Math.max(1, Math.min(cfg.nWorkers || 1, S.replicates));
      const counts = new Map(); let total = 0, asked = 0, left = Ws;
      const prog = new Array(Ws).fill(0);
      for (let i = 0; i < Ws; i++) {
        const w = spawn();
        const indices = []; for (let r = i; r < S.replicates; r += Ws) indices.push(r);
        w.onmessage = (m) => {
          if (stopped) return;
          if (m.type === 'sprogress') { prog[i] = m.p.replicate; const d = prog.reduce((a, b) => a + b, 0); emit(cfg.onProgress, { phase: S.method, replicate: d, of: S.replicates, fraction: d / S.replicates, workers: Ws }); }
          else if (m.type === 'support') {
            for (const [k, c] of m.counts) counts.set(k, (counts.get(k) || 0) + c);
            total += m.replicates; asked += m.asked; kill(w);
            if (--left === 0) {
              const freq = []; for (const [k, c] of counts) freq.push([k, c / Math.max(1, total)]);
              emit(cfg.onSupport, { freq, replicates: total, method: S.method, timedOut: total < asked });
              if (cfg.bremer && cfg.bremer.enabled) bremerPhase(); else done();
            }
          } else if (m.type === 'error') fail(m.message);
        };
        w.postMessage({ cmd: 'resample', indices, method: S.method, seed: opts.seed, timeLeft: T });
      }
    }
    function bremerPhase() {
      phase = 'bremer';
      const n = cfg.mx.taxa.length, B = cfg.bremer;
      const trees = finalTrees.trees.map((nt) => E.fromNested(n, nt));
      // clades of the strict consensus, largest first so the order is fixed
      const cnt = new Map();
      for (const t of trees) for (const c of new Set(E.clades(t))) cnt.set(c, (cnt.get(c) || 0) + 1);
      const clades = [...cnt].filter(([, v]) => v === trees.length).map(([k]) => k).sort();
      if (!clades.length) { emit(cfg.onBremer, { values: [], best: finalTrees.best, timedOut: false, skipped: [] }); done(); return; }
      const Wb = Math.max(1, Math.min(cfg.nWorkers || 1, clades.length));
      const bound = new Map(clades.map((k) => [k, Infinity]));
      const skipped = new Set(), finished = new Set(); let left = Wb, timedOut = false;
      const prog = new Array(Wb).fill(0);
      for (let i = 0; i < Wb; i++) {
        const w = spawn();
        const indices = []; for (let r = i; r < clades.length; r += Wb) indices.push(r);
        w.onmessage = (m) => {
          if (stopped) return;
          if (m.type === 'bprogress') { prog[i] = m.p.replicate; const d = prog.reduce((a, b) => a + b, 0); emit(cfg.onProgress, { phase: 'bremer', replicate: d, of: clades.length, fraction: d / clades.length, workers: Wb }); }
          else if (m.type === 'bremer') {
            for (const [k, v] of m.bounds) if (v !== null && v < bound.get(k)) bound.set(k, v);
            for (const ci of m.skipped) skipped.add(clades[ci]);
            for (const ci of m.done) finished.add(clades[ci]);
            if (m.timedOut) timedOut = true;
            kill(w);
            if (--left === 0) {
              // value = extra cost of the best tree found without the clade; searched = its own search ran
              const values = clades.map((k) => ({ key: k, bound: isFinite(bound.get(k)) ? bound.get(k) : null, searched: finished.has(k), enforced: skipped.has(k) }));
              let lowest = Infinity; for (const v of bound.values()) if (v < lowest) lowest = v;
              emit(cfg.onBremer, { values, best: finalTrees.best, lowest: isFinite(lowest) ? lowest : null, timedOut });
              done();
            }
          } else if (m.type === 'error') fail(m.message);
        };
        w.postMessage({ cmd: 'bremer', clades, indices, startTrees: finalTrees.trees.slice(0, 8), reps: B.reps || 2, ratchetIter: B.ratchetIter || 0, seed: opts.seed, timeLeft: B.timeLimit || T });
      }
    }
    function done() { phase = 'done'; killAll(); emit(cfg.onDone, {}); }

    // Best trees so far (used when the user stops the search early)
    function snapshot() {
      if (finalTrees) return finalTrees;
      let best = Infinity;
      for (const r of results) if (r && r.score < best) best = r.score;
      for (const s of snaps.values()) if (s.score < best) best = s.score;
      const seen = new Set(), trees = [];
      const add = (nt) => { const k = E.treeKey(E.rootOn(E.fromNested(cfg.mx.taxa.length, nt), 0)); if (!seen.has(k)) { seen.add(k); trees.push(nt); } };
      for (const r of results) if (r && r.score === best) r.trees.forEach(add);
      for (const s of snaps.values()) if (s.score === best) add(s.tree);
      return { best, trees };
    }

    for (let i = 0; i < W; i++) dispatch(spawn());

    return {
      stop: function () {
        const snap = snapshot();
        stopped = true;
        killAll();
        return snap;
      },
      nWorkers: W,
    };
  }
  const api = { runParallel };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.OCParallel = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);
