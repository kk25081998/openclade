// Character changes on branches and collapsing of unsupported branches, checked
// by brute force. For small random trees and characters every assignment of
// states to the internal nodes is enumerated, giving all most parsimonious
// reconstructions (MPRs). Then:
//   - "unambiguous" changes  == branches that change in every MPR;
//   - ACCTRAN and DELTRAN    are MPRs (their number of changes is the minimum),
//                              and of all MPRs ACCTRAN has its changes nearest the root
//                              and DELTRAN nearest the tips (sum over changes of the
//                              depth of the branch);
//   - a collapsed branch     == a branch that has, for every character, an MPR
//                              without a change on it.
const P = require('../src/parser.js');
const E = require('../src/engine.js');
let pass = 0, fail = 0;
const check = (n, c, x) => { if (c) pass++; else { fail++; if (fail < 12) console.log('  FAIL', n, x === undefined ? '' : JSON.stringify(x)); } };
const rng = new E.RNG(909);
const N = 7, NCHAR = 6, CASES = parseInt(process.env.CASES || '400', 10);

function randomMatrix(ordered) {
  const symbols = ['0', '1', '2', '3'];
  const cells = [];
  for (let t = 0; t < N; t++) {
    const row = [];
    for (let c = 0; c < NCHAR; c++) {
      const r = rng.next();
      if (r < 0.08) row.push('?');
      else if (r < 0.2) { const a = rng.int(3), b = a + 1 + rng.int(3 - a); row.push([symbols[a], symbols[b]]); }
      else row.push([symbols[rng.int(c % 2 ? 4 : 3)]]);
    }
    cells.push(row);
  }
  return { taxa: Array.from({ length: N }, (_, i) => 'T' + i), nchar: NCHAR, cells, symbols, charLabels: [], stateLabels: [],
    ordered: ordered ? Array.from({ length: NCHAR }, (_, c) => c % 3 === 0) : null, weights: null, excluded: null };
}

// All MPRs of character c on tree t. Returns { min, mprs: [Int8Array state per node] }
function bruteMPR(mx, t, c, ord) {
  const m = 2 * N - 1, S = mx.symbols.length;
  const allowed = [];
  for (let v = 0; v < m; v++) {
    if (v >= N) { allowed.push([0, 1, 2, 3]); continue; }
    const cell = mx.cells[v][c];
    if (cell === '?') allowed.push([0, 1, 2, 3]);
    else { const idx = cell.map((s) => mx.symbols.indexOf(s)); if (ord) { const a = []; for (let s = Math.min(...idx); s <= Math.max(...idx); s++) a.push(s); allowed.push(a); } else allowed.push(idx); }
  }
  const nodes = []; for (let v = 0; v < m; v++) if (v !== t.root) nodes.push(v);
  // dynamic programming for the minimum, then enumerate assignments of internal nodes and
  // let each tip take its best allowed state given its parent (ties: every best state)
  const internal = []; for (let v = N; v < m; v++) internal.push(v);
  const st = new Int8Array(m);
  let min = Infinity; const mprs = [];
  const d = (a, b) => (ord ? Math.abs(a - b) : (a === b ? 0 : 1));
  (function rec(i) {
    if (i === internal.length) {
      let cost = 0;
      for (const v of internal) if (v !== t.root) cost += d(st[t.parent[v]], st[v]);
      // tips
      const opts = [];
      for (let v = 0; v < N; v++) {
        let b = Infinity, list = [];
        for (const s of allowed[v]) { const x = d(st[t.parent[v]], s); if (x < b) { b = x; list = [s]; } else if (x === b) list.push(s); }
        cost += b; opts.push(list);
      }
      if (cost > min) return;
      if (cost < min) { min = cost; mprs.length = 0; }
      (function tips(v) {
        if (v === N) { mprs.push(Int8Array.from(st)); return; }
        for (const s of opts[v]) { st[v] = s; tips(v + 1); }
      })(0);
      return;
    }
    for (let s = 0; s < S; s++) { st[internal[i]] = s; rec(i + 1); }
  })(0);
  return { min, mprs };
}
function depths(t) { const m = 2 * N - 1, dep = new Int32Array(m); const st = [t.root]; while (st.length) { const v = st.pop(); if (v >= N) { for (const c of [t.left[v], t.right[v]]) { dep[c] = dep[v] + 1; st.push(c); } } } return dep; }

let cases = 0, branches = 0;
for (let k = 0; k < CASES; k++) {
  const ordered = k % 2 === 1;
  const mx = randomMatrix(ordered);
  const data = E.encode(mx, 'missing');
  const t = E.randomTree(N, rng);
  const sc = new E.Scorer(data);
  const un = E.branchChanges(data, t, 'unambiguous', sc), ac = E.branchChanges(data, t, 'acctran', sc), de = E.branchChanges(data, t, 'deltran', sc);
  const dep = depths(t), m = 2 * N - 1;
  const steps = sc.charSteps(t);
  const zeroPossible = new Uint8Array(m).fill(1); // branch v: every character has an MPR with no change on it
  cases++;
  for (let c = 0; c < NCHAR; c++) {
    const ord = ordered && c % 3 === 0 && data.charOrdered[c];
    const bm = bruteMPR(mx, t, c, ord);
    check('minimum over all reconstructions == engine steps', bm.min === steps[c], [bm.min, steps[c]]);
    let unOK = true;
    for (let v = 0; v < m; v++) {
      if (v === t.root) continue;
      const p = t.parent[v];
      const always = bm.mprs.every((a) => a[v] !== a[p]), never = bm.mprs.some((a) => a[v] === a[p]);
      if (!never) { /* changes in every MPR */ }
      const listed = un[v].some((x) => x.c === c);
      if (listed !== always) unOK = false;
      if (!never) zeroPossible[v] = 0;
      branches++;
    }
    check('unambiguous changes == changes in every MPR', unOK, { c, ord });
    // ACCTRAN / DELTRAN: cost and position
    const costOf = (list) => { let s = 0; for (let v = 0; v < m; v++) for (const x of list[v]) if (x.c === c) s += ord ? Math.abs(mx.symbols.indexOf(x.from[0]) - mx.symbols.indexOf(x.to[0])) : 1; return s; };
    check('ACCTRAN is a most parsimonious reconstruction', costOf(ac) === bm.min, [costOf(ac), bm.min]);
    check('DELTRAN is a most parsimonious reconstruction', costOf(de) === bm.min, [costOf(de), bm.min]);
    const depthSum = (a) => { let s = 0; for (let v = 0; v < m; v++) if (v !== t.root) { const p = t.parent[v]; s += (ord ? Math.abs(a[v] - a[p]) : (a[v] !== a[p] ? 1 : 0)) * dep[v]; } return s; };
    const listDepth = (list) => { let s = 0; for (let v = 0; v < m; v++) for (const x of list[v]) if (x.c === c) s += (ord ? Math.abs(mx.symbols.indexOf(x.from[0]) - mx.symbols.indexOf(x.to[0])) : 1) * dep[v]; return s; };
    const ds = bm.mprs.map(depthSum);
    check('ACCTRAN: changes as near the root as any MPR', listDepth(ac) === Math.min(...ds), [listDepth(ac), Math.min(...ds)]);
    check('DELTRAN: changes as near the tips as any MPR', listDepth(de) === Math.max(...ds), [listDepth(de), Math.max(...ds)]);
  }
  // collapsing
  const col = E.collapseTree(sc, t);
  const keys = E.nodeKeys(t), kept = new Set(col.keys);
  let colOK = true;
  for (let v = N; v < m; v++) { if (!keys[v] || keys[v][0] === 't') continue; if (kept.has(keys[v]) === !!zeroPossible[v]) colOK = false; }
  check('collapsed branches == branches with a zero-change MPR for every character', colOK);
}
console.log(`${cases} random trees of ${N} taxa x ${NCHAR} characters (unordered and ordered, with polymorphic and missing cells); ${branches} branch-character pairs`);
console.log(`${pass}/${pass + fail} reconstruction checks passed`);
process.exit(fail ? 1 : 0);
