// Ordered characters, weights, exclusions and implied weighting, checked against
// independent code: a Sankoff dynamic-programming scorer (unit costs for unordered
// characters, |i-j| costs for ordered ones) and a separate calculation of the
// implied-weights fit. Gaps are read as missing data so that every character
// has a well-defined Sankoff length.
const fs = require('fs');
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const D = process.env.DATASETS;
const NTREES = parseInt(process.env.NTREES || '12', 10);
let pass = 0, fail = 0;
const check = (n, c, x) => { if (c) pass++; else { fail++; if (fail < 12) console.log('  FAIL', n, x === undefined ? '' : JSON.stringify(x)); } };

// Steps of each character on a tree by Sankoff DP. Ordered: a cell with several
// states is the range from its lowest to its highest (as OpenClade and TNT read it).
function sankoffSteps(mx, tree, ordered) {
  const n = mx.taxa.length, S = mx.symbols.length;
  const sym = new Map(mx.symbols.map((s, i) => [s, i]));
  const post = E.postorder(tree, []);
  const out = new Int32Array(mx.nchar), cost = new Array(2 * n - 1);
  for (let c = 0; c < mx.nchar; c++) {
    const ord = ordered && ordered[c];
    for (let t = 0; t < n; t++) {
      const cell = mx.cells[t][c], v = new Float64Array(S);
      if (cell !== '?' && cell !== '-') {
        v.fill(Infinity);
        const idx = cell.map((s) => sym.get(s));
        if (ord) { for (let s = Math.min(...idx); s <= Math.max(...idx); s++) v[s] = 0; } else for (const s of idx) v[s] = 0;
      }
      cost[t] = v;
    }
    for (const v of post) {
      const a = cost[tree.left[v]], b = cost[tree.right[v]], o = new Float64Array(S);
      for (let s = 0; s < S; s++) {
        let ma = Infinity, mb = Infinity;
        for (let q = 0; q < S; q++) { const d = ord ? Math.abs(s - q) : (s === q ? 0 : 1); if (a[q] + d < ma) ma = a[q] + d; if (b[q] + d < mb) mb = b[q] + d; }
        o[s] = ma + mb;
      }
      cost[v] = o;
    }
    out[c] = Math.min(...cost[tree.root]);
  }
  return out;
}
// Fewest steps a character can have on any tree, by brute force over sets of states.
function minSteps(mx, c, ordered) {
  const sym = new Map(mx.symbols.map((s, i) => [s, i]));
  const cells = [];
  for (let t = 0; t < mx.taxa.length; t++) { const x = mx.cells[t][c]; if (x !== '?' && x !== '-') cells.push(x.map((s) => sym.get(s))); }
  if (!cells.length) return 0;
  if (ordered) return Math.max(0, Math.max(...cells.map((x) => Math.min(...x))) - Math.min(...cells.map((x) => Math.max(...x))));
  const states = [...new Set(cells.flat())];
  for (let k = 1; k <= states.length; k++) { // smallest set of states that touches every cell
    const pick = [];
    const go = (i) => {
      if (pick.length === k) return cells.every((x) => x.some((s) => pick.includes(s)));
      for (let j = i; j < states.length; j++) { pick.push(states[j]); if (go(j + 1)) return true; pick.pop(); }
      return false;
    };
    if (go(0)) return k - 1;
  }
  return states.length - 1;
}

const rng = new E.RNG(77);
const files = fs.readdirSync(D).filter((x) => x.endsWith('.nex')).sort();
let trees = 0, orderedChars = 0;
for (const f of files) {
  const mx = P.parseMatrix(fs.readFileSync(D + '/' + f, 'utf8'));
  const n = mx.taxa.length, nc = mx.nchar;
  // random settings: about a third of the multistate characters ordered, weights 0-3 and a few 0.5, some exclusions
  const nstates = (c) => { const s = new Set(); for (let t = 0; t < n; t++) { const x = mx.cells[t][c]; if (Array.isArray(x)) x.forEach((q) => s.add(q)); } return s.size; };
  mx.ordered = Array.from({ length: nc }, (_, c) => nstates(c) > 2 && rng.next() < 0.6);
  mx.weights = Array.from({ length: nc }, () => [1, 1, 1, 2, 3, 0.5][rng.int(6)]);
  mx.excluded = Array.from({ length: nc }, () => rng.next() < 0.1);
  orderedChars += mx.ordered.filter(Boolean).length;
  const data = E.encode(mx, 'missing');
  const sc = new E.Scorer(data);
  const mins = Array.from({ length: nc }, (_, c) => minSteps(mx, c, mx.ordered[c]));
  const st = E.charStats(data);
  let minOk = true;
  for (let c = 0; c < nc; c++) if (data.charCols[c].length && st.minC[c] !== mins[c]) minOk = false;
  check(f + ': minimum steps per character', minOk);
  const k = [3, 12, 7.5][rng.int(3)];
  const dIW = E.encode(mx, 'missing', { concavity: k }), scIW = new E.Scorer(dIW);
  const dPlain = E.encode(mx, 'missing', { applyOrdered: false, applyWeights: false, applyExcluded: false }), scPlain = new E.Scorer(dPlain);
  for (let i = 0; i < NTREES; i++) {
    const t = E.rootOn(E.randomTree(n, rng), rng.int(n));
    trees++;
    const ref = sankoffSteps(mx, t, mx.ordered);
    const got = sc.charSteps(t);
    let same = true, len = 0, fit = 0;
    for (let c = 0; c < nc; c++) {
      if (mx.excluded[c] || mx.weights[c] === 0) { if (got[c] !== 0) same = false; continue; }
      if (got[c] !== ref[c]) same = false;
      len += mx.weights[c] * ref[c];
      const e = ref[c] - mins[c];
      fit += mx.weights[c] * e / (e + k);
    }
    check(f + ': steps per character (ordered + unordered)', same);
    check(f + ': weighted length', Math.abs(sc.length(t) / data.scoreDiv - len) < 1e-9, [sc.length(t) / data.scoreDiv, len]);
    check(f + ': implied-weights fit, k=' + k, Math.abs(scIW.length(t) / dIW.scoreDiv - fit) < 1e-5 * nc, [scIW.length(t) / dIW.scoreDiv, fit]);
    // with the file's settings switched off the matrix is scored as before
    const plainRef = sankoffSteps(mx, t, null);
    let plain = 0; for (let c = 0; c < nc; c++) plain += plainRef[c];
    check(f + ': settings switched off', scPlain.length(t) === plain, [scPlain.length(t), plain]);
  }
}
console.log(`${trees} random trees on ${files.length} matrices, ${orderedChars} characters ordered at random, random weights and exclusions`);
console.log(`${pass}/${pass + fail} checks: engine == independent Sankoff / implied-weights calculation`);
process.exit(fail ? 1 : 0);
