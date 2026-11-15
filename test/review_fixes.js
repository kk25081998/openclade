// Regression tests for the defects found by the independent review of 0.2.0.
// Each check uses an answer worked out by hand or by enumerating every tree.
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const W = require('../src/worker.js');
let pass = 0, fail = 0;
const check = (n, c, x) => { if (c) pass++; else { fail++; console.log('  FAIL', n, x === undefined ? '' : JSON.stringify(x)); } };
const nex = (rows, extra = '', fmt = '') => `#NEXUS\nbegin data; dimensions ntax=${rows.length} nchar=${rows[0][1].replace(/[{(][^})]*[})]/g, 'x').length}; format ${fmt};\nmatrix\n${rows.map((r) => r.join(' ')).join('\n')}\n;\nend;\n${extra}`;

// every rooted binary tree on n tips (rooted on tip 0)
function allTrees(n, f) {
  const t = E.fromNested(n, [0, [1, 2]]);
  (function rec(kk) {
    if (kk === n) { f(t); return; }
    const w = n + kk - 1, edges = [];
    for (let v = 1; v < 2 * n - 1; v++) if (t.parent[v] >= 0) edges.push(v);
    for (const e of edges) {
      const p = t.parent[e];
      t.parent[w] = p; if (t.left[p] === e) t.left[p] = w; else t.right[p] = w;
      t.left[w] = e; t.right[w] = kk; t.parent[e] = w; t.parent[kk] = w;
      rec(kk + 1);
      if (t.left[p] === w) t.left[p] = e; else t.right[p] = e;
      t.parent[e] = p; t.parent[w] = -1; t.parent[kk] = -1; t.left[w] = -1; t.right[w] = -1;
    }
  })(3);
}
function rng(seed) { let s = seed >>> 0; return () => { s = (s + 0x6D2B79F5) >>> 0; let t = s; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }

// ---- parser: "state or gap" cells, MATCHCHAR, lower-case MISSING/GAP
{
  const m = P.parseMatrix(nex([['A', '0{0-}1'], ['B', '1(1-)1'], ['C', '0-0'], ['D', '1{01}0']]));
  check('{0-} keeps the gap', JSON.stringify(m.cells[0][1]) === '["-","0"]', m.cells[0]);
  check('{0-} is warned about', m.warnings.some((w) => /state or a gap/.test(w)), m.warnings);
  check('the gap is not a symbol', !m.symbols.includes('-'), m.symbols);
  const tree = E.fromNested(4, [0, [1, [2, 3]]]);
  // read as unknown: the character 0 ? 0 1 / ... ; with gaps as an extra state {0,-}
  const lenOf = (mode) => new E.Scorer(E.encode(m, mode)).length(tree);
  const asQ = P.parseMatrix(nex([['A', '0?1'], ['B', '1?1'], ['C', '0-0'], ['D', '1{01}0']]));
  check('{0-} read as unknown with gaps as inapplicable', lenOf('inapplicable') === new E.Scorer(E.encode(asQ, 'inapplicable')).length(tree));
  check('{0-} read as unknown with gaps as missing', lenOf('missing') === new E.Scorer(E.encode(asQ, 'missing')).length(tree));
  const asX = P.parseMatrix(nex([['A', '0{0X}1'], ['B', '1{1X}1'], ['C', '0X0'], ['D', '1{01}0']]));
  check('{0-} is {0, extra state} with gaps as an extra state', lenOf('extra') === new E.Scorer(E.encode(asX, 'missing')).length(tree));

  const mc = P.parseMatrix(nex([['A', '01011'], ['B', '..1..'], ['C', '1.0.?']], '', 'matchchar=.'));
  check('MATCHCHAR copies the first row', JSON.stringify(mc.cells.map((r) => r.map((c) => (Array.isArray(c) ? c.join('') : c)).join(''))) === '["01011","01111","110??"]'.replace('110??', '1101?'), mc.cells);
  check('MATCHCHAR is not a state', mc.symbols.join('') === '01', mc.symbols);
  const lc = P.parseMatrix(nex([['A', '01nx'], ['B', '1x0n'], ['C', '0101']], '', 'missing=n gap=x'));
  check('lower-case MISSING and GAP symbols', lc.symbols.join('') === '01' && lc.cells[0][2] === '?' && lc.cells[0][3] === '-', [lc.symbols, lc.cells[0]]);
}

// ---- ordered characters: distance by state value, not by position among the states present
{
  const m = P.parseMatrix(nex([['A', '0'], ['B', '0'], ['C', '3'], ['D', '3']], 'begin assumptions; ctype ord: 1; end;'));
  check('states 0-3 when only 0 and 3 occur', m.symbols.join('') === '0123', m.symbols);
  check('ordered 0 -> 3 is three steps', new E.Scorer(E.encode(m, 'inapplicable')).length(E.fromNested(4, [0, [1, [2, 3]]])) === 3);
  const dna = P.parseMatrix(nex([['A', 'AC'], ['B', 'GT'], ['C', 'AA']], '', 'datatype=dna'));
  check('letter states are not padded', dna.symbols.join('') === 'ACGT', dna.symbols);
}

// ---- TNT piwe before xread, cc without a space
{
  const a = P.parseMatrix('mxram 100;\npiwe=3;\nxread\n2 3\nA 01\nB 11\nC 10\n;\nproc/;');
  check('piwe before xread', a.concavity === 3, a.concavity);
  const b = P.parseMatrix('xread\n3 3\nA 012\nB 112\nC 100\n;\ncc+1 2;\nproc/;');
  check('cc+1 2; without a space', b.ordered[1] && b.ordered[2] && !b.ordered[0], b.ordered);
  const c = P.parseMatrix('piwe=5;\nxread\n2 3\nA 01\nB 11\nC 10\n;\npiwe-;\nproc/;');
  check('piwe- after the matrix turns it off', c.concavity === null, c.concavity);
}

// ---- unambiguous changes to and from the gap with gaps as an extra state
{
  const m = P.parseMatrix(nex([['A', '0'], ['B', '0'], ['C', '0'], ['D', '-'], ['E', '-'], ['F', '-']]));
  const data = E.encode(m, 'extra');
  const t = E.fromNested(6, [0, [1, [2, [3, [4, 5]]]]]);
  const ch = E.branchChanges(data, t, 'unambiguous').flat();
  check('0 -> - is listed as unambiguous', ch.some((x) => x.to.includes('-') || x.from.includes('-')), ch);
}

// ---- implied weights: trees with exactly the same fit stay in the tree set
{
  // two trees with extra steps {6,6,6} and {3,3,3,3} at k = 3 both have fit exactly 2
  const pool = new E.Pool(10, 7);
  pool.offer(E.fromNested(4, [0, [1, [2, 3]]]), 8589934593);
  pool.offer(E.fromNested(4, [0, [2, [1, 3]]]), 8589934592);
  check('pool keeps scores within the rounding tolerance', pool.map.size === 2, pool.map.size);
  check('stopping rule counts them as hits', E.repsNeeded([5, 4, 5, 4], 4, 10, 0, 1) === 4 && E.repsNeeded([5, 4, 5, 4], 4, 10, 0, 0) === 0);
  const m = P.parseMatrix(nex([['A', '00000'], ['B', '01101'], ['C', '11011'], ['D', '10110'], ['E', '11111']]));
  const data = E.encode(m, 'missing', { concavity: 3 });
  check('tolerance is set under implied weights only', data.tieTol > 0 && data.tieTol < data.scoreDiv / 1e6 && E.encode(m, 'missing').tieTol === 0, data.tieTol);
}

// ---- NEXUS TRANSLATE with quoted names containing apostrophes and commas
{
  const names = ["O'Brien", 'a,b', 'C', 'D'];
  const r = E.parseTrees("#NEXUS\nbegin trees;\n translate 1 'O''Brien', 2 'a,b', 3 C, 4 D;\n tree t = (1,(2,(3,4)));\nend;", names);
  check('translate with quoted apostrophe and comma', r.trees.length === 1 && !r.errors.length, r.errors);
}

// ---- Bremer and resampling keep the user's "must not" constraint (brute force, 7 taxa)
{
  const rnd = rng(11), n = 7;
  let values = 0, wrong = 0, forbidden = 0, reps = 0;
  for (let it = 0; it < 40 && values < 40; it++) {
    const rows = [];
    for (let t = 0; t < n; t++) { let s = ''; for (let c = 0; c < 12; c++) s += rnd() < 0.5 ? '0' : '1'; rows.push(['T' + t, s]); }
    const mx = P.parseMatrix(nex(rows));
    const data = E.encode(mx, 'missing'), sc = new E.Scorer(data);
    const all = []; allTrees(n, (t) => all.push({ score: sc.length(t), clades: new Set(E.clades(t)), nested: E.toNested(t) }));
    const best = Math.min(...all.map((x) => x.score)), opt = all.filter((x) => x.score === best);
    const strict = [...opt[0].clades].filter((c) => opt.every((x) => x.clades.has(c)));
    if (!strict.length) continue;
    const X = strict[0], allowed = all.filter((x) => !x.clades.has(X));
    const bestC = Math.min(...allowed.map((x) => x.score)), optC = allowed.filter((x) => x.score === bestC);
    const strictC = [...optC[0].clades].filter((c) => optC.every((x) => x.clades.has(c)));
    if (!strictC.length) continue;
    const constraint = { pos: [], neg: [E.cladeTips(X, n)] };
    const b = E.bremer(data, { clades: strictC, startTrees: optC.slice(0, 8).map((x) => x.nested), reps: 3, ratchetIter: 3, seed: 1, timeLimit: 60000, constraint });
    for (const [k, v] of b.bounds) {
      const truth = Math.min(...allowed.filter((x) => !x.clades.has(k)).map((x) => x.score)) - bestC;
      values++; if (v - bestC < truth) { wrong++; if (process.env.DEBUG) console.log('it', it, 'clade', E.cladeTips(k, n), 'got', v - bestC, 'truth', truth, 'X', E.cladeTips(X, n), 'done', b.done, 'skipped', b.skipped); }
    }
    const r = E.resample(data, { method: 'bootstrap', replicates: 5, seed: 2, timeLimit: 60000, constraint });
    reps += r.replicates; forbidden += r.counts.get(X) || 0;
  }
  check(`Bremer under a "not:" constraint never below the true value (${values} values)`, values >= 20 && wrong === 0, { values, wrong });
  check('bootstrap replicates never contain the forbidden clade', reps > 0 && forbidden === 0, { reps, forbidden });
}

// ---- searches cut by the time limit are not counted
{
  const rows = [];
  const rnd = rng(5);
  for (let t = 0; t < 40; t++) { let s = ''; for (let c = 0; c < 60; c++) s += rnd() < 0.5 ? '0' : rnd() < 0.5 ? '1' : '2'; rows.push(['T' + t, s]); }
  const mx = P.parseMatrix(nex(rows)), data = E.encode(mx, 'missing');
  const r = E.resample(data, { method: 'bootstrap', replicates: 3, seed: 1, timeLimit: 1 });
  check('a bootstrap replicate cut by the time limit is not counted', r.replicates === 0 && r.timedOut, { rep: r.replicates, t: r.timedOut });
  const S = new E.Searcher(data, { seed: 1 }), t = E.rootOn(S.climb(S.wagner(), null, 'spr').tree, 0);
  const b = E.bremer(data, { clades: E.clades(t).slice(0, 2), startTrees: [E.toNested(t)], reps: 2, ratchetIter: 2, seed: 1, timeLimit: 1 });
  check('a Bremer search cut by the time limit is not reported as done', b.done.length === 0 && b.timedOut, b.done);
}

// ---- re-rooting on several outgroup taxa keeps support values and vertex ids
{
  const nested = [[0, 1], [[2, 3], [4, 5]]];
  const t = E.rootOn(E.fromNested(6, nested), 0);
  const tree = E.nestedFromTree ? E.nestedFromTree(t) : null;
  if (tree) {
    (function mark(x) { if (x.children) { x.support = 0.5; x.children.forEach(mark); } })(tree);
    const rr = E.rerootOnOutgroup(tree, [0, 1], 6).node;
    let missing = 0;
    (function walk(x, isRoot) { if (x.children) { if (!isRoot && x.flipped && x.support === undefined) missing++; x.children.forEach((c) => walk(c, false)); } })(rr, true);
    check('flipped nodes keep a support value', missing === 0, missing);
  }
}

console.log(`${pass}/${pass + fail} review regression checks passed`);
process.exit(fail ? 1 : 0);
