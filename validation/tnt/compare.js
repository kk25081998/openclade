// Compares TNT's results (out/<name>.tre, written by the scripts from make_runs.js)
// with an OpenClade search of the same matrix, gaps as missing data, equal weights.
//   DATASETS=path/to/TreeSearch/inst/datasets node compare.js [seconds per matrix, default 300]
// For each matrix: the length of TNT's trees as scored by OpenClade, OpenClade's own
// best length, and whether the two strict consensus trees have the same clades.
// TNT collapses branches with no support, so consensus trees are compared after
// collapsing OpenClade's trees the same way.
// Set TREES=dir to read the tree files from another folder (used by the self-test).
const fs = require('fs');
const path = require('path');
const P = require('../../src/parser.js');
const E = require('../../src/engine.js');
const D = process.env.DATASETS;
const DIR = process.env.TREES || path.join(__dirname, 'out');
const LIMIT = parseInt(process.argv[2] || '300', 10) * 1000;
const cladeKeys = (nested, n) => { // clade keys of a tree that may have polytomies, rooted on tip 0
  const t = E.rootOn(E.fromNested(n, nested), 0); // resolves polytomies arbitrarily
  const all = new Set(E.clades(t));
  // keep only clades that are in the tree as given: a clade made by the arbitrary resolution is not
  const given = new Set();
  (function walk(x) {
    if (typeof x === 'number') return [x];
    const tips = x.flatMap(walk);
    const s = new Set(tips); let set = tips;
    if (s.has(0)) { set = []; for (let i = 0; i < n; i++) if (!s.has(i)) set.push(i); }
    if (set.length > 1 && set.length < n - 1) { const w = new Uint32Array(Math.ceil(n / 32)); for (const i of set) w[i >> 5] |= 1 << (i & 31); given.add(Array.from(w).join('.')); }
    return tips;
  })(nested);
  return [...given].filter((k) => all.has(k));
};
const strict = (sets) => { const first = sets[0]; return first.filter((k) => sets.every((s) => s.includes(k))).sort(); };
const rows = [];
let same = 0, shorterOC = 0, shorterTNT = 0, done = 0;
// SHARD=i/k runs every k-th matrix starting at i, to split the work across processes
const [shI, shK] = (process.env.SHARD || '0/1').split('/').map(Number);
const csv = process.env.SHARD ? `openclade_vs_tnt.shard${shI}.csv` : 'openclade_vs_tnt.csv';
const mats = fs.readdirSync(D).filter((x) => x.endsWith('.nex')).sort().filter((_, i) => i % shK === shI);
for (const f of mats) {
  const name = f.replace('.nex', '');
  // TNT's trees: from the random-addition run and, if present, the driven (xmult) run
  const files = [name + '.tre', name + '.xmult.tre'].map((x) => path.join(DIR, x)).filter((x) => fs.existsSync(x));
  if (!files.length) { console.log(name.padEnd(14), 'no TNT trees found (' + path.relative(process.cwd(), path.join(DIR, name + '.tre')) + ')'); continue; }
  const mx = P.parseMatrix(fs.readFileSync(path.join(D, f), 'utf8'));
  const n = mx.taxa.length, labels = mx.taxa.map((_, i) => 't' + (i + 1));
  const data = E.encode(mx, 'missing');
  const sc = new E.Scorer(data);
  const sets = [];
  for (const file of files) {
    const parsed = E.parseTrees(fs.readFileSync(file, 'utf8'), labels);
    if (!parsed.trees.length) { console.log(name.padEnd(14), 'could not read', path.basename(file) + ':', parsed.errors[0]); continue; }
    const lens = parsed.trees.map((t) => sc.length(E.rootOn(E.fromNested(n, t.nested), 0)));
    sets.push({ file: path.basename(file), trees: parsed.trees, lens, best: Math.min(...lens) });
  }
  if (!sets.length) continue;
  const tntBest = Math.min(...sets.map((x) => x.best));
  const lenOf = (suffix) => { const x = sets.find((y) => y.file === name + suffix); return x ? x.best : ''; };
  // every TNT tree at TNT's best length, from either run
  const tntTrees = []; for (const x of sets) x.trees.forEach((t, i) => { if (x.lens[i] === tntBest) tntTrees.push(t); });
  const poly = tntTrees.some((t) => !t.binary);
  const parsed = { trees: tntTrees };
  const t0 = Date.now();
  const res = E.analyse(data, { hits: 10, maxReps: 200, ratchetIter: 20, maxTrees: 5000, seed: 1, timeLimit: LIMIT });
  const secs = ((Date.now() - t0) / 1000).toFixed(1);
  const ocStrict = strict(E.collapseAll(data, res.trees).map((c) => c.keys));
  // both sides are collapsed by the same rule (minimum length zero) before the strict consensus
  const tntStrict = poly ? strict(parsed.trees.map((t) => cladeKeys(t.nested, n)))
    : strict(E.collapseAll(data, parsed.trees.map((t) => E.rootOn(E.fromNested(n, t.nested), 0))).map((c) => c.keys));
  const consSame = ocStrict.join('|') === tntStrict.join('|');
  done++;
  if (res.bestScore === tntBest) same++; else if (res.bestScore < tntBest) shorterOC++; else shorterTNT++;
  rows.push([name, n, lenOf('.tre'), lenOf('.xmult.tre'), tntBest, res.bestScore, res.hits + ' of ' + res.reps, secs, parsed.trees.length, res.trees.length, tntStrict.length, ocStrict.length, consSame ? 'yes' : 'no', res.timedOut ? 'time limit' : '', poly ? 'TNT trees have polytomies (length is an upper bound)' : ''].join(','));
  console.log(name.padEnd(14), 'TNT', String(tntBest).padStart(5), ' OpenClade', String(res.bestScore).padStart(5), res.bestScore === tntBest ? ' same length' : res.bestScore < tntBest ? ' OpenClade SHORTER' : ' TNT SHORTER', '| strict consensus', consSame ? 'identical' : `differs (${tntStrict.length} vs ${ocStrict.length} clades)`, res.timedOut ? '| time limit reached' : '');
}
fs.writeFileSync(path.join(__dirname, csv), 'matrix,taxa,tnt_mult,tnt_xmult,tnt_best,openclade_best,openclade_found_in,openclade_seconds,tnt_trees,openclade_trees,tnt_strict_clades,openclade_strict_clades,strict_identical,note,note2\n' + rows.join('\n') + '\n');
console.log(`\n${done} matrices compared: same best length on ${same}; OpenClade shorter on ${shorterOC}; TNT shorter on ${shorterTNT}`);
