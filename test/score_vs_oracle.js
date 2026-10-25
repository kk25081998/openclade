// Compares OpenClade tree lengths with MorphyLib (compiled C reference) on
// random trees for every dataset, under inapplicable and missing gap modes.
const fs = require('fs');
const { execFileSync } = require('child_process');
const P = require('../src/parser.js');
const E = require('../src/engine.js');

const D = process.env.DATASETS;
const ORACLE = process.env.ORACLE;
const NTREES = parseInt(process.env.NTREES || '200', 10);

function oracleInput(mx, trees, gapmode) {
  const n = mx.taxa.length;
  let mat = '';
  for (let t = 0; t < n; t++) for (const c of mx.cells[t]) {
    if (c === '?' || c === '-') mat += c;
    else if (c.length === 1) mat += c[0];
    else mat += '{' + c.join('') + '}';
  }
  mat += ';';
  let s = `${n} ${mx.nchar} ${trees.length} ${gapmode}\n${mat}\n`;
  for (const tr of trees) {
    // preorder renumbering of internal nodes
    const id = new Int32Array(2 * n - 1).fill(-1);
    let next = n;
    const order = [tr.root];
    while (order.length) {
      const v = order.shift();
      if (v < n) continue;
      id[v] = next++;
      order.push(tr.left[v], tr.right[v]);
    }
    for (let i = 0; i < n; i++) id[i] = i;
    const anc = new Array(2 * n - 1), left = new Array(n - 1), right = new Array(n - 1);
    for (let v = 0; v < 2 * n - 1; v++) {
      const nv = id[v];
      anc[nv] = tr.parent[v] < 0 ? nv : id[tr.parent[v]];
      if (v >= n) { left[nv - n] = id[tr.left[v]]; right[nv - n] = id[tr.right[v]]; }
    }
    s += anc.join(' ') + '\n' + left.join(' ') + '\n' + right.join(' ') + '\n';
  }
  return s;
}

let totalTrees = 0, mismatches = 0;
const files = fs.readdirSync(D).filter((f) => f.endsWith('.nex')).sort();
const rng = new E.RNG(12345);
for (const f of files) {
  const mx = P.parseMatrix(fs.readFileSync(D + '/' + f, 'utf8'));
  for (const [mode, gm] of [['inapplicable', 0], ['missing', 1]]) {
    const data = E.encode(mx, mode);
    const sc = new E.Scorer(data);
    const trees = [];
    for (let i = 0; i < NTREES; i++) {
      const t = E.randomTree(mx.taxa.length, rng);
      // also test unusual rootings for half the trees
      trees.push(i % 2 ? t : rerootRandom(t, rng));
    }
    const js = trees.map((t) => sc.length(t));
    const out = execFileSync(ORACLE, { input: oracleInput(mx, trees, gm), maxBuffer: 1 << 26 }).toString().trim().split('\n').map(Number);
    let bad = 0;
    for (let i = 0; i < trees.length; i++) if (js[i] !== out[i]) { bad++; if (bad <= 3) console.log('  MISMATCH', f, mode, 'tree', i, 'js', js[i], 'morphy', out[i]); }
    totalTrees += trees.length; mismatches += bad;
    console.log(f.padEnd(20), mode.padEnd(13), `${trees.length - bad}/${trees.length} match`, 'range', Math.min(...js), '-', Math.max(...js));
  }
}
console.log(`\nTOTAL: ${totalTrees - mismatches}/${totalTrees} trees scored identically to MorphyLib`);
process.exit(mismatches ? 1 : 0);

function rerootRandom(t, rng) {
  // root on a random tip (not necessarily tip 0)
  return E.rootOn(t, rng.int(t.n));
}
