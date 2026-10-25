// Verifies the search's rearrangement generator: the set of trees the engine
// visits from a starting tree must equal the full TBR neighbourhood, computed
// independently by bisecting every edge of the unrooted tree and reconnecting
// every pair of edges (brute force on an adjacency representation).
const E = require('../src/engine.js');

const { bruteTBR } = require('./tbr_brute.js');
const rng = new E.RNG(5);
let ok = 0, total = 0;
for (const n of [5, 6, 7, 8, 9, 10, 12]) {
  for (let rep = 0; rep < 6; rep++) {
    const t = E.randomTree(n, rng);
    const expected = bruteTBR(t);
    // engine neighbourhood: record every tree evaluated
    const data = E.encode({ taxa: [...Array(n).keys()].map(String), cells: Array.from({ length: n }, () => [['0']]), nchar: 1, symbols: ['0'] }, 'missing');
    const S = new E.Searcher(data, { seed: rep });
    const seen = new Set();
    // a constant character makes every tree score 0, so keep() sees every candidate
    S.rearrange(t.clone(), 0, null, 'tbr', 0, (tt) => seen.add(E.treeKey(E.rootOn(tt, 0))));
    seen.delete(E.treeKey(t));
    const missing = [...expected].filter((k) => !seen.has(k)).length;
    const extra = [...seen].filter((k) => !expected.has(k)).length;
    total++;
    if (!missing && !extra) ok++; else console.log(`  n=${n}: expected ${expected.size}, engine ${seen.size}, missing ${missing}, extra ${extra}`);
  }
}
console.log(`${ok}/${total} random trees: engine TBR neighbourhood == brute-force TBR neighbourhood`);
process.exit(ok === total ? 0 : 1);
