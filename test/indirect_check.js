// Verifies the fast indirect length calculation used during rearrangements:
// with verify mode on, every candidate tree is also scored from scratch and
// any disagreement throws. Runs SPR+TBR climbs with unit and random weights.
const fs = require('fs');
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const D = process.env.DATASETS;
const files = (process.env.FILES || 'Vinther2008,Longrich2010,Sansom2010,Sun2018,Aria2015,DeAssis2011,OLeary1999').split(',');
let candidates = 0;
for (const f of files) {
  const mx = P.parseMatrix(fs.readFileSync(D + '/' + f + '.nex', 'utf8'));
  for (const mode of ['inapplicable', 'missing', 'extra']) {
    const data = E.encode(mx, mode);
    const S = new E.Searcher(data, { seed: 11 });
    S.verify = true;
    const before = S.evaluations;
    const t = S.wagner();
    S.climb(t, null, 'tbr');
    const w = new Int32Array(data.nchar);
    for (let j = 0; j < data.nchar; j++) w[j] = S.rng.int(4); // includes zero weights
    S.climb(E.randomTree(data.ntax, S.rng), w, 'tbr');
    candidates += S.evaluations - before;
    console.log(f.padEnd(14), mode.padEnd(13), 'ok', S.evaluations - before, 'candidates verified');
  }
}
// Ordered characters, file weights, exclusions and implied weighting: the same check, plus
// the screening bound (verify mode throws if a bound exceeds the exact score).
let total2 = 0;
{
  const rng = new E.RNG(31);
  for (const f of (process.env.FILES2 || 'Vinther2008,Aria2015,Sun2018,Wills2012,Agnarsson2004').split(',')) {
    const mx = P.parseMatrix(fs.readFileSync(D + '/' + f + '.nex', 'utf8'));
    const nstates = (c) => { const st = new Set(); for (const r of mx.cells) if (Array.isArray(r[c])) r[c].forEach((q) => st.add(q)); return st.size; };
    mx.ordered = Array.from({ length: mx.nchar }, (_, c) => nstates(c) > 2 && rng.next() < 0.6);
    mx.weights = Array.from({ length: mx.nchar }, () => [1, 1, 2, 3, 0.5][rng.int(5)]);
    mx.excluded = Array.from({ length: mx.nchar }, () => rng.next() < 0.08);
    for (const mode of ['inapplicable', 'missing', 'extra']) for (const k of [null, 3, 12]) for (const regions of [undefined, 'decompose']) {
      if (regions && (mode !== 'inapplicable' || k === 3)) continue;
      const data = E.encode(mx, mode, { concavity: k, orderedRegions: regions });
      const S = new E.Searcher(data, { seed: 9 });
      S.verify = true;
      const before = S.evaluations;
      let r = S.climb(S.wagner(), null, 'spr', Date.now() + 4000);
      r = S.climb(r.tree, null, 'tbr', Date.now() + 4000);
      const mult = new Float64Array(data.nOrig); for (let c = 0; c < data.nOrig; c++) mult[c] = rng.next() < 0.25 ? 2 : 1;
      S.climb(r.tree.clone(), E.expandWeights(data, mult), 'spr', Date.now() + 2000);
      S.drift(r.tree.clone(), r.score, 6, Date.now() + 2000);
      if (new E.Scorer(data).length(r.tree) !== r.score) { console.log('FAIL: final score', f, mode, k); process.exit(1); }
      total2 += S.evaluations - before;
    }
    console.log(f.padEnd(14), 'ordered, weighted and excluded characters, equal and implied weights: ok');
  }
}
console.log(`\nAll ${candidates} rearrangement candidates: indirect score == full rescoring`);
console.log(`All ${total2} further candidates with ordered characters, weights, exclusions and implied weighting: indirect score == full rescoring, screening bound <= exact score`);
