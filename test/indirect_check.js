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
console.log(`\nAll ${candidates} rearrangement candidates: indirect score == full rescoring`);
