// Sectorial search (off by default in the app) must be exact: the tree it returns
// re-scores to the score it reports and is never worse than the tree it was given.
// With verify mode on, every rearrangement inside a sector is also checked against
// full rescoring, and for standard characters under equal weights the change of
// length inside the sector must equal the change of length of the whole tree.
const fs = require('fs');
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const D = process.env.DATASETS;
let runs = 0, improved = 0, fail = 0;
for (const f of ['Dikow2009', 'Wilson2003', 'Agnarsson2004', 'Conrad2008']) for (const mode of ['missing', 'inapplicable']) for (const k of [null, 8]) {
  const mx = P.parseMatrix(fs.readFileSync(D + '/' + f + '.nex', 'utf8'));
  const data = E.encode(mx, mode, k ? { concavity: k } : null);
  const S = new E.Searcher(data, { seed: 3 });
  S.verify = true;
  const cur = S.climb(S.wagner(), null, 'spr'), before = cur.score;
  const out = E.sectorial(S, cur, 6);
  runs++;
  const full = new E.Scorer(data).length(out.tree);
  if (full !== out.score || out.score > before) { fail++; console.log('  FAIL', f, mode, k, before, out.score, full); }
  if (out.score < before) improved++;
}
console.log(`${runs - fail}/${runs} sectorial-search runs exact (4 matrices x 2 gap treatments x equal and implied weights); ${improved} improved the starting tree`);
process.exit(fail ? 1 : 0);
