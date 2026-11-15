// Does a search with the Standard settings find the same best length from two seeds,
// and how long does it take? One core per search.
//   DATASETS=path/to/TreeSearch/inst/datasets node test/convergence.js [seed] [files...]
// Prints one JSON line per matrix: best length, replicates that found it / replicates
// used, the score of each replicate, trees found, seconds, and whether the time limit
// (LIMIT seconds, default 600) cut the search short.
const fs = require('fs');
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const D = process.env.DATASETS;
const seed = parseInt(process.argv[2] || '1', 10);
const LIMIT = parseInt(process.env.LIMIT || '600', 10);
const files = process.argv.length > 3 ? process.argv.slice(3) : fs.readdirSync(D).filter((x) => x.endsWith('.nex')).map((x) => x.replace('.nex', ''));
for (const f of files) {
  const mx = P.parseMatrix(fs.readFileSync(D + '/' + f + '.nex', 'utf8'));
  const data = E.encode(mx, process.env.GAP || 'inapplicable', process.env.K ? { concavity: +process.env.K } : null);
  const r = E.analyse(data, { hits: 5, minReps: 10, maxReps: 50, ratchetIter: 10, maxTrees: 100, seed, timeLimit: LIMIT * 1000 }); // the app's Standard effort
  console.log(JSON.stringify({ matrix: f, taxa: mx.taxa.length, characters: mx.nchar, inapplicableColumns: data.naChars.length, seed, best: r.bestScore / data.scoreDiv, foundIn: r.hits, replicates: r.reps, scores: r.repScores.map((x) => x / data.scoreDiv), trees: r.trees.length, seconds: +(r.elapsed / 1000).toFixed(1), timedOut: r.timedOut }));
}
