// Hand-worked inapplicable-data test cases (expected lengths transcribed from
// the MorphyLib test suite, testfitch.c, by Brazeau et al.). Each case is a
// single character on a 12-taxon tree; tested on both published rootings and
// on every possible rooting (the length should not depend on the root).
const P = require('../src/parser.js');
const E = require('../src/engine.js');

const cases = [
  ['23--1??--032', 5], ['1---1111---1', 2], ['1100----1100', 3], ['11-------100', 2],
  ['----1111---1', 1], ['01----010101', 5], ['01---1010101', 5], ['1??--??--100', 2],
  ['21--3??--032', 5], ['11--1??--111', 2], ['11--1000001-', 2], ['01------0101', 4],
  ['110--?---100', 3], ['11--1??--111', 2], ['210--100--21', 5], ['????----1???', 0],
  ['23--1----032', 5], ['1----1----1-', 2], ['-1-1-1--1-1-', 4], ['23--1??--032', 5],
  ['--------0101', 2], ['10101-----01', 4], ['011--?--0011', 3], ['110--??--100', 3],
  ['11--1000001-', 2], ['21--1----012', 5], ['11----111111', 1], ['10101-----01', 4],
  ['210210------', 4], ['----1111----', 0], ['230--??1--32', 5], ['023--??1--32', 5],
  ['023-???1--32', 4], ['23--1?1--023', 5], ['----1010----', 2], ['------11---1', 1],
  ['10----11---1', 3], ['320--??3--21', 5], ['-------1----', 0], ['0--11-111111', 1],
  ['1???---??011', 2],
];

const n = 12;
const caterpillar = [0, [1, [2, [3, [4, [5, [6, [7, [8, [9, [10, 11]]]]]]]]]]];
const balancedRoot = [[[[[[0, 1], 2], 3], 4], 5], [6, [7, [8, [9, [10, 11]]]]]];

let pass = 0, fail = 0, rootFail = 0;
cases.forEach(([pattern, expected], k) => {
  const taxa = [...Array(n).keys()].map((i) => 't' + i);
  const nex = '#NEXUS\nbegin data; dimensions ntax=12 nchar=1; format symbols="0123" missing=? gap=-; matrix\n' +
    taxa.map((t, i) => t + ' ' + pattern[i]).join('\n') + '\n;\nend;';
  const mx = P.parseMatrix(nex);
  const data = E.encode(mx, 'inapplicable');
  const sc = new E.Scorer(data);
  const a = sc.length(E.fromNested(n, caterpillar));
  const b = sc.length(E.fromNested(n, balancedRoot));
  const ok = a === expected && b === expected;
  if (ok) pass++; else { fail++; console.log(`  FAIL case ${k} ${pattern}: expected ${expected}, got ${a} (caterpillar) / ${b} (alt root)`); }
  // every rooting
  const lens = new Set();
  for (let r = 0; r < n; r++) lens.add(sc.length(E.rootOn(E.fromNested(n, caterpillar), r)));
  if (lens.size > 1) { rootFail++; console.log(`  root-dependence in case ${k} ${pattern}: lengths ${[...lens].join(',')}`); }
});
console.log(`${pass}/${cases.length} hand-worked inapplicable cases correct on both published rootings`);
console.log(`${cases.length - rootFail}/${cases.length} cases give the same length on all 12 rootings`);
process.exit(fail ? 1 : 0);
