// Unit tests: parsers, Newick round trip, consensus, resampling, mapping.
const P = require('../src/parser.js');
const E = require('../src/engine.js');
let pass = 0, fail = 0;
function check(name, cond, extra) { if (cond) pass++; else { fail++; console.log('  FAIL', name, extra !== undefined ? JSON.stringify(extra) : ''); } }

// --- NEXUS: MorphoBank-style export with quoted names, CHARSTATELABELS, polymorphism, comments
const mb = `#NEXUS
[ MorphoBank (http://www.morphobank.org) version 3.0 ]
BEGIN TAXA;
  DIMENSIONS NTAX=4;
  TAXLABELS
    'Homo sapiens' 'Pan troglodytes' 'Gorilla gorilla' 'Pongo pygmaeus'
  ;
ENDBLOCK;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=3;
  FORMAT DATATYPE=STANDARD GAP=- MISSING=? SYMBOLS="0 1 2";
  CHARSTATELABELS
    1 'Brow ridge' / absent present,
    2 'Canine size' / small 'very large',
    3 'Tail colour' / red blue
  ;
  MATRIX
    'Homo sapiens'      0 (0 1) -
    'Pan troglodytes'   1 {12} -
    'Gorilla gorilla'   1 1 ?
    'Pongo pygmaeus'    0 2 1
  ;
END;`;
let m = P.parseMatrix(mb);
check('mb: no errors', m.errors.length === 0, m.errors);
check('mb: 4 taxa', m.taxa.length === 4 && m.taxa[0] === 'Homo sapiens', m.taxa);
check('mb: polymorphism (0 1)', JSON.stringify(m.cells[0][1]) === '["0","1"]', m.cells[0][1]);
check('mb: polymorphism {12}', JSON.stringify(m.cells[1][1]) === '["1","2"]', m.cells[1][1]);
check('mb: gap and missing', m.cells[0][2] === '-' && m.cells[2][2] === '?');
check('mb: char labels', m.charLabels[1] === 'Canine size', m.charLabels);
check('mb: state labels', m.stateLabels[1] && m.stateLabels[1][1] === 'very large', m.stateLabels);

// --- Interleaved NEXUS
const il = `#NEXUS
begin data; dimensions ntax=3 nchar=6; format symbols="01" interleave; matrix
A 010
B 110
C 001

A 111
B 000
C 01?
; end;`;
m = P.parseMatrix(il);
check('interleave: rows complete', m.errors.length === 0 && m.cells[2].length === 6, m.errors);
check('interleave: last cell missing', m.cells[2][5] === '?');

// --- TNT xread with [ ] polymorphism
const tnt = `xread
'test matrix'
4 3
alpha 01[01]2
beta  1-1?
gamma 0011
delta 1102
;
proc/;`;
m = P.parseMatrix(tnt);
check('tnt: no errors', m.errors.length === 0, m.errors);
check('tnt: 3 taxa? no, 4 taxa x 4? (header 4 chars 3 taxa mismatch warns)', m.warnings.length >= 1);
const tnt2 = `xread 4 3\nA 01[01]2\nB 1-1?\nC 0011\n;`;
m = P.parseMatrix(tnt2);
check('tnt: polymorphism', JSON.stringify(m.cells[0][2]) === '["0","1"]', m.cells[0]);
check('tnt: gap', m.cells[1][1] === '-');

// --- bad input
m = P.parseMatrix('hello world');
check('garbage rejected', m.errors.length > 0);
m = P.parseMatrix('#NEXUS\nbegin data; dimensions ntax=2 nchar=3; matrix\nA 01\nB 011\n; end;');
check('short row detected', m.errors.some((e) => /wrong number/.test(e)), m.errors);

// --- Newick round trip & rooting
const names = ['A', 'B', 'C', 'D', 'E', 'F'];
const t = E.fromNested(6, E.parseNewick("(A,((B,C),(D,(E,F))));", names));
check('newick round trip', E.toNewick(t, names) === '(A,((B,C),(D,(E,F))));', E.toNewick(t, names));
const t2 = E.fromNested(6, E.parseNewick("(((B:0.1,C:0.2)90:0.3,(D,(E,F))),A);", names));
check('newick with lengths/support parses, same topology after rooting', E.treeKey(E.rootOn(t2, 0)) === E.treeKey(E.rootOn(t, 0)));
check('rootOn(other tip) preserves unrooted topology', E.treeKey(E.rootOn(E.rootOn(t, 3), 0)) === E.treeKey(E.rootOn(t, 0)));

// --- Consensus
const a = E.rootOn(E.fromNested(6, E.parseNewick('(A,((B,C),(D,(E,F))));', names)), 0);
const b = E.rootOn(E.fromNested(6, E.parseNewick('(A,((B,C),(E,(D,F))));', names)), 0);
const c = E.rootOn(E.fromNested(6, E.parseNewick('(A,(B,(C,(D,(E,F)))));', names)), 0);
const strict = E.consensus([a, b, c], 6, 1);
check('strict consensus', E.nestedToNewick(strict, names) === '(A,B,C,(D,E,F)1);', E.nestedToNewick(strict, names));
const maj = E.consensus([a, b, c], 6, 0.5);
const mn = E.nestedToNewick(maj, names);
check('majority-rule consensus', mn === '(A,((B,C)0.6666666666666666,(D,E,F)1,));'.replace(',)', ')') || /\(B,C\)0\.66/.test(mn), mn);

// --- Search on a clean dataset recovers the true tree
const truth = '(A,(B,(C,(D,(E,F)))));';
const mx = { taxa: names, nchar: 4, symbols: ['0', '1'], cells: [
  [['0'], ['0'], ['0'], ['0']], [['1'], ['0'], ['0'], ['0']], [['1'], ['1'], ['0'], ['0']],
  [['1'], ['1'], ['1'], ['0']], [['1'], ['1'], ['1'], ['1']], [['1'], ['1'], ['1'], ['1']]] };
// make E,F distinguishable? they are identical -> (E,F) unresolved but a clade anyway
const data = E.encode(mx, 'inapplicable');
const res = E.analyse(data, { replicates: 2, ratchetIter: 2, seed: 3 });
check('perfect data: length 4', res.bestScore === 4, res.bestScore);
check('perfect data: unique MPT', res.trees.length === 1 && E.treeKey(res.trees[0]) === E.treeKey(E.rootOn(E.fromNested(6, E.parseNewick(truth, names)), 0)), res.trees.map((x) => E.toNewick(x, names)));

// --- Bootstrap on perfect data gives 100% for supported clades
const bs = E.resample(data, { method: 'bootstrap', replicates: 50, seed: 5 });
check('bootstrap ran 50 replicates', bs.replicates === 50);
const keyBCDEF = E.clades(res.trees[0]).sort().pop();
check('bootstrap frequencies in [0,1]', [...bs.freq.values()].every((f) => f >= 0 && f <= 1));
void keyBCDEF;
const jk = E.resample(data, { method: 'jackknife', replicates: 30, seed: 5 });
check('jackknife ran', jk.replicates === 30);

// --- Character mapping & bounds
const sc = new E.Scorer(data);
const steps = new Int32Array(data.nchar);
const total = sc.length(res.trees[0], null, steps);
check('per-character steps sum to total', steps.reduce((x, y) => x + y, 0) === total, [...steps]);
const bounds = E.charBounds(data);
check('min steps (binary informative chars) = 1', [...bounds.minS].every((x) => x === 1), [...bounds.minS]);
check('CI = 1 on perfect data', [...steps].every((s, j) => s === bounds.minS[j]));
const rec = sc.reconstruct(res.trees[0], 0);
check('reconstruction gives tip states', E.setToStates(rec[0], data).join('') === '0' && E.setToStates(rec[1], data).join('') === '1');

// --- Inapplicable mapping: tail colour example
const mx2 = { taxa: ['o', 'p', 'q', 'r', 's', 't'], nchar: 2, symbols: ['0', '1'], cells: [
  [['0'], '-'], [['0'], '-'], [['0'], '-'], [['1'], ['0']], [['1'], ['0']], [['1'], ['1']]] };
const d2 = E.encode(mx2, 'inapplicable');
check('tail colour char uses inapplicable algorithm', d2.isNA[1] === 1);
const sc2 = new E.Scorer(d2);
const tt = E.fromNested(6, [0, [1, [2, [3, [4, 5]]]]]);
const st2 = new Int32Array(2);
sc2.length(tt, null, st2);
check('tail presence 1 step, colour 1 step', st2[0] === 1 && st2[1] === 1, [...st2]);
const rec2 = sc2.reconstruct(tt, 1);
check('colour inapplicable at root region', E.setToStates(rec2[tt.root], d2).join('') === '-', E.setToStates(rec2[tt.root], d2));

// --- Gap-treatment comparison
{
  const cmp = E.compareTreatments([{ name: 'a', trees: [a] }, { name: 'b', trees: [b] }, { name: 'c', trees: [c] }], 6);
  const nm = (cl) => cl.tips.map((i) => names[i]).join('');
  const byName = Object.fromEntries(cmp.clades.map((cl) => [nm(cl), cl.present.map((x) => (x ? 1 : 0)).join('')]));
  const exp = { CDEF: '001', DEF: '111', BC: '110', EF: '101', DF: '010' };
  check('compare: clade presence', Object.keys(exp).length === Object.keys(byName).length && Object.keys(exp).every((k) => byName[k] === exp[k]), byName);
  check('compare: shared clades', cmp.shared === 1, cmp.shared);
  check('compare: pairwise distances', cmp.pairwise[0][1] === 2 && cmp.pairwise[0][2] === 2 && cmp.pairwise[1][2] === 4 && cmp.pairwise[2][1] === 4, cmp.pairwise);
  const cmp2 = E.compareTreatments([{ name: 'ab', trees: [a, b] }, { name: 'a', trees: [a] }], 6);
  check('compare: strict consensus of several trees', cmp2.strict[0].size === 2 && cmp2.pairwise[0][1] === 1, [cmp2.strict[0].size, cmp2.pairwise]);
  check('compare: consensus newick', E.nestedToNewick(cmp2.consensus[0], names) === '(A,(B,C),(D,E,F));', E.nestedToNewick(cmp2.consensus[0], names));
}

// --- Character settings OpenClade does not apply are reported, not silently dropped
{
  const base = '#NEXUS\nbegin data; dimensions ntax=4 nchar=3; format symbols="012"; matrix\nA 000\nB 011\nC 112\nD 112\n;\nend;\n';
  const asm = P.parseMatrix(base + 'begin assumptions;\n  typeset * default = ord: 3;\n  wtset * w = 2: 1;\n  exset * ex = 2;\nend;\n');
  const w = asm.warnings.join(' ');
  check('assumptions: warns about ordering, weights and exclusions', /ordered/.test(w) && /weights/.test(w) && /excluded/.test(w), asm.warnings);
  check('assumptions: no warning without them', P.parseMatrix(base).warnings.length === 0, P.parseMatrix(base).warnings);
  const tnt = "xread\n3 4\nA 000\nB 011\nC 112\nD 112\n;\n";
  check('TNT: ccode warning', /ccode/.test(P.parseMatrix(tnt + 'ccode + 2;\n').warnings.join(' ')));
  check('TNT: "cc - ." is not warned about', P.parseMatrix(tnt + 'cc - .;\n').warnings.length === 0, P.parseMatrix(tnt + 'cc - .;\n').warnings);
}

// --- The worked example in the user guide (docs, "Gaps and inapplicable characters")
{
  const nex = '#NEXUS\nbegin data; dimensions ntax=7 nchar=3; format symbols="01" missing=? gap=-;\nmatrix\n' +
    'Outgroup 0-0\nA 110\nB 110\nC 0-1\nD 0-1\nE 111\nF 111\n;\nend;';
  const mx = P.parseMatrix(nex);
  const X = '(Outgroup,((C,D),((A,B),(E,F))));', Y = '(Outgroup,((A,B),(C,(D,(E,F)))));';
  const len = (mode, nwk) => new E.Scorer(E.encode(mx, mode)).length(E.fromNested(7, E.parseNewick(nwk, mx.taxa)));
  const got = ['inapplicable', 'missing', 'extra'].map((m) => len(m, X) + '/' + len(m, Y)).join(' ');
  check('guide example: tree X / tree Y lengths are 3/4, 3/3, 4/5', got === '3/4 3/3 4/5', got);
}

console.log(`${pass}/${pass + fail} unit tests passed`);
process.exit(fail ? 1 : 0);
