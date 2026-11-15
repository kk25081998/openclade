// Tree files, outgroup re-rooting, unsupported-branch collapsing on real trees,
// ordered characters with inapplicable tokens, and the inapplicable-coding check.
const fs = require('fs');
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const D = process.env.DATASETS;
let pass = 0, fail = 0;
const check = (n, c, x) => { if (c) pass++; else { fail++; console.log('  FAIL', n, x === undefined ? '' : JSON.stringify(x)); } };
const names = ['Alpha', 'Beta gamma', 'C', 'D', 'E', 'F'];
const key = (nested) => E.treeKey(E.rootOn(E.fromNested(6, nested), 0));
const ref = key([0, [[1, 2], [3, [4, 5]]]]);

// ---- tree files
{
  let r = E.parseTrees("(Alpha,(('Beta gamma',C),(D,(E,F))));\n((E,F),D,((C,Beta_gamma),Alpha));", names);
  check('newick: two trees', r.trees.length === 2 && !r.errors.length, r);
  check('newick: quoted and underscore names, same topology', r.trees.every((t) => key(t.nested) === ref));
  check('newick: basal trichotomy counts as resolved', r.trees[1].binary === true);
  r = E.parseTrees("#NEXUS\nbegin trees;\n translate 1 Alpha, 2 'Beta gamma', 3 C, 4 D, 5 E, 6 F;\n tree one = [&R] (1,((2,3),(4,(5,6))));\n tree * two = [&U] (1,(2,3),(4,5,6));\nend;", names);
  check('nexus: translate table', r.trees.length === 2 && key(r.trees[0].nested) === ref && r.trees[0].name === 'one', r);
  check('nexus: polytomy is flagged', r.trees[1].binary === false);
  r = E.parseTrees("tread 'two trees'\n(Alpha ((Beta_gamma C )(D (E F ))))*\n(0 ((1 2 )(3 (4 5 ))));\nproc-;", names);
  check('TNT: names and numbers', r.trees.length === 2 && r.trees.every((t) => key(t.nested) === ref), r);
  r = E.parseTrees('(Alpha,(Beta_gamma,C),(D,E));', names);
  check('a tree missing a taxon is refused with the reason', r.trees.length === 0 && /not in the tree/.test(r.errors[0]), r.errors);
  r = E.parseTrees('(Alpha,(Beta_gamma,C),(D,(E,(F,Zeta))));', names);
  check('an unknown taxon is named', /Zeta/.test(r.errors.join(' ')), r.errors);
  r = E.parseTrees('(Alpha,(Beta_gamma,C),(D,(E,(F,E))));', names);
  check('a repeated taxon is refused', r.trees.length === 0 && /more than once/.test(r.errors[0]), r.errors);
  check('no trees at all', E.parseTrees('hello', names).errors.length === 1);
}

// ---- outgroup of several taxa
{
  const t = E.fromNested(6, [0, [1, [2, [3, [4, 5]]]]]);
  const node = E.nestedFromTree(t);
  let r = E.rerootOnOutgroup(node, [0, 1], 6);
  const tips = (x) => (typeof x.tip === 'number' ? [x.tip] : x.children.flatMap(tips)).sort();
  check('outgroup {0,1}: root splits outgroup from ingroup', r.mono && r.node.children.length === 2 && tips(r.node.children[0]).join() === '0,1' && tips(r.node.children[1]).join() === '2,3,4,5', r.node.children.map(tips));
  r = E.rerootOnOutgroup(node, [0, 5], 6);
  check('outgroup {0,5}: ingroup is not a clade', r.mono === false && r.node === node);
  r = E.rerootOnOutgroup(node, [0], 6);
  check('single outgroup: unchanged', r.mono && r.node === node);
  // consensus with a polytomy at the base
  const cons = E.consensus([E.fromNested(6, [0, [1, [[2, 3], [4, 5]]]]), E.fromNested(6, [0, [[1, [2, 3]], [4, 5]]])], 6, 1);
  r = E.rerootOnOutgroup(cons, [0, 1], 6);
  check('outgroup on a consensus: ingroup must be a clade there too', r.mono === false);
  r = E.rerootOnOutgroup(cons, [0, 4, 5], 6);
  check('outgroup {0,4,5} on that consensus', r.mono === false || tips(r.node.children[1]).join() === '1,2,3');
}

// ---- ordered characters scored with the inapplicable algorithm
{
  // two tails that arose separately (A,B and E,F), with different lengths: 0 and 2
  const rows = (a, b) => `#NEXUS\nbegin data; dimensions ntax=9 nchar=1; format symbols="0123" gap=-; matrix\nO -\nA ${a}\nB ${a}\nC -\nD -\nE ${b}\nF ${b}\nG -\nH -\n;\nend;\nbegin assumptions; typeset * t = ord: 1; end;`;
  const tr = E.fromNested(9, E.parseNewick('(O,(((A,B),(C,D)),((E,F),(G,H))));', ['O', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']));
  const len = (text, co) => new E.Scorer(E.encode(P.parseMatrix(text), 'inapplicable', co)).length(tr);
  check('two separate regions, states 0 and 2: one step (the second region), no change counted between regions', len(rows(0, 2)) === 1, len(rows(0, 2)));
  check('the same with states renumbered 1 and 3', len(rows(1, 3)) === 1, len(rows(1, 3)));
  check('unordered gives the same here', len(rows(0, 2), { applyOrdered: false }) === 1);
  check("TreeSearch-compatible mode counts the region once per binary column (2 for states 0,2; 3 for states 1,3)", len(rows(0, 2), { orderedRegions: 'decompose' }) === 2 && len(rows(1, 3), { orderedRegions: 'decompose' }) === 3, [len(rows(0, 2), { orderedRegions: 'decompose' }), len(rows(1, 3), { orderedRegions: 'decompose' })]);
  // one region with a change 0 -> 2 inside it: two steps ordered, one unordered
  const one = `#NEXUS\nbegin data; dimensions ntax=9 nchar=1; format symbols="0123" gap=-; matrix\nO -\nA 0\nB 0\nC 2\nD 2\nE -\nF -\nG -\nH -\n;\nend;\nbegin assumptions; typeset * t = ord: 1; end;`;
  check('one region, 0 -> 2 inside it: 2 steps ordered, 1 unordered', len(one) === 2 && len(one, { applyOrdered: false }) === 1, [len(one), len(one, { applyOrdered: false })]);
  // a gap read as an extra state has no place in the order: the character is analysed unordered and reported
  const d = E.encode(P.parseMatrix(one), 'extra');
  check('ordered + gaps as extra state: analysed unordered and reported', d.orderedAsUnordered.length === 1 && d.charOrdered[0] === 0);
}

// ---- collapsing on real most parsimonious trees
if (D) {
  const mx = P.parseMatrix(fs.readFileSync(D + '/Vinther2008.nex', 'utf8'));
  const data = E.encode(mx, 'inapplicable');
  const res = E.analyse(data, { replicates: 4, ratchetIter: 4, maxTrees: 200, seed: 1, timeLimit: 60000 });
  const col = E.collapseAll(data, res.trees);
  check('Vinther: 52 binary trees', res.trees.length === 52, res.trees.length);
  check('collapsing gives fewer, distinct trees', col.length >= 1 && col.length < res.trees.length, col.length);
  check('every binary tree maps to one collapsed tree', col.reduce((a, c) => a + c.from.length, 0) === res.trees.length);
  const strictBin = E.consensus(res.trees, mx.taxa.length, 1), strictCol = E.consensusOfKeys(col.map((c) => c.keys), mx.taxa.length, 1);
  const count = (x) => (x.children ? x.children.reduce((a, c) => a + count(c), x.key ? 1 : 0) : 0);
  check('strict consensus of collapsed trees is no more resolved than of binary trees', count(strictCol) <= count(strictBin), [count(strictCol), count(strictBin)]);
  // every clade of the binary strict consensus has unambiguous or region support in at least one tree? not required; but
  // synapomorphies common to all trees must be changes of each single tree
  const common = E.commonChanges(data, res.trees, 'unambiguous');
  const first = E.branchChanges(data, res.trees[0], 'unambiguous'), k0 = E.nodeKeys(res.trees[0]);
  let ok = true, n = 0;
  for (const [k, list] of common) { const v = k0.indexOf(k); for (const ch of list) { n++; if (!first[v].some((x) => x.c === ch.c && x.from.join() === ch.from.join() && x.to.join() === ch.to.join())) ok = false; } }
  check('changes common to all trees are changes of each tree', ok && n > 0, n);
  console.log(`  Vinther2008: ${res.trees.length} binary trees -> ${col.length} after collapsing; ${n} changes common to all trees`);
}

// ---- inapplicable-coding check
{
  // character 1: tail absent/present; character 2: tail colour, with four coding slips
  const text = `#NEXUS\nbegin data; dimensions ntax=16 nchar=3; format symbols="01" gap=-; matrix
a 0-0\nb 0-1\nc 0-0\nd 0-1\ne 100\nf 111\ng 101\nh 0 1 0\ni 1 - 1\nj ? - 0\nk 0-1\nl 0-0\nm 0-1\nn 110\no 101\np 110\n;\nend;`;
  const mx = P.parseMatrix(text);
  const L = P.lintInapplicable(mx);
  const c2 = L.chars.find((x) => x.c === 1);
  check('lint: finds the parent character and its "absent" state', c2 && c2.parent && c2.parent.c === 0 && c2.parent.state === '0', c2 && c2.parent);
  const kinds = Object.fromEntries(c2.issues.map((i) => [mx.taxa[i.t], i.kind]));
  check('lint: scored although the parent is absent (h)', kinds.h === 'scored-parent-absent', kinds);
  check('lint: gap although the parent is present (i)', kinds.i === 'gap-parent-present', kinds);
  check('lint: gap although the parent is unknown (j)', kinds.j === 'gap-parent-unknown', kinds);
  check('lint: consistent cells are not flagged', Object.keys(kinds).length === 3, kinds);
  const miss = P.parseMatrix(text.replace('a 0-0', 'a 0?0'));
  const k2 = P.lintInapplicable(miss).chars.find((x) => x.c === 1).issues.find((i) => i.t === 0);
  check('lint: ? where the parent says absent', k2 && k2.kind === 'missing-parent-absent', k2);
  check('lint: a matrix without gaps has nothing to report', P.lintInapplicable(P.parseMatrix(text.replace(/-/g, '?'))).withGaps === 0);
  if (D) {
    let chars = 0, withParent = 0, issues = 0;
    for (const f of fs.readdirSync(D).filter((x) => x.endsWith('.nex'))) { const r = P.lintInapplicable(P.parseMatrix(fs.readFileSync(D + '/' + f, 'utf8'))); chars += r.withGaps; withParent += r.withParent; for (const k in r.counts) issues += r.counts[k]; }
    check('lint: runs on all 31 matrices', chars > 0 && withParent > 0);
    console.log(`  coding check on the 31 matrices: ${chars} characters with gaps, a controlling character found for ${withParent}, ${issues} cells flagged`);
  }
}
console.log(`${pass}/${pass + fail} feature checks passed`);
process.exit(fail ? 1 : 0);
