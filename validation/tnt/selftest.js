// Checks compare.js without TNT: OpenClade writes tree files in the format the TNT
// scripts ask for (NEXUS trees, taxa t1..tN, unsupported branches collapsed), and
// compare.js must then report the same length and an identical strict consensus.
//   DATASETS=path/to/TreeSearch/inst/datasets node selftest.js
const fs = require('fs');
const os = require('os');
const path = require('path');
const { execFileSync } = require('child_process');
const P = require('../../src/parser.js');
const E = require('../../src/engine.js');
const D = process.env.DATASETS;
const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'oc-tnt-'));
const pick = ['Vinther2008', 'Longrich2010', 'DeAssis2011'];
for (const name of pick) {
  const mx = P.parseMatrix(fs.readFileSync(path.join(D, name + '.nex'), 'utf8'));
  const n = mx.taxa.length, labels = mx.taxa.map((_, i) => 't' + (i + 1));
  const data = E.encode(mx, 'missing');
  const res = E.analyse(data, { hits: 5, maxReps: 30, ratchetIter: 8, maxTrees: 5000, seed: 7, timeLimit: 120000 });
  const col = E.collapseAll(data, res.trees);
  const strip = (x) => { if (x.children) { x.support = null; x.children.forEach(strip); } return x; };
  fs.writeFileSync(path.join(tmp, name + '.tre'), '#NEXUS\nbegin trees;\n' + col.map((c, i) => `tree tnt_${i} = [&U] ${E.nestedToNewick(strip(c.node), labels)}`).join('\n') + '\nend;\n');
}
const out = execFileSync('node', [path.join(__dirname, 'compare.js'), '120'], { env: Object.assign({}, process.env, { TREES: tmp, CSVOUT: path.join(tmp, 'selftest.csv') }), encoding: 'utf8' });
const lines = out.split('\n').filter((l) => pick.some((p) => l.startsWith(p)));
let ok = lines.length === pick.length;
for (const l of lines) { console.log(l); if (!/same length/.test(l) || !/strict consensus identical/.test(l)) ok = false; }
console.log(ok ? 'compare.js self-test passed (3 matrices; no TNT involved)' : 'SELF-TEST FAILED');
process.exit(ok ? 0 : 1);
