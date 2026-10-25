// Checks the user guide: every link from the app into the guide, and every link
// inside the guide, lands on a section that exists; every image it uses was built;
// and the guide's worked example matches what the engine calculates.
const fs = require('fs');
const path = require('path');
const { JSDOM } = require('jsdom');
const P = require('../src/parser.js');
const E = require('../src/engine.js');
const DIST = path.join(__dirname, '..', 'dist');

let pass = 0, fail = 0;
const check = (n, c, x) => { if (c) pass++; else { fail++; console.log('  FAIL', n, x === undefined ? '' : JSON.stringify(x)); } };

const docsHtml = fs.readFileSync(path.join(DIST, 'docs', 'index.html'), 'utf8');
const doc = new JSDOM(docsHtml).window.document;
const ids = new Set([...doc.querySelectorAll('[id]')].map((e) => e.id));

// 1. Links inside the guide
const inner = [...doc.querySelectorAll('a[href^="#"]')].map((a) => a.getAttribute('href').slice(1));
const brokenInner = inner.filter((id) => !ids.has(id));
check('guide: every in-page link has a target', brokenInner.length === 0, brokenInner);
check('guide: has the sections the app links to', ['start', 'matrix', 'gaps', 'outgroup', 'search', 'more-settings', 'support', 'tree', 'characters', 'compare', 'export', 'cite', 'trouble'].every((id) => ids.has(id)));

// 2. Links from the app (both builds) into the guide
for (const [file, prefix] of [['index.html', 'docs/'], ['openclade.html', 'https://kartikkhosa.com/openclade/docs/']]) {
  const app = new JSDOM(fs.readFileSync(path.join(DIST, file), 'utf8')).window.document;
  const links = [...app.querySelectorAll('a[href]')].map((a) => a.getAttribute('href')).filter((h) => h.startsWith(prefix));
  const targets = links.map((h) => h.slice(prefix.length)).filter((h) => h.startsWith('#')).map((h) => h.slice(1));
  const broken = targets.filter((id) => !ids.has(id));
  check(`${file}: links to the guide`, links.length >= 10, links.length);
  check(`${file}: every guide link has a target`, broken.length === 0, broken);
  check(`${file}: guide links open in a new tab`, [...app.querySelectorAll('a[href]')].filter((a) => a.getAttribute('href').startsWith(prefix)).every((a) => a.target === '_blank'));
  check(`${file}: no unfilled placeholders`, !/\{\{(APP|DOCS)\}\}/.test(app.documentElement.outerHTML));
}
// The app's script builds a few guide links at run time from the top-bar link
const appSrc = fs.readFileSync(path.join(__dirname, '..', 'src', 'app.html'), 'utf8');
const jsTargets = [...appSrc.matchAll(/\$\{DOCS\}#([a-z-]+)/g)].map((m) => m[1]);
check('app script: guide links have targets', jsTargets.length > 0 && jsTargets.every((id) => ids.has(id)), jsTargets);

// 3. Images
const imgs = [...doc.querySelectorAll('img[src], source[srcset]')].map((e) => e.getAttribute('src') || e.getAttribute('srcset'));
const missing = imgs.filter((s) => !fs.existsSync(path.join(DIST, 'docs', s)));
check('guide: every image was built', imgs.length >= 3 && missing.length === 0, missing);
check('guide: every image has alt text', [...doc.querySelectorAll('img')].every((i) => (i.getAttribute('alt') || '').length > 20));
check('guide: no unfilled placeholders', !/\{\{(APP|DOCS)\}\}/.test(docsHtml));

// 4. The worked example: the lengths printed in the guide are what the engine gives
const rows = [...doc.querySelectorAll('#gaps ~ * table, section[aria-labelledby="gaps"] table')].pop().querySelectorAll('tbody tr');
const printed = [...rows].map((tr) => [...tr.children].slice(1, 3).map((td) => td.textContent.trim()).join('/')).join(' ');
const mx = P.parseMatrix('#NEXUS\nbegin data; dimensions ntax=7 nchar=3; format symbols="01" missing=? gap=-;\nmatrix\nOutgroup 0-0\nA 110\nB 110\nC 0-1\nD 0-1\nE 111\nF 111\n;\nend;');
const X = '(Outgroup,((C,D),((A,B),(E,F))));', Y = '(Outgroup,((A,B),(C,(D,(E,F)))));';
const len = (mode, nwk) => new E.Scorer(E.encode(mx, mode)).length(E.fromNested(7, E.parseNewick(nwk, mx.taxa)));
const computed = ['inapplicable', 'missing', 'extra'].map((m) => len(m, X) + '/' + len(m, Y)).join(' ');
check('guide example: printed lengths match the engine', printed === computed, [printed, computed]);
const coded = [...doc.querySelectorAll('.example td.c')].map((td) => td.textContent.replace(/\s/g, ''));
check('guide example: coded rows match the test matrix', coded.join(',') === '0-0,110,0-1,111', coded);

console.log(`${pass}/${pass + fail} guide checks passed`);
process.exit(fail ? 1 : 0);
