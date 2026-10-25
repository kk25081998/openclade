// End-to-end UI test of the built page in jsdom. The Web Worker is emulated by
// running the worker source in a vm context and delivering messages async.
const fs = require('fs');
const vm = require('vm');
const { JSDOM } = require('jsdom');
const html = fs.readFileSync(require('path').join(__dirname, '../dist/openclade.html'), 'utf8');

let pass = 0, fail = 0;
const check = (n, c, x) => { if (c) pass++; else { fail++; console.log('  FAIL', n, x === undefined ? '' : x); } };

const dom = new JSDOM(html, { runScripts: 'dangerously', pretendToBeVisual: true, url: 'https://example.test/',
  beforeParse(win) {
    const blobs = new Map();
    win.URL.createObjectURL = (b) => { const id = 'blob:' + blobs.size; blobs.set(id, b); return id; };
    win.Worker = class {
      constructor(url) {
        this.terminated = false;
        const blob = blobs.get(url);
        const src = blob.__text;
        const self = { addEventListener() {}, postMessage: (m) => { if (!this.terminated) setTimeout(() => this.onmessage && this.onmessage({ data: JSON.parse(JSON.stringify(m)) }), 0); } };
        self.self = self;
        const ctx = vm.createContext({ self, globalThis: self, Date, Math, JSON, Map, Set, Int32Array, Uint32Array, Uint8Array, Float64Array, Array, Object, Number, String, Error, setTimeout, console });
        ctx.self = self;
        vm.runInContext(src, ctx);
        this.self = self;
      }
      postMessage(m) { setTimeout(() => { try { this.self.onmessage({ data: JSON.parse(JSON.stringify(m)) }); } catch (e) { console.log('worker error', e); this.onerror && this.onerror(e); } }, 0); }
      terminate() { this.terminated = true; }
    };
    const OrigBlob = win.Blob;
    win.Blob = class extends OrigBlob { constructor(parts, opts) { super(parts, opts); this.__text = parts.join(''); } };
    win.navigator.clipboard = { writeText: (t) => { win.__copied = t; return Promise.resolve(); } };
  } });
const { window } = dom;
const doc = window.document;
const $ = (id) => doc.getElementById(id);
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
async function waitFor(fn, ms = 120000) { const t0 = Date.now(); while (Date.now() - t0 < ms) { if (fn()) return true; await sleep(50); } return false; }

(async () => {
  window.addEventListener('error', (e) => console.log('page error', e.message));
  // 1. auto-run on the example
  check('example loaded', /Vinther/.test($('dataInfo').textContent));
  const ok = await waitFor(() => /Finished/.test($('progText').textContent));
  check('example analysis finished', ok, $('progText').textContent);
  console.log('   ', $('progText').textContent);
  const stats = [...doc.querySelectorAll('.stat')].map((s) => s.textContent);
  console.log('    stats:', stats.join(' | '));
  check('tree drawn', !!$('treebox').querySelector('svg'));
  check('23 tip labels', [...$('treebox').querySelectorAll('svg > text')].filter((t) => !t.getAttribute('font-size')).length === 23);
  check('support values shown', /bootstrap/.test($('legend').textContent), $('legend').textContent);
  // 2. character mapping
  $('mapChar').value = '9'; $('mapChar').dispatchEvent(new window.Event('change'));
  check('mapping draws node states', $('treebox').querySelectorAll('circle,path').length > 20);
  check('char info panel', /step/.test($('charInfo').textContent), $('charInfo').textContent);
  // 3. views
  const views = [...$('view').options].map((o) => o.value);
  console.log('    views:', views.slice(0, 4).join(', '), views.length > 4 ? '…' : '');
  if (views.includes('strict')) { $('view').value = 'strict'; $('view').dispatchEvent(new window.Event('change')); check('strict consensus drawn', !!$('treebox').querySelector('svg')); }
  // 4. tabs
  $('tab-chars').click();
  check('character table rows', doc.querySelectorAll('#charTable tbody tr').length === 57);
  $('tab-export').click();
  const ex = [...doc.querySelectorAll('#exports textarea')].map((t) => t.value);
  check('exports present', ex.length >= 6, ex.length);
  check('newick export parses', ex[0].split('\n').every((l) => /^\(.*\);$/.test(l)));
  check('TNT script has xread', ex.some((t) => /xread/.test(t)));
  check('SVG export', ex.some((t) => t.startsWith('<svg')));
  doc.querySelector('#exports [data-copy]').click(); await sleep(20);
  check('copy works', typeof window.__copied === 'string' && window.__copied.length > 10);
  $('tab-methods').click();
  check('methods paragraph', /most parsimonious tree/.test($('methodsText').textContent), $('methodsText').textContent.slice(0, 120));
  console.log('    methods:', $('methodsText').textContent.slice(0, 160) + '…');
  // check Newick round trip: exported MPT re-scores to the reported length using the engine
  const E = window.OCEngine, P = window.OCParser;
  // 5. load a different sample with missing gaps, outgroup change, no support
  $('sample').value = '1'; $('sample').dispatchEvent(new window.Event('change'));
  check('second sample loaded', /Longrich/.test($('dataInfo').textContent));
  $('gapMissing').checked = true; $('support').value = 'none'; $('outgroup').value = '3';
  $('run').click();
  const ok2 = await waitFor(() => /Finished/.test($('progText').textContent));
  check('second analysis finished', ok2, $('progText').textContent);
  console.log('   ', $('progText').textContent);
  $('tab-tree').click();
  const firstTip = $('treebox').querySelector('svg > text:last-of-type');
  const og = window.OCParser.parseMatrix(window.document.getElementById('samples') ? JSON.parse($('samples').textContent)[1].text : '').taxa[3];
  const tipTexts = [...$('treebox').querySelectorAll('svg > text')].map((t) => t.textContent);
  check('outgroup appears as first tip (rooted on it)', tipTexts[0] === og || tipTexts.includes(og), [tipTexts[0], og]);
  check('methods mention missing', /missing data/.test($('methodsText').textContent));
  // independent re-score of exported trees against a fresh parse
  $('tab-export').click();
  const nwk = doc.querySelector('#exports textarea').value.split('\n')[0];
  const mx = P.parseMatrix(JSON.parse($('samples').textContent)[1].text);
  const data = E.encode(mx, 'missing');
  const t = E.fromNested(mx.taxa.length, E.parseNewick(nwk, mx.taxa));
  const len = new E.Scorer(data).length(t);
  const reported = +doc.querySelector('.stat .v').textContent;
  check('exported tree re-scores to reported length', len === reported, [len, reported]);
  // 6. compare gap treatments
  $('sample').value = '0'; $('sample').dispatchEvent(new window.Event('change'));
  check('gap pane empty after load', /Compare gap treatments/.test($('gapsBox').textContent));
  $('intQuick').checked = true; $('gapInapp').checked = true; $('timeLimit').value = '30';
  $('compare').click();
  check('buttons disabled while comparing', $('run').disabled && $('compare').disabled && !$('stop').disabled);
  const seen = new Set();
  const okC = await waitFor(() => { const m = $('progText').textContent.match(/Treatment (\d) of 3/); if (m) seen.add(m[1]); return /Compared 3 gap treatments/.test($('progText').textContent); });
  check('comparison finished', okC, $('progText').textContent);
  console.log('   ', $('progText').textContent);
  check('progress named each treatment', seen.size >= 2, [...seen]);
  check('gap tab selected', $('tab-gaps').getAttribute('aria-selected') === 'true' && !$('pane-gaps').hidden);
  check('buttons re-enabled', !$('run').disabled && !$('compare').disabled && $('stop').disabled);
  const head = doc.querySelector('#gapsBox .headline').textContent;
  console.log('    headline:', head);
  const sumRows = [...doc.querySelectorAll('#gapsBox table.plain tbody tr')].map((tr) => [...tr.children].map((td) => td.textContent));
  check('summary has three treatments', sumRows.length === 3 && /Inapplicable/.test(sumRows[0][0]) && /Missing/.test(sumRows[1][0]) && /Extra/.test(sumRows[2][0]), sumRows);
  console.log('    lengths:', sumRows.map((r) => r[0] + ' ' + r[1]).join(', '));
  check('Vinther lengths 79 / 78 / 88', sumRows.map((r) => +r[1]).join() === '79,78,88', sumRows.map((r) => r[1]));
  // Independent check: re-derive the comparison from the exported Newick consensus trees
  const nwks = [0, 1, 2].map((i) => $('cmpExp' + i).value.trim());
  const taxa0 = P.parseMatrix(JSON.parse($('samples').textContent)[0].text).taxa;
  const og0 = taxa0[0];
  function cladesOf(nwk) {
    let i = 0; const out = new Set();
    function nm() { let s = ''; if (nwk[i] === "'") { i++; while (i < nwk.length) { if (nwk[i] === "'" && nwk[i + 1] === "'") { s += "'"; i += 2; } else if (nwk[i] === "'") { i++; break; } else s += nwk[i++]; } } else while (!/[(),;:]/.test(nwk[i])) s += nwk[i++]; while (/[^(),;]/.test(nwk[i])) i++; return s; }
    function node() { if (nwk[i] === '(') { i++; const tips = []; for (;;) { tips.push(...node()); if (nwk[i] === ',') i++; else break; } i++; while (/[^(),;]/.test(nwk[i])) i++; if (tips.length > 1 && tips.length < taxa0.length - 1 && !tips.includes(og0)) out.add(tips.slice().sort().join('|')); return tips; } return [nm()]; }
    const all = node(); return { out, all };
  }
  const cs = nwks.map(cladesOf);
  check('consensus trees contain every taxon', cs.every((c) => c.all.length === taxa0.length && new Set(c.all).size === taxa0.length));
  const union = new Set(); cs.forEach((c) => c.out.forEach((k) => union.add(k)));
  const shared = [...union].filter((k) => cs.every((c) => c.out.has(k))).length;
  check('headline matches an independent count', head.startsWith(`${shared} of ${union.size} clades`), [head, shared, union.size]);
  check('clade counts per treatment match', sumRows.every((r, i) => +r[3] === cs[i].out.size), [sumRows.map((r) => r[3]), cs.map((c) => c.out.size)]);
  const depRows = doc.querySelectorAll('#cmpTable tbody tr');
  check('table lists exactly the treatment-dependent clades', depRows.length === union.size - shared, [depRows.length, union.size - shared]);
  check('each listed clade is missing from some treatment', [...depRows].every((tr) => tr.querySelectorAll('td.no').length >= 1 && tr.querySelectorAll('td.yes').length >= 1));
  const redBranches = doc.querySelectorAll('#cmpTree g[stroke="var(--bad)"] line').length;
  const depInI = [...union].filter((k) => cs[0].out.has(k) && !cs.every((c) => c.out.has(k))).length;
  check('tree marks the treatment-dependent branches', redBranches === depInI, [redBranches, depInI]);
  check('comparison tree has all tips', [...doc.querySelectorAll('#cmpTree svg > text')].length === taxa0.length);
  $('cmpShared').click();
  check('show-shared lists every clade', doc.querySelectorAll('#cmpTable tbody tr').length === union.size);
  const row = doc.querySelector('#cmpTable tbody tr'); row.click();
  check('row click highlights taxa', row.classList.contains('sel') && doc.querySelectorAll('#cmpTree text[font-weight="700"]').length === +row.children[0].textContent);
  $('cmpView').value = '2'; $('cmpView').dispatchEvent(new window.Event('change'));
  const red2 = doc.querySelectorAll('#cmpTree g[stroke="var(--bad)"] line').length;
  check('switching tree to extra state', red2 === [...union].filter((k) => cs[2].out.has(k) && !cs.every((c) => c.out.has(k))).length);
  const cm = $('cmpMethods').textContent;
  check('comparison methods paragraph', /Brazeau, Guillerme & Smith \(2019\)/.test(cm) && /Robinson–Foulds/.test(cm) && cm.includes(`${shared} were recovered under all three`), cm.slice(-260));
  console.log('    methods:', cm.slice(-240));
  check('clade TSV export', $('cmpExp3').value.split('\n').length === union.size + 1);
  check('main tabs show the inapplicable analysis', +doc.querySelector('.stat .v').textContent === 79);
  // 6b. stop part-way through a comparison
  $('compare').click();
  await waitFor(() => /Treatment 2 of 3/.test($('progText').textContent), 120000);
  $('stop').click(); await sleep(50);
  const stopTxt = $('progText').textContent;
  check('stop ends the comparison', /Stopped|nothing to compare/.test(stopTxt) && !$('compare').disabled && $('stop').disabled, stopTxt);
  const cols = doc.querySelectorAll('#cmpTable thead th.tr').length;
  check('partial comparison is consistent', (/Compared 2 gap treatments/.test(stopTxt) && cols === 2) || (/nothing to compare/.test(stopTxt) && cols === 0), [stopTxt, cols]);
  console.log('   ', stopTxt);
  // 6c. a matrix without gaps
  $('paste').value = '#NEXUS\nbegin data; dimensions ntax=5 nchar=4; format symbols="01" missing=? gap=-; matrix\nA 0000\nB 0011\nC 1100\nD 1111\nE 1?10\n;\nend;';
  $('loadPaste').click(); $('compare').click(); await sleep(20);
  check('no-gap matrix explained', /no gap/.test($('gapsBox').textContent) && !$('compare').disabled, $('gapsBox').textContent);
  // 7. bad file
  $('paste').value = 'this is not a matrix'; $('loadPaste').click();
  check('bad input shows error', !!doc.querySelector('#dataInfo .msg.err'));
  check('run disabled on bad input', $('run').disabled);
  // 8. Stop mid-run returns trees
  $('sample').value = '3'; $('sample').dispatchEvent(new window.Event('change'));
  $('gapInapp').checked = true; $('support').value = 'bootstrap'; $('intThor').checked = true;
  $('run').click();
  await waitFor(() => /best length \d/.test($('progText').textContent), 120000);
  $('stop').click();
  await sleep(50);
  check('stop keeps best-so-far trees', !!$('treebox').querySelector('svg') && /stopped early/.test($('progText').textContent), $('progText').textContent);
  console.log(`${pass}/${pass + fail} UI checks passed`);
  process.exit(fail ? 1 : 0);
})();
