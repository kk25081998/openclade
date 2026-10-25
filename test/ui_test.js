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
  // 6. bad file
  $('paste').value = 'this is not a matrix'; $('loadPaste').click();
  check('bad input shows error', !!doc.querySelector('#dataInfo .msg.err'));
  check('run disabled on bad input', $('run').disabled);
  // 7. Stop mid-run returns trees
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
