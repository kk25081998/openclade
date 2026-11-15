// End-to-end test of the features added in v0.2 (implied weights, Bremer support, constraints,
// outgroups, synapomorphies, data check, your own trees, saved analyses, robustness, PDF figure),
// on the built page in jsdom. The Web Worker is emulated by
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
    win.URL.revokeObjectURL = () => {};
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
    // Downloads: record the file name and the contents of the blob it points to
    win.__downloads = [];
    win.HTMLAnchorElement.prototype.click = function () { const b = blobs.get(this.href); win.__downloads.push({ name: this.download, text: b && b.__text, type: b && b.type }); };
  } });
const { window } = dom;
const doc = window.document;
const $ = (id) => doc.getElementById(id);
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
async function waitFor(fn, ms = 120000) { const t0 = Date.now(); while (Date.now() - t0 < ms) { if (fn()) return true; await sleep(50); } return false; }

(async () => {
  window.addEventListener('error', (e) => console.log('page error', e.message));
  const E = window.OCEngine, P = window.OCParser;
  const stat = (label) => { const el = [...doc.querySelectorAll('.stat')].find((x) => x.querySelector('.k').textContent.startsWith(label)); return el ? el.querySelector('.v').textContent : null; };
  const run = async () => { $('run').click(); const ok = await waitFor(() => /Finished|Not started|error/.test($('progText').textContent)); await sleep(60); return ok; };
  const fire = (id, type) => $(id).dispatchEvent(new window.Event(type || 'change'));
  const exportsText = () => { $('tab-export').click(); return Object.fromEntries([...doc.querySelectorAll('#exports .export')].map((x) => [x.querySelector('.lbl').textContent, x.querySelector('textarea').value])); };
  await waitFor(() => /Finished/.test($('progText').textContent));
  const taxa = P.parseMatrix(JSON.parse($('samples').textContent)[0].text).taxa;

  // ---- 1. implied weighting + Bremer support
  $('wtImplied').checked = true; fire('wtImplied');
  check('k field appears for implied weights', !$('kRow').hidden && $('concavity').value === '12');
  $('support').value = 'none'; $('bremer').checked = true; $('intQuick').checked = true;
  check('implied-weights run finished', await run(), $('progText').textContent);
  await waitFor(() => /Bremer/.test($('legend').textContent));
  const fit = stat('Fit (implied weights');
  check('fit is shown with 5 decimals', /^\d+\.\d{5}$/.test(fit || ''), fit);
  check('tree length in steps is shown too', +stat('Tree length') >= 79, stat('Tree length'));
  check('replicate count shown', /^\d+ of \d+$/.test(stat('Found in replicates') || ''), stat('Found in replicates'));
  // independent re-score: exported first tree under implied weights, k = 12
  let ex = exportsText();
  const mx0 = P.parseMatrix(JSON.parse($('samples').textContent)[0].text);
  const dIW = E.encode(mx0, 'inapplicable', { concavity: 12 });
  const firstTree = E.fromNested(mx0.taxa.length, E.parseNewick(ex['Shortest trees (Newick)'].split('\n')[0].replace(/;$/, ''), mx0.taxa));
  check('exported tree re-scores to the reported fit', (new E.Scorer(dIW).length(firstTree) / dIW.scoreDiv).toFixed(5) === fit, [new E.Scorer(dIW).length(firstTree) / dIW.scoreDiv, fit]);
  check('methods name implied weighting and Bremer', /implied weighting \(Goloboff, 1993\) with concavity constant k = 12/.test($('methodsText').textContent) && /Bremer support \(Bremer, 1994\)/.test($('methodsText').textContent));
  check('Bremer export has a value on every clade', /\)\d+\.\d\d/.test(ex['Strict consensus with Bremer support (Newick)'] || ''), (ex['Strict consensus with Bremer support (Newick)'] || '').slice(0, 120));
  $('tab-tree').click(); $('view').value = 'strict'; fire('view');
  check('Bremer values drawn below branches', [...$('treebox').querySelectorAll('text')].some((t) => /^\d+\.\d\d$/.test(t.textContent)) && /Bremer support \(difference in fit\)/.test($('legend').textContent), $('legend').textContent);
  $('tab-chars').click();
  check('character table has a cost column under implied weights', /Cost/.test(doc.querySelector('#charTable thead').textContent));
  $('wtEqual').checked = true; fire('wtEqual'); $('bremer').checked = false;

  // ---- 2. equal weights with Bremer, exact check on integers
  $('bremer').checked = true;
  check('equal-weights run finished', await run());
  await waitFor(() => /Bremer/.test($('legend').textContent) || /Bremer/.test($('resWhat').textContent));
  ex = exportsText();
  const bre = ex['Strict consensus with Bremer support (Newick)'];
  check('Bremer values are whole steps of at least 1', /\)\d+[,)]/.test(bre) && !/\)0[,)]/.test(bre), bre.slice(0, 160));
  $('bremer').checked = false;

  // ---- 3. synapomorphies
  $('tab-syn').click();
  check('opening the tab turns node numbers on', $('nodeNums').checked);
  const synRows = () => doc.querySelectorAll('#synTable tbody tr').length;
  check('strict consensus: changes shared by all trees', synRows() > 5 && /shared by all 52 trees/.test($('synNote').textContent), $('synNote').textContent);
  const sharedCount = +$('synNote').textContent.match(/(\d+) changes?/)[1];
  $('tab-tree').click(); $('view').value = 't0'; fire('view'); $('tab-syn').click();
  const t0count = +$('synNote').textContent.match(/(\d+) changes?/)[1];
  check('single tree lists at least as many changes as are shared by all', t0count >= sharedCount && /on tree 1/.test($('synNote').textContent), [t0count, sharedCount]);
  $('synMode').value = 'acctran'; fire('synMode');
  const acc = +$('synNote').textContent.match(/(\d+) changes?/)[1];
  check('ACCTRAN lists more changes than unambiguous only', acc > t0count && /ACCTRAN/.test($('synNote').textContent), [acc, t0count]);
  $('synTips').checked = true; fire('synTips');
  check('terminal branches can be included', synRows() > 23);
  check('synapomorphy TSV export', $('synTsv').value.split('\n').length === synRows() + 1 && /^node\ttaxa/.test($('synTsv').value));
  $('tab-tree').click();
  check('node numbers drawn on the tree', $('treebox').querySelectorAll('text[font-weight="700"]').length >= 20);
  $('synMode').value = 'unambiguous'; $('synTips').checked = false;

  // ---- 4. data check
  $('tab-check').click();
  check('data check reports the coding of gaps', /Coding of gaps/.test($('checkBox').textContent) && /characters have gaps/.test($('checkBox').textContent), $('checkBox').textContent.slice(0, 200));
  check('data check lists file settings', /Settings in the file/.test($('checkBox').textContent));

  // ---- 5. constraints
  const d0 = E.encode(mx0, 'inapplicable');
  const a = taxa[5], b = taxa[20];
  $('constraints').value = `${a}, No such taxon`;
  await run();
  check('unknown taxon in a constraint stops the run', /Not started/.test($('progText').textContent) && /No such taxon/.test($('conMsg').textContent), $('progText').textContent);
  $('constraints').value = `${a}, ${b}`;
  check('constrained run finished', await run() && /Finished/.test($('progText').textContent), $('progText').textContent);
  ex = exportsText();
  const conTrees = ex['Shortest trees (Newick)'].split('\n').map((l) => E.rootOn(E.fromNested(23, E.parseNewick(l.replace(/;$/, ''), taxa)), 0));
  const kw = new Uint32Array(1); kw[0] = (1 << 5) | (1 << 20);
  check('every tree has the constrained group', conTrees.every((t) => E.clades(t).includes(String(kw[0]))), conTrees.length);
  check('the constraint costs steps', +stat('Tree length') > 79, stat('Tree length'));
  check('constrained trees re-score to the reported length', conTrees.every((t) => new E.Scorer(d0).length(t) === +stat('Tree length')));
  check('methods mention the constraint', /1 topological constraint/.test($('methodsText').textContent));
  check('MrBayes file carries the constraint', /constraint group1 = /.test(ex['MrBayes file for a Bayesian analysis of the same matrix (Mk model)']) && /prset topologypr = constraints\(group1\)/.test(ex['MrBayes file for a Bayesian analysis of the same matrix (Mk model)']));
  // ordered characters together with a constraint: files for TreeSearch and MrBayes (run for real in validation)
  const multi = []; for (let c = 0; c < mx0.nchar && multi.length < 3; c++) { const st = new Set(); let gaps = 0; for (const r of mx0.cells) { if (Array.isArray(r[c])) r[c].forEach((x) => st.add(x)); else if (r[c] === '-') gaps++; } if (st.size >= 3 && gaps < 3) multi.push(c + 1); }
  $('ordered').value = multi.join(' '); fire('ordered', 'input');
  await run();
  ex = exportsText();
  check('ordered characters reach the exports', multi.length >= 2 && new RegExp('ordered <- c\\(' + multi.join(', ') + '\\)').test(ex['R script to re-score these trees in TreeSearch']) && new RegExp('ctype ordered: ' + multi.join(' ') + ';').test(ex['MrBayes file for a Bayesian analysis of the same matrix (Mk model)']) && new RegExp('cc \\+ ' + multi.map((x) => x - 1).join(' ') + ';').test(ex['TNT script to repeat the search in TNT']), multi);
  fs.writeFileSync(require('path').join(__dirname, '../dist/example-ordered-treesearch.R'), ex['R script to re-score these trees in TreeSearch']);
  fs.writeFileSync(require('path').join(__dirname, '../dist/example-ordered-mrbayes.nex'), ex['MrBayes file for a Bayesian analysis of the same matrix (Mk model)']);
  fs.writeFileSync(require('path').join(__dirname, '../dist/example-ordered.tnt'), ex['TNT script to repeat the search in TNT']);
  $('ordered').value = ''; fire('ordered', 'input');
  $('constraints').value = `not: ${taxa[1]}, ${taxa[2]}`;
  await run();
  check('converse constraint accepted', /Finished/.test($('progText').textContent) && /required not to be monophyletic/.test($('methodsText').textContent), $('progText').textContent);
  $('constraints').value = '';

  // ---- 6. several outgroup taxa
  $('ogDetails').open = true;
  for (const o of $('outgroupMore').options) o.selected = o.textContent === taxa[1] || o.textContent === taxa[2];
  $('ingroupMono').checked = true;
  check('multi-outgroup run finished', await run() && /Finished/.test($('progText').textContent), $('progText').textContent);
  check('header names the outgroup taxa', /rooted on 3 outgroup taxa/.test($('resWhat').textContent), $('resWhat').textContent);
  check('methods describe the outgroup and the ingroup constraint', /rooted on the outgroup \(/.test($('methodsText').textContent) && /ingroup was constrained to be monophyletic/.test($('methodsText').textContent));
  $('tab-tree').click(); $('view').value = 't0'; fire('view');
  const tips = [...$('treebox').querySelectorAll('svg > text')].filter((t) => !t.getAttribute('font-size')).map((t) => t.textContent);
  check('the three outgroup taxa are drawn together at the base', new Set(tips.slice(0, 3)).size === 3 && [taxa[0], taxa[1], taxa[2]].every((x) => tips.slice(0, 3).includes(x)), tips.slice(0, 4));
  check('no warning that the outgroup is split', !/do not form a group/.test($('legend').textContent));
  $('ogDetails').open = false; $('ingroupMono').checked = false;

  // ---- 7. collapsing
  $('collapse').checked = true;
  await run();
  check('collapsed trees are listed', [...$('view').options].some((o) => /collapsed/.test(o.textContent)) && +stat('After collapsing') < +stat('Shortest trees found'), [stat('After collapsing'), stat('Shortest trees found')]);
  $('tab-tree').click(); $('view').value = 'c0'; fire('view');
  check('a collapsed tree is drawn with a note', /no character supports are collapsed/.test($('legend').textContent));
  $('mapChar').value = '3'; fire('mapChar');
  check('mapping a character on a collapsed tree draws its resolved tree', $('treebox').querySelectorAll('circle,path').length > 20);
  $('mapChar').value = '-1'; fire('mapChar');
  $('collapse').checked = false;

  // ---- 8. PDF figure, R script, MrBayes file
  window.__downloads.length = 0;
  $('figPdf').click();
  const pdf = (window.__downloads[0] || {}).text || '';
  check('PDF downloads', /\.pdf$/.test((window.__downloads[0] || {}).name || '') && pdf.startsWith('%PDF-1.4') && /%%EOF\n$/.test(pdf));
  const xrefAt = +pdf.match(/startxref\n(\d+)/)[1];
  const offs = [...pdf.slice(xrefAt).matchAll(/(\d{10}) 00000 n /g)].map((m) => +m[1]);
  check('PDF cross-reference table points at every object', pdf.slice(xrefAt, xrefAt + 4) === 'xref' && offs.length === 7 && offs.every((o, i) => pdf.slice(o, o + 7) === `${i + 1} 0 obj`), offs);
  const len = +pdf.match(/\/Length (\d+)/)[1], st = pdf.indexOf('stream\n') + 7;
  check('PDF stream length is right', pdf.slice(st + len, st + len + 10) === '\nendstream', pdf.slice(st + len, st + len + 12));
  check('PDF is plain ASCII with every taxon', !/[^\x00-\x7f]/.test(pdf) && taxa.every((t) => pdf.includes('(' + t.replace(/[()\\]/g, '\\$&') + ') Tj')));
  $('figSvg').click();
  check('SVG downloads', /\.svg$/.test(window.__downloads[1].name) && window.__downloads[1].text.startsWith('<svg'));
  await run();
  ex = exportsText();
  const R = ex['R script to re-score these trees in TreeSearch'];
  check('R script re-scores the trees', /TreeLength\(tr, dat, concavity = Inf\)/.test(R) && /openclade <- c\(79(, 79)*\)/.test(R) && !/NOTE:/.test(R), R.slice(0, 300));
  fs.writeFileSync(require('path').join(__dirname, '../dist/example-treesearch.R'), R);
  const MB = ex['MrBayes file for a Bayesian analysis of the same matrix (Mk model)'];
  check('MrBayes file', /begin mrbayes;/.test(MB) && /lset coding=variable rates=gamma;/.test(MB) && /dimensions ntax=23 nchar=57;/.test(MB));
  check('MrBayes file is plain ASCII (MrBayes stops at other characters in comments)', !/[^\x00-\x7f]/.test(MB));
  fs.writeFileSync(require('path').join(__dirname, '../dist/example-mrbayes.nex'), MB);

  // ---- 9. saved analysis: open it again and reproduce it
  const saved = ex['Analysis file: matrix, settings and result, to open in OpenClade and run again'];
  const rec = JSON.parse(saved);
  check('analysis file holds matrix, settings and result', rec.format === 'openclade-analysis' && /#NEXUS/i.test(rec.matrix.text) && rec.form.seed === '1' && rec.result.score === 79 && rec.result.trees === 52 && rec.result.treesNewick.length === 52, [rec.result.score, rec.result.trees]);
  $('sample').value = '1'; fire('sample');
  $('seed').value = '9'; $('gapMissing').checked = true;
  $('paste').value = saved; $('loadPaste').click();
  check('opening the file restores the matrix and settings', /Vinther/.test($('dataInfo').textContent) && $('seed').value === '1' && $('gapInapp').checked && /Saved analysis/.test($('savedInfo').textContent) && /52 trees, length 79/.test($('savedInfo').textContent), $('savedInfo').textContent);
  $('rerun').click();
  await waitFor(() => !!$('savedVerdict'));
  check('running it again reproduces it', /Reproduced: the same length \(79\) and the same 52 trees/.test($('savedVerdict').textContent), $('savedVerdict') && $('savedVerdict').textContent);
  // the same file on a different number of cores
  $('cores').value = $('cores').value === '1' ? '2' : '1';
  $('rerun').click(); await sleep(100);
  await waitFor(() => /Finished/.test($('progText').textContent) && !!$('savedVerdict'));
  check('and again on a different number of cores', /Reproduced/.test($('savedVerdict').textContent), $('savedVerdict').textContent);
  // a tampered file is noticed
  const bad = JSON.parse(saved); bad.result.score = 78; bad.result.scoreShown = '78';
  $('paste').value = JSON.stringify(bad); $('loadPaste').click(); $('rerun').click();
  await waitFor(() => !!$('savedVerdict'));
  check('a file whose result does not match is reported', /Not reproduced: length 79 now, 78 in the file/.test($('savedVerdict').textContent), $('savedVerdict').textContent);
  $('paste').value = '{"format":"openclade-analysis"}'; $('loadPaste').click();
  check('a broken analysis file is refused', /does not look like an OpenClade analysis file/.test($('dataInfo').textContent));

  // ---- 10. your own trees
  $('sample').value = '0'; fire('sample');
  const loadTrees = async (text, name) => {
    Object.defineProperty($('treeFile'), 'files', { value: [new window.File([text], name)], configurable: true });
    fire('treeFile'); await sleep(300);
  };
  const strictNwk = rec.result.strictNewick, two = rec.result.treesNewick.slice(0, 2).join('\n');
  await loadTrees(two + '\n' + strictNwk, 'mine.tre');
  check('own trees are scored', /Scored 3 trees from mine\.tre: shortest length 79/.test($('progText').textContent) && /1 tree had polytomies/.test($('progText').textContent), $('progText').textContent);
  check('own trees appear in the view list with their lengths', [...$('view').options].filter((o) => /length \d+/.test(o.textContent)).length === 3 && stat('Trees scored') === '3', [...$('view').options].map((o) => o.textContent));
  $('tab-syn').click();
  check('synapomorphies work on own trees', doc.querySelectorAll('#synTable tbody tr').length > 5);
  check('methods describe scoring of the trees', /The 3 trees in mine\.tre were scored/.test($('methodsText').textContent));
  await loadTrees('(A,B,(C,D));', 'wrong.tre');
  check('trees with other taxa are refused with a reason', /No tree could be read from wrong\.tre/.test($('progText').textContent), $('progText').textContent);

  // ---- 11. robustness: gap treatment x weighting
  $('tab-gaps').click();
  $('planGapE').checked = false; $('planK').value = '3, 12'; fire('planK', 'input');
  check('plan counts the analyses', /6 analyses/.test($('planCount').textContent), $('planCount').textContent);
  $('planGapM').checked = false; fire('planGapM');
  check('plan updates', /3 analyses/.test($('planCount').textContent), $('planCount').textContent);
  $('timeLimit').value = '60'; $('intQuick').checked = true;
  $('planRun').click();
  const okP = await waitFor(() => /Compared 3 analyses/.test($('progText').textContent));
  check('robustness comparison finished', okP, $('progText').textContent);
  const rows = [...doc.querySelectorAll('#gapsBox table.plain tbody tr')].map((tr) => [...tr.children].map((td) => td.textContent));
  check('one row per analysis: equal weights, k = 3, k = 12', rows.length === 3 && /equal weights/.test(rows[0][0]) && /k = 3/.test(rows[1][0]) && /k = 12/.test(rows[2][0]), rows.map((r) => r[0]));
  check('scores: steps for equal weights, fits for implied', rows[0][1] === '79' && /^\d+\.\d{5}$/.test(rows[1][1]) && /^\d+\.\d{5}$/.test(rows[2][1]), rows.map((r) => r[1]));
  check('clade table has a column per analysis and a count', doc.querySelectorAll('#cmpTable thead th.tr').length === 3 && /Found in/.test(doc.querySelector('#cmpTable thead').textContent));
  const head = doc.querySelector('#gapsBox .headline').textContent;
  check('headline counts clades over the three analyses', /clades are found under all three analyses|no difference/.test(head), head);
  check('robustness methods paragraph', /robustness of the results to analytical choices/.test($('cmpMethods').textContent) && /implied weights \(k = 3\)/.test($('cmpMethods').textContent), $('cmpMethods').textContent.slice(0, 200));
  const rob = [...doc.querySelectorAll('#gapsBox .export')].find((x) => /number of analyses/.test(x.querySelector('.lbl').textContent));
  check('consensus labelled with the number of analyses', rob && /\)[123][,)]/.test(rob.querySelector('textarea').value), rob && rob.querySelector('textarea').value.slice(0, 100));
  $('planK').value = ''; fire('planK', 'input'); $('planGapI').checked = true; $('planGapM').checked = false; $('planGapE').checked = false; fire('planGapE');
  $('planRun').click(); await sleep(30);
  check('a plan with one analysis is not run', /at least two/.test($('planCount').textContent) && !$('run').disabled, $('planCount').textContent);

  console.log(`${pass}/${pass + fail} checks of the new features passed`);
  process.exit(fail ? 1 : 0);
})();
