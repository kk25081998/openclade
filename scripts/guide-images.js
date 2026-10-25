// Regenerates the screenshots used in the user guide (docs/img) from the built app.
//   npm run build && npm i --no-save playwright sharp && node scripts/guide-images.js
// Writes WebP images, which build.js copies to dist/docs/img.
// FONT_DIR: optional folder of @fontsource packages, used when Google Fonts
// cannot be reached (for example in a sandbox without internet access).
const http = require('http');
const fs = require('fs');
const path = require('path');
const { chromium } = require('playwright');
const sharp = require('sharp');

const ROOT = path.join(__dirname, '..');
const DIST = path.join(ROOT, 'dist');
const OUT = path.join(ROOT, 'docs', 'img');
const FONT_DIR = process.env.FONT_DIR;
const TYPES = { '.html': 'text/html; charset=utf-8', '.webp': 'image/webp', '.png': 'image/png', '.woff2': 'font/woff2', '.json': 'application/json' };

function fontCSS() {
  const f = (pkg, file) => `/__fonts/${pkg}/files/${file}`;
  return [
    `@font-face{font-family:"Bricolage Grotesque";font-style:normal;font-weight:200 800;font-stretch:75% 100%;src:url(${f('@fontsource-variable/bricolage-grotesque', 'bricolage-grotesque-latin-standard-normal.woff2')}) format("woff2-variations")}`,
    ...[400, 500, 600].map((w) => `@font-face{font-family:"Instrument Sans";font-style:normal;font-weight:${w};src:url(${f('@fontsource/instrument-sans', `instrument-sans-latin-${w}-normal.woff2`)}) format("woff2")}`),
    ...[400, 500].map((w) => `@font-face{font-family:"JetBrains Mono";font-style:normal;font-weight:${w};src:url(${f('@fontsource/jetbrains-mono', `jetbrains-mono-latin-${w}-normal.woff2`)}) format("woff2")}`),
  ].join('\n');
}

const server = http.createServer((req, res) => {
  let p = decodeURIComponent(req.url.split('?')[0]);
  let file;
  if (p.startsWith('/__fonts/') && FONT_DIR) file = path.join(FONT_DIR, p.slice('/__fonts/'.length));
  else { if (p.endsWith('/')) p += 'index.html'; file = path.join(DIST, path.normalize(p).replace(/^(\.\.[/\\])+/, '')); }
  fs.readFile(file, (err, buf) => {
    if (err) { res.writeHead(404); res.end(); return; }
    res.writeHead(200, { 'content-type': TYPES[path.extname(file)] || 'application/octet-stream' });
    res.end(buf);
  });
});

async function newPage(browser, port, scheme, viewport) {
  const page = await browser.newPage({ viewport, deviceScaleFactor: 2, colorScheme: scheme });
  page.on('pageerror', (e) => console.log('page error:', e.message));
  if (FONT_DIR) {
    await page.route('https://fonts.googleapis.com/**', (r) => r.fulfill({ contentType: 'text/css', body: fontCSS().split('/__fonts/').join(`http://127.0.0.1:${port}/__fonts/`) }));
    await page.route('https://fonts.gstatic.com/**', (r) => r.abort());
  }
  // A first-time visitor would see the example notice; the screenshots leave it out.
  await page.addInitScript(() => { try { localStorage.setItem('openclade.welcomeDismissed', '1'); } catch (e) { /* ignore */ } });
  await page.goto(`http://127.0.0.1:${port}/index.html`);
  await page.evaluate(() => document.fonts.ready);
  return page;
}
const finished = (page, re) => page.waitForFunction((src) => new RegExp(src).test(document.getElementById('progText').textContent), re.source, { timeout: 240000 });

async function save(buf, name) {
  await sharp(buf).webp({ quality: 88 }).toFile(path.join(OUT, name + '.webp'));
  console.log('wrote', name + '.webp');
}

server.listen(0, '127.0.0.1', async () => {
  const port = server.address().port;
  fs.mkdirSync(OUT, { recursive: true });
  const browser = await chromium.launch();
  try {
    // 1. Overview of the whole app after the example analysis, light and dark.
    for (const scheme of ['light', 'dark']) {
      const page = await newPage(browser, port, scheme, { width: 1440, height: 900 });
      await finished(page, /^Finished/);
      await page.waitForTimeout(300);
      await save(await page.screenshot(), 'app-overview' + (scheme === 'dark' ? '-dark' : ''));
      await page.close();
    }
    // 2. A character mapped on a single tree (Vinther 2008): the first ordinary character with 3 steps.
    {
      const page = await newPage(browser, port, 'light', { width: 1200, height: 900 });
      await finished(page, /^Finished/);
      const j = await page.evaluate(() => {
        const rows = [...document.querySelectorAll('#charTable tbody tr')];
        const r = rows.find((tr) => /Fitch/.test(tr.children[2].textContent) && +tr.children[3].textContent === 3) || rows[0];
        return r.dataset.j;
      });
      await page.selectOption('#mapChar', j);
      await page.waitForTimeout(200);
      const top = await page.locator('#pane-tree .toolbar').boundingBox();
      const tree = await page.locator('#treebox').boundingBox();
      const legend = await page.locator('#legend').boundingBox();
      const clip = { x: tree.x - 8, y: top.y - 8, width: tree.width + 16, height: Math.min(legend.y + legend.height + 8, tree.y + 640) - top.y + 8 };
      await save(await page.screenshot({ clip, fullPage: true }), 'tree-mapping');
      await page.close();
    }
    // 3. Gap-treatment comparison on Vinther et al. (2008), short enough to show the whole tree.
    {
      const page = await newPage(browser, port, 'light', { width: 1200, height: 900 });
      await finished(page, /^Finished/);
      await page.check('#intStd', { force: true });
      await page.evaluate(() => { document.querySelector('.more').open = true; document.getElementById('timeLimit').value = '90'; });
      await page.click('#compare');
      await finished(page, /Compared 3 gap treatments/);
      await page.waitForTimeout(300);
      const head = await page.locator('#gapsBox .headline').boundingBox();
      const legend = await page.locator('#cmpLegend').boundingBox();
      const box = await page.locator('#gapsBox').boundingBox();
      const clip = { x: box.x - 8, y: head.y - 10, width: box.width + 16, height: legend.y + legend.height - head.y + 20 };
      await save(await page.screenshot({ clip, fullPage: true }), 'gap-comparison');
      console.log('comparison:', await page.textContent('#progText'));
      await page.close();
    }
  } finally {
    await browser.close();
    server.close();
  }
});
