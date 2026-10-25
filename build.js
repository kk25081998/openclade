// Builds the single-file app: inlines parser, engine and example matrices.
const fs = require('fs');
const D = process.env.DATASETS;
const samples = [
  ['Vinther2008', 'Vinther et al. 2008 · machaeridians (23 taxa)', 'Vinther, J., Van Roy, P. & Briggs, D. E. G. 2008. Machaeridians are Palaeozoic armoured annelids. Nature 451: 185–188.'],
  ['Longrich2010', 'Longrich et al. 2010 · pachycephalosaurs (20 taxa)', 'Longrich, N. R., Sankey, J. & Tanke, D. 2010. Texacephale langstoni, a new genus of pachycephalosaurid from the upper Campanian Aguja Formation, Texas. Cretaceous Research 31: 274–284.'],
  ['Sansom2010', 'Sansom et al. 2010 · jawless vertebrates (23 taxa)', 'Sansom, R. S., Freedman, K., Gabbott, S. E., Aldridge, R. J. & Purnell, M. A. 2010. Taphonomy and affinity of an enigmatic Silurian vertebrate, Jamoytius kerwoodi. Palaeontology 53: 1393–1409.'],
  ['OLeary1999', "O'Leary & Geisler 1999 · Cetacea in Mammalia (40 taxa)", "O'Leary, M. A. & Geisler, J. H. 1999. The position of Cetacea within Mammalia: phylogenetic analysis of morphological data from extinct and extant taxa. Systematic Biology 48: 455–490."],
  ['Aria2015', 'Aria et al. 2015 · stem arthropods (35 taxa)', 'Aria, C., Caron, J.-B. & Gaines, R. 2015. A large new leanchoiliid from the Burgess Shale and the influence of inapplicable states on stem arthropod phylogeny. Palaeontology 58: 629–660.'],
].map(([f, label, citation]) => ({ label, citation, text: fs.readFileSync(D + '/' + f + '.nex', 'utf8') }));
let html = fs.readFileSync('src/app.html', 'utf8');
const safe = (s) => s.replace(/<\/script/gi, '<\\/script');
html = html.split('{{PARSER}}').join(safe(fs.readFileSync('src/parser.js', 'utf8')))
  .split('{{ENGINE}}').join(safe(fs.readFileSync('src/engine.js', 'utf8')))
  .split('{{WORKER}}').join(safe(fs.readFileSync('src/worker.js', 'utf8')))
  .split('{{PARALLEL}}').join(safe(fs.readFileSync('src/parallel.js', 'utf8')))
  .split('{{SAMPLES}}').join(JSON.stringify(samples).replace(/</g, '\\u003c'));
fs.mkdirSync('dist', { recursive: true });
fs.writeFileSync('dist/openclade.html', html);
// Standalone page for static hosting (GitHub Pages): full document skeleton
const title = (html.match(/<title>[\s\S]*?<\/title>/) || [''])[0];
const body = html.replace(title, '');
const standalone = '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n' + title +
  '\n<meta name="description" content="OpenClade: parsimony analysis of morphological character matrices in the browser, with correct handling of inapplicable characters.">\n<meta name="author" content="Kartik Khosa">\n<style>body{margin:0}img{max-width:100%}[hidden]{display:none!important}</style>\n</head>\n<body>\n' + body + '\n</body>\n</html>\n';
fs.writeFileSync('dist/index.html', standalone);
console.log('dist/index.html', (standalone.length / 1024).toFixed(0) + ' KB');
console.log('dist/openclade.html', (html.length / 1024).toFixed(0) + ' KB');
