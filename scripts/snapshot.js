// Keeps the built app of this version at its own address, so that an analysis
// file made with it can always be opened in the version that made it.
//   npm run build && npm run snapshot
// Copies dist/openclade.html (the single-file build) to releases/v<version>/index.html.
// The Pages workflow publishes the releases folder as /v/<version>/.
const fs = require('fs');
const path = require('path');
const root = path.join(__dirname, '..');
const version = require(path.join(root, 'package.json')).version;
const src = path.join(root, 'dist', 'openclade.html');
if (!fs.existsSync(src)) { console.error('Run "npm run build" first.'); process.exit(1); }
const html = fs.readFileSync(src, 'utf8');
if (!html.includes('OpenClade v' + version)) { console.error('dist/openclade.html was not built from version ' + version + '; run "npm run build".'); process.exit(1); }
const dir = path.join(root, 'releases', 'v' + version);
fs.mkdirSync(dir, { recursive: true });
fs.writeFileSync(path.join(dir, 'index.html'), html);
console.log('releases/v' + version + '/index.html', (html.length / 1024).toFixed(0), 'KB');
