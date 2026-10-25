/* OpenClade — matrix parser (NEXUS and TNT).
 * Returns a normalised matrix object:
 *   { taxa: [names], nchar, cells: [[cell...] per taxon], symbols: [..],
 *     charLabels: [..], stateLabels: [[..]], format, warnings: [..], errors: [..] }
 * A cell is '?' (missing), '-' (gap / inapplicable) or an array of symbols
 * (length 1 = single state, >1 = polymorphism/uncertainty).
 * License: GPL-3.0-or-later
 */
(function (root) {
  'use strict';

  function stripNexusComments(text) {
    let out = '', depth = 0, inQuote = false;
    for (let i = 0; i < text.length; i++) {
      const ch = text[i];
      if (inQuote) {
        out += ch;
        if (ch === "'") {
          if (text[i + 1] === "'") { out += "'"; i++; }
          else inQuote = false;
        }
        continue;
      }
      if (depth === 0 && ch === "'") { inQuote = true; out += ch; continue; }
      if (ch === '[') { depth++; continue; }
      if (ch === ']' && depth > 0) { depth--; continue; }
      if (depth === 0) out += ch;
    }
    return out;
  }

  // Split text into NEXUS blocks: [{name, body}]
  function nexusBlocks(text) {
    const blocks = [];
    const re = /begin\s+([A-Za-z_]+)\s*;([\s\S]*?)\bend(?:block)?\s*;/gi;
    let m;
    while ((m = re.exec(text))) blocks.push({ name: m[1].toLowerCase(), body: m[2] });
    return blocks;
  }

  // Split a block body into commands (top-level ';' outside quotes)
  function commands(body) {
    const out = [];
    let cur = '', inQuote = false;
    for (let i = 0; i < body.length; i++) {
      const ch = body[i];
      if (ch === "'") {
        if (inQuote && body[i + 1] === "'") { cur += "''"; i++; continue; }
        inQuote = !inQuote;
      }
      if (ch === ';' && !inQuote) { out.push(cur.trim()); cur = ''; continue; }
      cur += ch;
    }
    if (cur.trim()) out.push(cur.trim());
    return out.filter(Boolean);
  }

  function getOpt(cmd, key) {
    const re = new RegExp('\\b' + key + '\\s*=\\s*("[^"]*"|\'[^\']*\'|[^\\s;]+)', 'i');
    const m = cmd.match(re);
    if (!m) return null;
    let v = m[1];
    if ((v[0] === '"' && v.endsWith('"')) || (v[0] === "'" && v.endsWith("'"))) v = v.slice(1, -1);
    return v;
  }

  // Tokenise a name: quoted 'x y' or bare word. Returns [name, nextIndex]
  function readName(s, i) {
    while (i < s.length && /\s/.test(s[i])) i++;
    if (i >= s.length) return [null, i];
    if (s[i] === "'" || s[i] === '"') {
      const q = s[i];
      let name = '';
      i++;
      while (i < s.length) {
        if (s[i] === q) {
          if (s[i + 1] === q) { name += q; i += 2; continue; }
          i++; break;
        }
        name += s[i++];
      }
      return [name, i];
    }
    let name = '';
    while (i < s.length && !/\s/.test(s[i])) name += s[i++];
    return [name.replace(/_/g, ' ') === name ? name : name, i];
  }

  const OPEN = { '(': ')', '{': '}', '[': ']' };

  // Read one cell starting at i (skipping whitespace). Returns [cell, nextIndex]
  function readCell(s, i, ctx) {
    while (i < s.length && /\s/.test(s[i])) i++;
    if (i >= s.length) return [null, i];
    const ch = s[i];
    if (OPEN[ch] && ch !== ctx.noOpen) {
      const close = OPEN[ch];
      const j = s.indexOf(close, i + 1);
      if (j < 0) throw new Error('Unclosed polymorphism starting at "' + s.slice(i, i + 20) + '"');
      const inner = s.slice(i + 1, j).replace(/[\s,&\/]/g, '');
      const states = [];
      for (const c of inner) {
        const sym = ctx.norm(c);
        if (sym === ctx.missing || sym === '?') continue;
        if (sym === ctx.gap || sym === '-') continue;
        if (!states.includes(sym)) states.push(sym);
      }
      return [states.length ? states.sort() : '?', j + 1];
    }
    const sym = ctx.norm(ch);
    if (sym === ctx.missing || sym === '?') return ['?', i + 1];
    if (sym === ctx.gap || sym === '-') return ['-', i + 1];
    if (ctx.equate && ctx.equate[sym]) return [ctx.equate[sym].slice(), i + 1];
    return [[sym], i + 1];
  }

  function makeCtx(fmt) {
    const respectCase = fmt.respectCase;
    return {
      missing: fmt.missing || '?',
      gap: fmt.gap || '-',
      norm: respectCase ? (c) => c : (c) => c.toUpperCase(),
      equate: fmt.equate || null,
      noOpen: fmt.noOpen || null,
    };
  }

  const IUPAC = {
    R: ['A', 'G'], Y: ['C', 'T'], S: ['C', 'G'], W: ['A', 'T'], K: ['G', 'T'], M: ['A', 'C'],
    B: ['C', 'G', 'T'], D: ['A', 'G', 'T'], H: ['A', 'C', 'T'], V: ['A', 'C', 'G'], N: null,
  };

  function parseMatrixRows(body, ntax, nchar, fmt, interleave, warnings) {
    const ctx = makeCtx(fmt);
    const taxa = [], cells = [];
    const index = new Map();
    if (!interleave) {
      let i = 0;
      while (i < body.length) {
        const [name, j] = readName(body, i);
        if (name === null || name === '') break;
        i = j;
        const row = [];
        while (row.length < nchar) {
          const [cell, k] = readCell(body, i, ctx);
          if (cell === null) break;
          row.push(cell);
          i = k;
        }
        if (index.has(name)) warnings.push('Duplicate taxon name "' + name + '"');
        index.set(name, taxa.length);
        taxa.push(name);
        cells.push(row);
      }
    } else {
      const lines = body.split(/\r?\n/);
      for (const line of lines) {
        if (!line.trim()) continue;
        const [name, j] = readName(line, 0);
        if (!name) continue;
        let t = index.get(name);
        if (t === undefined) { t = taxa.length; index.set(name, t); taxa.push(name); cells.push([]); }
        let i = j;
        while (true) {
          const [cell, k] = readCell(line, i, ctx);
          if (cell === null) break;
          cells[t].push(cell);
          i = k;
        }
      }
    }
    return { taxa, cells };
  }

  function parseNexus(text) {
    const warnings = [], errors = [];
    const clean = stripNexusComments(text);
    const blocks = nexusBlocks(clean);
    let taxaFromBlock = null;
    let result = null;
    let pendingLabels = null;
    for (const b of blocks) {
      if (b.name === 'taxa') {
        for (const cmd of commands(b.body)) {
          if (/^taxlabels\b/i.test(cmd)) {
            const s = cmd.replace(/^taxlabels/i, '');
            const names = [];
            let i = 0;
            while (true) { const [n, j] = readName(s, i); if (!n) break; names.push(n); i = j; }
            taxaFromBlock = names;
          }
        }
      }
      if (b.name === 'data' || b.name === 'characters') {
        let ntax = taxaFromBlock ? taxaFromBlock.length : null, nchar = null;
        const fmt = { missing: '?', gap: '-', respectCase: false, symbols: null, datatype: 'standard' };
        let interleave = false;
        let charLabels = null, stateLabels = null, matrixBody = null;
        for (const cmd of commands(b.body)) {
          const head = cmd.split(/\s+/)[0].toLowerCase();
          if (head === 'dimensions') {
            const nt = getOpt(cmd, 'ntax'); if (nt) ntax = parseInt(nt, 10);
            const nc = getOpt(cmd, 'nchar'); if (nc) nchar = parseInt(nc, 10);
          } else if (head === 'format') {
            const miss = getOpt(cmd, 'missing'); if (miss) fmt.missing = miss;
            const gap = getOpt(cmd, 'gap'); if (gap) fmt.gap = gap;
            const dt = getOpt(cmd, 'datatype'); if (dt) fmt.datatype = dt.toLowerCase();
            const sy = cmd.match(/symbols\s*=\s*"([^"]*)"/i) || cmd.match(/symbols\s*=\s*'([^']*)'/i);
            if (sy) fmt.symbols = sy[1].replace(/\s+/g, '').split('');
            if (/\brespectcase\b/i.test(cmd)) fmt.respectCase = true;
            const il = cmd.match(/\binterleave(?:\s*=\s*(\w+))?/i);
            if (il) interleave = !il[1] || /^(yes|true)$/i.test(il[1]);
            const eq = getOpt(cmd, 'equate');
            if (eq) warnings.push('EQUATE macros are not supported and were ignored');
          } else if (head === 'charstatelabels') {
            const parsed = parseCharStateLabels(cmd.replace(/^charstatelabels/i, ''));
            charLabels = parsed.charLabels; stateLabels = parsed.stateLabels;
          } else if (head === 'charlabels') {
            const s = cmd.replace(/^charlabels/i, '');
            const names = []; let i = 0;
            while (true) { const [n, j] = readName(s, i); if (!n) break; names.push(n); i = j; }
            charLabels = names;
          } else if (head === 'statelabels') {
            const parts = cmd.replace(/^statelabels/i, '').split(',');
            stateLabels = stateLabels || [];
            for (const p of parts) {
              const m = p.trim().match(/^(\d+)\s*([\s\S]*)$/);
              if (!m) continue;
              const st = []; let i = 0;
              while (true) { const [n, j] = readName(m[2], i); if (!n) break; st.push(n); i = j; }
              stateLabels[parseInt(m[1], 10) - 1] = st;
            }
          } else if (head === 'matrix') {
            matrixBody = cmd.replace(/^matrix/i, '');
          }
        }
        if (matrixBody === null && (charLabels || stateLabels)) {
          // label-only block (e.g. Mesquite/MorphoBank exports): keep for later
          pendingLabels = { charLabels, stateLabels };
          continue;
        }
        if (!nchar) { errors.push('NEXUS DIMENSIONS NCHAR is missing'); continue; }
        if (matrixBody === null) { errors.push('No MATRIX command found'); continue; }
        if (/dna|rna|nucleotide/.test(fmt.datatype)) {
          fmt.equate = Object.assign({}, IUPAC, { U: ['T'] });
          for (const k of Object.keys(fmt.equate)) if (fmt.equate[k] === null) delete fmt.equate[k];
          fmt.equate.N = null; delete fmt.equate.N;
          fmt.missingExtra = 'N';
        }
        let rows;
        try {
          rows = parseMatrixRows(matrixBody, ntax, nchar, fmt, interleave, warnings);
        } catch (e) { errors.push(e.message); continue; }
        if (fmt.missingExtra) {
          for (const r of rows.cells) for (let c = 0; c < r.length; c++) if (Array.isArray(r[c]) && r[c][0] === 'N') r[c] = '?';
        }
        result = finalise(rows.taxa, rows.cells, nchar, ntax, fmt, charLabels, stateLabels, warnings, errors, 'NEXUS');
      }
    }
    if (result && pendingLabels) {
      if (pendingLabels.charLabels && !result.charLabels.length) result.charLabels = pendingLabels.charLabels;
      if (pendingLabels.stateLabels && !result.stateLabels.length) result.stateLabels = pendingLabels.stateLabels;
    }
    if (!result) {
      if (!errors.length) errors.push('No DATA or CHARACTERS block found');
      return { taxa: [], cells: [], nchar: 0, symbols: [], warnings, errors, format: 'NEXUS' };
    }
    return result;
  }

  function parseCharStateLabels(s) {
    // Format: 1 name / state0 state1, 2 name / s0 s1, ...
    const charLabels = [], stateLabels = [];
    const parts = [];
    let cur = '', inQuote = false;
    for (let i = 0; i < s.length; i++) {
      const ch = s[i];
      if (ch === "'") {
        if (inQuote && s[i + 1] === "'") { cur += "''"; i++; continue; }
        inQuote = !inQuote;
      }
      if (ch === ',' && !inQuote) { parts.push(cur); cur = ''; continue; }
      cur += ch;
    }
    if (cur.trim()) parts.push(cur);
    for (const p of parts) {
      const m = p.trim().match(/^(\d+)\s*([\s\S]*)$/);
      if (!m) continue;
      const idx = parseInt(m[1], 10) - 1;
      let rest = m[2];
      let name = '', states = [];
      const slash = findTopLevel(rest, '/');
      const namePart = slash >= 0 ? rest.slice(0, slash) : rest;
      const statePart = slash >= 0 ? rest.slice(slash + 1) : '';
      const [n] = readName(namePart, 0);
      name = n || '';
      let i = 0;
      while (true) { const [st, j] = readName(statePart, i); if (!st) break; states.push(st); i = j; }
      charLabels[idx] = name;
      stateLabels[idx] = states;
    }
    return { charLabels, stateLabels };
  }

  function findTopLevel(s, ch) {
    let inQuote = false;
    for (let i = 0; i < s.length; i++) {
      if (s[i] === "'") inQuote = !inQuote;
      else if (s[i] === ch && !inQuote) return i;
    }
    return -1;
  }

  function parseTNT(text) {
    const warnings = [], errors = [];
    // remove TNT single-quoted titles after xread and comments
    const m = text.match(/\bxread\b([\s\S]*)/i);
    if (!m) return { taxa: [], cells: [], nchar: 0, symbols: [], warnings, errors: ['No xread command found'], format: 'TNT' };
    let s = m[1];
    s = s.replace(/^\s*'[^']*'/, '');
    const dm = s.match(/^\s*(\d+)\s+(\d+)/);
    if (!dm) return { taxa: [], cells: [], nchar: 0, symbols: [], warnings, errors: ['Could not read nchar/ntax after xread'], format: 'TNT' };
    const nchar = parseInt(dm[1], 10), ntax = parseInt(dm[2], 10);
    s = s.slice(dm[0].length);
    const end = s.indexOf(';');
    const body = end >= 0 ? s.slice(0, end) : s;
    const fmt = { missing: '?', gap: '-', respectCase: false, noOpen: '(' };
    const ctx = makeCtx(fmt);
    const taxa = [], cells = [];
    const index = new Map();
    // TNT data may be in blocks introduced by &[...]
    const chunks = body.split(/&\s*\[[^\]]*\]|&/);
    const multiBlock = chunks.length > 1;
    for (const chunk of chunks) {
      let i = 0;
      while (i < chunk.length) {
        const [name, j] = readName(chunk, i);
        if (!name) break;
        i = j;
        let t = index.get(name);
        if (t === undefined) { t = taxa.length; index.set(name, t); taxa.push(name); cells.push([]); }
        const row = cells[t];
        // read until end of line for multi-block; otherwise until nchar
        if (multiBlock) {
          const eol = chunk.indexOf('\n', i);
          const lineEnd = eol < 0 ? chunk.length : eol;
          while (i < lineEnd) {
            const [cell, k] = readCell(chunk.slice(0, lineEnd), i, ctx);
            if (cell === null) break;
            row.push(cell); i = k;
          }
        } else {
          while (row.length < nchar) {
            const [cell, k] = readCell(chunk, i, ctx);
            if (cell === null) break;
            row.push(cell); i = k;
          }
        }
      }
    }
    return finalise(taxa, cells, nchar, ntax, fmt, null, null, warnings, errors, 'TNT');
  }

  function finalise(taxa, cells, nchar, ntax, fmt, charLabels, stateLabels, warnings, errors, format) {
    if (ntax && taxa.length !== ntax) warnings.push('Header says ' + ntax + ' taxa but ' + taxa.length + ' rows were read');
    const bad = [];
    for (let t = 0; t < taxa.length; t++) {
      if (cells[t].length !== nchar) bad.push(taxa[t] + ' (' + cells[t].length + ')');
    }
    if (bad.length) errors.push('Rows with the wrong number of characters (expected ' + nchar + '): ' + bad.slice(0, 5).join(', ') + (bad.length > 5 ? ' …' : ''));
    const symSet = new Set();
    for (const row of cells) for (const c of row) if (Array.isArray(c)) for (const s of c) symSet.add(s);
    let symbols;
    if (fmt.symbols && fmt.symbols.length) {
      symbols = fmt.symbols.map((c) => (fmt.respectCase ? c : c.toUpperCase()));
      const extra = [...symSet].filter((x) => !symbols.includes(x)).sort();
      if (extra.length) { warnings.push('Symbols not declared in FORMAT SYMBOLS: ' + extra.join(' ')); symbols = symbols.concat(extra); }
    } else {
      symbols = [...symSet].sort();
    }
    if (symbols.length > 30) errors.push('Too many distinct state symbols (' + symbols.length + '); maximum is 30');
    // empty taxa
    for (let t = 0; t < taxa.length; t++) {
      if (cells[t].length && cells[t].every((c) => c === '?' || c === '-')) warnings.push('Taxon "' + taxa[t] + '" has no scored characters');
    }
    return { taxa, cells, nchar, symbols, charLabels: charLabels || [], stateLabels: stateLabels || [], format, warnings, errors };
  }

  function parseMatrix(text) {
    if (/^\s*#nexus/i.test(text) || /begin\s+(data|characters)\s*;/i.test(text)) return parseNexus(text);
    if (/\bxread\b/i.test(text)) return parseTNT(text);
    return { taxa: [], cells: [], nchar: 0, symbols: [], warnings: [], errors: ['Unrecognised file: expected a NEXUS (#NEXUS … MATRIX) or TNT (xread) matrix'], format: '?' };
  }

  // Summary statistics for the data check panel
  function summarise(mx) {
    const n = mx.taxa.length, nc = mx.nchar;
    let missing = 0, gaps = 0, poly = 0, informative = 0, invariant = 0, inapplicChars = 0;
    for (let c = 0; c < nc; c++) {
      const counts = new Map();
      let g = 0;
      for (let t = 0; t < n; t++) {
        const cell = mx.cells[t][c];
        if (cell === '?') missing++;
        else if (cell === '-') { gaps++; g++; }
        else {
          if (cell.length > 1) poly++;
          else counts.set(cell[0], (counts.get(cell[0]) || 0) + 1);
        }
      }
      if (g > 0) inapplicChars++;
      const vals = [...counts.values()];
      const multi = vals.filter((v) => v >= 2).length;
      if (multi >= 2) informative++;
      if (counts.size <= 1) invariant++;
    }
    const total = n * nc || 1;
    return {
      ntax: n, nchar: nc, informative, invariant, inapplicChars,
      pctMissing: (100 * missing) / total, pctGap: (100 * gaps) / total, polymorphic: poly,
    };
  }

  const api = { parseMatrix, parseNexus, parseTNT, summarise };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.OCParser = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);
