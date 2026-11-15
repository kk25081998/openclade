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
      let unknown = false, gap = false;
      for (const c of inner) {
        const sym = ctx.norm(c);
        if (sym === ctx.missing || sym === '?') { unknown = true; continue; }
        if (sym === ctx.gap || sym === '-') { gap = true; continue; }
        if (ctx.match && sym === ctx.match) throw new Error('The MATCHCHAR symbol "' + c + '" cannot be used inside a polymorphic cell');
        if (!states.includes(sym)) states.push(sym);
      }
      // "{0?}" (state 0 or unknown) is unknown. "{0-}" (state 0 or inapplicable) keeps
      // the gap as an element of the set, before the states: ['-', '0'].
      if (unknown) return ['?', j + 1];
      if (!states.length) return [gap ? '-' : '?', j + 1];
      states.sort();
      if (gap) states.unshift('-');
      return [states, j + 1];
    }
    const sym = ctx.norm(ch);
    if (sym === ctx.missing || sym === '?') return ['?', i + 1];
    if (sym === ctx.gap || sym === '-') return ['-', i + 1];
    if (ctx.match && sym === ctx.match) return [MATCH, i + 1];
    if (ctx.equate && ctx.equate[sym]) return [ctx.equate[sym].slice(), i + 1];
    return [[sym], i + 1];
  }

  const MATCH = '\u0001match'; // placeholder for a MATCHCHAR cell until the first row is known

  function makeCtx(fmt) {
    const respectCase = fmt.respectCase;
    const norm = respectCase ? (c) => c : (c) => c.toUpperCase();
    return {
      // the declared symbols are compared with cells after the same case folding
      missing: norm(fmt.missing || '?'),
      gap: norm(fmt.gap || '-'),
      match: fmt.matchchar ? norm(fmt.matchchar) : null,
      norm,
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
    if (ctx.match) {
      for (let t = 0; t < cells.length; t++) for (let c = 0; c < cells[t].length; c++) {
        if (cells[t][c] !== MATCH) continue;
        const first = cells[0][c];
        if (t === 0 || first === undefined || first === MATCH) throw new Error('The MATCHCHAR symbol "' + fmt.matchchar + '" is used where there is no first-row state to copy (taxon "' + taxa[t] + '", character ' + (c + 1) + ')');
        cells[t][c] = Array.isArray(first) ? first.slice() : first;
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
            const mc = getOpt(cmd, 'matchchar'); if (mc) fmt.matchchar = mc;
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
    if (result) applyNexusSettings(blocks, result);
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
    const res = finalise(taxa, cells, nchar, ntax, fmt, null, null, warnings, errors, 'TNT');
    applyTntSettings(text, res);
    return res;
  }

  function finalise(taxa, cells, nchar, ntax, fmt, charLabels, stateLabels, warnings, errors, format) {
    if (ntax && taxa.length !== ntax) warnings.push('Header says ' + ntax + ' taxa but ' + taxa.length + ' rows were read');
    const bad = [];
    for (let t = 0; t < taxa.length; t++) {
      if (cells[t].length !== nchar) bad.push(taxa[t] + ' (' + cells[t].length + ')');
    }
    if (bad.length) errors.push('Rows with the wrong number of characters (expected ' + nchar + '): ' + bad.slice(0, 5).join(', ') + (bad.length > 5 ? ' …' : ''));
    const symSet = new Set();
    for (const row of cells) for (const c of row) if (Array.isArray(c)) for (const s of c) if (s !== '-') symSet.add(s);
    let symbols;
    if (fmt.symbols && fmt.symbols.length) {
      symbols = fmt.symbols.map((c) => (fmt.respectCase ? c : c.toUpperCase()));
      const extra = [...symSet].filter((x) => !symbols.includes(x)).sort();
      if (extra.length) { warnings.push('Symbols not declared in FORMAT SYMBOLS: ' + extra.join(' ')); symbols = symbols.concat(extra); }
    } else {
      symbols = [...symSet];
      // No SYMBOLS list: the states are 0-9 then A, B, ... (as in TNT). States that no
      // taxon shows still count in the order, so that 0 -> 3 is three steps of an
      // ordered character even if no taxon has state 2.
      const seq = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ';
      if (symbols.length && symbols.every((c) => seq.includes(c))) {
        const top = Math.max(...symbols.map((c) => seq.indexOf(c)));
        const hasDigit = symbols.some((c) => c <= '9');
        // letters alone (A C G T ...) are not read as states 10 and up
        if (hasDigit && top < 30) symbols = seq.slice(0, top + 1).split('');
      }
    }
    // Digits then letters: this is the order of the states of an ordered character.
    const rank = (c) => (/[0-9]/.test(c) ? c.charCodeAt(0) - 48 : /[A-Za-z]/.test(c) ? 10 + c.toUpperCase().charCodeAt(0) - 65 + (c >= 'a' ? 0.5 : 0) : 100 + c.charCodeAt(0));
    symbols = [...new Set(symbols)].sort((a, b) => rank(a) - rank(b));
    if (symbols.length > 30) errors.push('Too many distinct state symbols (' + symbols.length + '); maximum is 30');
    let stateOrGap = 0;
    for (const row of cells) for (const c of row) if (Array.isArray(c) && c[0] === '-') stateOrGap++;
    if (stateOrGap) warnings.push(stateOrGap + ' cell' + (stateOrGap === 1 ? ' is' : 's are') + ' coded as "a state or a gap", such as {0-}. With gaps as an extra state the gap stays one of the possibilities; otherwise such a cell is read as unknown (?). TreeSearch scores these cells differently, so lengths can differ from it.');
    // empty taxa
    for (let t = 0; t < taxa.length; t++) {
      if (cells[t].length && cells[t].every((c) => c === '?' || c === '-')) warnings.push('Taxon "' + taxa[t] + '" has no scored characters');
    }
    return { taxa, cells, nchar, symbols, charLabels: charLabels || [], stateLabels: stateLabels || [], format, warnings, errors,
      ordered: null, weights: null, excluded: null, concavity: null, settings: { ordered: 0, weighted: 0, excluded: 0, unsupported: [] } };
  }


  /* ------------------------------------------------------------------ */
  /* Character settings stored in the file: ordering, weights, exclusions */
  /* ------------------------------------------------------------------ */
  // "1-5 8 10-.\3 name all" -> zero-based character indices. sets: Map of CHARSET names.
  function charList(str, nchar, sets, unknown) {
    const out = [];
    const toks = str.replace(/\s*-\s*/g, '-').replace(/\s*\\\s*/g, '\\').trim().split(/[\s,]+/).filter(Boolean);
    for (const tk of toks) {
      const low = tk.toLowerCase().replace(/^'|'$/g, '');
      if (low === 'all' || tk === '.') { for (let i = 0; i < nchar; i++) out.push(i); continue; }
      if (sets && sets.has(low)) { out.push(...sets.get(low)); continue; }
      const m = tk.match(/^(\d+)(?:-(\d+|\.))?(?:\\(\d+))?$/);
      if (!m) { if (unknown) unknown.push(tk); continue; }
      const a = +m[1], b = m[2] === undefined ? a : (m[2] === '.' ? nchar : +m[2]), st = m[3] ? Math.max(1, +m[3]) : 1;
      for (let i = a; i <= b; i += st) if (i >= 1 && i <= nchar) out.push(i - 1);
    }
    return out;
  }
  function summariseSettings(mx) {
    const st = mx.settings;
    st.ordered = mx.ordered ? mx.ordered.filter(Boolean).length : 0;
    st.weighted = mx.weights ? mx.weights.filter((w) => w !== 1).length : 0;
    st.excluded = mx.excluded ? mx.excluded.filter(Boolean).length : 0;
    if (!st.ordered) mx.ordered = null;
    if (!st.weighted) mx.weights = null;
    if (!st.excluded) mx.excluded = null;
    if (st.unsupported.length) mx.warnings.push('This file uses ' + [...new Set(st.unsupported)].join(', ') + ', which OpenClade does not support. Those characters are analysed as unordered, so tree lengths can differ from the original analysis.');
  }
  // NEXUS: ASSUMPTIONS (TYPESET, WTSET, EXSET, OPTIONS DEFTYPE), SETS (CHARSET) and
  // PAUP (CTYPE, WEIGHTS, EXCLUDE, INCLUDE) blocks.
  function applyNexusSettings(blocks, mx) {
    const nc = mx.nchar;
    const ordered = new Array(nc).fill(false), weights = new Array(nc).fill(1), excluded = new Array(nc).fill(false);
    const sets = new Map(), unknown = [], st = mx.settings;
    const typeOf = (name) => {
      const t = name.toLowerCase();
      if (/^ord(ered)?$/.test(t)) return 'ord';
      if (/^unord(ered)?$/.test(t)) return 'unord';
      return t;
    };
    const setType = (type, list) => {
      const t = typeOf(type);
      if (t !== 'ord' && t !== 'unord') st.unsupported.push(/^irrev|^dollo/.test(t) ? t.replace(/\..*/, '') + ' characters' : 'user-defined character type "' + type + '"');
      for (const i of list) ordered[i] = t === 'ord';
    };
    // "name: list, name: list" or a (VECTOR) of one value per character
    const pairs = (body, each, vector) => {
      if (vector) {
        const v = body.trim().split(/\s+/).filter(Boolean);
        if (v.length === nc) v.forEach((x, i) => each(x, [i]));
        else mx.warnings.push('A (VECTOR) character setting has ' + v.length + ' values for ' + nc + ' characters and was skipped');
        return;
      }
      for (const part of body.split(',')) {
        const k = part.indexOf(':');
        if (k < 0) continue;
        each(part.slice(0, k).trim(), charList(part.slice(k + 1), nc, sets, unknown));
      }
    };
    for (const b of blocks) if (b.name === 'sets' || b.name === 'assumptions' || b.name === 'paup') {
      for (const cmd of commands(b.body)) {
        const m = cmd.match(/^charset\s+\*?\s*(\S+)\s*(?:\([^)]*\))?\s*=\s*([\s\S]*)$/i);
        if (m) sets.set(m[1].toLowerCase().replace(/^'|'$/g, ''), charList(m[2], nc, sets, unknown));
      }
    }
    for (const b of blocks) {
      if (b.name !== 'assumptions' && b.name !== 'paup') continue;
      const cmds = commands(b.body);
      // a starred set is the one in force; without a star the last one is used
      const pick = (re) => { const all = cmds.filter((c) => re.test(c)); return all.find((c) => /^\w+\s+\*/.test(c)) || all[all.length - 1]; };
      for (const cmd of cmds) {
        if (/^(options|set)\b/i.test(cmd)) {
          const dt = getOpt(cmd, 'deftype');
          if (dt) { if (typeOf(dt) === 'ord') ordered.fill(true); else if (typeOf(dt) !== 'unord') st.unsupported.push('default character type "' + dt + '"'); }
        }
      }
      const head = (cmd) => cmd.replace(/^\w+\s+\*?\s*/, '');
      const ts = pick(/^typeset\b/i);
      if (ts) { const m = ts.match(/^typeset\s+\*?\s*\S+\s*(\([^)]*\))?\s*=\s*([\s\S]*)$/i); if (m) pairs(m[2], setType, m[1] && /vector/i.test(m[1])); }
      const ws = pick(/^wtset\b/i);
      const setW = (w, list) => { const x = parseFloat(w); if (isFinite(x) && x >= 0) for (const i of list) weights[i] = x; };
      if (ws) { const m = ws.match(/^wtset\s+\*?\s*\S+\s*(\([^)]*\))?\s*=\s*([\s\S]*)$/i); if (m) pairs(m[2], setW, m[1] && /vector/i.test(m[1])); }
      const es = pick(/^exset\b/i);
      if (es) { const m = es.match(/^exset\s+\*?\s*\S+\s*(?:\([^)]*\))?\s*=\s*([\s\S]*)$/i); if (m) for (const i of charList(m[1], nc, sets, unknown)) excluded[i] = true; }
      for (const cmd of cmds) { // PAUP commands, in order
        if (/^ctype\b/i.test(cmd)) pairs(head(cmd), setType, false);
        else if (/^weights\b/i.test(cmd)) pairs(head(cmd), setW, false);
        else if (/^exclude\b/i.test(cmd)) { for (const i of charList(head(cmd).replace(/\/.*$/, ''), nc, sets, unknown)) excluded[i] = true; }
        else if (/^include\b/i.test(cmd)) { for (const i of charList(head(cmd).replace(/\/.*$/, ''), nc, sets, unknown)) excluded[i] = false; }
        else if (/^usertype\b/i.test(cmd)) st.unsupported.push('step matrices (USERTYPE)');
      }
    }
    if (unknown.length) mx.warnings.push('Character lists in the file contain names OpenClade could not resolve and skipped: ' + [...new Set(unknown)].slice(0, 5).join(', '));
    mx.ordered = ordered; mx.weights = weights; mx.excluded = excluded;
    summariseSettings(mx);
  }
  // TNT: ccode (+ additive, - non-additive, [ active, ] inactive, /N weight; characters
  // are numbered from 0, "a.b" is a range, "." is every character) and piwe (implied weighting).
  function applyTntSettings(text, mx) {
    const nc = mx.nchar, st = mx.settings;
    const ordered = new Array(nc).fill(false), weights = new Array(nc).fill(1), excluded = new Array(nc).fill(false);
    const xr = text.search(/\bxread\b/i);
    const tail = text.slice(xr);
    const afterMatrix = tail.slice(tail.indexOf(';') + 1);
    const beforeMatrix = text.slice(0, xr);
    // "ccode + 1 2;" and "cc+1 2;" are both accepted
    for (const m of afterMatrix.matchAll(/(?:^|[\s;])(?:ccode|cc)(?=[\s+\-\[\]\/(])([^;]*);/gi)) {
      const s = m[1];
      let spec = {}, fresh = true, i = 0;
      const apply = (list) => { for (const c of list) { if (spec.ord !== undefined) ordered[c] = spec.ord; if (spec.active !== undefined) excluded[c] = !spec.active; if (spec.weight !== undefined) weights[c] = spec.weight; } fresh = false; };
      const start = () => { if (!fresh) { spec = {}; fresh = true; } };
      while (i < s.length) {
        const ch = s[i];
        if (/\s/.test(ch)) { i++; continue; }
        if (ch === '+') { start(); spec.ord = true; i++; }
        else if (ch === '-') { start(); spec.ord = false; i++; }
        else if (ch === '[') { start(); spec.active = true; i++; }
        else if (ch === ']') { start(); spec.active = false; i++; }
        else if (ch === '(') { start(); st.unsupported.push('step-matrix (Sankoff) characters'); i++; }
        else if (ch === ')') { start(); i++; }
        else if (ch === '/') { start(); const w = s.slice(i + 1).match(/^\s*(\d+(?:\.\d+)?)/); if (w) { spec.weight = parseFloat(w[1]); i += 1 + w[0].length; } else i++; }
        else {
          const tk = s.slice(i).match(/^(\d+)?(\.)?(\d+)?/);
          if (!tk || !tk[0]) { i++; continue; }
          i += tk[0].length;
          let a, b;
          if (tk[2]) { a = tk[1] !== undefined ? +tk[1] : 0; b = tk[3] !== undefined ? +tk[3] : nc - 1; } else { a = +tk[1]; b = a; }
          const list = []; for (let c = a; c <= b; c++) if (c >= 0 && c < nc) list.push(c);
          apply(list);
        }
      }
    }
    // TNT wants piwe before xread; accept it on either side of the matrix. The last one counts.
    let pi = null;
    for (const part of [beforeMatrix, afterMatrix]) for (const m of part.matchAll(/(?:^|[\s;])piwe\s*([=\-])\s*(\d+(?:\.\d+)?)?\s*;/gi)) pi = m;
    if (pi) mx.concavity = pi[1] === '-' ? null : pi[2] ? parseFloat(pi[2]) : 3; // "piwe-" turns it off; TNT's default concavity is 3
    mx.ordered = ordered; mx.weights = weights; mx.excluded = excluded;
    summariseSettings(mx);
  }

  /* ------------------------------------------------------------------ */
  /* Check of inapplicable coding                                        */
  /* ------------------------------------------------------------------ */
  // For each character with gaps, finds the character and state that best explain
  // where the gaps are (for example "tail: absent" explains the gaps of "tail
  // colour"), then lists the cells that contradict that reading:
  //   scored-parent-absent : a state is scored although the parent says the structure is absent
  //   gap-parent-present   : a gap although the parent says the structure is present
  //   gap-parent-unknown   : a gap although the parent is unknown (?), where ? is usual
  //   missing-parent-absent: ? although the parent says absent, where a gap carries more information
  // Characters with gaps and no explaining character are reported with parent null.
  // opts.minAgree (default 0.75): share of taxa that must fit "gap exactly where the parent has that state".
  function lintInapplicable(mx, opts) {
    const minAgree = (opts && opts.minAgree) || 0.75;
    const n = mx.taxa.length, nc = mx.nchar;
    const single = (cell) => (Array.isArray(cell) && cell.length === 1 ? cell[0] : null);
    const out = [];
    const counts = { 'scored-parent-absent': 0, 'gap-parent-present': 0, 'gap-parent-unknown': 0, 'missing-parent-absent': 0 };
    let noParent = 0;
    for (let c = 0; c < nc; c++) {
      let gaps = 0;
      for (let t = 0; t < n; t++) if (mx.cells[t][c] === '-') gaps++;
      if (!gaps) continue;
      let best = null;
      if (gaps >= 2) for (let p = 0; p < nc; p++) {
        if (p === c) continue;
        const tally = new Map(); // state of p -> [gap here, scored here]
        let gapElsewhere = 0;
        for (let t = 0; t < n; t++) {
          const x = mx.cells[t][c];
          if (x === '?') continue;
          const s = single(mx.cells[t][p]);
          if (s === null) continue;
          let e = tally.get(s); if (!e) { e = [0, 0]; tally.set(s, e); }
          if (x === '-') { e[0]++; gapElsewhere++; } else e[1]++;
        }
        if (tally.size < 2) continue; // a parent must also have a "present" state
        for (const [s, e] of tally) {
          if (e[0] < 2) continue;
          const agree = e[0] / (gapElsewhere + e[1]); // gaps explained / (all gaps + scored cells under that state)
          // ties: the character that explains more gaps, then the nearer one (dependent characters usually follow their parent)
          if (agree >= minAgree && (!best || agree > best.agree || (agree === best.agree && (e[0] > best.explained || (e[0] === best.explained && Math.abs(p - c) < Math.abs(best.c - c)))))) best = { c: p, state: s, agree, explained: e[0] };
        }
      }
      const rec = { c, gaps, parent: best, issues: [] };
      if (best) {
        for (let t = 0; t < n; t++) {
          const x = mx.cells[t][c], y = mx.cells[t][best.c], s = single(y);
          let kind = null;
          if (x === '-') { if (y === '?') kind = 'gap-parent-unknown'; else if (s !== null && s !== best.state) kind = 'gap-parent-present'; }
          else if (x === '?') { if (s === best.state) kind = 'missing-parent-absent'; }
          else if (s === best.state) kind = 'scored-parent-absent';
          if (kind) { rec.issues.push({ t, kind }); counts[kind]++; }
        }
      } else noParent++;
      out.push(rec);
    }
    return { chars: out, counts, noParent, withGaps: out.length, withParent: out.length - noParent };
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

  const api = { parseMatrix, parseNexus, parseTNT, summarise, lintInapplicable, charList };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.OCParser = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);
