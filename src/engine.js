/* OpenClade — parsimony engine.
 *
 * Scoring follows MorphyLib (Brazeau, Smith & Guillerme) as used in the
 * TreeSearch R package: characters with more than two inapplicable tokens are
 * scored with the inapplicable-aware algorithm of Brazeau, Guillerme & Smith
 * (2019, Systematic Biology 68:619-631); other characters use Fitch (1971)
 * parsimony with the gap token treated as missing data.
 *
 * License: GPL-3.0-or-later (algorithmic port of GPL-3 MorphyLib logic).
 */
(function (root) {
  'use strict';

  const NA = 1;            // inapplicable token bit
  const MISSING = -1;      // all bits set
  const ISAPPLIC = ~1;     // every bit except NA
  const NACUTOFF = 2;      // > NACUTOFF gaps => inapplicable algorithm
  const EMPTY = new Int32Array(0);

  /* ------------------------------------------------------------------ */
  /* Random numbers (seedable so analyses are reproducible)              */
  /* ------------------------------------------------------------------ */
  function RNG(seed) {
    let s = (seed >>> 0) || 0x9e3779b9;
    this.next = function () { // mulberry32
      s = (s + 0x6d2b79f5) | 0;
      let t = s;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
    this.int = (n) => Math.floor(this.next() * n);
    this.shuffle = (a) => {
      for (let i = a.length - 1; i > 0; i--) { const j = this.int(i + 1); const x = a[i]; a[i] = a[j]; a[j] = x; }
      return a;
    };
  }

  /* ------------------------------------------------------------------ */
  /* Data encoding                                                       */
  /* ------------------------------------------------------------------ */
  // A character of the matrix becomes one or more *columns*:
  //   - an unordered character is one column;
  //   - an ordered (additive) character with states lo..hi is hi-lo binary
  //     columns, one per threshold ("state >= th"). The Fitch lengths of the
  //     columns sum to the Farris (1970) length of the ordered character.
  // Excluded characters and characters of weight 0 have no columns.
  //
  // gapMode: 'inapplicable' (default) | 'missing' | 'extra'
  // co: { applyOrdered, applyWeights, applyExcluded (default true: use mx.ordered,
  //       mx.weights, mx.excluded when the matrix has them), concavity (k of
  //       implied weighting; null or Infinity = equal weights) }
  // An ordered character scored with the inapplicable algorithm: every column
  // counts its changes between applicable states, but only the first column
  // counts extra applicable regions, so a region costs one step however many
  // states the character has and however they are numbered.
  // co.orderedRegions = 'decompose' instead reproduces TreeSearch with
  // TreeTools::Decompose(): one column per state above the lowest symbol of the
  // matrix, gaps kept in every column, so a region costs one step per column.
  function encode(mx, gapMode, co) {
    gapMode = gapMode || 'inapplicable';
    co = co || {};
    const ntax = mx.taxa.length, nOrig = mx.nchar;
    const symIndex = new Map();
    mx.symbols.forEach((s, i) => symIndex.set(s, i));
    const extraBit = mx.symbols.length; // used for 'extra' mode
    const push = gapMode === 'missing' ? 0 : 1; // keep bit 0 for NA
    const ordIn = co.applyOrdered !== false && mx.ordered ? mx.ordered : null;
    const wIn = co.applyWeights !== false && mx.weights ? mx.weights : null;
    const exIn = co.applyExcluded !== false && mx.excluded ? mx.excluded : null;
    const nGaps = new Int32Array(nOrig), charNA = new Uint8Array(nOrig), charOrdered = new Uint8Array(nOrig), charIncluded = new Uint8Array(nOrig);
    const charWeight = new Float64Array(nOrig).fill(1), charLo = new Int32Array(nOrig);
    let wScale = 1;
    for (let c = 0; c < nOrig; c++) {
      let w = wIn ? Number(wIn[c]) : 1;
      if (wIn && (wIn[c] == null || !isFinite(w))) w = 1;
      if (w < 0) w = 0;
      charWeight[c] = w;
      if (Math.abs(w - Math.round(w)) > 1e-9) wScale = 1000;
    }
    const cols = [], charCols = [], orderedAsUnordered = [];
    for (let c = 0; c < nOrig; c++) {
      let g = 0, lo = Infinity, hi = -Infinity;
      for (let t = 0; t < ntax; t++) {
        const cell = mx.cells[t][c];
        if (cell === '-') g++;
        else if (cell !== '?') for (const s of cell) { const i = symIndex.get(s); if (i < lo) lo = i; if (i > hi) hi = i; }
      }
      nGaps[c] = g;
      charCols.push([]);
      const inc = !(exIn && exIn[c]) && charWeight[c] > 0;
      charIncluded[c] = inc ? 1 : 0;
      if (!inc) continue;
      const na = gapMode === 'inapplicable' && g > NACUTOFF;
      charNA[c] = na ? 1 : 0;
      let ord = !!(ordIn && ordIn[c]) && hi - lo > 1; // two-state characters are the same ordered or not
      // A gap read as one more state has no place in the order of the others.
      if (ord && gapMode === 'extra' && g > 0) { ord = false; orderedAsUnordered.push(c); }
      charOrdered[c] = ord ? 1 : 0;
      charLo[c] = lo === Infinity ? 0 : lo;
      if (!ord) { charCols[c].push(cols.length); cols.push({ c, th: -1, na, rf: 0 }); }
      else if (co.orderedRegions === 'decompose') { charLo[c] = 0; for (let th = 1; th <= hi; th++) { charCols[c].push(cols.length); cols.push({ c, th, na, rf: 0 }); } }
      else for (let th = lo + 1; th <= hi; th++) { charCols[c].push(cols.length); cols.push({ c, th, na, rf: na && th > lo + 1 ? 1 : 0 }); }
    }
    const nchar = cols.length;
    const tips = new Int32Array(ntax * nchar), screenTips = new Int32Array(ntax * nchar);
    const isNA = new Uint8Array(nchar), regionFree = new Uint8Array(nchar), colChar = new Int32Array(nchar), colTh = new Int32Array(nchar);
    const baseW = new Float64Array(nchar);
    const LOW = 1 << push, HIGH = 1 << (push + 1);
    for (let j = 0; j < nchar; j++) {
      const col = cols[j], c = col.c;
      isNA[j] = col.na ? 1 : 0; regionFree[j] = col.rf; colChar[j] = c; colTh[j] = col.th;
      baseW[j] = Math.round(charWeight[c] * wScale);
      for (let t = 0; t < ntax; t++) {
        const cell = mx.cells[t][c];
        let v;
        if (cell === '?') v = MISSING;
        else if (cell === '-') {
          if (gapMode === 'extra') v = 1 << (extraBit + push);
          else v = col.na ? NA : MISSING;
        } else if (col.th < 0) {
          v = 0;
          for (const s of cell) v |= 1 << (symIndex.get(s) + push);
        } else { // ordered: a cell with several states is the range from its lowest to its highest
          let mn = Infinity, mx2 = -Infinity;
          for (const s of cell) { const i = symIndex.get(s); if (i < mn) mn = i; if (i > mx2) mx2 = i; }
          v = (mn < col.th ? LOW : 0) | (mx2 >= col.th ? HIGH : 0);
        }
        tips[t * nchar + j] = v;
        // Screening copy: gaps as missing data gives a lower bound on the length of an
        // inapplicable-type column (checked in test/screen_bound.js). Columns that do not
        // count regions have no such bound, so they never count in the screen.
        screenTips[t * nchar + j] = col.na ? (col.rf ? MISSING : (v === NA ? MISSING : v)) : v;
      }
    }
    const naChars = [], stdChars = [];
    for (let j = 0; j < nchar; j++) (isNA[j] ? naChars : stdChars).push(j);
    // Work order: standard columns first, then inapplicable-type columns. Groups = characters.
    const work = Int32Array.from(stdChars.concat(naChars));
    const pos = new Int32Array(nchar); for (let i = 0; i < nchar; i++) pos[work[i]] = i;
    const gStart = [], gLen = [], gCol = [], gChar = [], gIsNA = [];
    for (const pass of [0, 1]) for (let c = 0; c < nOrig; c++) {
      if (!charCols[c].length || charNA[c] !== pass) continue;
      gStart.push(pos[charCols[c][0]]); gLen.push(charCols[c].length); gCol.push(charCols[c][0]); gChar.push(c); gIsNA.push(pass);
    }
    let nStdGroups = 0; for (const x of gIsNA) if (!x) nStdGroups++;
    const data = {
      ntax, nchar, nOrig, tips, screenTips, isNA, regionFree, nGaps, gapMode,
      naChars: Int32Array.from(naChars), stdChars: Int32Array.from(stdChars), work,
      symbols: mx.symbols.slice(), push, colChar, colTh, charCols, charNA, charOrdered, charIncluded, charWeight, charLo,
      baseW, wScale, scoreDiv: wScale, orderedAsUnordered,
      gStart: Int32Array.from(gStart), gLen: Int32Array.from(gLen), gCol: Int32Array.from(gCol), gChar: Int32Array.from(gChar),
      gIsNA: Uint8Array.from(gIsNA), nGroups: gStart.length, nStdGroups, iw: null,
    };
    const k = co.concavity;
    if (k != null && isFinite(k) && k > 0) setImplied(data, k);
    return data;
  }

  // Implied weighting (Goloboff 1993): a character with e extra steps costs e / (e + k).
  // Costs are held as integers (fit * FIX, rounded per character) so that sums are exact
  // and equally good trees compare as equal whatever the order of addition.
  // FIX is the largest power of two (up to 2^32) that keeps every total below 2^52.
  function setImplied(data, k) {
    let wmax = 1, inc = 0;
    for (let j = 0; j < data.nchar; j++) if (data.baseW[j] > wmax) wmax = data.baseW[j];
    for (let c = 0; c < data.nOrig; c++) if (data.charCols[c].length) inc++;
    // totals can be doubled by ratchet weights; a bootstrap replicate cannot exceed inc * wmax either
    const FIX = Math.pow(2, Math.max(16, Math.min(32, Math.floor(Math.log2(Math.pow(2, 52) / (4 * Math.max(1, inc) * wmax))))));
    const b = charBounds(data);
    const min = new Int32Array(data.nGroups);
    let maxLen = 1;
    for (let g = 0; g < data.nGroups; g++) {
      let s = 0;
      for (let q = 0; q < data.gLen[g]; q++) s += b.minS[data.gCol[g] + q];
      min[g] = s;
      if (data.gLen[g] > maxLen) maxLen = data.gLen[g];
    }
    const size = maxLen * data.ntax + 2;
    const F = new Float64Array(size);
    for (let e = 0; e < size; e++) F[e] = Math.round((FIX * e) / (e + k));
    data.iw = { k, F, min, FIX };
    data.scoreDiv = data.wScale * FIX;
    return data;
  }

  // Cost of a tree from its raw steps per column: weighted steps, or implied-weights fit.
  function cost(d, cs, w, naOnly) {
    let s = 0;
    if (!d.iw) {
      if (naOnly) { const nac = d.naChars; for (let i = 0; i < nac.length; i++) { const j = nac[i]; s += w[j] * cs[j]; } }
      else for (let j = 0; j < d.nchar; j++) s += w[j] * cs[j];
      return s;
    }
    const F = d.iw.F, min = d.iw.min, gCol = d.gCol, gLen = d.gLen;
    for (let g = naOnly ? d.nStdGroups : 0; g < d.nGroups; g++) {
      const c0 = gCol[g];
      let e = -min[g];
      for (let q = 0; q < gLen[g]; q++) e += cs[c0 + q];
      if (e > 0) s += w[c0] * F[e];
    }
    return s;
  }

  // Per-column weights from per-character multipliers (ratchet, bootstrap, jackknife).
  function expandWeights(data, mult) {
    const w = new Float64Array(data.nchar);
    for (let j = 0; j < data.nchar; j++) w[j] = data.baseW[j] * mult[data.colChar[j]];
    return w;
  }

  // Seeds for independent replicates: a function of (seed, index, purpose) only, so a
  // replicate gives the same result whichever worker runs it.
  function mixSeed(seed, a, b) {
    let h = ((seed >>> 0) ^ 0x9e3779b9) >>> 0;
    h = Math.imul(h ^ ((a + 0x85ebca6b) >>> 0), 0xc2b2ae35) >>> 0; h ^= h >>> 15;
    h = Math.imul(h ^ ((Math.imul(b + 1, 0x27d4eb2f)) >>> 0), 0x165667b1) >>> 0; h ^= h >>> 13;
    h = Math.imul(h, 0x9e3779b1) >>> 0;
    return ((h ^ (h >>> 16)) >>> 0) || 1;
  }

  /* ------------------------------------------------------------------ */
  /* Trees: rooted binary, tips 0..n-1, internal n..2n-2                 */
  /* ------------------------------------------------------------------ */
  function Tree(n) {
    this.n = n;
    const m = 2 * n - 1;
    this.parent = new Int32Array(m).fill(-1);
    this.left = new Int32Array(m).fill(-1);
    this.right = new Int32Array(m).fill(-1);
    this.root = -1;
  }
  Tree.prototype.clone = function () {
    const t = new Tree(this.n);
    t.parent.set(this.parent); t.left.set(this.left); t.right.set(this.right); t.root = this.root;
    return t;
  };

  // postorder list of internal nodes (children before parents)
  function postorder(tree, out) {
    const res = out || [];
    res.length = 0;
    const stack = [tree.root];
    const tmp = [];
    while (stack.length) {
      const v = stack.pop();
      if (v < tree.n) continue;
      tmp.push(v);
      stack.push(tree.left[v], tree.right[v]);
    }
    for (let i = tmp.length - 1; i >= 0; i--) res.push(tmp[i]);
    return res;
  }

  // Build a tree from nested arrays of tip indices, e.g. [0,[1,[2,3]]]; resolves polytomies left-to-right
  function fromNested(n, nested) {
    const t = new Tree(n);
    let next = n;
    function build(x) {
      if (typeof x === 'number') return x;
      if (x.length === 1) return build(x[0]);
      let acc = build(x[0]);
      for (let i = 1; i < x.length; i++) {
        const r = build(x[i]);
        const v = next++;
        t.left[v] = acc; t.right[v] = r; t.parent[acc] = v; t.parent[r] = v;
        acc = v;
      }
      return acc;
    }
    t.root = build(nested);
    t.parent[t.root] = -1;
    return t;
  }

  function toNested(tree, v) {
    if (v === undefined) v = tree.root;
    if (v < tree.n) return v;
    return [toNested(tree, tree.left[v]), toNested(tree, tree.right[v])];
  }

  // Reroot so that tip `out` is sister to the rest (root = parent of out)
  function rootOn(tree, out) {
    const nested = toNestedUnrooted(tree);
    return fromNested(tree.n, rerootNested(nested, out));
  }

  // Convert to adjacency then build rooted nested on tip `out`
  function adjacency(tree) {
    const m = 2 * tree.n - 1;
    const adj = Array.from({ length: m }, () => []);
    for (let v = 0; v < m; v++) {
      const p = tree.parent[v];
      if (p >= 0) { adj[v].push(p); adj[p].push(v); }
    }
    // suppress root of degree 2
    const r = tree.root;
    if (adj[r].length === 2) {
      const [a, b] = adj[r];
      adj[a] = adj[a].filter((x) => x !== r); adj[b] = adj[b].filter((x) => x !== r);
      adj[a].push(b); adj[b].push(a); adj[r] = [];
    }
    return adj;
  }
  function toNestedUnrooted(tree) { return adjacency(tree); }
  function rerootNested(adj, out) {
    const nb = adj[out][0];
    function rec(v, from) {
      if (adj[v].length <= 1 && v !== nb) return v;
      const kids = adj[v].filter((x) => x !== from).map((x) => rec(x, v));
      return kids.length === 1 ? kids[0] : kids;
    }
    // tips are numbered < n; leaf detection via degree
    const rest = rec(nb, out);
    return [out, rest];
  }

  function randomTree(n, rng) {
    const order = rng.shuffle([...Array(n).keys()]);
    const t = new Tree(n);
    const r = n, u = n + 1;
    const [a, b, c] = order;
    t.root = r; t.left[r] = a; t.right[r] = u; t.parent[a] = r; t.parent[u] = r;
    t.left[u] = b; t.right[u] = c; t.parent[b] = u; t.parent[c] = u;
    const present = [a, b, c, u];
    let next = n + 2;
    for (let i = 3; i < n; i++) {
      const e = present[rng.int(present.length)];
      const w = next++;
      insertTip(t, order[i], e, w);
      present.push(order[i], w);
    }
    return rootOn(t, 0);
  }

  // Insert tip `tip` on the edge above node `e`, using new internal node `w`
  function insertTip(t, tip, e, w) {
    const p = t.parent[e];
    t.parent[w] = p;
    if (p >= 0) { if (t.left[p] === e) t.left[p] = w; else t.right[p] = w; } else t.root = w;
    t.left[w] = e; t.right[w] = tip; t.parent[e] = w; t.parent[tip] = w;
  }

  /* ------------------------------------------------------------------ */
  /* Newick                                                              */
  /* ------------------------------------------------------------------ */
  function quoteName(s) {
    if (/^[A-Za-z0-9_.\-]+$/.test(s)) return s;
    if (/^[^'()\[\]:;,]+$/.test(s) && !/\s/.test(s)) return s;
    return "'" + s.replace(/'/g, "''") + "'";
  }
  function toNewick(tree, names, labels) {
    function rec(v) {
      if (v < tree.n) return quoteName(names[v]);
      const s = '(' + rec(tree.left[v]) + ',' + rec(tree.right[v]) + ')';
      return labels && labels[v] != null ? s + labels[v] : s;
    }
    return rec(tree.root) + ';';
  }
  // Nested (possibly multifurcating) consensus structure to Newick
  function nestedToNewick(node, names) {
    function rec(x) {
      if (typeof x.tip === 'number') return quoteName(names[x.tip]);
      const s = '(' + x.children.map(rec).join(',') + ')';
      return x.support != null ? s + x.support : s;
    }
    return rec(node) + ';';
  }

  // Parse Newick into nested arrays of tip indices (names must match taxa)
  function parseNewick(str, names, resolve) {
    const idx = new Map(names.map((n, i) => [n, i]));
    const idxLoose = new Map(names.map((n, i) => [n.replace(/[\s_]+/g, '_').toLowerCase(), i]));
    let i = 0;
    const s = str.trim();
    function name() {
      let out = '';
      if (s[i] === "'") {
        i++;
        while (i < s.length) { if (s[i] === "'") { if (s[i + 1] === "'") { out += "'"; i += 2; continue; } i++; break; } out += s[i++]; }
      } else {
        while (i < s.length && !/[(),:;\[]/.test(s[i])) out += s[i++];
      }
      return out.trim();
    }
    function skipMeta() {
      if (s[i] === '[') { while (i < s.length && s[i] !== ']') i++; i++; }
      if (s[i] === ':') { i++; while (i < s.length && /[0-9eE.+\-]/.test(s[i])) i++; }
      if (s[i] === '[') { while (i < s.length && s[i] !== ']') i++; i++; }
    }
    function node() {
      while (/\s/.test(s[i])) i++;
      if (s[i] === '(') {
        i++;
        const kids = [node()];
        while (s[i] === ',') { i++; kids.push(node()); }
        if (s[i] !== ')') throw new Error('Newick: expected ) at ' + i);
        i++;
        name(); skipMeta();
        return kids;
      }
      const nm = name();
      skipMeta();
      let t = resolve ? resolve(nm) : undefined;
      if (t === undefined) t = idx.get(nm);
      if (t === undefined) t = idxLoose.get(nm.replace(/[\s_]+/g, '_').toLowerCase());
      if (t === undefined) throw new Error('Newick: unknown taxon "' + nm + '"');
      return t;
    }
    return node();
  }

  /* ------------------------------------------------------------------ */
  /* Scorer                                                              */
  /* ------------------------------------------------------------------ */
  function Scorer(data) {
    const n = data.ntax, nc = data.nchar, m = 2 * n - 1;
    this.data = data;
    this.down1 = new Int32Array(m * nc);
    this.up1 = new Int32Array(m * nc);
    this.down2 = new Int32Array(m * nc);
    this.up2 = new Int32Array(m * nc);
    this.act = new Int32Array(m * nc);
    this.weights = data.baseW;
    this.cs = new Int32Array(nc);
    this.post = [];
    // tip rows
    for (let t = 0; t < n; t++) this.down1.set(data.tips.subarray(t * nc, t * nc + nc), t * nc);
  }

  // Raw (unweighted) steps per column, written to this.cs and returned.
  // In the second downpass a step is either a change between applicable states
  // or one more applicable region; columns flagged regionFree count only the
  // first kind (see encode, orderedRegions).
  Scorer.prototype.colSteps = function (tree, naOnly) {
    const d = this.data, nc = d.nchar, n = d.ntax;
    const D1 = this.down1, U1 = this.up1, D2 = this.down2, ACT = this.act, cs = this.cs, rf = d.regionFree;
    const std = naOnly ? EMPTY : d.stdChars, nac = d.naChars;
    const post = postorder(tree, this.post);
    const L = tree.left, R = tree.right, P = tree.parent;
    cs.fill(0);
    // First downpass
    for (let k = 0; k < post.length; k++) {
      const v = post[k], l = L[v] * nc, r = R[v] * nc, o = v * nc;
      for (let i = 0; i < std.length; i++) {
        const j = std[i];
        let x = D1[l + j] & D1[r + j];
        if (x === 0) { x = D1[l + j] | D1[r + j]; cs[j]++; }
        D1[o + j] = x;
      }
      for (let i = 0; i < nac.length; i++) {
        const j = nac[i];
        const a = D1[l + j], b = D1[r + j];
        let x = a & b;
        if (x === 0) {
          x = a | b;
          if ((a & ISAPPLIC) && (b & ISAPPLIC)) x &= ISAPPLIC;
        } else if (x === NA) {
          if ((a & ISAPPLIC) && (b & ISAPPLIC)) x = a | b;
        }
        D1[o + j] = x;
      }
    }
    if (nac.length === 0) return cs;
    // Root update (root treated as its own ancestor)
    const ro = tree.root * nc;
    for (let i = 0; i < nac.length; i++) {
      const j = nac[i];
      const x = D1[ro + j];
      const y = (x & ISAPPLIC) ? (x & ISAPPLIC) : NA;
      D1[ro + j] = y; U1[ro + j] = y;
    }
    // First uppass, internal nodes in preorder
    for (let k = post.length - 1; k >= 0; k--) {
      const v = post[k], o = v * nc, l = L[v] * nc, r = R[v] * nc;
      const a = (P[v] >= 0 ? P[v] : v) * nc;
      for (let i = 0; i < nac.length; i++) {
        const j = nac[i];
        const npre = D1[o + j];
        const anc = U1[a + j];
        let res;
        if (npre & NA) {
          if (npre & ISAPPLIC) {
            res = anc === NA ? NA : (npre & ISAPPLIC);
          } else if (anc === NA) {
            res = NA;
          } else {
            const u = (D1[l + j] | D1[r + j]) & ISAPPLIC;
            res = u ? u : NA;
          }
        } else res = npre;
        U1[o + j] = res;
      }
    }
    // Tip update
    for (let t = 0; t < n; t++) {
      if (P[t] < 0) continue; // tip not (yet) attached
      const o = t * nc, a = P[t] * nc;
      for (let i = 0; i < nac.length; i++) {
        const j = nac[i];
        const tp = D1[o + j], as = U1[a + j];
        ACT[o + j] = (tp & as) ? (tp & as & ISAPPLIC) : (tp & ISAPPLIC);
        let t2 = tp;
        if ((t2 & as) && (as & ISAPPLIC)) t2 &= ISAPPLIC;
        U1[o + j] = t2; D2[o + j] = t2;
      }
    }
    // Second downpass (counts steps)
    for (let k = 0; k < post.length; k++) {
      const v = post[k], o = v * nc, l = L[v] * nc, r = R[v] * nc;
      for (let i = 0; i < nac.length; i++) {
        const j = nac[i];
        const nif = U1[o + j], a = D2[l + j], b = D2[r + j];
        const la = ACT[l + j], ra = ACT[r + j];
        let npre;
        if (nif & ISAPPLIC) {
          const tmp = a & b;
          if (tmp) {
            npre = (tmp & ISAPPLIC) ? (tmp & ISAPPLIC) : tmp;
          } else {
            npre = (a | b) & ISAPPLIC;
            if ((a & ISAPPLIC) && (b & ISAPPLIC)) cs[j]++;
            else if (la && ra && !rf[j]) cs[j]++;
          }
        } else {
          npre = nif;
          if (la && ra && !rf[j]) cs[j]++;
        }
        D2[o + j] = npre;
        ACT[o + j] = (la | ra) & ISAPPLIC;
      }
    }
    return cs;
  };

  // Score of a tree: weighted length, or the implied-weights fit when data.iw is set
  // (both as integers; divide by data.scoreDiv to display). If `steps` (one entry per
  // column) is supplied, weighted steps per column are written to it.
  Scorer.prototype.length = function (tree, weights, steps, naOnly) {
    const cs = this.colSteps(tree, naOnly), w = weights || this.weights;
    if (steps) for (let j = 0; j < cs.length; j++) steps[j] = w[j] * cs[j];
    return cost(this.data, cs, w, naOnly);
  };

  // Unweighted steps per character of the matrix (columns of an ordered character summed).
  Scorer.prototype.charSteps = function (tree) {
    const d = this.data, cs = this.colSteps(tree), out = new Int32Array(d.nOrig);
    for (let j = 0; j < d.nchar; j++) out[d.colChar[j]] += cs[j];
    return out;
  };

  // Final (most parsimonious reconstruction) state sets for every node and column:
  // Int32Array[(2n-1) * ncolumns]. Standard columns: Fitch final sets. Inapplicable-type
  // columns: the final pass of the Brazeau et al. algorithm.
  Scorer.prototype.finalSets = function (tree) {
    const d = this.data, nc = d.nchar, n = d.ntax, m = 2 * n - 1;
    this.colSteps(tree); // fills D1 (and U1/D2 for inapplicable-type columns)
    const D1 = this.down1, D2 = this.down2, isNA = d.isNA;
    const L = tree.left, R = tree.right, P = tree.parent;
    const post = postorder(tree, []);
    const fin = new Int32Array(m * nc);
    const ro = tree.root * nc;
    for (let j = 0; j < nc; j++) fin[ro + j] = isNA[j] ? D2[ro + j] : D1[ro + j];
    for (let k = post.length - 1; k >= 0; k--) {
      const v = post[k];
      const o = v * nc, lo = L[v] * nc, r2 = R[v] * nc, ao = (v === tree.root ? v : P[v]) * nc;
      for (let j = 0; j < nc; j++) {
        let x;
        if (!isNA[j]) {
          if (v === tree.root) continue;
          const npre = D1[o + j], anc = fin[ao + j], l = D1[lo + j], r = D1[r2 + j];
          x = anc & npre;
          if (x !== anc) x = (l & r) ? (npre | (anc & (l | r))) : (npre | anc);
        } else {
          const npre = D2[o + j], anc = fin[ao + j], l = D2[lo + j], r = D2[r2 + j];
          if (npre & ISAPPLIC) {
            if (anc & ISAPPLIC) {
              if ((anc & npre) === anc) x = anc & npre;
              else if (l & r) x = npre | (anc & (l | r));
              else if ((l | r) & NA) x = ((l | r) & anc) ? anc : ((l | r | anc) & ISAPPLIC);
              else x = npre | anc;
            } else x = npre;
          } else x = npre;
        }
        fin[o + j] = x;
      }
    }
    for (let t = 0; t < n; t++) {
      if (P[t] < 0) continue;
      const o = t * nc, ao = P[t] * nc;
      for (let j = 0; j < nc; j++) { const tp = D1[o + j], a = fin[ao + j]; fin[o + j] = (tp & a) ? (tp & a) : tp; }
    }
    return fin;
  };

  // Ancestral state sets of one column, for character mapping.
  Scorer.prototype.reconstruct = function (tree, j) {
    const d = this.data, nc = d.nchar, m = 2 * d.ntax - 1;
    const all = this.finalSets(tree), fin = new Int32Array(m);
    for (let v = 0; v < m; v++) fin[v] = all[v * nc + j];
    return fin;
  };

  // State symbols of one character of the matrix at node v, from final sets `fin`
  // (as returned by finalSets). Returns '?' or an array of symbols, '-' = inapplicable.
  function charStatesAt(data, fin, v, c) {
    const cols = data.charCols[c], nc = data.nchar;
    if (!cols.length) return '?';
    if (!data.charOrdered[c]) return setToStates(fin[v * nc + cols[0]], data);
    // ordered: states consistent with every threshold column
    const LOW = 1 << data.push, HIGH = 1 << (data.push + 1), lo = data.charLo[c];
    let miss = true, na = false, appl = false;
    for (const j of cols) { const x = fin[v * nc + j]; if (x !== MISSING) miss = false; if (x & NA && data.isNA[j]) na = true; if (x & (LOW | HIGH)) appl = true; }
    if (miss) return '?';
    const out = [];
    if (na) out.push('-');
    if (appl) {
      const before = out.length;
      for (let s = lo; s <= lo + cols.length; s++) {
        let ok = true;
        for (let q = 0; q < cols.length && ok; q++) { const x = fin[v * nc + cols[q]]; if (!(x & (s >= lo + q + 1 ? HIGH : LOW))) ok = false; }
        if (ok) out.push(data.symbols[s]);
      }
      if (out.length === before) { // thresholds disagree (not expected): fall back to the count of "high" columns
        let hi = 0; for (const j of cols) if (fin[v * nc + j] & HIGH) hi++;
        out.push(data.symbols[lo + hi]);
      }
    }
    return out;
  }

  /* ------------------------------------------------------------------ */
  /* Incremental scorer for inapplicable-type columns                    */
  /* ------------------------------------------------------------------ */
  // Holds the four MorphyLib pass arrays for a "base" tree. evalDelta(tree)
  // recomputes only nodes whose inputs changed relative to the base (dirty
  // propagation through the four passes), journalling old values so revert()
  // restores the base exactly; commit() makes the evaluated tree the new base.
  // Results are identical to Scorer.length(tree, w, null, true).
  function NAInc(data) {
    const n = data.ntax, m = 2 * n - 1, k = data.naChars.length, nc = data.nchar;
    this.n = n; this.m = m; this.k = k; this.data = data;
    this.A = [new Int32Array(m * k), new Int32Array(m * k), new Int32Array(m * k), new Int32Array(m * k), new Int32Array(m * k)];
    [this.D1, this.U1, this.D2, this.ACT, this.SC] = this.A; // SC = steps (0/1) per cell
    this.w = new Float64Array(k);
    this.rf = new Uint8Array(k);
    for (let i = 0; i < k; i++) this.rf[i] = data.regionFree[data.naChars[i]];
    // implied weighting: steps per character (group of columns)
    const ng = data.nGroups - data.nStdGroups;
    this.grp = new Int32Array(k); this.gs = new Int32Array(ng); this.gw = new Float64Array(ng); this.gmin = new Int32Array(ng);
    for (let g = 0; g < ng; g++) {
      const G = data.nStdGroups + g, s = data.gStart[G] - data.stdChars.length;
      for (let q = 0; q < data.gLen[G]; q++) this.grp[s + q] = g;
    }
    this.bL = new Int32Array(m); this.bR = new Int32Array(m); this.bP = new Int32Array(m); this.bRoot = -1;
    // per-node lists of columns whose value changed in this evaluation
    this.L1 = new Int32Array(m * k); this.C1 = new Int32Array(m);
    this.LU = new Int32Array(m * k); this.CU = new Int32Array(m);
    this.L2 = new Int32Array(m * k); this.C2 = new Int32Array(m);
    this.sD = new Uint8Array(m);
    this.mark = new Int32Array(k); this.epoch = 0; this.cand = new Int32Array(k);
    this.touched = []; this.post = [];
    this.jA = new Uint8Array(4096); this.jI = new Int32Array(4096); this.jO = new Int32Array(4096); this.jn = 0;
    this.total = 0;
    for (let t = 0; t < n; t++) for (let i = 0; i < k; i++) this.D1[t * k + i] = data.tips[t * nc + data.naChars[i]];
  }
  NAInc.prototype.full = function (tree, weights) {
    const d = this.data, k = this.k, n = this.n, nac = d.naChars, iw = d.iw;
    const W = weights || d.baseW;
    for (let i = 0; i < k; i++) this.w[i] = W[nac[i]];
    const post = postorder(tree, this.post);
    const L = tree.left, R = tree.right, P = tree.parent, root = tree.root;
    const D1 = this.D1, U1 = this.U1, D2 = this.D2, ACT = this.ACT, SC = this.SC, w = this.w, rf = this.rf;
    // cells of nodes that are not in this tree (a partial tree during stepwise addition) must not count
    SC.fill(0);
    for (const v of post) { const o = v * k, lo = L[v] * k, ro = R[v] * k; for (let i = 0; i < k; i++) D1[o + i] = d1cell(D1[lo + i], D1[ro + i], v === root); }
    for (let q = post.length - 1; q >= 0; q--) {
      const v = post[q], o = v * k, lo = L[v] * k, ro = R[v] * k, ao = (P[v] >= 0 ? P[v] : v) * k;
      for (let i = 0; i < k; i++) U1[o + i] = u1cell(D1[o + i], v === root ? D1[o + i] : U1[ao + i], D1[lo + i], D1[ro + i]);
    }
    for (let t = 0; t < n; t++) {
      if (P[t] < 0) continue;
      const o = t * k, ao = P[t] * k;
      for (let i = 0; i < k; i++) { const tp = D1[o + i], as = U1[ao + i]; ACT[o + i] = tipAct(tp, as); D2[o + i] = tipD2(tp, as); }
    }
    let total = 0;
    if (iw) this.gs.fill(0);
    for (const v of post) {
      const o = v * k, lo = L[v] * k, ro = R[v] * k;
      for (let i = 0; i < k; i++) {
        d2cell(U1[o + i], D2[lo + i], D2[ro + i], ACT[lo + i], ACT[ro + i]);
        const sc = STEPOUT === 1 || (STEPOUT === 2 && !rf[i]) ? 1 : 0;
        D2[o + i] = D2OUT; ACT[o + i] = ACTOUT; SC[o + i] = sc;
        if (sc) { if (iw) this.gs[this.grp[i]]++; else total += w[i]; }
      }
    }
    if (iw) {
      const F = iw.F, ns = d.nStdGroups;
      for (let g = 0; g < this.gs.length; g++) {
        this.gw[g] = W[d.gCol[ns + g]]; this.gmin[g] = iw.min[ns + g];
        const e = this.gs[g] - this.gmin[g];
        if (e > 0) total += this.gw[g] * F[e];
      }
    }
    this.total = total;
    this.bL.set(L); this.bR.set(R); this.bP.set(P); this.bRoot = root;
    return total;
  };
  // Cell formulas (MorphyLib passes); d2cell writes results to module globals.
  // STEPOUT: 0 = no step, 1 = change between applicable states, 2 = one more applicable region.
  let D2OUT = 0, ACTOUT = 0, STEPOUT = 0;
  function d1cell(a, b, isRoot) {
    let x = a & b;
    if (x === 0) { x = a | b; if ((a & ISAPPLIC) && (b & ISAPPLIC)) x &= ISAPPLIC; }
    else if (x === NA) { if ((a & ISAPPLIC) && (b & ISAPPLIC)) x = a | b; }
    if (isRoot) x = (x & ISAPPLIC) ? (x & ISAPPLIC) : NA;
    return x;
  }
  function u1cell(npre, anc, l, r) {
    if (npre & NA) {
      if (npre & ISAPPLIC) return anc === NA ? NA : (npre & ISAPPLIC);
      if (anc === NA) return NA;
      const u = (l | r) & ISAPPLIC; return u ? u : NA;
    }
    return npre;
  }
  function tipAct(tp, as) { return (tp & as) ? (tp & as & ISAPPLIC) : (tp & ISAPPLIC); }
  function tipD2(tp, as) { return ((tp & as) && (as & ISAPPLIC)) ? (tp & ISAPPLIC) : tp; }
  function d2cell(nif, a, b, la, ra) {
    let npre, st = 0;
    if (nif & ISAPPLIC) {
      const tmp = a & b;
      if (tmp) npre = (tmp & ISAPPLIC) ? (tmp & ISAPPLIC) : tmp;
      else { npre = (a | b) & ISAPPLIC; if ((a & ISAPPLIC) && (b & ISAPPLIC)) st = 1; else if (la && ra) st = 2; }
    } else { npre = nif; if (la && ra) st = 2; }
    D2OUT = npre; ACTOUT = (la | ra) & ISAPPLIC; STEPOUT = st;
    return 0;
  }
  NAInc.prototype.set = function (ai, idx, val) {
    const arr = this.A[ai];
    if (arr[idx] === val) return false;
    if (this.jn === this.jA.length) {
      const N = this.jn * 2;
      const a = new Uint8Array(N); a.set(this.jA); this.jA = a;
      const b = new Int32Array(N); b.set(this.jI); this.jI = b;
      const c = new Int32Array(N); c.set(this.jO); this.jO = c;
    }
    this.jA[this.jn] = ai; this.jI[this.jn] = idx; this.jO[this.jn] = arr[idx]; this.jn++;
    arr[idx] = val;
    return true;
  };
  // Union of changed-column lists (list, count, node triples) into this.cand
  NAInc.prototype.gather = function (A1, C1, v1, A2, C2, v2, A3, C3, v3, A4, C4, v4) {
    const k = this.k, mark = this.mark, cand = this.cand;
    const ep = ++this.epoch;
    let c = 0, cnt, o;
    cnt = C1[v1]; o = v1 * k; for (let z = 0; z < cnt; z++) { const i = A1[o + z]; if (mark[i] !== ep) { mark[i] = ep; cand[c++] = i; } }
    if (A2) { cnt = C2[v2]; o = v2 * k; for (let z = 0; z < cnt; z++) { const i = A2[o + z]; if (mark[i] !== ep) { mark[i] = ep; cand[c++] = i; } } }
    if (A3) { cnt = C3[v3]; o = v3 * k; for (let z = 0; z < cnt; z++) { const i = A3[o + z]; if (mark[i] !== ep) { mark[i] = ep; cand[c++] = i; } } }
    if (A4) { cnt = C4[v4]; o = v4 * k; for (let z = 0; z < cnt; z++) { const i = A4[o + z]; if (mark[i] !== ep) { mark[i] = ep; cand[c++] = i; } } }
    return c;
  };
  NAInc.prototype.evalDelta = function (tree) {
    const n = this.n, k = this.k, L = tree.left, R = tree.right, P = tree.parent, root = tree.root;
    if (root !== this.bRoot) throw new Error('NAInc: root changed');
    const bL = this.bL, bR = this.bR, bP = this.bP, sD = this.sD;
    const D1 = this.D1, U1 = this.U1, D2 = this.D2, ACT = this.ACT, SC = this.SC, w = this.w, rf = this.rf;
    const L1 = this.L1, C1 = this.C1, LU = this.LU, CU = this.CU, L2 = this.L2, C2 = this.C2;
    const cand = this.cand, touched = this.touched;
    const iw = this.data.iw, F = iw ? iw.F : null, grp = this.grp, gs = this.gs, gw = this.gw, gmin = this.gmin;
    const post = postorder(tree, this.post);
    this.savedTotal = this.total;
    let total = this.total;
    this.jn = 0;
    // child order is irrelevant (all formulas are symmetric in left/right)
    for (const v of post) { const l = L[v], r = R[v]; if (!((l === bL[v] && r === bR[v]) || (l === bR[v] && r === bL[v]))) { sD[v] = 1; touched.push(v); } }
    const allChars = () => { for (let i = 0; i < k; i++) cand[i] = i; return k; };
    // first downpass
    for (const v of post) {
      const l = L[v], r = R[v];
      if (!sD[v] && !C1[l] && !C1[r]) continue;
      const c = sD[v] ? allChars() : this.gather(L1, C1, l, L1, C1, r);
      if (!c) continue;
      const o = v * k, lo = l * k, ro = r * k, isRoot = v === root;
      let cnt = 0;
      for (let z = 0; z < c; z++) { const i = cand[z]; if (this.set(0, o + i, d1cell(D1[lo + i], D1[ro + i], isRoot))) L1[o + cnt++] = i; }
      if (cnt) { C1[v] = cnt; touched.push(v); }
    }
    // first uppass, root first
    for (let q = post.length - 1; q >= 0; q--) {
      const v = post[q], p = P[v], l = L[v], r = R[v];
      let c;
      const pc = p >= 0 && p !== bP[v];
      if (!sD[v] && !pc && !C1[v] && !C1[l] && !C1[r] && !(p >= 0 && CU[p])) continue;
      if (sD[v] || pc) c = allChars();
      else c = p >= 0 ? this.gather(L1, C1, v, L1, C1, l, L1, C1, r, LU, CU, p) : this.gather(L1, C1, v, L1, C1, l, L1, C1, r);
      if (!c) continue;
      const o = v * k, lo = l * k, ro = r * k, ao = (p >= 0 ? p : v) * k;
      let cnt = 0;
      for (let z = 0; z < c; z++) {
        const i = cand[z], npre = D1[o + i];
        if (this.set(1, o + i, u1cell(npre, p >= 0 ? U1[ao + i] : npre, D1[lo + i], D1[ro + i]))) LU[o + cnt++] = i;
      }
      if (cnt) { CU[v] = cnt; touched.push(v); }
    }
    // tips
    for (let t = 0; t < n; t++) {
      const p = P[t];
      if (p < 0) continue;
      let c;
      if (p !== bP[t]) c = allChars(); else { if (!CU[p]) continue; c = this.gather(LU, CU, p); }
      const o = t * k, ao = p * k;
      let cnt = 0;
      for (let z = 0; z < c; z++) {
        const i = cand[z], tp = D1[o + i], as = U1[ao + i];
        const ch1 = this.set(3, o + i, tipAct(tp, as)), ch2 = this.set(2, o + i, tipD2(tp, as));
        if (ch1 || ch2) L2[o + cnt++] = i;
      }
      if (cnt) { C2[t] = cnt; touched.push(t); }
    }
    // second downpass
    for (const v of post) {
      const l = L[v], r = R[v];
      if (!sD[v] && !CU[v] && !C2[l] && !C2[r]) continue;
      const c = sD[v] ? allChars() : this.gather(LU, CU, v, L2, C2, l, L2, C2, r);
      if (!c) continue;
      const o = v * k, lo = l * k, ro = r * k;
      let cnt = 0;
      for (let z = 0; z < c; z++) {
        const i = cand[z];
        d2cell(U1[o + i], D2[lo + i], D2[ro + i], ACT[lo + i], ACT[ro + i]);
        const ch1 = this.set(2, o + i, D2OUT), ch2 = this.set(3, o + i, ACTOUT);
        const sc = STEPOUT === 1 || (STEPOUT === 2 && !rf[i]) ? 1 : 0;
        if (SC[o + i] !== sc) {
          const dl = sc - SC[o + i];
          if (iw) {
            const g = grp[i], e0 = gs[g] - gmin[g], e1 = e0 + dl;
            gs[g] += dl;
            total += gw[g] * ((e1 > 0 ? F[e1] : 0) - (e0 > 0 ? F[e0] : 0));
          } else total += dl * w[i];
          this.set(4, o + i, sc);
        }
        if (ch1 || ch2) L2[o + cnt++] = i;
      }
      if (cnt) { C2[v] = cnt; touched.push(v); }
    }
    this.total = total;
    return total;
  };
  NAInc.prototype.clearFlags = function () {
    for (const v of this.touched) { this.C1[v] = 0; this.CU[v] = 0; this.C2[v] = 0; this.sD[v] = 0; }
    this.touched.length = 0;
  };
  NAInc.prototype.revert = function () {
    const A = this.A, jA = this.jA, jI = this.jI, jO = this.jO, k = this.k;
    const iw = this.data.iw, SC = this.SC;
    for (let q = this.jn - 1; q >= 0; q--) {
      if (iw && jA[q] === 4) this.gs[this.grp[jI[q] % k]] += jO[q] - SC[jI[q]];
      A[jA[q]][jI[q]] = jO[q];
    }
    this.jn = 0;
    this.total = this.savedTotal;
    this.clearFlags();
  };
  NAInc.prototype.commit = function (tree) {
    this.bL.set(tree.left); this.bR.set(tree.right); this.bP.set(tree.parent); this.bRoot = tree.root;
    this.jn = 0;
    this.clearFlags();
  };

  /* ------------------------------------------------------------------ */
  /* Character statistics                                                */
  /* ------------------------------------------------------------------ */
  // Minimum and maximum (star tree) Fitch steps for standard characters.
  function charBounds(data) {
    const n = data.ntax, nc = data.nchar;
    const minS = new Int32Array(nc), maxS = new Int32Array(nc);
    for (let j = 0; j < nc; j++) {
      const sets = [];
      for (let t = 0; t < n; t++) {
        let v = data.tips[t * nc + j];
        if (v === MISSING) continue;
        if (data.push) v &= ISAPPLIC; // bit 0 is the inapplicable token only when states start at bit 1
        if (v === 0) continue; // inapplicable
        sets.push(v >>> 0);
      }
      // min = minimum hitting set size - 1
      const uniq = [...new Set(sets)];
      let best = 0;
      if (uniq.length) {
        let all = 0; for (const s of uniq) all |= s;
        const bits = []; for (let b = 0; b < 31; b++) if (all & (1 << b)) bits.push(1 << b);
        best = bits.length;
        const nb = bits.length;
        if (nb <= 16) {
          for (let mask = 1; mask < (1 << nb); mask++) {
            let pc = 0, cov = 0;
            for (let b = 0; b < nb; b++) if (mask & (1 << b)) { pc++; cov |= bits[b]; }
            if (pc >= best) continue;
            if (uniq.every((s) => s & cov)) best = pc;
          }
        }
      }
      minS[j] = Math.max(0, best - 1);
      // max on a bush: choose root state hitting most tips
      let bestHit = 0, all2 = 0;
      for (const s of sets) all2 |= s;
      for (let b = 0; b < 31; b++) {
        if (!(all2 & (1 << b))) continue;
        let h = 0; for (const s of sets) if (s & (1 << b)) h++;
        if (h > bestHit) bestHit = h;
      }
      maxS[j] = sets.length - bestHit;
    }
    return { minS, maxS };
  }

  // Per character of the matrix: minimum and maximum steps (columns summed), and
  // whether it is scored with the inapplicable algorithm.
  function charStats(data) {
    const b = charBounds(data), minC = new Int32Array(data.nOrig), maxC = new Int32Array(data.nOrig);
    for (let j = 0; j < data.nchar; j++) { minC[data.colChar[j]] += b.minS[j]; maxC[data.colChar[j]] += b.maxS[j]; }
    return { minC, maxC };
  }

  /* ------------------------------------------------------------------ */
  /* Search                                                              */
  /* ------------------------------------------------------------------ */
  // Remove subtree rooted at s (s not root, not child of root). Returns {sib, p}
  // where p is the freed internal node and sib the former sibling.
  function prune(t, s) {
    const p = t.parent[s];
    const sib = t.left[p] === s ? t.right[p] : t.left[p];
    const g = t.parent[p];
    if (g >= 0) { if (t.left[g] === p) t.left[g] = sib; else t.right[g] = sib; }
    else t.root = sib;
    t.parent[sib] = g;
    t.parent[p] = -1;
    // p keeps s as one child; other slot freed
    t.left[p] = s; t.right[p] = -1;
    return { sib, p };
  }
  // Regraft subtree (held at freed node p with left child s) onto edge above e
  function regraft(t, p, e) {
    const g = t.parent[e];
    t.parent[p] = g;
    if (g >= 0) { if (t.left[g] === e) t.left[g] = p; else t.right[g] = p; } else t.root = p;
    t.right[p] = e; t.parent[e] = p;
  }

  function collectNodes(t, v, out) { // nodes in subtree of v
    const st = [v];
    while (st.length) { const x = st.pop(); out.push(x); if (x >= t.n) { st.push(t.left[x], t.right[x]); } }
    return out;
  }

  // Re-root the subtree whose root is s so that its root edge lies above node
  // `e` (a node inside the subtree, not s, not a child of s). Node s is reused
  // as the new subtree root; other node ids keep their identity.
  function rerootSubtree(t, s, e, fGiven) {
    const nodes = collectNodes(t, s, []);
    const adj = new Map();
    for (const v of nodes) adj.set(v, []);
    const c1 = t.left[s], c2 = t.right[s];
    for (const v of nodes) {
      if (v === s) continue;
      const p = t.parent[v];
      if (p === s) continue;
      adj.get(v).push(p); adj.get(p).push(v);
    }
    adj.get(c1).push(c2); adj.get(c2).push(c1); // join across suppressed s
    // e's neighbour on the far side of the new root edge
    let f = fGiven;
    if (f === undefined) {
      f = t.parent[e];
      if (f === s) f = e === c1 ? c2 : c1;
    }
    function build(v, from, parentId) {
      t.parent[v] = parentId;
      if (v < t.n) return;
      const kids = adj.get(v).filter((x) => x !== from);
      t.left[v] = kids[0]; t.right[v] = kids[1];
      build(kids[0], v, v); build(kids[1], v, v);
    }
    t.left[s] = e; t.right[s] = f;
    build(e, f, s); build(f, e, s);
    t.parent[s] = -1;
    return s;
  }

  function Searcher(data, opts) {
    this.data = data;
    this.n = data.ntax;
    this.scorer = new Scorer(data);
    this.opts = Object.assign({ maxTrees: 100, seed: 1 }, opts || {});
    this.rng = new RNG(this.opts.seed);
    this.evaluations = 0; // candidate trees considered
    this.naEvals = 0;     // of which needed the inapplicable algorithm
    this.constraint = null;
    const n = this.n, m = 2 * n - 1, nw = data.nchar, wk = data.work;
    this.nw = nw;
    this.Dm = new Int32Array(m * nw); this.Um = new Int32Array(m * nw); this.X = new Int32Array(m * nw);
    this.bs = new Int32Array(nw); this.jw = new Float64Array(nw);
    this.pre = []; this.stk = [];
    // tip rows of the directional sets never change: screening tips in work order
    for (let t = 0; t < n; t++) for (let i = 0; i < nw; i++) this.Dm[t * nw + i] = data.screenTips[t * nw + wk[i]];
    this.mg = []; // implied weighting: ordered standard characters (several columns)
    this.baseStd = 0; this.baseScr = 0;
  }

  Searcher.prototype.score = function (t, w) { this.evaluations++; return this.scorer.length(t, w); };

  /* ---- topological constraints ---- */
  // constraint = { pos: [tip-index arrays], neg: [tip-index arrays] }. A set that
  // contains tip 0 (the root tip) is replaced by its complement, which is the
  // same bipartition. pos: the set must be a clade; neg: it must not be.
  Searcher.prototype.setConstraint = function (c) {
    const n = this.n;
    const norm = (tips) => {
      const s = new Uint8Array(n); let k = 0;
      for (const x of tips) if (!s[x]) { s[x] = 1; k++; }
      if (s[0]) { for (let i = 0; i < n; i++) s[i] = s[i] ? 0 : 1; k = n - k; }
      return { mask: s, size: k };
    };
    const pos = ((c && c.pos) || []).map(norm).filter((x) => x.size > 1 && x.size < n - 1);
    const neg = ((c && c.neg) || []).map(norm).filter((x) => x.size > 1 && x.size < n - 1);
    this.constraint = pos.length || neg.length ? { pos, neg } : null;
    this.cnt = new Int32Array(2 * n - 1); this.siz = new Int32Array(2 * n - 1);
    return this.constraint;
  };
  // Is the tip set a clade of t? With `present` (Uint8Array), only tips in the
  // tree so far are considered (stepwise addition).
  Searcher.prototype.isClade = function (t, c, present) {
    const n = this.n, cnt = this.cnt, siz = this.siz, mask = c.mask;
    let need = 0, inTree = 0;
    for (let i = 0; i < n; i++) {
      if (present && !present[i]) { cnt[i] = 0; siz[i] = 0; continue; }
      inTree++; siz[i] = 1; cnt[i] = mask[i]; need += mask[i];
    }
    if (need <= 1 || need >= inTree - 1) return true; // trivially a clade among the tips present (tip 0 is outside)
    const post = postorder(t, this.pre);
    for (let q = 0; q < post.length; q++) {
      const v = post[q], l = t.left[v], r = t.right[v];
      cnt[v] = cnt[l] + cnt[r]; siz[v] = siz[l] + siz[r];
      if (cnt[v] === need) return siz[v] === need;
    }
    return false;
  };
  Searcher.prototype.ok = function (t, present) {
    const C = this.constraint;
    if (!C) return true;
    for (const c of C.pos) if (!this.isClade(t, c, present)) return false;
    if (!present) for (const c of C.neg) if (this.isClade(t, c)) return false;
    return true;
  };
  // Move one tip of a forbidden clade elsewhere, so that the tree starts out
  // without the clade. Returns false if no allowed move was found.
  Searcher.prototype.breakClades = function (t) {
    const C = this.constraint, n = this.n, rng = this.rng;
    if (!C) return true;
    for (const c of C.neg) {
      let tries = 0;
      while (this.isClade(t, c) && tries++ < 200) {
        const inside = []; for (let i = 1; i < n; i++) if (c.mask[i] && t.parent[i] !== t.root) inside.push(i);
        if (!inside.length) break;
        const x = inside[rng.int(inside.length)];
        const { sib, p } = prune(t, x);
        const targets = collectNodes(t, t.root, []).filter((e) => e !== t.root && e !== 0 && e !== sib);
        let done = false;
        rng.shuffle(targets);
        for (const e of targets) {
          regraft(t, p, e);
          if (!this.isClade(t, c) && C.pos.every((q) => this.isClade(t, q))) { done = true; break; }
          const g = t.parent[p];
          if (t.left[g] === p) t.left[g] = e; else t.right[g] = e;
          t.parent[e] = g; t.parent[p] = -1; t.right[p] = -1;
        }
        if (!done) regraft(t, p, sib);
      }
      if (this.isClade(t, c)) return false;
    }
    return true;
  };

  // Directional Fitch sets on the tree hanging from `top` (top has two children),
  // for every column in work order, using the screening tips. D[v] = set of the
  // subtree below v; U[v] = set of everything on the other side of v's parent
  // edge. Steps are added to bs (one counter per column).
  Searcher.prototype.dirSets = function (t, top, D, U, bs) {
    const nw = this.nw, n = this.n;
    const pre = this.pre, st = this.stk;
    pre.length = 0; st.length = 0; st.push(top);
    while (st.length) { const v = st.pop(); pre.push(v); if (v >= n) st.push(t.left[v], t.right[v]); }
    for (let k = pre.length - 1; k >= 0; k--) {
      const v = pre[k];
      if (v < n) continue;
      const o = v * nw, l = t.left[v] * nw, r = t.right[v] * nw;
      for (let i = 0; i < nw; i++) {
        let x = D[l + i] & D[r + i];
        if (x === 0) { x = D[l + i] | D[r + i]; bs[i]++; }
        D[o + i] = x;
      }
    }
    for (let k = 1; k < pre.length; k++) {
      const v = pre[k], pp = t.parent[v];
      const sib = t.left[pp] === v ? t.right[pp] : t.left[pp];
      const o = v * nw, so = sib * nw;
      if (pp === top) { for (let i = 0; i < nw; i++) U[o + i] = D[so + i]; }
      else {
        const po = pp * nw;
        for (let i = 0; i < nw; i++) { let x = D[so + i] & U[po + i]; if (x === 0) x = D[so + i] | U[po + i]; U[o + i] = x; }
      }
    }
  };

  // From the steps of the two parts (bs): the exact cost of the standard columns
  // before joining (baseStd), the screening cost of the inapplicable-type columns
  // (baseScr, a lower bound), and jw[i] = what one more step in column i adds.
  Searcher.prototype.clipCosts = function (W) {
    const d = this.data, wk = d.work, ns = d.stdChars.length, nw = this.nw, bs = this.bs, jw = this.jw;
    let baseStd = 0, baseScr = 0;
    this.mg.length = 0;
    if (!d.iw) {
      for (let i = 0; i < ns; i++) { const wi = W[wk[i]]; jw[i] = wi; baseStd += wi * bs[i]; }
      for (let i = ns; i < nw; i++) { const wi = W[wk[i]]; jw[i] = wi; baseScr += wi * bs[i]; }
    } else {
      const F = d.iw.F, min = d.iw.min;
      for (let g = 0; g < d.nGroups; g++) {
        const s = d.gStart[g], len = d.gLen[g], wg = W[d.gCol[g]];
        let e = -min[g];
        if (!d.gIsNA[g]) {
          for (let q = 0; q < len; q++) e += bs[s + q];
          baseStd += wg * (e > 0 ? F[e] : 0);
          if (len === 1) jw[s] = wg * ((e + 1 > 0 ? F[e + 1] : 0) - (e > 0 ? F[e] : 0));
          else { for (let q = 0; q < len; q++) jw[s + q] = 0; this.mg.push(s, len, e, wg); }
        } else { // only the first column of the character is screened
          e += bs[s];
          baseScr += wg * (e > 0 ? F[e] : 0);
          jw[s] = wg * ((e + 1 > 0 ? F[e + 1] : 0) - (e > 0 ? F[e] : 0));
          for (let q = 1; q < len; q++) jw[s + q] = 0;
        }
      }
    }
    this.baseStd = baseStd; this.baseScr = baseScr;
  };
  // Cost added by joining the edges with sets X[a..] and X[b..]: standard columns only.
  Searcher.prototype.joinStd = function (a, b) {
    const X = this.X, jw = this.jw, ns = this.data.stdChars.length, mg = this.mg;
    let sc = 0;
    for (let i = 0; i < ns; i++) if ((X[a + i] & X[b + i]) === 0) sc += jw[i];
    if (mg.length) {
      const F = this.data.iw.F;
      for (let z = 0; z < mg.length; z += 4) {
        const s = mg[z], len = mg[z + 1], e = mg[z + 2];
        let c = 0;
        for (let q = 0; q < len; q++) if ((X[a + s + q] & X[b + s + q]) === 0) c++;
        if (c) sc += mg[z + 3] * ((e + c > 0 ? F[e + c] : 0) - (e > 0 ? F[e] : 0));
      }
    }
    return sc;
  };
  Searcher.prototype.joinScr = function (a, b) {
    const X = this.X, jw = this.jw, nw = this.nw;
    let sc = 0;
    for (let i = this.data.stdChars.length; i < nw; i++) if ((X[a + i] & X[b + i]) === 0) sc += jw[i];
    return sc;
  };
  Searcher.prototype.edgeSets = function (list) {
    const nw = this.nw, Dm = this.Dm, Um = this.Um, X = this.X;
    for (const e of list) {
      const o = e * nw;
      for (let i = 0; i < nw; i++) { let x = Dm[o + i] & Um[o + i]; if (x === 0) x = Dm[o + i] | Um[o + i]; X[o + i] = x; }
    }
  };

  // Random addition sequence: taxa are added in random order, each at the best
  // place found (ties broken at random). Under topological constraints only
  // places that keep the constraints are tried.
  Searcher.prototype.wagner = function (w) {
    const n = this.n, rng = this.rng, d = this.data, nw = this.nw;
    const W = w || d.baseW;
    const hasNA = d.naChars.length > 0, C = this.constraint;
    const order = [0].concat(rng.shuffle([...Array(n - 1).keys()].map((x) => x + 1)));
    // Initial 3-taxon tree; other tips absent (postorder only visits attached nodes)
    const t = new Tree(n);
    const a = order[0], b = order[1], c = order[2];
    let next = n;
    const r = next++, u = next++;
    t.root = r; t.left[r] = a; t.right[r] = u; t.parent[a] = r; t.parent[u] = r;
    t.left[u] = b; t.right[u] = c; t.parent[b] = u; t.parent[c] = u;
    const present = [a, b, c, r, u];
    const inTree = new Uint8Array(n); inTree[a] = inTree[b] = inTree[c] = 1;
    let NAS = null;
    if (hasNA) { if (!this.na) this.na = new NAInc(d); NAS = this.na; }
    const Dm = this.Dm, X = this.X;
    for (let i = 3; i < n; i++) {
      const tip = order[i];
      const w2 = next++;
      const edges = present.filter((v) => v !== t.root && v !== a);
      this.bs.fill(0);
      this.dirSets(t, t.root, Dm, this.Um, this.bs);
      this.clipCosts(W);
      this.edgeSets(edges);
      if (NAS) NAS.full(t, W);
      const to = tip * nw;
      for (let q = 0; q < nw; q++) X[to + q] = Dm[to + q];
      inTree[tip] = 1;
      let best = Infinity, bestE = [];
      for (const e of edges) {
        this.evaluations++;
        let sc = this.baseStd + this.joinStd(to, e * nw);
        if (sc > best) continue;
        if (NAS && sc + this.baseScr + this.joinScr(to, e * nw) > best) continue;
        if (!NAS && !C) { if (sc < best) { best = sc; bestE = [e]; } else bestE.push(e); continue; }
        insertTip(t, tip, e, w2);
        if (!C || this.ok(t, inTree)) {
          if (NAS) { this.naEvals++; sc += NAS.evalDelta(t); NAS.revert(); }
          if (sc < best) { best = sc; bestE = [e]; } else if (sc === best) bestE.push(e);
        }
        const g = t.parent[w2];
        if (t.left[g] === w2) t.left[g] = e; else t.right[g] = e;
        t.parent[e] = g; t.parent[w2] = -1; t.parent[tip] = -1; t.left[w2] = -1; t.right[w2] = -1;
      }
      if (!bestE.length) throw new Error('The constraints cannot all be satisfied together.');
      insertTip(t, tip, bestE[rng.int(bestE.length)], w2);
      present.push(tip, w2);
    }
    if (C && C.neg.length && !this.breakClades(t)) throw new Error('No starting tree satisfies the constraints.');
    return t;
  };

  // One sweep of SPR or TBR rearrangements in random order, returning the first
  // better tree found (null if there is none). Carrying on with the sweep after
  // an improvement was tried and was no faster to the same local optima. Standard columns are scored exactly by indirect calculation
  // (cost of the two parts + one step per column whose Fitch sets at the two
  // joined edges are disjoint). Inapplicable-type columns are first screened
  // with a lower bound (gaps read as missing), and only candidates that could
  // still be good enough are scored with the inapplicable algorithm.
  // `keep(tree, score)` is called for every equally good tree met.
  // With `accept(score, curScore)` (tree drifting) a worse tree may be taken.
  // With `scan` = { best, tree, test(tree) } nothing is accepted: the whole
  // neighbourhood is visited and the best tree that passes test() is recorded.
  Searcher.prototype.rearrange = function (t, curScore, w, mode, deadline, keep, accept, scan) {
    const n = this.n, m = 2 * n - 1, rng = this.rng, d = this.data, nw = this.nw;
    const W = w || d.baseW;
    const hasNA = d.naChars.length > 0, C = this.constraint;
    const Dm = this.Dm, Um = this.Um, X = this.X;
    const verify = this.verify;
    let NAS = null;
    if (hasNA && !this.noIncremental) {
      if (!this.na) this.na = new NAInc(d);
      NAS = this.na;
      NAS.full(t, W);
    }
    const nodes = [];
    for (let v = 0; v < m; v++) if (v !== t.root && t.parent[v] !== t.root) nodes.push(v);
    rng.shuffle(nodes);
    let tick = 0;
    for (const s of nodes) {
      if (deadline && (++tick & 7) === 0 && Date.now() > deadline) return null;
      if (t.parent[s] < 0 || s === t.root || t.parent[s] === t.root) continue;
      const { sib, p } = prune(t, s);
      const main = collectNodes(t, t.root, []).filter((e) => e !== t.root && e !== 0);
      const sub = mode === 'tbr' && s >= n ? collectNodes(t, s, []) : [s];
      const c1 = t.left[s], c2 = t.right[s];
      const origPar = sub.map((v) => t.parent[v]);
      this.bs.fill(0);
      this.dirSets(t, t.root, Dm, Um, this.bs);
      this.dirSets(t, s, Dm, Um, this.bs);
      this.clipCosts(W);
      this.edgeSets(main);
      {
        const o = s * nw;
        for (let i = 0; i < nw; i++) X[o + i] = Dm[o + i];
        const inner = [];
        for (let q = 0; q < sub.length; q++) if (sub[q] !== s && origPar[q] !== s) inner.push(sub[q]);
        this.edgeSets(inner);
      }
      const baseStd = this.baseStd, baseScr = this.baseScr;
      let curQ = 0; // index in `sub` of the current physical rooting (sub[0] === s)
      const ensure = (q) => {
        if (q === curQ) return;
        if (sub[q] === s) rerootSubtree(t, s, c1, c2); else rerootSubtree(t, s, sub[q], origPar[q]);
        t.parent[s] = p; curQ = q;
      };
      t.parent[s] = p;
      for (let q = 0; q < sub.length; q++) {
        const rs = sub[q];
        if (rs !== s && origPar[q] === s) continue;
        rng.shuffle(main);
        const a = rs * nw;
        for (const e of main) {
          if (e === sib && rs === s) continue;
          this.evaluations++;
          const b = e * nw;
          let sc = baseStd + this.joinStd(a, b);
          if (scan) curScore = scan.best;
          // candidates that cannot beat (or, when collecting, equal) the current tree are dropped
          if (!verify && !accept && (keep ? sc > curScore : sc >= curScore)) continue;
          let lb = sc;
          if (hasNA) {
            lb += baseScr;
            if (!verify && !accept && (keep ? lb > curScore : lb >= curScore)) continue;
            lb += this.joinScr(a, b);
            if (!verify && !accept && (keep ? lb > curScore : lb >= curScore)) continue;
          }
          if (accept && !verify && lb >= curScore && !accept(lb, curScore, true)) continue; // too much worse to be taken
          ensure(q);
          regraft(t, p, e);
          if (hasNA) { this.naEvals++; sc += NAS ? NAS.evalDelta(t) : this.scorer.length(t, W, null, true); }
          if (verify) {
            const full = this.scorer.length(t, W);
            if (full !== sc) throw new Error('indirect score ' + sc + ' != full score ' + full);
            if (lb > sc) throw new Error('screening bound ' + lb + ' > exact score ' + sc);
          }
          if (scan) { if (sc < scan.best && scan.test(t)) { scan.best = sc; scan.tree = t.clone(); } }
          else if (accept) {
            if ((sc < curScore || accept(sc, curScore, false)) && (!C || this.ok(t))) { if (NAS) NAS.commit(t); return { tree: t, score: sc }; }
          } else if (sc < curScore && (!C || this.ok(t))) {
            if (NAS) NAS.commit(t);
            return { tree: t, score: sc };
          }
          if (NAS) NAS.revert();
          if (keep && sc === curScore && (!C || this.ok(t))) keep(t, sc);
          const g = t.parent[p];
          if (t.left[g] === p) t.left[g] = e; else t.right[g] = e;
          t.parent[e] = g; t.parent[p] = -1; t.right[p] = -1;
        }
      }
      ensure(0);
      regraft(t, p, sib);
    }
    // Moves of the outgroup (tip 0) itself: equivalent to re-rooting the
    // ingroup clade C on each of its edges while tip 0 stays at the root.
    const r = this.outgroupMoves(t, curScore, w, deadline, keep, scan);
    if (r) return r;
    return null;
  };

  Searcher.prototype.outgroupMoves = function (t, curScore, w, deadline, keep, scan) {
    const Cn = t.left[t.root] === 0 ? t.right[t.root] : t.left[t.root];
    if (Cn < this.n) return null;
    const c1 = t.left[Cn], c2 = t.right[Cn];
    const inner = collectNodes(t, Cn, []).filter((v) => v !== Cn && t.parent[v] !== Cn);
    const par = inner.map((v) => t.parent[v]);
    const order = this.rng.shuffle([...inner.keys()]);
    let moved = false;
    for (const q of order) {
      if (deadline && Date.now() > deadline) break;
      rerootSubtree(t, Cn, inner[q], par[q]);
      t.parent[Cn] = t.root;
      moved = true;
      const sc = this.score(t, w);
      if (scan) { if (sc < scan.best && scan.test(t)) { scan.best = sc; scan.tree = t.clone(); } continue; }
      if (sc < curScore && this.ok(t)) return { tree: t, score: sc };
      if (keep && sc === curScore && this.ok(t)) keep(t, sc);
    }
    if (moved) { rerootSubtree(t, Cn, c1, c2); t.parent[Cn] = t.root; }
    return null;
  };

  // Hill-climb to a local optimum. Returns {tree, score}
  Searcher.prototype.climb = function (t, w, mode, deadline, keep) {
    let score = this.score(t, w);
    while (true) {
      const r = this.rearrange(t, score, w, mode, deadline, keep);
      if (!r) break;
      t = r.tree; score = r.score;
      if (!w && this.onImprove) this.onImprove(t, score);
      if (deadline && Date.now() > deadline) break;
    }
    return { tree: t, score };
  };

  // Tree drifting (Goloboff 1999, simplified): a short walk of SPR moves in which
  // equally good trees are always taken and slightly worse ones sometimes, to
  // leave a local optimum; the caller then climbs back with TBR.
  Searcher.prototype.drift = function (t, score, moves, deadline) {
    const rng = this.rng, d = this.data;
    // typical cost of one step, to scale "slightly worse"
    let unit = 0, cnt = 0;
    for (let j = 0; j < d.nchar; j++) if (d.baseW[j] > 0) { unit += d.baseW[j]; cnt++; }
    unit = cnt ? unit / cnt : 1;
    if (d.iw) unit *= d.iw.F[1] / 4;
    const accept = (sc, cur, boundOnly) => {
      const worse = (sc - cur) / unit;
      if (worse <= 0) return true;
      if (worse > 3) return false;
      if (boundOnly) return true; // could still be acceptable: needs the exact score
      return rng.next() < 1 / (2 + 3 * worse);
    };
    for (let i = 0; i < moves; i++) {
      if (deadline && Date.now() > deadline) break;
      const r = this.rearrange(t, score, null, 'spr', deadline, null, accept);
      if (!r) break;
      t = r.tree; score = r.score;
    }
    return { tree: t, score };
  };

  /* ------------------------------------------------------------------ */
  /* Tree identity, splits, consensus                                    */
  /* ------------------------------------------------------------------ */
  // Clades (as sorted tip-index strings) of a tree rooted on tip 0; excludes trivial clades
  function clades(t) {
    const n = t.n;
    const words = Math.ceil(n / 32);
    const sets = new Map();
    const res = [];
    const post = postorder(t, []);
    const bits = new Array(2 * n - 1);
    for (let i = 0; i < n; i++) { const b = new Uint32Array(words); b[i >> 5] |= 1 << (i & 31); bits[i] = b; }
    for (const v of post) {
      const b = new Uint32Array(words);
      const l = bits[t.left[v]], r = bits[t.right[v]];
      for (let k = 0; k < words; k++) b[k] = l[k] | r[k];
      bits[v] = b;
      let cnt = 0; for (let k = 0; k < words; k++) cnt += popcnt(b[k]);
      if (cnt > 1 && cnt < n - 1 && !(b[0] & 1)) res.push(Array.from(b).join('.'));
    }
    void sets;
    return res;
  }
  function popcnt(x) { x = x - ((x >>> 1) & 0x55555555); x = (x & 0x33333333) + ((x >>> 2) & 0x33333333); return (((x + (x >>> 4)) & 0x0f0f0f0f) * 0x01010101) >>> 24; }
  function treeKey(t) {
    // canonical: sorted clade list; tree assumed rooted on tip 0
    return clades(t).sort().join('|');
  }

  function cladeTips(key, n) {
    const w = key.split('.').map(Number);
    const tips = [];
    for (let i = 0; i < n; i++) if (w[i >> 5] & (1 << (i & 31))) tips.push(i);
    return tips;
  }

  // Build a (multifurcating) tree from compatible clades. counts: Map key->frequency (0..1)
  function buildFromClades(n, keys, freq) {
    const list = keys.map((k) => ({ key: k, tips: cladeTips(k, n) }))
      .sort((a, b) => b.tips.length - a.tips.length);
    const rootNode = { children: [], tips: [...Array(n).keys()] };
    const owner = new Array(n).fill(rootNode);
    const nodes = [];
    for (const c of list) {
      // all tips must share the same current owner
      const o = owner[c.tips[0]];
      if (!c.tips.every((x) => owner[x] === o)) continue; // incompatible (shouldn't happen)
      const node = { children: [], tips: c.tips, support: freq ? freq.get(c.key) : null, key: c.key };
      for (const x of c.tips) owner[x] = node;
      node.parent = o;
      nodes.push(node);
    }
    for (const node of nodes) node.parent.children.push(node);
    for (let i = 0; i < n; i++) owner[i].children.push({ tip: i });
    // order children: tips & clades sorted by smallest tip index for stable display
    function minTip(x) { return typeof x.tip === 'number' ? x.tip : Math.min(...x.tips); }
    function sortRec(x) { if (x.children) { x.children.sort((a, b) => minTip(a) - minTip(b)); x.children.forEach(sortRec); } }
    sortRec(rootNode);
    return rootNode;
  }

  function consensus(trees, n, level) {
    // level: 1 = strict, 0.5 = majority-rule
    const counts = new Map();
    for (const t of trees) for (const c of new Set(clades(t))) counts.set(c, (counts.get(c) || 0) + 1);
    const keep = [], freq = new Map();
    for (const [k, v] of counts) {
      const f = v / trees.length;
      if ((level >= 1 && v === trees.length) || (level < 1 && f > level)) { keep.push(k); freq.set(k, f); }
    }
    return buildFromClades(n, keep, freq);
  }

  function nestedFromTree(t) { // binary tree -> consensus-style node structure
    function rec(v) {
      if (v < t.n) return { tip: v };
      return { children: [rec(t.left[v]), rec(t.right[v])], node: v };
    }
    return rec(t.root);
  }

  /* ------------------------------------------------------------------ */
  /* High-level analysis driver                                          */
  /* ------------------------------------------------------------------ */
  // A search is a series of independent replicates, numbered 0, 1, 2, ...
  // Replicate i depends only on (seed, i): random addition sequence, SPR then
  // TBR, then rounds of the parsimony ratchet (Nixon 1999) alternating with
  // tree drifting. The search ends at the first replicate by which the best
  // score has been found `hits` times (or after maxReps replicates), so the
  // result does not depend on the speed of the computer or on how many
  // workers shared the replicates. The time limit only cuts a search short,
  // and a search cut short is reported as such.
  function searchOpts(opts) {
    const o = Object.assign({ ratchetIter: 10, mode: 'tbr', maxTrees: 100, seed: 1, timeLimit: 60000, fuse: true, drift: true, keepPerRep: 20 }, opts || {});
    if (o.hits == null || o.hits <= 0) { o.hits = 0; o.maxReps = o.maxReps || o.replicates || 5; }
    else if (!o.maxReps) o.maxReps = Math.max(o.hits * 10, 20);
    return o;
  }

  // Number of replicates (in index order) that make up the search, or 0 if that
  // cannot be told yet. scores[i] is undefined while replicate i is not finished.
  function repsNeeded(scores, hits, maxReps) {
    let best = Infinity, count = 0;
    for (let i = 0; i < maxReps; i++) {
      const s = scores[i];
      if (s === undefined || s === null) return 0;
      if (s < best) { best = s; count = 1; } else if (s === best) count++;
      if (hits > 0 && count >= hits) return i + 1;
    }
    return maxReps;
  }

  /* ---- sectorial search ---- */
  // Random sectorial search (after Goloboff 1999). A sector is a clade of
  // moderate size. It is cut out as a small data set: its own tips plus one
  // extra tip that stands for the rest of the tree. For standard columns that
  // tip carries the Fitch set of the rest of the tree, so a change of length in
  // the sector is exactly the change of length of the whole tree; for
  // inapplicable-type columns it carries the reconstruction just outside the
  // sector, which is only an approximation. New trees for the sector are built
  // from scratch and by drifting, and a better one is put back only if the
  // whole tree, scored in full, is better.
  function sectorData(S, tree, v, tipList, fin) {
    const d = S.data, nc = d.nchar, nw = S.nw, k = tipList.length, wk = d.work;
    const tips = new Int32Array((k + 1) * nc), screen = new Int32Array((k + 1) * nc);
    const p = tree.parent[v];
    for (let i = 0; i < nw; i++) {
      const j = wk[i];
      let x;
      if (!d.isNA[j]) x = S.Um[v * nw + i];
      else x = fin[p * nc + j];
      tips[j] = x;
      screen[j] = d.isNA[j] ? ((x & NA) || d.regionFree[j] ? MISSING : x) : x;
    }
    for (let q = 0; q < k; q++) {
      const t = tipList[q];
      tips.set(d.tips.subarray(t * nc, t * nc + nc), (q + 1) * nc);
      screen.set(d.screenTips.subarray(t * nc, t * nc + nc), (q + 1) * nc);
    }
    return Object.assign({}, d, { ntax: k + 1, tips, screenTips: screen, iw: d.iw ? { k: d.iw.k, F: d.iw.F, FIX: d.iw.FIX, min: d.iw.min.slice() } : null });
  }
  // Tries `count` sectors on cur = {tree, score}. Returns the best tree found (cur if none).
  function sectorial(S, cur, count, deadline) {
    if (S.constraint) return cur;
    const n = S.n, d = S.data, rng = S.rng;
    const lo = Math.max(6, Math.round(0.15 * n)), hi = Math.min(Math.max(lo, Math.round(0.5 * n)), 45);
    if (n < 2 * lo + 2) return cur;
    for (let it = 0; it < count; it++) {
      if (deadline && Date.now() > deadline) break;
      const t = cur.tree, size = new Int32Array(2 * n - 1);
      for (let i = 0; i < n; i++) size[i] = 1;
      const cands = [];
      for (const v of postorder(t, [])) { size[v] = size[t.left[v]] + size[t.right[v]]; if (v !== t.root && t.parent[v] !== t.root && size[v] >= lo && size[v] <= hi) cands.push(v); }
      if (!cands.length) break;
      const v = cands[rng.int(cands.length)];
      const tipList = collectNodes(t, v, []).filter((x) => x < n).sort((a, b) => a - b);
      const map = new Map(tipList.map((x, i) => [x, i + 1]));
      S.bs.fill(0);
      S.dirSets(t, t.root, S.Dm, S.Um, S.bs);
      const fin = d.naChars.length ? S.scorer.finalSets(t) : null;
      const red = sectorData(S, t, v, tipList, fin);
      const toRed = (x) => (x < n ? map.get(x) : [toRed(t.left[x]), toRed(t.right[x])]);
      const start = fromNested(tipList.length + 1, [0, toRed(v)]);
      if (red.iw) { // steps outside the sector still count towards each character's homoplasy
        const full = Int32Array.from(S.scorer.colSteps(t)), part = new Scorer(red).colSteps(start);
        for (let g = 0; g < d.nGroups; g++) { let rest = 0; for (let q = 0; q < d.gLen[g]; q++) rest += Math.max(0, full[d.gCol[g] + q] - part[d.gCol[g] + q]); red.iw.min[g] -= rest; }
      }
      const R = new Searcher(red, { seed: 1 });
      R.rng = new RNG(rng.int(0x7fffffff) + 1);
      R.verify = S.verify;
      const s0 = R.score(start);
      let best = null, bestScore = s0;
      const tryTree = (c) => { if (c.score < bestScore) { bestScore = c.score; best = c.tree; } };
      for (let r = 0; r < 2; r++) {
        let c = R.climb(R.wagner(), null, 'spr', deadline);
        c = R.climb(c.tree, null, 'tbr', deadline);
        tryTree(c);
      }
      {
        const dr = R.drift(start.clone(), s0, Math.max(6, tipList.length >> 2), deadline);
        let c = R.climb(dr.tree, null, 'spr', deadline);
        c = R.climb(c.tree, null, 'tbr', deadline);
        tryTree(c);
      }
      S.evaluations += R.evaluations; S.naEvals += R.naEvals;
      if (!best) continue;
      const rt = rootOn(best, 0);
      const inner = rt.left[rt.root] === 0 ? rt.right[rt.root] : rt.left[rt.root];
      const back = (x) => (x <= tipList.length ? tipList[x - 1] : [back(rt.left[x]), back(rt.right[x])]);
      const cand = fromNested(n, replaceNested(t, v, back(inner)));
      const sc = S.score(cand);
      if (S.verify && !d.naChars.length && !d.iw && sc - cur.score !== bestScore - s0) throw new Error('sector: change in the sector ' + (bestScore - s0) + ' != change in the whole tree ' + (sc - cur.score));
      if (sc < cur.score) cur = { tree: cand, score: sc, improved: true };
    }
    return cur;
  }

  function replicate(S, index, opts, deadline, onStep) {
    S.rng = new RNG(mixSeed(opts.seed, index, 1));
    const d = S.data, mode = opts.mode || 'tbr', rng = S.rng;
    const pool = new Pool(opts.keepPerRep || 20);
    const keeper = (t, sc) => pool.offer(t, sc);
    let r = S.climb(S.wagner(), null, 'spr', deadline);
    if (mode === 'tbr') r = S.climb(r.tree, null, 'tbr', deadline);
    pool.offer(r.tree, r.score);
    let cur = r;
    // Sectorial search is off unless asked for: on the benchmark matrices (up to 88 taxa)
    // it did not raise the share of replicates that reach the best score (VALIDATION.md).
    const sectors = opts.sectors || 0;
    const sect = () => {
      if (!sectors) return;
      const x = sectorial(S, cur, sectors, deadline);
      if (x !== cur) { cur = mode === 'tbr' ? S.climb(x.tree, null, 'tbr', deadline) : x; pool.offer(cur.tree, cur.score); }
    };
    sect();
    if (onStep) onStep({ phase: 'search', replicate: index + 1, iteration: 0, best: pool.best });
    const mult = new Float64Array(d.nOrig);
    const driftMoves = Math.max(8, Math.round(d.ntax / 4));
    for (let it = 0; it < (opts.ratchetIter || 0); it++) {
      if (deadline && Date.now() > deadline) break;
      let pert;
      if (opts.drift && (it & 1)) pert = S.drift(cur.tree.clone(), cur.score, driftMoves, deadline);
      else {
        for (let c = 0; c < d.nOrig; c++) mult[c] = rng.next() < 0.25 ? 2 : 1; // upweight about a quarter of the characters
        pert = S.climb(cur.tree.clone(), expandWeights(d, mult), 'spr', deadline);
      }
      const back = S.climb(pert.tree, null, mode, deadline, keeper);
      if (back.score <= cur.score) cur = back;
      pool.offer(back.tree, back.score);
      if (it % 3 === 2) sect();
      if (onStep) onStep({ phase: 'ratchet', replicate: index + 1, iteration: it + 1, best: pool.best });
    }
    return { index, score: pool.best, trees: pool.list(), timedOut: !!(deadline && Date.now() > deadline) };
  }

  // Single-thread search. opts: {hits, maxReps | replicates, ratchetIter, mode, maxTrees,
  // seed, timeLimit (ms), constraint, drift, fuse, onProgress, searchOnly}
  function analyse(data, opts) {
    opts = searchOpts(opts);
    const S = new Searcher(data, { seed: opts.seed, maxTrees: opts.maxTrees });
    if (opts.constraint) S.setConstraint(opts.constraint);
    const start = Date.now(), deadline = start + opts.timeLimit;
    const results = [];
    let best = Infinity;
    const progress = (o) => { if (o.best < best) best = o.best; if (opts.onProgress) opts.onProgress(Object.assign(o, { best, elapsed: Date.now() - start, evaluations: S.evaluations, fraction: Math.min(1, results.length / Math.max(1, opts.hits ? opts.hits : opts.maxReps)) })); };
    let used = 0;
    for (let i = 0; i < opts.maxReps; i++) {
      if (Date.now() > deadline) break;
      results.push(replicate(S, i, opts, deadline, progress));
      used = repsNeeded(results.map((r) => r.score), opts.hits, opts.maxReps);
      if (used) break;
    }
    const use = used ? results.slice(0, used) : results;
    const cut = !used || use.some((r) => r.timedOut);
    if (opts.searchOnly) return { results: use, evaluations: S.evaluations, naEvals: S.naEvals, elapsed: Date.now() - start, timedOut: cut };
    const fin = finishTrees(S, use, opts, deadline, progress);
    return {
      bestScore: fin.bestScore, trees: fin.trees, hits: fin.hits, reps: use.length, repScores: use.map((r) => r.score),
      evaluations: S.evaluations, naEvals: S.naEvals, elapsed: Date.now() - start, timedOut: cut || fin.timedOut,
    };
  }

  // Pool of distinct equally-short trees (rooted on tip 0)
  function Pool(max) { this.max = max; this.best = Infinity; this.map = new Map(); }
  Pool.prototype.offer = function (t, sc) {
    if (sc > this.best) return false;
    if (sc < this.best) { this.best = sc; this.map = new Map(); }
    if (this.map.size >= this.max) return false;
    const rt = rootOn(t, 0);
    const k = treeKey(rt);
    if (this.map.has(k)) return false;
    this.map.set(k, rt);
    return true;
  };
  Pool.prototype.list = function () { return [...this.map.values()]; };
  Pool.prototype.size = function () { return this.map.size; };

  // After the replicates: fuse their results (Goloboff 1999), then swap on every
  // optimal tree to collect all equally good trees. `results` are replicate
  // results in index order (trees as Tree objects or nested arrays).
  function finishTrees(S, results, opts, deadline, progress) {
    const n = S.n;
    S.rng = new RNG(mixSeed(opts.seed, 0, 2));
    const P = new Pool(opts.maxTrees || 100);
    const asTree = (t) => (Array.isArray(t) ? fromNested(n, t) : t);
    const repBests = [];
    for (const r of results) {
      const ts = r.trees.map(asTree);
      for (const t of ts) P.offer(t, r.score);
      if (ts.length) repBests.push(ts[0]);
    }
    finish(S, P, repBests, opts, deadline, progress);
    let hits = 0;
    for (const r of results) if (r.score === P.best) hits++;
    return { bestScore: P.best, trees: P.list(), hits, timedOut: !!(deadline && Date.now() > deadline) };
  }

  function finish(S, P, repBests, opts, deadline, progress) {
    if (opts.fuse !== false && Date.now() < deadline) {
      const donors = repBests.concat(P.list().slice(0, 20));
      if (donors.length > 1) {
        const f = fuse(S, donors, deadline);
        if (f.score <= P.best) {
          const c = S.climb(f.tree, null, opts.mode || 'tbr', deadline, (t, sc) => P.offer(t, sc));
          P.offer(c.tree, c.score);
        }
        if (progress) progress({ phase: 'fuse', best: P.best, trees: P.size() });
      }
    }
    let queue = P.list();
    while (queue.length && P.size() < P.max && Date.now() < deadline) {
      const t = queue.shift();
      const before = P.best;
      const better = collectEqual(S, t, P.best, P.map, P.max, queue, deadline);
      if (better) {
        const c = S.climb(better.tree, null, opts.mode || 'tbr', deadline);
        P.offer(c.tree, c.score);
        queue = P.list();
      } else if (P.best !== before) queue = P.list();
      if (progress) progress({ phase: 'collect', best: P.best, trees: P.size() });
    }
  }

  // Merge results of independent searches and finish (kept for scripts that used it).
  // results: [{trees: [nested], best | score}]
  function mergeAndFinish(data, results, opts) {
    opts = searchOpts(opts);
    const S = new Searcher(data, { seed: opts.seed });
    if (opts.constraint) S.setConstraint(opts.constraint);
    const start = Date.now(), deadline = start + opts.timeLimit;
    const rs = results.map((r, i) => {
      const trees = r.trees.map((t) => (Array.isArray(t) ? fromNested(data.ntax, t) : t));
      return { index: i, trees, score: r.score != null ? r.score : (r.best != null ? r.best : (trees.length ? S.score(trees[0]) : Infinity)) };
    });
    const fin = finishTrees(S, rs, opts, deadline, opts.onProgress);
    return { bestScore: fin.bestScore, trees: fin.trees, hits: fin.hits, evaluations: S.evaluations, elapsed: Date.now() - start, timedOut: fin.timedOut };
  }

  // Tree fusing (Goloboff 1999): exchange clades with identical taxon sets
  // between trees whenever the exchange shortens the recipient.
  function fuse(S, trees, deadline) {
    const n = S.n;
    let scored = trees.map((t) => ({ tree: rootOn(t, 0), score: S.score(t) })).sort((a, b) => a.score - b.score);
    let best = scored[0];
    let improved = true, rounds = 0;
    while (improved && rounds < 20 && Date.now() < deadline) {
      improved = false; rounds++;
      const bInfo = subtreeInfo(best.tree);
      for (let d = 0; d < scored.length && !improved; d++) {
        const donor = scored[d].tree;
        if (donor === best.tree) continue;
        const dInfo = subtreeInfo(donor);
        for (const [key, node] of bInfo.byKey) {
          const dn = dInfo.byKey.get(key);
          if (dn === undefined) continue;
          if (bInfo.size.get(key) < 4) continue;
          if (bInfo.shape.get(node) === dInfo.shape.get(dn)) continue;
          const nested = replaceNested(best.tree, node, toNested(donor, dn));
          const cand = fromNested(n, nested);
          const sc = S.score(cand);
          if (sc < best.score && S.ok(cand)) { best = { tree: cand, score: sc }; improved = true; break; }
        }
      }
    }
    return best;
  }
  // For each internal node: clade key, size, and a canonical shape string
  function subtreeInfo(t) {
    const n = t.n, byKey = new Map(), size = new Map(), shape = new Map();
    const words = Math.ceil(n / 32);
    const bits = new Array(2 * n - 1), sh = new Array(2 * n - 1), cnt = new Int32Array(2 * n - 1);
    for (let i = 0; i < n; i++) { const b = new Uint32Array(words); b[i >> 5] |= 1 << (i & 31); bits[i] = b; sh[i] = String(i); cnt[i] = 1; }
    for (const v of postorder(t, [])) {
      const l = t.left[v], r = t.right[v], b = new Uint32Array(words);
      for (let q = 0; q < words; q++) b[q] = bits[l][q] | bits[r][q];
      bits[v] = b; cnt[v] = cnt[l] + cnt[r];
      const a = sh[l], c = sh[r];
      sh[v] = a < c ? '(' + a + ',' + c + ')' : '(' + c + ',' + a + ')';
      if (v !== t.root) { const k = Array.from(b).join('.'); byKey.set(k, v); size.set(k, cnt[v]); shape.set(v, sh[v]); }
    }
    return { byKey, size, shape, bits };
  }
  function replaceNested(t, node, repl) {
    function rec(v) { if (v === node) return repl; if (v < t.n) return v; return [rec(t.left[v]), rec(t.right[v])]; }
    return rec(t.root);
  }

  // Swap (full TBR + outgroup moves) on tree t0, adding every equally short
  // tree to the pool. Returns a strictly shorter tree if one is met, else null.
  function collectEqual(S, t0, best, pool, maxTrees, queue, deadline) {
    const keep = (tt, sc) => {
      if (sc !== best || pool.size >= maxTrees) return;
      const rt = rootOn(tt, 0);
      const k = treeKey(rt);
      if (!pool.has(k)) { pool.set(k, rt); queue.push(rt); }
    };
    const r = S.rearrange(t0.clone(), best, null, 'tbr', deadline, keep);
    return r && r.score < best ? r : null;
  }

  /* ---- resampling ---- */
  // One bootstrap or jackknife replicate (characters of the matrix are resampled;
  // an ordered character keeps all its columns together). Returns the clades of its tree.
  function resampleRep(S, index, opts, deadline) {
    const d = S.data;
    S.rng = new RNG(mixSeed(opts.seed, index, 3));
    const rng = S.rng, mult = new Float64Array(d.nOrig);
    const inc = []; for (let c = 0; c < d.nOrig; c++) if (d.charCols[c].length) inc.push(c);
    if (opts.method === 'jackknife') { for (const c of inc) mult[c] = rng.next() < Math.exp(-1) ? 0 : 1; }
    else for (let q = 0; q < inc.length; q++) mult[inc[rng.int(inc.length)]]++;
    const w = expandWeights(d, mult);
    const r = S.climb(S.wagner(w), w, opts.mode || 'spr', deadline);
    return clades(rootOn(r.tree, 0));
  }
  // Resampling support. opts.indices: which replicates to run (default 0..replicates-1).
  // Returns Map cladeKey -> frequency.
  function resample(data, opts) {
    opts = Object.assign({ method: 'bootstrap', replicates: 100, seed: 7, mode: 'spr', timeLimit: 60000 }, opts || {});
    const S = new Searcher(data, { seed: opts.seed });
    if (opts.constraint) S.setConstraint({ pos: opts.constraint.pos });
    const deadline = Date.now() + opts.timeLimit;
    const counts = new Map();
    const idx = opts.indices || [...Array(opts.replicates).keys()];
    let done = 0;
    for (const i of idx) {
      if (Date.now() > deadline) break;
      for (const c of resampleRep(S, i, opts, deadline)) counts.set(c, (counts.get(c) || 0) + 1);
      done++;
      if (opts.onProgress) opts.onProgress({ phase: opts.method, replicate: done, of: idx.length, fraction: done / idx.length });
    }
    const freq = new Map();
    for (const [k, v] of counts) freq.set(k, v / Math.max(1, done));
    return { freq, counts, replicates: done, timedOut: done < idx.length };
  }

  /* ---- Bremer support ---- */
  // Decay index by converse constraints: for each clade, the best tree found that
  // does not contain it. Every tree met in any of these searches also updates the
  // bound of each other clade it lacks. The searches are heuristic, so a value is
  // an upper bound on the true Bremer support.
  // opts: { clades: [keys], indices (which to search; default all), startTrees: [nested MPTs],
  //         reps, ratchetIter, seed, timeLimit, constraint, onProgress }
  // Returns { bounds: [[key, bestScoreWithoutClade]], done: [indices], skipped: [indices], timedOut }
  function bremer(data, opts) {
    opts = Object.assign({ reps: 2, ratchetIter: 2, seed: 1, timeLimit: 60000 }, opts || {});
    const S = new Searcher(data, { seed: opts.seed });
    const n = data.ntax, all = opts.clades;
    const deadline = Date.now() + opts.timeLimit;
    const bound = new Map(all.map((k) => [k, Infinity]));
    const note = (tree, score) => {
      const has = new Set(clades(rootOn(tree, 0)));
      for (const k of all) if (!has.has(k) && score < bound.get(k)) bound.set(k, score);
    };
    const userPos = (opts.constraint && opts.constraint.pos) || [];
    const mult = new Float64Array(data.nOrig);
    const done = [], skipped = [];
    const idx = opts.indices || [...all.keys()];
    for (const ci of idx) {
      if (Date.now() > deadline) break;
      const tips = cladeTips(all[ci], n);
      S.setConstraint({ pos: userPos, neg: [tips] });
      let any = false;
      for (let rep = 0; rep < opts.reps; rep++) {
        if (Date.now() > deadline) break;
        S.rng = new RNG(mixSeed(opts.seed, ci * 64 + rep, 4));
        let t = null;
        try {
          if (rep === 0 && opts.startTrees && opts.startTrees.length) {
            // the best tree without the clade among all TBR neighbours of an optimal tree
            const t0 = rootOn(fromNested(n, opts.startTrees[ci % opts.startTrees.length]), 0);
            const scan = { best: Infinity, tree: null, test: (x) => S.ok(x) };
            S.rearrange(t0, Infinity, null, 'tbr', deadline, null, null, scan);
            t = scan.tree;
          }
          if (!t) t = S.wagner();
        } catch (e) { t = null; } // the clade is enforced by a constraint: no tree lacks it
        if (!t) break;
        any = true;
        let cur = S.climb(t, null, 'spr', deadline);
        cur = S.climb(cur.tree, null, 'tbr', deadline);
        note(cur.tree, cur.score);
        for (let it = 0; it < opts.ratchetIter; it++) {
          if (Date.now() > deadline) break;
          for (let c = 0; c < data.nOrig; c++) mult[c] = S.rng.next() < 0.25 ? 2 : 1;
          const pert = S.climb(cur.tree.clone(), expandWeights(data, mult), 'spr', deadline);
          const back = S.climb(pert.tree, null, 'tbr', deadline);
          note(back.tree, back.score);
          if (back.score <= cur.score) cur = back;
        }
      }
      (any ? done : skipped).push(ci);
      if (opts.onProgress) opts.onProgress({ phase: 'bremer', replicate: done.length + skipped.length, of: idx.length, fraction: (done.length + skipped.length) / idx.length });
    }
    return { bounds: [...bound], done, skipped, timedOut: done.length + skipped.length < idx.length };
  }

  // Compare analyses of the same taxa (e.g. the three gap treatments).
  // sets: [{name, trees: [Tree rooted on tip 0]}]. For each analysis the strict
  // consensus of its trees is taken; every clade found in any consensus is
  // reported with the analyses that recover it. pairwise[i][j] = number of
  // clades found in exactly one of consensus i and consensus j (Robinson-Foulds
  // distance between the consensus trees, counting clades).
  function compareTreatments(sets, n) {
    const strict = sets.map((s) => {
      const counts = new Map();
      for (const t of s.trees) for (const c of new Set(clades(t))) counts.set(c, (counts.get(c) || 0) + 1);
      return new Set([...counts].filter(([, v]) => v === s.trees.length).map(([k]) => k));
    });
    const all = new Set();
    strict.forEach((ks) => ks.forEach((c) => all.add(c)));
    const out = [...all].map((key) => ({ key, tips: cladeTips(key, n), present: strict.map((ks) => ks.has(key)) }))
      .sort((a, b) => b.tips.length - a.tips.length || a.tips[0] - b.tips[0]);
    const pairwise = sets.map(() => sets.map(() => 0));
    for (let i = 0; i < sets.length; i++) {
      for (let j = 0; j < sets.length; j++) {
        if (i === j) continue;
        let d = 0;
        for (const c of strict[i]) if (!strict[j].has(c)) d++;
        for (const c of strict[j]) if (!strict[i].has(c)) d++;
        pairwise[i][j] = d;
      }
    }
    return {
      names: sets.map((s) => s.name), clades: out, strict,
      shared: out.filter((c) => c.present.every(Boolean)).length,
      pairwise, consensus: sets.map((s, i) => buildFromClades(n, [...strict[i]], null)),
    };
  }

  // Bits -> state symbols for display
  function setToStates(v, data) {
    if (v === MISSING) return '?';
    const out = [];
    if (data.gapMode === 'inapplicable' && (v & NA)) out.push('-');
    for (let i = 0; i < data.symbols.length; i++) if (v & (1 << (i + data.push))) out.push(data.symbols[i]);
    if (data.gapMode === 'extra' && (v & (1 << (data.symbols.length + data.push)))) out.push('-');
    return out;
  }

  /* ------------------------------------------------------------------ */
  /* Unsupported branches                                                */
  /* ------------------------------------------------------------------ */
  // A branch is "unsupported" when every character can be optimised with no
  // change along it: the final state sets at its two ends share a state in
  // every column (TNT's collapsing rule 1, minimum length zero). For
  // inapplicable-type columns only applicable states are compared.
  // Returns { node, keys }: the tree with those branches collapsed, as a nested
  // {children|tip} structure, and its non-trivial clade keys.
  function collapseTree(scorer, tree) {
    const d = scorer.data, nc = d.nchar, n = d.ntax, isNA = d.isNA, w = d.baseW;
    const fin = scorer.finalSets(tree);
    const P = tree.parent, root = tree.root;
    const supported = new Uint8Array(2 * n - 1);
    for (let v = n; v < 2 * n - 1; v++) {
      const p = P[v];
      if (p < 0) continue;
      if (p === root) { supported[v] = 1; continue; } // this is the branch to the root tip
      const o = v * nc, po = p * nc;
      for (let j = 0; j < nc; j++) {
        if (!(w[j] > 0)) continue;
        let a = fin[o + j], b = fin[po + j];
        if (isNA[j]) { a &= ISAPPLIC; b &= ISAPPLIC; if (!a || !b) continue; }
        if ((a & b) === 0) { supported[v] = 1; break; }
      }
    }
    const keys = [], words = Math.ceil(n / 32);
    function rec(v) {
      if (v < n) { const b = new Uint32Array(words); b[v >> 5] |= 1 << (v & 31); return { node: { tip: v }, bits: b, size: 1 }; }
      const kids = [], bits = new Uint32Array(words);
      let size = 0;
      for (const c of [tree.left[v], tree.right[v]]) {
        const r = rec(c);
        for (let q = 0; q < words; q++) bits[q] |= r.bits[q];
        size += r.size;
        if (c >= n && !supported[c]) kids.push(...r.node.children); else kids.push(r.node);
      }
      const node = { children: kids, node: v };
      if (v !== root && supported[v] && size > 1 && size < n - 1 && !(bits[0] & 1)) { node.key = Array.from(bits).join('.'); keys.push(node.key); }
      return { node, bits, size };
    }
    return { node: rec(root).node, keys };
  }
  // Distinct trees after collapsing: [{ node, keys, from: [indices of the binary trees] }]
  function collapseAll(data, trees) {
    const sc = new Scorer(data), seen = new Map(), list = [];
    trees.forEach((t, i) => {
      const c = collapseTree(sc, t), id = c.keys.slice().sort().join('|');
      if (seen.has(id)) list[seen.get(id)].from.push(i);
      else { seen.set(id, list.length); list.push({ node: c.node, keys: c.keys, from: [i] }); }
    });
    return list;
  }
  // Consensus from sets of clade keys (e.g. collapsed trees). level 1 = strict, 0.5 = majority rule.
  function consensusOfKeys(keySets, n, level) {
    const counts = new Map();
    for (const ks of keySets) for (const c of new Set(ks)) counts.set(c, (counts.get(c) || 0) + 1);
    const keep = [], freq = new Map();
    for (const [k, v] of counts) {
      const f = v / keySets.length;
      if ((level >= 1 && v === keySets.length) || (level < 1 && f > level)) { keep.push(k); freq.set(k, f); }
    }
    return buildFromClades(n, keep, freq);
  }

  /* ------------------------------------------------------------------ */
  /* Character changes on branches (synapomorphies)                      */
  /* ------------------------------------------------------------------ */
  // Key of every node of a tree rooted on tip 0: clade key for internal nodes
  // (null for the root and for the all-but-root-tip node), 't<i>' for tips.
  function nodeKeys(tree) {
    const n = tree.n, words = Math.ceil(n / 32), bits = new Array(2 * n - 1), keys = new Array(2 * n - 1).fill(null);
    for (let i = 0; i < n; i++) { const b = new Uint32Array(words); b[i >> 5] |= 1 << (i & 31); bits[i] = b; keys[i] = 't' + i; }
    for (const v of postorder(tree, [])) {
      const b = new Uint32Array(words), l = bits[tree.left[v]], r = bits[tree.right[v]];
      let cnt = 0;
      for (let q = 0; q < words; q++) { b[q] = l[q] | r[q]; cnt += popcnt(b[q]); }
      bits[v] = b;
      if (cnt > 1 && cnt < n - 1 && !(b[0] & 1)) keys[v] = Array.from(b).join('.');
    }
    return keys;
  }

  // Changes on each branch of a binary tree: Array(2n-1), entry v = changes on the
  // branch from parent(v) to v, each { c (character index), from: [symbols], to: [symbols] }.
  //  mode 'unambiguous': changes required by every most parsimonious reconstruction
  //    (the final state sets at the two ends share no state).
  //  mode 'acctran' / 'deltran': one most parsimonious reconstruction (MPR) per
  //    character: of all MPRs, the one whose changes are nearest the root (ACCTRAN)
  //    or nearest the tips (DELTRAN), measured as the sum over changes of the depth
  //    of the branch; ties go to the lower state. Used for standard characters;
  //    characters scored with the inapplicable algorithm are listed as in 'unambiguous'.
  function branchChanges(data, tree, mode, scorer) {
    const sc = scorer || new Scorer(data);
    const n = data.ntax, m = 2 * n - 1, nc = data.nchar;
    const fin = sc.finalSets(tree);
    const P = tree.parent, L = tree.left, R = tree.right, root = tree.root;
    const pre = [];
    { const st = [root]; while (st.length) { const v = st.pop(); pre.push(v); if (v >= n) st.push(L[v], R[v]); } }
    const one = mode === 'acctran' || mode === 'deltran';
    const out = new Array(m);
    for (let v = 0; v < m; v++) out[v] = [];
    let depth = null;
    if (one) { depth = new Int32Array(m); for (let q = 1; q < pre.length; q++) depth[pre[q]] = depth[P[pre[q]]] + 1; }
    const sign = mode === 'deltran' ? -1 : 1;
    for (let c = 0; c < data.nOrig; c++) {
      const cols = data.charCols[c];
      if (!cols.length) continue;
      if (one && !data.charNA[c]) {
        const st = oneMPR(data, tree, c, pre, depth, sign);
        if (!st) continue;
        for (let q = 1; q < pre.length; q++) {
          const v = pre[q];
          if (st.state[v] !== st.state[P[v]]) out[v].push({ c, from: [st.sym[st.state[P[v]]]], to: [st.sym[st.state[v]]] });
        }
        continue;
      }
      for (let q = 1; q < pre.length; q++) {
        const v = pre[q], p = P[v];
        const A = charStatesAt(data, fin, p, c), B = charStatesAt(data, fin, v, c);
        if (A === '?' || B === '?') continue;
        const pa = A.filter((x) => x !== '-'), ca = B.filter((x) => x !== '-');
        if (pa.length && ca.length && !pa.some((x) => ca.includes(x))) out[v].push({ c, from: pa, to: ca });
      }
    }
    return out;
  }
  // One MPR of a standard character by Sankoff dynamic programming: fewest steps,
  // then least (sign = 1) or greatest (sign = -1) total depth of the changes.
  // Returns { state: Int32Array(2n-1) of state numbers, sym: symbol of each state }.
  function oneMPR(data, tree, c, pre, depth, sign) {
    const n = data.ntax, m = 2 * n - 1, nc = data.nchar, cols = data.charCols[c], ord = data.charOrdered[c];
    const L = tree.left, R = tree.right, P = tree.parent;
    const LOW = 1 << data.push, HIGH = 1 << (data.push + 1);
    // states and, for each tip, which of them it may have
    let S, sym, allowed;
    if (!ord) {
      let all = 0;
      for (let t = 0; t < n; t++) { const x = data.tips[t * nc + cols[0]]; if (x !== MISSING) all |= x; }
      const bits = []; for (let b = data.push; b < 32; b++) if (all & (1 << b)) bits.push(b);
      S = bits.length; if (!S) return null;
      sym = bits.map((b) => (b - data.push < data.symbols.length ? data.symbols[b - data.push] : '-'));
      allowed = (t, s) => (data.tips[t * nc + cols[0]] & (1 << bits[s])) !== 0;
    } else {
      S = cols.length + 1;
      sym = Array.from({ length: S }, (_, s) => data.symbols[data.charLo[c] + s]);
      allowed = (t, s) => { for (let q = 0; q < cols.length; q++) if (!(data.tips[t * nc + cols[q]] & (s >= q + 1 ? HIGH : LOW))) return false; return true; };
    }
    const d = ord ? (a, b) => Math.abs(a - b) : (a, b) => (a === b ? 0 : 1);
    const st = new Float64Array(m * S), dp = new Float64Array(m * S); // steps, depth-weighted changes
    for (let q = pre.length - 1; q >= 0; q--) {
      const v = pre[q], o = v * S;
      if (v < n) { for (let s = 0; s < S; s++) { st[o + s] = allowed(v, s) ? 0 : Infinity; dp[o + s] = 0; } continue; }
      for (let s = 0; s < S; s++) {
        let ts = 0, td = 0;
        for (const ch of [L[v], R[v]]) {
          const co = ch * S; let bs = Infinity, bd = Infinity;
          for (let x = 0; x < S; x++) { const a = st[co + x] + d(s, x), e = dp[co + x] + sign * d(s, x) * depth[ch]; if (a < bs || (a === bs && e < bd)) { bs = a; bd = e; } }
          ts += bs; td += bd;
        }
        st[o + s] = ts; dp[o + s] = td;
      }
    }
    const state = new Int32Array(m);
    const ro = tree.root * S; let best = 0;
    for (let s = 1; s < S; s++) if (st[ro + s] < st[ro + best] || (st[ro + s] === st[ro + best] && dp[ro + s] < dp[ro + best])) best = s;
    state[tree.root] = best;
    for (let q = 1; q < pre.length; q++) {
      const v = pre[q], a = state[P[v]], o = v * S;
      let bx = -1, bs = Infinity, bd = Infinity;
      for (let x = 0; x < S; x++) { const y = st[o + x] + d(a, x), e = dp[o + x] + sign * d(a, x) * depth[v]; if (y < bs || (y === bs && e < bd)) { bs = y; bd = e; bx = x; } }
      state[v] = bx;
    }
    return { state, sym };
  }
  // Changes found on the same branch in every tree: Map nodeKey -> [{c, from, to}].
  // Keys as in nodeKeys(); only branches present in all trees are listed.
  function commonChanges(data, trees, mode) {
    const sc = new Scorer(data);
    let acc = null;
    for (const t of trees) {
      const keys = nodeKeys(t), bc = branchChanges(data, t, mode, sc), cur = new Map();
      for (let v = 0; v < keys.length; v++) {
        if (!keys[v]) continue;
        const mset = new Map();
        for (const ch of bc[v]) mset.set(ch.c + '|' + ch.from.join('') + '>' + ch.to.join(''), ch);
        cur.set(keys[v], mset);
      }
      if (!acc) acc = cur;
      else for (const [k, mset] of acc) {
        const other = cur.get(k);
        if (!other) { acc.delete(k); continue; }
        for (const id of [...mset.keys()]) if (!other.has(id)) mset.delete(id);
      }
    }
    const out = new Map();
    if (acc) for (const [k, mset] of acc) out.set(k, [...mset.values()].sort((a, b) => a.c - b.c));
    return out;
  }

  /* ------------------------------------------------------------------ */
  /* Trees from files                                                    */
  /* ------------------------------------------------------------------ */
  // Reads trees from a NEXUS trees block (with or without TRANSLATE), plain
  // Newick, or TNT parenthetical notation (tread; names or 0-based numbers).
  // names: taxa in matrix order. Returns { trees: [{name, nested, binary}], errors: [] }
  // where nested holds tip indices. A tree must contain every taxon exactly once.
  function parseTrees(text, names) {
    const out = { trees: [], errors: [] };
    const n = names.length;
    const loose = (s) => s.replace(/^'+|'+$/g, '').replace(/[\s_]+/g, '_').toLowerCase();
    const byLoose = new Map(names.map((x, i) => [loose(x), i]));
    const strip = (s) => { let o = '', depth = 0, q = false; for (const ch of s) { if (ch === "'" && depth === 0) q = !q; if (!q && ch === '[') { depth++; continue; } if (!q && ch === ']' && depth) { depth--; continue; } if (!depth) o += ch; } return o; };
    const add = (name, nested, src) => {
      const seen = new Uint8Array(n);
      let dup = null, count = 0, binary = true;
      (function walk(x, top) {
        if (typeof x === 'number') { if (seen[x]) dup = names[x]; seen[x] = 1; count++; return; }
        if (x.length === 1) return walk(x[0], top);
        if (x.length > (top ? 3 : 2)) binary = false;
        for (const k of x) walk(k, false);
      })(nested, true);
      if (dup) { out.errors.push(`${name}: taxon "${dup}" appears more than once`); return; }
      if (count !== n) { const miss = names.filter((_, i) => !seen[i]); out.errors.push(`${name}: ${n - count} taxa of the matrix are not in the tree (${miss.slice(0, 3).join(', ')}${miss.length > 3 ? ', …' : ''})`); return; }
      out.trees.push({ name, nested, binary, source: src });
    };
    const tryNewick = (name, str, resolve) => {
      try { add(name, parseNewick(str, names, resolve), 'newick'); } catch (e) { out.errors.push(name + ': ' + e.message.replace(/^Newick: /, '')); }
    };
    if (/\btread\b/i.test(text) && !/begin\s+trees/i.test(text)) {
      let body = text.replace(/^[\s\S]*?\btread\b/i, '').replace(/^\s*'[^']*'/, '');
      const end = body.indexOf(';'); if (end >= 0) body = body.slice(0, end);
      const numeric = !/[A-Za-z]/.test(body.replace(/[()\s*]/g, ' ').trim().split(/\s+/).filter((t) => !byLoose.has(loose(t))).join(''));
      body.split('*').map((s) => s.trim()).filter(Boolean).forEach((s, k) => {
        let i = 0;
        const name = 'Tree ' + (k + 1);
        try {
          const node = () => {
            while (i < s.length && /\s/.test(s[i])) i++;
            if (s[i] === '(') {
              i++;
              const kids = [];
              while (true) { while (i < s.length && /\s/.test(s[i])) i++; if (i >= s.length) throw new Error('unbalanced parentheses'); if (s[i] === ')') { i++; break; } kids.push(node()); }
              return kids;
            }
            let tok = '';
            while (i < s.length && !/[\s()]/.test(s[i])) tok += s[i++];
            let t = byLoose.get(loose(tok));
            if (t === undefined && numeric && /^\d+$/.test(tok) && +tok < n) t = +tok;
            if (t === undefined) throw new Error('unknown taxon "' + tok + '"');
            return t;
          };
          add(name, node(), 'tnt');
        } catch (e) { out.errors.push(name + ': ' + e.message); }
      });
      return out;
    }
    const clean = strip(text);
    if (/begin\s+trees\s*;/i.test(clean)) {
      const blk = clean.match(/begin\s+trees\s*;([\s\S]*?)\bend(?:block)?\s*;/i);
      const body = blk ? blk[1] : clean;
      const tr = new Map();
      const tm = body.match(/\btranslate\b([\s\S]*?);/i);
      if (tm) for (const part of tm[1].split(',')) { const mm = part.trim().match(/^(\S+)\s+([\s\S]+)$/); if (mm) { const t = byLoose.get(loose(mm[2].trim())); if (t !== undefined) tr.set(mm[1], t); } }
      const re = /\b[uU]?tree\s+\*?\s*('[^']*'|[^\s=]+)\s*=\s*([^;]+);/gi;
      let mm;
      while ((mm = re.exec(body))) tryNewick(mm[1].replace(/^'|'$/g, ''), mm[2].trim(), tr.size ? (nm) => tr.get(nm) : null);
      if (!out.trees.length && !out.errors.length) out.errors.push('No TREE commands found in the trees block');
      return out;
    }
    clean.split(';').map((s) => s.trim()).filter((s) => s.includes('(')).forEach((s, k) => tryNewick('Tree ' + (k + 1), s.slice(s.indexOf('(')), null));
    if (!out.trees.length && !out.errors.length) out.errors.push('No trees found: expected Newick, a NEXUS trees block or TNT tread notation');
    return out;
  }

  /* ------------------------------------------------------------------ */
  /* Drawing with several outgroup taxa                                  */
  /* ------------------------------------------------------------------ */
  // node: nested {children|tip} structure rooted on tip 0 (as from consensus() or
  // nestedFromTree()). outTips: tip indices of the outgroup (includes 0). If the
  // remaining taxa form a clade, returns the tree re-rooted between the two
  // groups ({node, mono: true}); otherwise {node: unchanged, mono: false}.
  function rerootOnOutgroup(node, outTips, n) {
    const out = new Uint8Array(n); let k = 0;
    for (const t of outTips) if (!out[t]) { out[t] = 1; k++; }
    if (k <= 1 || k >= n - 1) return { node, mono: true };
    const par = new Map();
    let target = null;
    (function walk(x) {
      if (typeof x.tip === 'number') { x._in = out[x.tip] ? 0 : 1; x._sz = 1; return; }
      x._in = 0; x._sz = 0;
      for (const c of x.children) { par.set(c, x); walk(c); x._in += c._in; x._sz += c._sz; }
      if (!target && x._in === n - k && x._sz === n - k) target = x;
    })(node);
    if (!target || !par.has(target)) return { node, mono: false };
    // rebuild the rest of the tree as seen from the parent of the ingroup node
    function up(x, from) {
      const kids = x.children.filter((c) => c !== from);
      const p = par.get(x);
      if (p) kids.push(up(p, x));
      if (kids.length === 1) return kids[0];
      return { children: kids, flipped: true };
    }
    const rest = up(par.get(target), target);
    return { node: { children: [rest, target] }, mono: true };
  }

  const api = {
    NA, MISSING, ISAPPLIC, RNG, encode, setImplied, cost, expandWeights, mixSeed, Tree, postorder, fromNested, toNested, rootOn, randomTree,
    toNewick, nestedToNewick, parseNewick, parseTrees, Scorer, NAInc, Searcher, charBounds, charStats, charStatesAt, clades, treeKey, consensus, consensusOfKeys,
    nestedFromTree, analyse, searchOpts, repsNeeded, replicate, sectorial, finishTrees, mergeAndFinish, fuse, Pool, resample, resampleRep, bremer,
    compareTreatments, compareAnalyses: compareTreatments, setToStates, prune, regraft, rerootSubtree, collectNodes, cladeTips,
    collapseTree, collapseAll, nodeKeys, branchChanges, commonChanges, rerootOnOutgroup,
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.OCEngine = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);
