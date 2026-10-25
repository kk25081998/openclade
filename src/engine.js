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
  // gapMode: 'inapplicable' (default) | 'missing' | 'extra'
  function encode(mx, gapMode) {
    gapMode = gapMode || 'inapplicable';
    const ntax = mx.taxa.length, nchar = mx.nchar;
    const symIndex = new Map();
    mx.symbols.forEach((s, i) => symIndex.set(s, i));
    const extraBit = mx.symbols.length; // used for 'extra' mode
    const push = gapMode === 'missing' ? 0 : 1; // keep bit 0 for NA
    const tips = new Int32Array(ntax * nchar); // [taxon*nchar + char]
    const isNA = new Uint8Array(nchar);
    const nGaps = new Int32Array(nchar);
    for (let c = 0; c < nchar; c++) {
      let g = 0;
      for (let t = 0; t < ntax; t++) if (mx.cells[t][c] === '-') g++;
      nGaps[c] = g;
      isNA[c] = gapMode === 'inapplicable' && g > NACUTOFF ? 1 : 0;
    }
    for (let t = 0; t < ntax; t++) {
      for (let c = 0; c < nchar; c++) {
        const cell = mx.cells[t][c];
        let v;
        if (cell === '?') v = MISSING;
        else if (cell === '-') {
          if (gapMode === 'extra') v = 1 << (extraBit + push);
          else v = isNA[c] ? NA : MISSING;
        } else {
          v = 0;
          for (const s of cell) v |= 1 << (symIndex.get(s) + push);
        }
        tips[t * nchar + c] = v;
      }
    }
    const naChars = [], stdChars = [];
    for (let c = 0; c < nchar; c++) (isNA[c] ? naChars : stdChars).push(c);
    return {
      ntax, nchar, tips, isNA, nGaps, gapMode,
      naChars: Int32Array.from(naChars), stdChars: Int32Array.from(stdChars),
      symbols: mx.symbols.slice(), push,
    };
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
  function parseNewick(str, names) {
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
      let t = idx.get(nm);
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
    this.weights = new Int32Array(nc).fill(1);
    this.post = [];
    // tip rows
    for (let t = 0; t < n; t++) this.down1.set(data.tips.subarray(t * nc, t * nc + nc), t * nc);
  }

  // Full length; if `steps` (Int32Array nchar) supplied, per-character steps are written.
  Scorer.prototype.length = function (tree, weights, steps, naOnly) {
    const d = this.data, nc = d.nchar, n = d.ntax;
    const w = weights || this.weights;
    const D1 = this.down1, U1 = this.up1, D2 = this.down2, ACT = this.act;
    const std = naOnly ? EMPTY : d.stdChars, nac = d.naChars;
    const post = postorder(tree, this.post);
    const L = tree.left, R = tree.right, P = tree.parent;
    let score = 0;
    if (steps) steps.fill(0);
    // First downpass
    for (let k = 0; k < post.length; k++) {
      const v = post[k], l = L[v] * nc, r = R[v] * nc, o = v * nc;
      for (let i = 0; i < std.length; i++) {
        const j = std[i];
        let x = D1[l + j] & D1[r + j];
        if (x === 0) { x = D1[l + j] | D1[r + j]; score += w[j]; if (steps) steps[j] += w[j]; }
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
    if (nac.length === 0) return score;
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
            if (((a & ISAPPLIC) && (b & ISAPPLIC)) || (la && ra)) { score += w[j]; if (steps) steps[j] += w[j]; }
          }
        } else {
          npre = nif;
          if (la && ra) { score += w[j]; if (steps) steps[j] += w[j]; }
        }
        D2[o + j] = npre;
        ACT[o + j] = (la | ra) & ISAPPLIC;
      }
    }
    return score;
  };

  // Ancestral state sets for character mapping. Returns Int32Array(2n-1) of final sets for char j.
  Scorer.prototype.reconstruct = function (tree, j) {
    const d = this.data, nc = d.nchar, n = d.ntax, m = 2 * n - 1;
    this.length(tree); // fills D1/U1/D2/ACT for NA chars and D1 for std chars
    const D1 = this.down1, U1 = this.up1, D2 = this.down2;
    const L = tree.left, R = tree.right, P = tree.parent;
    const post = postorder(tree, []);
    const fin = new Int32Array(m);
    if (!d.isNA[j]) {
      fin[tree.root] = D1[tree.root * nc + j];
      for (let k = post.length - 1; k >= 0; k--) {
        const v = post[k];
        if (v === tree.root) continue;
        const npre = D1[v * nc + j], anc = fin[P[v]];
        const l = D1[L[v] * nc + j], r = D1[R[v] * nc + j];
        let x = anc & npre;
        if (x !== anc) x = (l & r) ? (npre | (anc & (l | r))) : (npre | anc);
        fin[v] = x;
      }
      for (let t = 0; t < n; t++) {
        const tp = D1[t * nc + j], a = fin[P[t]];
        fin[t] = (tp & a) ? (tp & a) : tp;
      }
      return fin;
    }
    // inapplicable character: second uppass
    fin[tree.root] = D2[tree.root * nc + j];
    for (let k = post.length - 1; k >= 0; k--) {
      const v = post[k];
      const npre = D2[v * nc + j];
      const anc = v === tree.root ? fin[v] : fin[P[v]];
      const l = D2[L[v] * nc + j], r = D2[R[v] * nc + j];
      let x;
      if (npre & ISAPPLIC) {
        if (anc & ISAPPLIC) {
          if ((anc & npre) === anc) x = anc & npre;
          else if (l & r) x = npre | (anc & (l | r));
          else if ((l | r) & NA) x = ((l | r) & anc) ? anc : ((l | r | anc) & ISAPPLIC);
          else x = npre | anc;
        } else x = npre;
      } else x = npre;
      fin[v] = x;
    }
    for (let t = 0; t < n; t++) {
      const tp = D1[t * nc + j], a = fin[P[t]];
      fin[t] = (tp & a) ? (tp & a) : tp;
    }
    return fin;
  };

  /* ------------------------------------------------------------------ */
  /* Incremental scorer for inapplicable-type characters                 */
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
    [this.D1, this.U1, this.D2, this.ACT, this.SC] = this.A; // SC = per-cell weighted steps
    this.w = new Int32Array(k);
    this.bL = new Int32Array(m); this.bR = new Int32Array(m); this.bP = new Int32Array(m); this.bRoot = -1;
    // per-node lists of characters whose value changed in this evaluation
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
    const k = this.k, n = this.n, nac = this.data.naChars;
    for (let i = 0; i < k; i++) this.w[i] = weights ? weights[nac[i]] : 1;
    const post = postorder(tree, this.post);
    const L = tree.left, R = tree.right, P = tree.parent, root = tree.root;
    const D1 = this.D1, U1 = this.U1, D2 = this.D2, ACT = this.ACT, SC = this.SC, w = this.w;
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
    for (const v of post) {
      const o = v * k, lo = L[v] * k, ro = R[v] * k;
      for (let i = 0; i < k; i++) {
        d2cell(U1[o + i], D2[lo + i], D2[ro + i], ACT[lo + i], ACT[ro + i]);
        D2[o + i] = D2OUT; ACT[o + i] = ACTOUT; SC[o + i] = STEPOUT ? w[i] : 0; total += SC[o + i];
      }
    }
    this.total = total;
    this.bL.set(L); this.bR.set(R); this.bP.set(P); this.bRoot = root;
    return total;
  };
  // Cell formulas (MorphyLib passes); d2cell writes results to module globals
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
      else { npre = (a | b) & ISAPPLIC; if (((a & ISAPPLIC) && (b & ISAPPLIC)) || (la && ra)) st = 1; }
    } else { npre = nif; if (la && ra) st = 1; }
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
  // Union of changed-character lists (list, count, node triples) into this.cand
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
    const D1 = this.D1, U1 = this.U1, D2 = this.D2, ACT = this.ACT, SC = this.SC, w = this.w;
    const L1 = this.L1, C1 = this.C1, LU = this.LU, CU = this.CU, L2 = this.L2, C2 = this.C2;
    const cand = this.cand, touched = this.touched;
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
        const sc = STEPOUT ? w[i] : 0;
        if (SC[o + i] !== sc) { total += sc - SC[o + i]; this.set(4, o + i, sc); }
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
    const A = this.A, jA = this.jA, jI = this.jI, jO = this.jO;
    for (let q = this.jn - 1; q >= 0; q--) A[jA[q]][jI[q]] = jO[q];
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
        v &= ISAPPLIC;
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
    this.evaluations = 0;
  }

  Searcher.prototype.score = function (t, w) { this.evaluations++; return this.scorer.length(t, w); };

  // Random addition sequence / stepwise addition (greedy)
  Searcher.prototype.wagner = function (w) {
    const n = this.n, rng = this.rng;
    const order = [0].concat(rng.shuffle([...Array(n - 1).keys()].map((x) => x + 1)));
    // Build initial 3-taxon tree using only the first three taxa; other tips absent
    const t = new Tree(n);
    const a = order[0], b = order[1], c = order[2];
    let next = n;
    const r = next++, u = next++;
    t.root = r; t.left[r] = a; t.right[r] = u; t.parent[a] = r; t.parent[u] = r;
    t.left[u] = b; t.right[u] = c; t.parent[b] = u; t.parent[c] = u;
    // For scoring partial trees we need a scorer that ignores absent tips;
    // since postorder only visits attached nodes this works naturally.
    const present = [a, b, c, r, u];
    for (let i = 3; i < n; i++) {
      const tip = order[i];
      const w2 = next++;
      let best = Infinity, bestE = [];
      const edges = present.filter((v) => v !== t.root && v !== a);
      for (const e of edges) {
        insertTip(t, tip, e, w2);
        const s = this.score(t, w);
        if (s < best) { best = s; bestE = [e]; } else if (s === best) bestE.push(e);
        // undo
        const g = t.parent[w2];
        if (t.left[g] === w2) t.left[g] = e; else t.right[g] = e;
        t.parent[e] = g; t.parent[w2] = -1; t.parent[tip] = -1; t.left[w2] = -1; t.right[w2] = -1;
      }
      insertTip(t, tip, bestE[rng.int(bestE.length)], w2);
      present.push(tip, w2);
    }
    return t;
  };

  // One pass of rearrangements (SPR or TBR); returns improved tree or null.
  // Uses first-improvement hill climbing with randomised order.
  // Directional Fitch sets for the standard (non-inapplicable) characters on
  // the tree hanging from `top` (top has two children). D[v] = set of the
  // subtree below v; U[v] = set of everything on the other side of v's
  // parent edge. Rows are indexed by node id (ns ints each). Returns the
  // weighted Fitch length of that tree, which is independent of its rooting.
  Searcher.prototype.dirSets = function (t, top, w, D, U) {
    const d = this.data, std = d.stdChars, ns = std.length, nc = d.nchar, n = this.n, tips = d.tips;
    const pre = [];
    const st = [top];
    while (st.length) { const v = st.pop(); pre.push(v); if (v >= n) st.push(t.left[v], t.right[v]); }
    let L = 0;
    for (let k = pre.length - 1; k >= 0; k--) {
      const v = pre[k], o = v * ns;
      if (v < n) { const to = v * nc; for (let i = 0; i < ns; i++) D[o + i] = tips[to + std[i]]; continue; }
      const l = t.left[v] * ns, r = t.right[v] * ns;
      for (let i = 0; i < ns; i++) {
        let x = D[l + i] & D[r + i];
        if (x === 0) { x = D[l + i] | D[r + i]; L += w[std[i]]; }
        D[o + i] = x;
      }
    }
    for (let k = 1; k < pre.length; k++) {
      const v = pre[k], pp = t.parent[v];
      const sib = t.left[pp] === v ? t.right[pp] : t.left[pp];
      const o = v * ns, so = sib * ns;
      if (pp === top) { for (let i = 0; i < ns; i++) U[o + i] = D[so + i]; }
      else {
        const po = pp * ns;
        for (let i = 0; i < ns; i++) { let x = D[so + i] & U[po + i]; if (x === 0) x = D[so + i] | U[po + i]; U[o + i] = x; }
      }
    }
    return L;
  };

  // One sweep of SPR or TBR rearrangements with first-improvement hill
  // climbing. Standard characters are scored exactly by indirect calculation
  // (length(S) + length(main) + one step per character whose Fitch sets at
  // the two joined edges are disjoint); inapplicable-type characters are
  // scored by a full pass on the rearranged tree.
  Searcher.prototype.rearrange = function (t, curScore, w, mode, deadline, keep) {
    const n = this.n, m = 2 * n - 1, rng = this.rng, d = this.data;
    const W = w || this.scorer.weights;
    const std = d.stdChars, ns = std.length, hasNA = d.naChars.length > 0;
    if (!this.Dm || this.Dm.length < m * ns) {
      this.Dm = new Int32Array(m * ns); this.Um = new Int32Array(m * ns); this.X = new Int32Array(m * ns);
    }
    const Dm = this.Dm, Um = this.Um, X = this.X;
    const verify = this.verify;
    let NAS = null;
    if (hasNA && !this.noIncremental) {
      if (!this.na) this.na = new NAInc(d);
      NAS = this.na;
      NAS.full(t, w);
    }
    const nodes = [];
    for (let v = 0; v < m; v++) if (v !== t.root && t.parent[v] !== t.root) nodes.push(v);
    rng.shuffle(nodes);
    for (const s of nodes) {
      if (deadline && Date.now() > deadline) return null;
      if (t.parent[s] < 0 || s === t.root || t.parent[s] === t.root) continue;
      const { sib, p } = prune(t, s);
      const main = collectNodes(t, t.root, []).filter((e) => e !== t.root && e !== 0);
      const sub = mode === 'tbr' && s >= n ? collectNodes(t, s, []) : [s];
      const c1 = t.left[s], c2 = t.right[s];
      const origPar = sub.map((v) => t.parent[v]);
      let base = 0;
      if (ns) {
        base = this.dirSets(t, t.root, W, Dm, Um) + this.dirSets(t, s, W, Dm, Um);
        for (const e of main) {
          const o = e * ns;
          for (let i = 0; i < ns; i++) { let x = Dm[o + i] & Um[o + i]; if (x === 0) x = Dm[o + i] | Um[o + i]; X[o + i] = x; }
        }
        for (let q = 0; q < sub.length; q++) {
          const rs = sub[q], o = rs * ns;
          if (rs === s) { for (let i = 0; i < ns; i++) X[o + i] = Dm[o + i]; continue; }
          if (origPar[q] === s) continue;
          for (let i = 0; i < ns; i++) { let x = Dm[o + i] & Um[o + i]; if (x === 0) x = Dm[o + i] | Um[o + i]; X[o + i] = x; }
        }
      }
      let curQ = 0; // index in `sub` of the current physical rooting (sub[0] === s)
      const ensure = (q) => {
        if (q === curQ) return;
        if (sub[q] === s) rerootSubtree(t, s, c1, c2); else rerootSubtree(t, s, sub[q], origPar[q]);
        t.parent[s] = p; curQ = q;
      };
      t.parent[s] = p;
      const thresh = keep ? curScore : curScore - 1; // candidates above this are useless
      for (let q = 0; q < sub.length; q++) {
        const rs = sub[q];
        if (rs !== s && origPar[q] === s) continue;
        rng.shuffle(main);
        const a = rs * ns;
        for (const e of main) {
          if (e === sib && rs === s) continue;
          this.evaluations++;
          let sc = base;
          if (ns) { const b = e * ns; for (let i = 0; i < ns; i++) if ((X[a + i] & X[b + i]) === 0) sc += W[std[i]]; }
          if (!hasNA && !verify && sc > thresh) continue;
          if (hasNA && !verify && sc > thresh) continue; // NA part >= 0
          ensure(q);
          regraft(t, p, e);
          if (hasNA) sc += NAS ? NAS.evalDelta(t) : this.scorer.length(t, W, null, true);
          if (verify) {
            const full = this.scorer.length(t, W);
            if (full !== sc) throw new Error('indirect score ' + sc + ' != full score ' + full);
          }
          if (sc < curScore) { if (NAS) NAS.commit(t); return { tree: t, score: sc }; }
          if (NAS) NAS.revert();
          if (keep && sc === curScore) keep(t, sc);
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
    const r = this.outgroupMoves(t, curScore, w, deadline, keep);
    if (r) return r;
    return null;
  };

  Searcher.prototype.outgroupMoves = function (t, curScore, w, deadline, keep) {
    const C = t.left[t.root] === 0 ? t.right[t.root] : t.left[t.root];
    if (C < this.n) return null;
    const c1 = t.left[C], c2 = t.right[C];
    const inner = collectNodes(t, C, []).filter((v) => v !== C && t.parent[v] !== C);
    const par = inner.map((v) => t.parent[v]);
    const order = this.rng.shuffle([...inner.keys()]);
    let moved = false;
    for (const q of order) {
      if (deadline && Date.now() > deadline) break;
      rerootSubtree(t, C, inner[q], par[q]);
      t.parent[C] = t.root;
      moved = true;
      const sc = this.score(t, w);
      if (sc < curScore) return { tree: t, score: sc };
      if (keep && sc === curScore) keep(t, sc);
    }
    if (moved) { rerootSubtree(t, C, c1, c2); t.parent[C] = t.root; }
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
  // opts: {replicates, ratchetIter, mode:'tbr'|'spr', maxTrees, seed, timeLimit (ms), onProgress}
  // Search driver. Phases: (1) random-addition replicates, each climbed with
  // SPR then TBR and improved by the parsimony ratchet; (2) tree fusing of the
  // replicate results; (3) TBR swapping on every optimal tree to collect all
  // equally parsimonious trees. With opts.searchOnly, phases 2-3 are skipped
  // (used when several workers search in parallel and results are merged).
  function analyse(data, opts) {
    opts = Object.assign({ replicates: 5, ratchetIter: 10, mode: 'tbr', maxTrees: 100, seed: 1, timeLimit: 60000, fuse: true }, opts || {});
    const S = new Searcher(data, { seed: opts.seed, maxTrees: opts.maxTrees });
    const start = Date.now();
    const deadline = start + opts.timeLimit;
    const nc = data.nchar;
    const P = new Pool(opts.maxTrees);
    const repBests = [];
    const progress = (o) => { if (opts.onProgress) opts.onProgress(Object.assign(o, { snapshot: () => ({ best: P.best, trees: P.list(), repBests }) })); };
    const keeper = (t, sc) => P.offer(t, sc);
    // report intermediate improvements (throttled) so an early stop has trees
    let lastImp = 0;
    S.onImprove = (t, sc) => {
      const now = Date.now();
      if (sc < P.best && now - lastImp > 500) {
        lastImp = now; P.offer(t, sc);
        progress({ phase: 'climb', best: P.best, trees: P.size(), fraction: done / totalSteps, evaluations: S.evaluations, elapsed: now - start });
      }
    };
    let done = 0;
    const totalSteps = opts.replicates * (1 + opts.ratchetIter);
    // Each replicate gets an equal share of ~85% of the time budget, so every
    // replicate runs (diversity for fusing) and time remains for phases 2-3.
    const share = (0.85 * opts.timeLimit) / opts.replicates;
    for (let rep = 0; rep < opts.replicates; rep++) {
      if (Date.now() > deadline) break;
      const repDeadline = Math.min(deadline, start + share * (rep + 1));
      let t = S.wagner();
      let r = S.climb(t, null, 'spr', deadline);
      if (opts.mode === 'tbr') r = S.climb(r.tree, null, 'tbr', Math.max(repDeadline, Math.min(deadline, Date.now() + share / 2)));
      P.offer(r.tree, r.score);
      done++;
      progress({ phase: 'search', replicate: rep + 1, best: P.best, trees: P.size(), fraction: done / totalSteps, evaluations: S.evaluations, elapsed: Date.now() - start });
      let cur = r;
      for (let it = 0; it < opts.ratchetIter; it++) {
        if (Date.now() > repDeadline) break;
        const w = new Int32Array(nc);
        for (let j = 0; j < nc; j++) w[j] = S.rng.next() < 0.25 ? 2 : 1; // upweight ~25% of characters
        const pert = S.climb(cur.tree.clone(), w, 'spr', repDeadline);
        const back = S.climb(pert.tree, null, opts.mode, deadline, keeper);
        if (back.score <= cur.score) cur = back;
        P.offer(back.tree, back.score);
        done++;
        progress({ phase: 'ratchet', replicate: rep + 1, iteration: it + 1, best: P.best, trees: P.size(), fraction: done / totalSteps, evaluations: S.evaluations, elapsed: Date.now() - start });
      }
      repBests.push(rootOn(cur.tree, 0));
    }
    if (!opts.searchOnly) finish(S, P, repBests, opts, deadline, (o) => progress(Object.assign(o, { elapsed: Date.now() - start })));
    return {
      bestScore: P.best, trees: P.list(), repBests, evaluations: S.evaluations,
      elapsed: Date.now() - start, timedOut: Date.now() > deadline,
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

  // Phases 2-3: fuse replicate results, then collect equally short trees.
  function finish(S, P, repBests, opts, deadline, progress) {
    if (opts.fuse !== false && Date.now() < deadline) {
      const donors = repBests.concat(P.list().slice(0, 20));
      if (donors.length > 1) {
        const f = fuse(S, donors, deadline);
        if (f.score <= P.best) {
          const c = S.climb(f.tree, null, opts.mode || 'tbr', deadline, (t, sc) => P.offer(t, sc));
          P.offer(c.tree, c.score);
        }
        if (progress) progress({ phase: 'fuse', best: P.best, trees: P.size(), fraction: 1, evaluations: S.evaluations });
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
      if (progress) progress({ phase: 'collect', best: P.best, trees: P.size(), fraction: 1, evaluations: S.evaluations });
    }
  }

  // Merge results of independent searches (e.g. parallel workers) and finish.
  // results: [{trees: [nested], repBests: [nested]}]
  function mergeAndFinish(data, results, opts) {
    opts = Object.assign({ maxTrees: 100, mode: 'tbr', seed: 1, timeLimit: 60000 }, opts || {});
    const S = new Searcher(data, { seed: opts.seed });
    const n = data.ntax, start = Date.now(), deadline = start + opts.timeLimit;
    const P = new Pool(opts.maxTrees);
    const repBests = [];
    for (const r of results) {
      for (const nt of r.trees) { const t = fromNested(n, nt); P.offer(t, S.score(t)); }
      for (const nt of r.repBests || []) repBests.push(fromNested(n, nt));
    }
    finish(S, P, repBests, opts, deadline, opts.onProgress);
    return { bestScore: P.best, trees: P.list(), evaluations: S.evaluations, elapsed: Date.now() - start, timedOut: Date.now() > deadline };
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
          if (sc < best.score) { best = { tree: cand, score: sc }; improved = true; break; }
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
    return { byKey, size, shape };
  }
  function replaceNested(t, node, repl) {
    function rec(v) { if (v === node) return repl; if (v < t.n) return v; return [rec(t.left[v]), rec(t.right[v])]; }
    return rec(t.root);
  }

  // Swap on tree t, adding any equally-short trees to pool (TBR neighbourhood)
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

  // Resampling support (bootstrap or jackknife). Returns Map cladeKey -> frequency
  function resample(data, opts) {
    opts = Object.assign({ method: 'bootstrap', replicates: 100, seed: 7, mode: 'spr', timeLimit: 60000 }, opts || {});
    const S = new Searcher(data, { seed: opts.seed });
    const nc = data.nchar, rng = S.rng;
    const deadline = Date.now() + opts.timeLimit;
    const counts = new Map();
    let done = 0;
    for (let rep = 0; rep < opts.replicates; rep++) {
      if (Date.now() > deadline) break;
      const w = new Int32Array(nc);
      if (opts.method === 'jackknife') {
        for (let j = 0; j < nc; j++) w[j] = rng.next() < Math.exp(-1) ? 0 : 1;
      } else {
        for (let j = 0; j < nc; j++) w[rng.int(nc)]++;
      }
      const t = S.wagner(w);
      const r = S.climb(t, w, opts.mode, deadline);
      const rt = rootOn(r.tree, 0);
      for (const c of clades(rt)) counts.set(c, (counts.get(c) || 0) + 1);
      done++;
      if (opts.onProgress) opts.onProgress({ phase: opts.method, replicate: done, of: opts.replicates, fraction: done / opts.replicates });
    }
    const freq = new Map();
    for (const [k, v] of counts) freq.set(k, v / Math.max(1, done));
    return { freq, replicates: done };
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

  const api = {
    NA, MISSING, ISAPPLIC, RNG, encode, Tree, postorder, fromNested, toNested, rootOn, randomTree,
    toNewick, nestedToNewick, parseNewick, Scorer, Searcher, charBounds, clades, treeKey, consensus,
    nestedFromTree, analyse, mergeAndFinish, fuse, Pool, resample, setToStates, prune, regraft, rerootSubtree, collectNodes, cladeTips,
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.OCEngine = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);
