const E = require('../src/engine.js');
function toAdj(t) {
  const n = t.n, m = 2 * n - 1;
  const adj = Array.from({ length: m }, () => new Set());
  for (let v = 0; v < m; v++) { const p = t.parent[v]; if (p >= 0) { adj[v].add(p); adj[p].add(v); } }
  const r = t.root; // suppress degree-2 root
  const [a, b] = [...adj[r]];
  adj[a].delete(r); adj[b].delete(r); adj[a].add(b); adj[b].add(a); adj[r].clear();
  return adj;
}
function keyFromAdj(adj, n) {
  // root on tip 0 and compute clade key via a temporary Tree
  const nb = [...adj[0]][0];
  function rec(v, from) { if (v < n) return v; return [...adj[v]].filter((x) => x !== from).map((x) => rec(x, v)); }
  return E.treeKey(E.fromNested(n, [0, rec(nb, 0)]));
}
function edges(adj) { const out = []; adj.forEach((s, v) => s.forEach((u) => { if (v < u) out.push([v, u]); })); return out; }
function component(adj, start, blocked) {
  const seen = new Set([start]); const st = [start];
  while (st.length) { const v = st.pop(); for (const u of adj[v]) if (!(v === blocked[0] && u === blocked[1]) && !(v === blocked[1] && u === blocked[0]) && !seen.has(u)) { seen.add(u); st.push(u); } }
  return seen;
}
function bruteTBR(t) {
  const n = t.n;
  const adj0 = toAdj(t);
  const startKey = keyFromAdj(adj0, n);
  const out = new Set();
  for (const [a, b] of edges(adj0)) {
    const A = component(adj0, a, [a, b]);
    // remove edge a-b and suppress degree-2 nodes a and/or b (if internal)
    const adj = adj0.map((s) => new Set(s));
    adj[a].delete(b); adj[b].delete(a);
    const freeNodes = [];
    for (const x of [a, b]) {
      if (x >= n) { const [p, q] = [...adj[x]]; adj[p].delete(x); adj[q].delete(x); adj[p].add(q); adj[q].add(p); adj[x].clear(); freeNodes.push(x); }
    }
    const eA = edges(adj).filter(([u, v]) => A.has(u) && A.has(v));
    const eB = edges(adj).filter(([u, v]) => !A.has(u) && !A.has(v) && adj[u].size && adj[v].size);
    const optsA = eA.length ? eA : [[a, null]]; // single tip side: attach the tip directly
    const optsB = eB.length ? eB : [[b, null]];
    for (const ea of optsA) for (const eb of optsB) {
      const g = adj.map((s) => new Set(s));
      let next = 0;
      const pick = () => freeNodes[next++];
      function subdivide(e) {
        if (e[1] === null) return e[0];
        const w = pick();
        g[e[0]].delete(e[1]); g[e[1]].delete(e[0]); g[w] = new Set([e[0], e[1]]); g[e[0]].add(w); g[e[1]].add(w);
        return w;
      }
      const x = subdivide(ea), y = subdivide(eb);
      g[x].add(y); g[y].add(x);
      const k = keyFromAdj(g, n);
      if (k !== startKey) out.add(k);
    }
  }
  return out;
}

module.exports = { bruteTBR };
