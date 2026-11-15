// Experiment: full four-pass inapplicable scoring of k columns, compiled to WebAssembly.
#define NA 1
#define ISAPPLIC (~1)
typedef int i32;
__attribute__((export_name("na_full")))
i32 na_full(i32 n, i32 k, const i32* post, i32 npost, const i32* L, const i32* R, const i32* P, i32 root,
            i32* D1, i32* U1, i32* D2, i32* ACT, const i32* w) {
  i32 total = 0;
  for (i32 q = 0; q < npost; q++) {
    i32 v = post[q]; i32* o = D1 + v * k; const i32* a = D1 + L[v] * k; const i32* b = D1 + R[v] * k;
    i32 isRoot = v == root;
    for (i32 i = 0; i < k; i++) {
      i32 x = a[i] & b[i];
      i32 both = ((a[i] & ISAPPLIC) != 0) & ((b[i] & ISAPPLIC) != 0);
      if (x == 0) { x = a[i] | b[i]; if (both) x &= ISAPPLIC; }
      else if (x == NA) { if (both) x = a[i] | b[i]; }
      if (isRoot) x = (x & ISAPPLIC) ? (x & ISAPPLIC) : NA;
      o[i] = x;
    }
  }
  for (i32 q = npost - 1; q >= 0; q--) {
    i32 v = post[q]; i32 p = P[v] >= 0 ? P[v] : v;
    i32* o = U1 + v * k; const i32* d = D1 + v * k; const i32* anc = (v == root ? D1 : U1) + p * k; const i32* l = D1 + L[v] * k; const i32* r = D1 + R[v] * k;
    for (i32 i = 0; i < k; i++) {
      i32 npre = d[i], res;
      if (npre & NA) {
        if (npre & ISAPPLIC) res = anc[i] == NA ? NA : (npre & ISAPPLIC);
        else if (anc[i] == NA) res = NA;
        else { i32 u = (l[i] | r[i]) & ISAPPLIC; res = u ? u : NA; }
      } else res = npre;
      o[i] = res;
    }
  }
  for (i32 t = 0; t < n; t++) {
    if (P[t] < 0) continue;
    const i32* tp = D1 + t * k; const i32* as = U1 + P[t] * k; i32* act = ACT + t * k; i32* d2 = D2 + t * k;
    for (i32 i = 0; i < k; i++) {
      i32 x = tp[i] & as[i];
      act[i] = x ? (x & ISAPPLIC) : (tp[i] & ISAPPLIC);
      d2[i] = (x && (as[i] & ISAPPLIC)) ? (tp[i] & ISAPPLIC) : tp[i];
    }
  }
  for (i32 q = 0; q < npost; q++) {
    i32 v = post[q];
    const i32* nif = U1 + v * k; const i32* a = D2 + L[v] * k; const i32* b = D2 + R[v] * k; const i32* la = ACT + L[v] * k; const i32* ra = ACT + R[v] * k;
    i32* o = D2 + v * k; i32* act = ACT + v * k;
    for (i32 i = 0; i < k; i++) {
      i32 npre, st = 0;
      if (nif[i] & ISAPPLIC) {
        i32 tmp = a[i] & b[i];
        if (tmp) npre = (tmp & ISAPPLIC) ? (tmp & ISAPPLIC) : tmp;
        else { npre = (a[i] | b[i]) & ISAPPLIC; if (((a[i] & ISAPPLIC) && (b[i] & ISAPPLIC)) || (la[i] && ra[i])) st = 1; }
      } else { npre = nif[i]; if (la[i] && ra[i]) st = 1; }
      o[i] = npre; act[i] = (la[i] | ra[i]) & ISAPPLIC; total += st * w[i];
    }
  }
  return total;
}
