/* ec3d-math.js -- the mathematical core of the EllipticCurve3D web app.
 *
 * Exact part (BigInt rationals): parsing an affine equation, recognising Weierstrass and
 * hyperelliptic-type equations, a-invariants, discriminant, j-invariant, genus.
 * Numerical part (doubles): the period lattice by the AGM (a port of Sage's
 * PeriodLattice_ell._compute_periods_real and normalise_periods), and the Weierstrass
 * function by its q-series in the normalised basis, evaluated on a whole grid.
 *
 * Plain script so it can be inlined into a single HTML file; also loadable from Node
 * (module.exports) for the tests in web/test.
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.EC3D = factory();
})(typeof self !== 'undefined' ? self : this, function () {
'use strict';

// ------------------------------------------------------------------ rationals
function babs(a) { return a < 0n ? -a : a; }
function bgcd(a, b) { a = babs(a); b = babs(b); while (b) { const t = a % b; a = b; b = t; } return a; }

class Q {
  constructor(n, d = 1n) {
    n = BigInt(n); d = BigInt(d);
    if (d === 0n) throw new Error('division by zero');
    if (d < 0n) { n = -n; d = -d; }
    const g = bgcd(n, d) || 1n;
    this.n = n / g; this.d = d / g;
  }
  static of(x) {
    if (x instanceof Q) return x;
    if (typeof x === 'bigint') return new Q(x);
    if (typeof x === 'number' && Number.isInteger(x)) return new Q(BigInt(x));
    if (typeof x === 'string') return Q.parse(x);
    throw new Error('not a rational number: ' + x);
  }
  static parse(s) {
    s = String(s).trim();
    let m = s.match(/^([+-]?\d+)\s*\/\s*(\d+)$/);
    if (m) return new Q(BigInt(m[1]), BigInt(m[2]));
    m = s.match(/^([+-]?)(\d*)\.(\d+)$/);
    if (m) { const sign = m[1] === '-' ? -1n : 1n; return new Q(sign * BigInt((m[2] || '0') + m[3]), 10n ** BigInt(m[3].length)); }
    if (/^[+-]?\d+$/.test(s)) return new Q(BigInt(s));
    throw new Error('bad number: ' + s);
  }
  add(o) { return new Q(this.n * o.d + o.n * this.d, this.d * o.d); }
  sub(o) { return new Q(this.n * o.d - o.n * this.d, this.d * o.d); }
  mul(o) { return new Q(this.n * o.n, this.d * o.d); }
  div(o) { return new Q(this.n * o.d, this.d * o.n); }
  neg() { return new Q(-this.n, this.d); }
  pow(k) { let r = Q.ONE, b = this; for (; k > 0; k >>= 1) { if (k & 1) r = r.mul(b); b = b.mul(b); } return r; }
  isZero() { return this.n === 0n; }
  isOne() { return this.n === 1n && this.d === 1n; }
  isInteger() { return this.d === 1n; }
  sign() { return this.n === 0n ? 0 : (this.n < 0n ? -1 : 1); }
  eq(o) { return this.n === o.n && this.d === o.d; }
  cmp(o) { const l = this.n * o.d, r = o.n * this.d; return l < r ? -1 : (l > r ? 1 : 0); }
  toNumber() {
    const an = babs(this.n); if (an === 0n) return 0;
    if (an < 2n ** 53n && this.d < 2n ** 53n) return Number(this.n) / Number(this.d);
    const k = this.d.toString(2).length - an.toString(2).length + 64;          // scale so the quotient has ~64 bits
    const v = k >= 0 ? Number(this.n * (1n << BigInt(k)) / this.d) : Number(this.n / (this.d * (1n << BigInt(-k))));
    return v * Math.pow(2, -k / 2) * Math.pow(2, -k / 2);
  }
  toString() { return this.d === 1n ? this.n.toString() : this.n + '/' + this.d; }
}
Q.ZERO = new Q(0n); Q.ONE = new Q(1n);

// --------------------------------------------------- univariate polynomials over Q
// Represented as arrays of Q, index = degree, trailing zeros trimmed (the zero polynomial is []).
function ptrim(p) { let i = p.length; while (i > 0 && p[i - 1].isZero()) i--; return p.slice(0, i); }
function pdeg(p) { return p.length - 1; }               // -1 for the zero polynomial
function padd(p, q) { const r = []; for (let i = 0; i < Math.max(p.length, q.length); i++) r.push((p[i] || Q.ZERO).add(q[i] || Q.ZERO)); return ptrim(r); }
function psub(p, q) { return padd(p, pneg(q)); }
function pneg(p) { return p.map(c => c.neg()); }
function pscale(p, c) { return c.isZero() ? [] : p.map(a => a.mul(c)); }
function pmul(p, q) {
  if (!p.length || !q.length) return [];
  const r = new Array(p.length + q.length - 1).fill(Q.ZERO);
  for (let i = 0; i < p.length; i++) for (let j = 0; j < q.length; j++) r[i + j] = r[i + j].add(p[i].mul(q[j]));
  return ptrim(r);
}
function pdivmod(p, q) {
  if (!q.length) throw new Error('polynomial division by zero');
  let r = p.slice(); const dq = pdeg(q), lc = q[dq]; const quot = new Array(Math.max(0, p.length - dq)).fill(Q.ZERO);
  while (r.length && pdeg(r) >= dq) {
    const k = pdeg(r) - dq, c = r[pdeg(r)].div(lc);
    quot[k] = c;
    for (let i = 0; i <= dq; i++) r[i + k] = r[i + k].sub(c.mul(q[i]));
    r = ptrim(r);
  }
  return [ptrim(quot), r];
}
function pderiv(p) { const r = []; for (let i = 1; i < p.length; i++) r.push(p[i].mul(new Q(BigInt(i)))); return ptrim(r); }
function pmonic(p) { return p.length ? pscale(p, Q.ONE.div(p[pdeg(p)])) : p; }
function pgcd(p, q) {
  p = ptrim(p); q = ptrim(q);
  while (q.length) { const [, r] = pdivmod(p, q); p = q; q = r; }
  return pmonic(p);
}
function pevalNum(p, x) { let r = 0; for (let i = p.length - 1; i >= 0; i--) r = r * x + p[i].toNumber(); return r; }
function pnums(p) { return p.map(c => c.toNumber()); }
function pconst(c) { return c.isZero() ? [] : [c]; }
// squarefree decomposition g = prod g_i^i  (Yun);  returns {a, b} with g = a^2 * b, b squarefree, both monic up to lc(g)
function squarefreeParts(g) {
  g = ptrim(g);
  if (!g.length) return { a: [], b: [] };
  const lc = g[pdeg(g)];
  let f = pmonic(g), a = [Q.ONE], b = [Q.ONE];
  let i = 1, c = pgcd(f, pderiv(f)), w = pdivmod(f, c)[0];
  while (pdeg(w) > 0) {
    const y = pgcd(w, c), z = pdivmod(w, y)[0];          // z = g_i
    for (let k = 0; k < Math.floor(i / 2); k++) a = pmul(a, z);
    if (i % 2 === 1) b = pmul(b, z);
    w = y; c = pdivmod(c, y)[0]; i++;
  }
  return { a, b: pscale(b, lc) };
}
function pformat(p, v = 'x') {
  if (!p.length) return '0';
  const terms = [];
  for (let i = pdeg(p); i >= 0; i--) {
    const c = p[i]; if (c.isZero()) continue;
    const mag = new Q(babs(c.n), c.d), s = c.sign() < 0 ? '-' : '+';
    let t = (mag.isOne() && i > 0) ? '' : mag.toString();
    if (i > 0) t += (t && !mag.isInteger() ? '·' : '') + v + (i > 1 ? sup(i) : '');
    terms.push([s, t]);
  }
  return terms.map(([s, t], k) => (k === 0 ? (s === '-' ? '-' : '') : ' ' + s + ' ') + t).join('');
}
function sup(n) { return String(n).replace(/[0-9]/g, d => '⁰¹²³⁴⁵⁶⁷⁸⁹'[d]); }
// TeX versions of the formatters, for KaTeX in the web page
function qtex(c) {                       // a rational (Q) or a number
  if (c instanceof Q) return c.d === 1n ? c.n.toString() : (c.n < 0n ? '-' : '') + `\\frac{${babs(c.n)}}{${c.d}}`;
  return Number.isInteger(c) ? String(c) : Number(c.toPrecision(6)).toString();
}
function ptex(p, v = 'x') {
  if (!p.length) return '0';
  const terms = [];
  for (let i = pdeg(p); i >= 0; i--) {
    const c = p[i]; if (c.isZero()) continue;
    const mag = new Q(babs(c.n), c.d), s = c.sign() < 0 ? '-' : '+';
    let t = (mag.isOne() && i > 0) ? '' : qtex(mag);
    if (i > 0) t += v + (i > 1 ? `^{${i}}` : '');
    terms.push([s, t]);
  }
  return terms.map(([s, t], k) => (k === 0 ? (s === '-' ? '-' : '') : ' ' + s + ' ') + t).join('');
}
function btex(p) {
  if (!p.size) return '0';
  const keys = [...p.keys()].map(k => k.split(',').map(Number)).sort((a, b) => (b[1] - a[1]) || (b[0] - a[0]));
  const terms = keys.map(([i, j]) => {
    const c = p.get(bkey(i, j)), mag = new Q(babs(c.n), c.d);
    let t = (mag.isOne() && (i || j)) ? '' : qtex(mag);
    if (i) t += 'x' + (i > 1 ? `^{${i}}` : '');
    if (j) t += 'y' + (j > 1 ? `^{${j}}` : '');
    return [c.sign() < 0 ? '-' : '+', t];
  });
  return terms.map(([s, t], k) => (k === 0 ? (s === '-' ? '-' : '') : ' ' + s + ' ') + t).join('');
}
function formatWeierstrassTeX(a) {      // a: Q[] or numbers -> "y^2 + xy + y = x^3 - x^2 - 10x - 20"
  const term = (c, mono) => {
    const q = c instanceof Q;
    if (q ? c.isZero() : c === 0) return '';
    const neg = q ? c.sign() < 0 : c < 0, mag = q ? new Q(babs(c.n), c.d) : Math.abs(c), one = q ? mag.isOne() : mag === 1;
    return (neg ? ' - ' : ' + ') + (one && mono ? '' : qtex(mag)) + mono;
  };
  const [a1, a2, a3, a4, a6] = a;
  return 'y^2' + term(a1, 'xy') + term(a3, 'y') + ' = x^3' + term(a2, 'x^2') + term(a4, 'x') + term(a6, '');
}

// ------------------------------------------------- bivariate polynomials (Map 'i,j' -> Q)
function bkey(i, j) { return i + ',' + j; }
function bconst(c) { const m = new Map(); if (!c.isZero()) m.set(bkey(0, 0), c); return m; }
function bvar(name) { const m = new Map(); m.set(name === 'x' ? bkey(1, 0) : bkey(0, 1), Q.ONE); return m; }
function badd(p, q) {
  const r = new Map(p);
  for (const [k, c] of q) { const s = (r.get(k) || Q.ZERO).add(c); if (s.isZero()) r.delete(k); else r.set(k, s); }
  return r;
}
function bneg(p) { const r = new Map(); for (const [k, c] of p) r.set(k, c.neg()); return r; }
function bsub(p, q) { return badd(p, bneg(q)); }
function bscale(p, c) { const r = new Map(); if (c.isZero()) return r; for (const [k, v] of p) r.set(k, v.mul(c)); return r; }
function bmul(p, q) {
  const r = new Map();
  for (const [k1, c1] of p) { const [i1, j1] = k1.split(',').map(Number);
    for (const [k2, c2] of q) { const [i2, j2] = k2.split(',').map(Number); const k = bkey(i1 + i2, j1 + j2);
      const s = (r.get(k) || Q.ZERO).add(c1.mul(c2)); if (s.isZero()) r.delete(k); else r.set(k, s); } }
  return r;
}
function bpow(p, e) { let r = bconst(Q.ONE); for (let i = 0; i < e; i++) r = bmul(r, p); return r; }
function bisConst(p) { for (const k of p.keys()) if (k !== '0,0') return false; return true; }
function bconstValue(p) { return p.get('0,0') || Q.ZERO; }
function bdegY(p) { let d = -1; for (const k of p.keys()) d = Math.max(d, Number(k.split(',')[1])); return d; }
function bdegX(p) { let d = -1; for (const k of p.keys()) d = Math.max(d, Number(k.split(',')[0])); return d; }
function btotalDeg(p) { let d = -1; for (const k of p.keys()) { const [i, j] = k.split(',').map(Number); d = Math.max(d, i + j); } return d; }
// coefficient of y^j as a univariate polynomial in x
function bcoeffY(p, j) { const r = []; for (const [k, c] of p) { const [i, jj] = k.split(',').map(Number); if (jj === j) r[i] = c; } for (let i = 0; i < r.length; i++) if (!r[i]) r[i] = Q.ZERO; return ptrim(r); }
function bformat(p) {
  if (!p.size) return '0';
  const keys = [...p.keys()].map(k => k.split(',').map(Number)).sort((a, b) => (b[1] - a[1]) || (b[0] - a[0]));
  const terms = keys.map(([i, j]) => {
    const c = p.get(bkey(i, j)), mag = new Q(babs(c.n), c.d);
    let t = (mag.isOne() && (i || j)) ? '' : mag.toString();
    if (i) t += (t && !mag.isInteger() ? '·' : '') + 'x' + (i > 1 ? sup(i) : '');
    if (j) t += 'y' + (j > 1 ? sup(j) : '');
    return [c.sign() < 0 ? '-' : '+', t];
  });
  return terms.map(([s, t], k) => (k === 0 ? (s === '-' ? '-' : '') : ' ' + s + ' ') + t).join('');
}

// ------------------------------------------------------------------ the parser
// Accepts e.g. "y^2 + y = x^3 - x^2", "y² = x³ − 3x + 2", "2*y**2 = x^3 - x", "(y+1)^2 = x(x-1)(x+1)".
function normalizeInput(s) {
  return s.replace(/\*\*/g, '^').replace(/[−–]/g, '-').replace(/[·×⋅]/g, '*')
          .replace(/[⁰¹²³⁴⁵⁶⁷⁸⁹]+/g, m => '^' + [...m].map(c => '⁰¹²³⁴⁵⁶⁷⁸⁹'.indexOf(c)).join(''))
          .replace(/[XY]/g, c => c.toLowerCase()).replace(/\s+/g, ' ').trim();
}
function tokenize(s) {
  const toks = []; let i = 0;
  while (i < s.length) {
    const c = s[i];
    if (c === ' ') { i++; continue; }
    if (/[0-9.]/.test(c)) { let j = i; while (j < s.length && /[0-9.]/.test(s[j])) j++; toks.push({ t: 'num', v: s.slice(i, j) }); i = j; continue; }
    if (c === 'x' || c === 'y') { toks.push({ t: 'var', v: c }); i++; continue; }
    if ('+-*/^()='.includes(c)) { toks.push({ t: c }); i++; continue; }
    throw new Error(`unexpected character "${c}" at position ${i + 1}`);
  }
  return toks;
}
function parsePolynomial(text) {
  const src = normalizeInput(text);
  const sides = src.split('=');
  if (sides.length > 2) throw new Error('more than one "=" sign');
  const parseSide = (str) => {
    const toks = tokenize(str); let pos = 0;
    const peek = () => toks[pos], next = () => toks[pos++];
    const startsFactor = (tk) => tk && (tk.t === 'num' || tk.t === 'var' || tk.t === '(');
    function expr() {
      let v = term();
      while (peek() && (peek().t === '+' || peek().t === '-')) { const op = next().t; const w = term(); v = op === '+' ? badd(v, w) : bsub(v, w); }
      return v;
    }
    function term() {
      let v = unary();
      for (;;) {
        const tk = peek();
        if (tk && tk.t === '*') { next(); v = bmul(v, unary()); }
        else if (tk && tk.t === '/') { next(); const w = unary(); if (!bisConst(w) || bconstValue(w).isZero()) throw new Error('can only divide by a nonzero constant'); v = bscale(v, Q.ONE.div(bconstValue(w))); }
        else if (startsFactor(tk)) v = bmul(v, power());        // implicit multiplication: 2x, x y, 3(x+1)
        else return v;
      }
    }
    function unary() {
      const tk = peek();
      if (tk && tk.t === '-') { next(); return bneg(unary()); }
      if (tk && tk.t === '+') { next(); return unary(); }
      return power();
    }
    function power() {
      const base = atom();
      if (peek() && peek().t === '^') {
        next();
        let e;
        if (peek() && peek().t === '(') { next(); const inner = expr(); if (!peek() || next().t !== ')') throw new Error('missing ")" in exponent'); e = inner; }
        else { const tk = next(); if (!tk || tk.t !== 'num') throw new Error('exponent must be a nonnegative integer'); e = bconst(Q.parse(tk.v)); }
        if (!bisConst(e) || !bconstValue(e).isInteger() || bconstValue(e).sign() < 0) throw new Error('exponent must be a nonnegative integer');
        return bpow(base, Number(bconstValue(e).n));
      }
      return base;
    }
    function atom() {
      const tk = next();
      if (!tk) throw new Error('unexpected end of input');
      if (tk.t === 'num') return bconst(Q.parse(tk.v));
      if (tk.t === 'var') return bvar(tk.v);
      if (tk.t === '(') { const v = expr(); const cl = next(); if (!cl || cl.t !== ')') throw new Error('missing ")"'); return v; }
      throw new Error('unexpected "' + (tk.v || tk.t) + '"');
    }
    const v = expr();
    if (pos < toks.length) throw new Error('unexpected "' + (toks[pos].v || toks[pos].t) + '"');
    return v;
  };
  const lhs = parseSide(sides[0]);
  const F = sides.length === 2 ? bsub(lhs, parseSide(sides[1])) : lhs;
  if (!F.size) throw new Error('the equation is trivial (0 = 0)');
  return F;
}

// ------------------------------------------------------------- Weierstrass invariants
function invariantsQ(a) {              // a = [a1,a2,a3,a4,a6] as Q
  const [a1, a2, a3, a4, a6] = a;
  const b2 = a1.mul(a1).add(a2.mul(Q.of(4)));
  const b4 = a4.mul(Q.of(2)).add(a1.mul(a3));
  const b6 = a3.mul(a3).add(a6.mul(Q.of(4)));
  const b8 = a1.mul(a1).mul(a6).add(a2.mul(a6).mul(Q.of(4))).sub(a1.mul(a3).mul(a4)).add(a2.mul(a3).mul(a3)).sub(a4.mul(a4));
  const c4 = b2.mul(b2).sub(b4.mul(Q.of(24)));
  const c6 = b2.pow(3).neg().add(b2.mul(b4).mul(Q.of(36))).sub(b6.mul(Q.of(216)));
  const disc = b2.mul(b2).mul(b8).neg().sub(b4.pow(3).mul(Q.of(8))).sub(b6.mul(b6).mul(Q.of(27))).add(b2.mul(b4).mul(b6).mul(Q.of(9)));
  const j = disc.isZero() ? null : c4.pow(3).div(disc);
  return { b2, b4, b6, b8, c4, c6, disc, j };
}
function invariantsNum(a) {            // same with doubles (for models with irrational coefficients)
  const [a1, a2, a3, a4, a6] = a;
  const b2 = a1 * a1 + 4 * a2, b4 = 2 * a4 + a1 * a3, b6 = a3 * a3 + 4 * a6;
  const b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4;
  const c4 = b2 * b2 - 24 * b4, c6 = -(b2 ** 3) + 36 * b2 * b4 - 216 * b6;
  const disc = -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6;
  return { b2, b4, b6, b8, c4, c6, disc, j: disc === 0 ? null : c4 ** 3 / disc };
}
function formatWeierstrass(a) {        // a: Q[] or numbers -> "y² + xy + y = x³ - x² - 10x - 20"
  const q = a.map(v => (v instanceof Q) ? v : null);
  if (q.some(v => v === null)) return null;
  const [a1, a2, a3, a4, a6] = q;
  const lhs = pformatY([Q.ZERO, a3, Q.ONE], a1);
  const rhs = pformat([a6, a4, a2, Q.ONE]);
  return lhs + ' = ' + rhs;
  function pformatY(cy, a1) {          // y^2 + a1 x y + a3 y
    let s = 'y²';
    if (!a1.isZero()) s += (a1.sign() < 0 ? ' - ' : ' + ') + (babs(a1.n) === 1n && a1.d === 1n ? '' : new Q(babs(a1.n), a1.d)) + 'xy';
    if (!cy[1].isZero()) s += (cy[1].sign() < 0 ? ' - ' : ' + ') + (babs(cy[1].n) === 1n && cy[1].d === 1n ? '' : new Q(babs(cy[1].n), cy[1].d)) + 'y';
    return s;
  }
}

// ---------------------------------------------------------------- complex helpers
const cadd = (a, b) => [a[0] + b[0], a[1] + b[1]];
const csub = (a, b) => [a[0] - b[0], a[1] - b[1]];
const cmul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
const cdiv = (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; };
const cabs = (a) => Math.hypot(a[0], a[1]);
const csqrt = (a) => {                                   // principal branch, stable near the negative real axis (Kahan)
  const r = cabs(a); if (r === 0) return [0, 0];
  const t = Math.sqrt((Math.abs(a[0]) + r) / 2);
  if (a[0] >= 0) return [t, a[1] / (2 * t)];
  return [Math.abs(a[1]) / (2 * t), a[1] < 0 ? -t : t];
};
function cpolyEval(p, z) { let r = [0, 0]; for (let i = p.length - 1; i >= 0; i--) r = cadd(cmul(r, z), [p[i], 0]); return r; }

// real AGM
function agm(a, b) { for (let i = 0; i < 60; i++) { const c = (a + b) / 2, d = Math.sqrt(a * b); if (Math.abs(a - b) <= 1e-16 * Math.abs(c)) return c; a = c; b = d; } return (a + b) / 2; }

// roots of a real cubic  c3 x^3 + c2 x^2 + c1 x + c0  (returns three complex numbers).
// One real root from the closed formula polished by Newton, then deflation to a quadratic for the other two.
function cubicRoots(c3, c2, c1, c0) {
  const a = c2 / c3, b = c1 / c3, c = c0 / c3;                    // monic: x^3 + a x^2 + b x + c
  const p = b - a * a / 3, q = 2 * a ** 3 / 27 - a * b / 3 + c;    // depressed: t^3 + p t + q
  const D = q * q / 4 + p ** 3 / 27;
  let x0;
  if (D <= 0) {                                                      // three real roots: take the largest
    const r = Math.sqrt(Math.max(0, -p / 3));
    const phi = Math.acos(Math.max(-1, Math.min(1, r === 0 ? 0 : -q / (2 * r ** 3))));
    x0 = 2 * r * Math.cos(phi / 3) - a / 3;
  } else {
    const s = Math.sqrt(D);
    x0 = Math.cbrt(-q / 2 + s) + Math.cbrt(-q / 2 - s) - a / 3;
  }
  const f = x => ((x + a) * x + b) * x + c, df = x => (3 * x + 2 * a) * x + b;
  for (let i = 0; i < 50; i++) { const d = df(x0); if (!d) break; const step = f(x0) / d; x0 -= step; if (Math.abs(step) <= 1e-17 * Math.max(1, Math.abs(x0))) break; }
  // deflate: x^3 + a x^2 + b x + c = (x - x0)(x^2 + p1 x + p0)
  const p1 = a + x0, p0 = b + x0 * p1;
  const disc = p1 * p1 - 4 * p0;
  let r1, r2;
  if (disc >= 0) { const sq = Math.sqrt(disc); const t = p1 >= 0 ? -(p1 + sq) / 2 : -(p1 - sq) / 2; r1 = [t, 0]; r2 = [t !== 0 ? p0 / t : -p1 - t, 0]; }
  else { const sq = Math.sqrt(-disc); r1 = [-p1 / 2, sq / 2]; r2 = [-p1 / 2, -sq / 2]; }
  // polish all three on the original polynomial
  const F = z => cpolyEval([c0, c1, c2, c3], z), dF = z => cpolyEval([c1, 2 * c2, 3 * c3], z);
  return [[x0, 0], r1, r2].map(z => { for (let i = 0; i < 6; i++) { const d = dF(z); if (cabs(d) === 0) break; const st = cdiv(F(z), d); z = csub(z, st); if (cabs(st) <= 1e-17 * Math.max(1, cabs(z))) break; } return z; });
}

// exact rational from a double (every double is a dyadic rational)
function qFromDouble(x) {
  if (!isFinite(x)) throw new Error('not finite');
  let k = 0; while (!Number.isInteger(x) && k < 1100) { x *= 2; k++; }
  return new Q(BigInt(x), 1n << BigInt(k));
}
// Newton refinement of a root of a polynomial with Q coefficients, with EXACT evaluation, so that
// nearby roots are still found to full precision (the differences e_i - e_j then come out exact).
function qRound(q, bits = 320) {                // keep the rationals small: round to a dyadic with 2^bits denominator
  const d = 1n << BigInt(bits); let n = (q.n * d * 2n) / q.d; n = (n + (n < 0n ? -1n : 1n)) / 2n; return new Q(n, d);
}
function refineRootExact(pQ, z, maxSteps = 8) {
  const dQ = pderiv(pQ);
  let zr = qFromDouble(z[0]), zi = qFromDouble(z[1]);
  const isReal = z[1] === 0;
  for (let it = 0; it < maxSteps; it++) {
    const [fr, fi] = cqEval(pQ, zr, zi), [gr, gi] = cqEval(dQ, zr, zi);
    const den = gr.mul(gr).add(gi.mul(gi));
    if (den.isZero()) break;
    const sr = fr.mul(gr).add(fi.mul(gi)).div(den), si = fi.mul(gr).sub(fr.mul(gi)).div(den);
    zr = qRound(zr.sub(sr)); zi = isReal ? Q.ZERO : qRound(zi.sub(si));
    if (Math.hypot(sr.toNumber(), si.toNumber()) <= 1e-60 * (1 + Math.hypot(zr.toNumber(), zi.toNumber()))) break;
  }
  return [zr, zi];
}
// Roots of a cubic with exact coefficients pQ = [c0, c1, c2, c3] to full double precision even when two roots are
// close: one well-separated root is refined by exact Newton steps, the cubic is deflated exactly, and the other two
// roots come from the quadratic.  Returns {roots: [[re, im] x 3], diff(i, j)}: accurate roots[i] - roots[j].
function cubicRootsExact(pQ) {
  const c = pQ.map(v => v.toNumber());
  const approx = cubicRoots(c[3], c[2], c[1], c[0]);
  let best = 0, bestd = -1;
  for (let i = 0; i < 3; i++) { let d = Infinity; for (let j = 0; j < 3; j++) if (j !== i) d = Math.min(d, cabs(csub(approx[i], approx[j]))); if (d > bestd) { bestd = d; best = i; } }
  let z0 = approx[best]; if (Math.abs(z0[1]) <= 1e-8 * (1 + Math.abs(z0[0]))) z0 = [z0[0], 0];
  const [zr, zi] = refineRootExact(pQ, z0);
  // monic cubic x^3 + A x^2 + B x + C = (x - z)(x^2 + p1 x + p0):  p1 = A + z,  p0 = B + z p1   (complex rationals)
  const A = pQ[2].div(pQ[3]), B = pQ[1].div(pQ[3]);
  const p1r = A.add(zr), p1i = zi;
  const p0r = B.add(zr.mul(p1r).sub(zi.mul(p1i))), p0i = zr.mul(p1i).add(zi.mul(p1r));
  const four = Q.of(4), two = Q.of(2);
  const dr = p1r.mul(p1r).sub(p1i.mul(p1i)).sub(p0r.mul(four)), di = p1r.mul(p1i).mul(two).sub(p0i.mul(four));   // p1^2 - 4 p0, exact
  const sq = csqrt([dr.toNumber(), di.toNumber()]), hs = [sq[0] / 2, sq[1] / 2];
  const mp1 = [p1r.neg().toNumber() / 2, p1i.neg().toNumber() / 2];
  const z = [zr.toNumber(), zi.toNumber()], r1 = cadd(mp1, hs), r2 = csub(mp1, hs);
  const zp = [zr.add(p1r.div(two)).toNumber(), zi.add(p1i.div(two)).toNumber()];        // z + p1/2, exact then rounded
  const D = [[[0, 0], csub(zp, hs), cadd(zp, hs)], [null, [0, 0], sq], [null, null, [0, 0]]];
  const diff = (i, j) => i === j ? [0, 0] : (i < j ? D[i][j] : [-D[j][i][0], -D[j][i][1]]);
  return { roots: [z, r1, r2], diff };
}
function cqEval(pQ, zr, zi) {
  let rr = Q.ZERO, ri = Q.ZERO;
  for (let i = pQ.length - 1; i >= 0; i--) { const nr = rr.mul(zr).sub(ri.mul(zi)).add(pQ[i]); const ni = rr.mul(zi).add(ri.mul(zr)); rr = nr; ri = ni; }
  return [rr, ri];
}

// ---------------------------------------------------------------- the period lattice
// Port of Sage's PeriodLattice_ell for a curve with real a-invariants (numbers).
// Returns {w1 (real), w2 ([re,im]), discSign, e: roots of the 2-division polynomial, normalised: {n1, n2, mat, tau}}.
function periodLattice(a, discSign, aQ) {
  const inv = invariantsNum(a);
  if (discSign === undefined) discSign = Math.sign(inv.disc);
  if (discSign === 0) throw new Error('singular curve: discriminant is zero');
  // 2-division polynomial 4x^3 + b2 x^2 + 2 b4 x + b6
  let roots = cubicRoots(4, inv.b2, 2 * inv.b4, inv.b6);
  let w1, w2, ei, diffs;
  const PI = Math.PI;
  // with exact a-invariants, refine the roots exactly so that close roots keep full precision
  let ex = null;
  if (aQ) { const iq = invariantsQ(aQ); ex = cubicRootsExact([iq.b6, iq.b4.mul(Q.of(2)), iq.b2, Q.of(4)]); roots = ex.roots; }
  if (discSign > 0) {
    const order = [0, 1, 2].sort((i, j) => roots[i][0] - roots[j][0]);          // e1 < e2 < e3
    ei = order.map(i => roots[i][0]);
    if (ex) diffs = [ex.diff(order[2], order[0])[0], ex.diff(order[2], order[1])[0], ex.diff(order[1], order[0])[0]];
    else diffs = [ei[2] - ei[0], ei[2] - ei[1], ei[1] - ei[0]];
    const A = Math.sqrt(diffs[0]), B = Math.sqrt(diffs[1]), C = Math.sqrt(diffs[2]);
    w1 = PI / agm(A, B); w2 = [0, PI / agm(A, C)];
  } else {
    const order = [0, 1, 2].sort((i, j) => roots[i][1] - roots[j][1]);          // by imaginary part: e1 (Im<0), e3 (real), e2 (Im>0)
    const e1 = roots[order[0]], e3 = [roots[order[1]][0], 0], e2 = roots[order[2]];
    ei = [e1, e2, e3];
    const d31 = ex ? ex.diff(order[1], order[0]) : csub(e3, e1);
    const A = csqrt([d31[0], d31[1]]);
    const x = Math.abs(A[0]), y = Math.abs(A[1]), r = cabs(A);
    w1 = PI / agm(r, x);
    const w2i = PI / agm(r, y);
    w2 = [w1 / 2, w2i / 2];
  }
  const { basis, mat } = normalisePeriods([w1, 0], w2);
  const tau = cdiv(basis[0], basis[1]);
  return { w1, w2, discSign, e: ei, normalised: { n1: basis[0], n2: basis[1], mat, tau }, inv };
}
// Sage's reduce_tau / normalise_periods, verbatim in spirit: tau = w1/w2 is moved to the
// fundamental domain, and the SL2(Z) matrix [a,b,c,d] with (n1, n2) = (a w1 + b w2, c w1 + d w2) is returned.
const roundHalfAway = x => Math.sign(x) * Math.round(Math.abs(x));      // Sage's RealNumber.round()
function reduceTau(tau) {
  let a = 1, b = 0, c = 0, d = 1;
  let k = roundHalfAway(tau[0]); tau = [tau[0] - k, tau[1]]; a -= k * c; b -= k * d;
  let guard = 0;
  while (cabs(tau) < 0.999 && guard++ < 200) {
    tau = cdiv([-1, 0], tau);
    [a, b, c, d] = [c, d, -a, -b];
    k = roundHalfAway(tau[0]); tau = [tau[0] - k, tau[1]]; a -= k * c; b -= k * d;
  }
  return { tau, mat: [a, b, c, d] };
}
function normalisePeriods(w1, w2) {
  let tau = cdiv(w1, w2), s = 1;
  if (tau[1] < 0) { w2 = [-w2[0], -w2[1]]; tau = [-tau[0], -tau[1]]; s = -1; }
  const { mat } = reduceTau(tau);
  const [a, b, c, d] = mat;
  // the new basis uses the (possibly negated) w2; the reported matrix refers to the original w2, as in Sage
  const n1 = cadd(cmul([a, 0], w1), cmul([b, 0], w2)), n2 = cadd(cmul([c, 0], w1), cmul([d, 0], w2));
  return { basis: [n1, n2], mat: s < 0 ? [a, -b, c, -d] : [a, b, c, d] };
}

// ---------------------------------------------------------------- Weierstrass wp by q-series
// makeWp(lattice) returns wp(m, n): the pair (wp(z), wp'(z)) for z = m*n1 + n*n2 in the normalised basis,
// or null at a lattice point.  Uses  wp(zeta; tau) = (2 pi i)^2 [1/12 + sum_{k in Z} q^k u/(1-q^k u)^2 - 2 sum_{k>=1} q^k/(1-q^k)^2 ]
// with u = e^{2 pi i zeta}, q = e^{2 pi i tau}, zeta = z/n2, tau = n1/n2 in the fundamental domain (so |q| < 0.0044).
function makeWp(lat, nterms = 10) {
  const { n2, tau } = lat.normalised;
  const tr = tau[0], ti = tau[1];
  const qr0 = Math.exp(-2 * Math.PI * ti), qa = 2 * Math.PI * tr;
  const qn = []; let c0r = 1 / 12, c0i = 0;
  for (let k = 1; k <= nterms; k++) {
    const r = qr0 ** k, ang = qa * k, qkr = r * Math.cos(ang), qki = r * Math.sin(ang);
    qn.push([qkr, qki]);
    // -2 q^k/(1-q^k)^2
    const d = [1 - qkr, -qki], t = cdiv([qkr, qki], cmul(d, d));
    c0r -= 2 * t[0]; c0i -= 2 * t[1];
  }
  const n2sq = cmul(n2, n2), n2cu = cmul(n2sq, n2);
  const K2 = cdiv([-4 * Math.PI * Math.PI, 0], n2sq);          // (2 pi i)^2 / n2^2
  const K3 = cdiv([0, -8 * Math.PI ** 3], n2cu);               // (2 pi i)^3 / n2^3
  const TWO_PI = 2 * Math.PI;
  return function wp(m, n) {
    m -= Math.round(m);                                         // translate by n1 so that |Im zeta| <= Im tau / 2
    const r = Math.exp(-TWO_PI * m * ti), ang = TWO_PI * (m * tr + n);
    const ur = r * Math.cos(ang), ui = r * Math.sin(ang);
    // u/(1-u)^2 and u(1+u)/(1-u)^3
    let dr = 1 - ur, di = -ui;
    if (dr * dr + di * di < 1e-24) return null;                 // the pole
    let Pr, Pi, Dr, Di;
    { const d2 = cmul([dr, di], [dr, di]); const t = cdiv([ur, ui], d2); Pr = c0r + t[0]; Pi = c0i + t[1];
      const s = cdiv(cmul([ur, ui], [1 + ur, ui]), cmul(d2, [dr, di])); Dr = s[0]; Di = s[1]; }
    for (let k = 0; k < nterms; k++) {
      const [qr, qi] = qn[k];
      const v = [qr * ur - qi * ui, qr * ui + qi * ur];          // q^k u
      const w = cdiv([qr, qi], [ur, ui]);                        // q^k / u
      const dv = [1 - v[0], -v[1]], dw = [1 - w[0], -w[1]];
      const dv2 = cmul(dv, dv), dw2 = cmul(dw, dw);
      const tv = cdiv(v, dv2), tw = cdiv(w, dw2);
      Pr += tv[0] + tw[0]; Pi += tv[1] + tw[1];
      const sv = cdiv(cmul(v, [1 + v[0], v[1]]), cmul(dv2, dv)), sw = cdiv(cmul(w, [1 + w[0], w[1]]), cmul(dw2, dw));
      Dr += sv[0] - sw[0]; Di += sv[1] - sw[1];
    }
    return [cmul(K2, [Pr, Pi]), cmul(K3, [Dr, Di])];
  };
}
// (s, t) lattice coordinates w.r.t. the real basis (w1, w2)  ->  (m, n) w.r.t. the normalised basis
function toNormalisedCoords(lat, s, t) {
  const [a, b, c, d] = lat.normalised.mat;               // n1 = a w1 + b w2, n2 = c w1 + d w2, det = 1
  return [s * d - t * c, t * a - s * b];                 // z = s w1 + t w2 = m n1 + n n2
}

// ---------------------------------------------------------------- models and the map back
// A "model" is a Weierstrass equation with real coefficients together with the map from its
// coordinates (X, Y) back to the coordinates (x, y) of the equation the user typed.
//   kind 'direct':  x = X/k, y = Y/k                       (generalised Weierstrass input, exact a-invariants)
//   kind 'square':  x = X/k, y' = Y/k, y = a(x) y' - h(x)/2 (input y^2 + h(x) y = f(x) whose completed square is a(x)^2 b(x), b cubic)
//   kind 'quartic': x = r + 1/t, t = T/k, y' = (Y/k)/t^2 * (...)  (b quartic with a real root r; see makeQuarticModel)
function modelFromAinvs(aQ, extra = {}) {
  const a = aQ.map(v => v.toNumber());
  const inv = invariantsQ(aQ);
  if (inv.disc.isZero()) return { singular: true, ainvs: a, ainvsQ: aQ, inv, message: 'singular Weierstrass cubic (discriminant 0)' };
  return Object.assign({ kind: 'direct', ainvs: a, ainvsQ: aQ, inv, k: 1, aPoly: [1], hPoly: [], exact: true,
    toOriginal: (X, Y) => [X, Y] }, extra);
}
function scaledCubicModel(bQ, aPolyQ, hPolyQ) {
  // model for  y'^2 = b(x)  with b = k x^3 + b2 x^2 + b1 x + b0 (Q), plus the map y = a(x) y' - h(x)/2
  const k = bQ[3], b2 = bQ[2], b1 = bQ[1] || Q.ZERO, b0 = bQ[0] || Q.ZERO;
  const aQ = [Q.ZERO, b2, Q.ZERO, k.mul(b1), k.mul(k).mul(b0)];          // (X, Y) = (k x, k y')
  const kn = k.toNumber(), aP = pnums(aPolyQ), hP = pnums(hPolyQ);
  const inv = invariantsQ(aQ);
  return { kind: 'square', ainvs: aQ.map(v => v.toNumber()), ainvsQ: aQ, inv, k: kn, aPoly: aP, hPoly: hP, exact: true,
    singular: inv.disc.isZero(),
    toOriginal: (X, Y) => {
      const x = [X[0] / kn, X[1] / kn], yp = [Y[0] / kn, Y[1] / kn];
      const ax = cpolyEval(aP, x), hx = cpolyEval(hP, x);
      return [x, csub(cmul(ax, yp), [hx[0] / 2, hx[1] / 2])];
    } };
}
function quarticModel(bQ, aPolyQ, hPolyQ) {
  // y'^2 = b(x), b a squarefree quartic with a real root r:  x = r + 1/t,  y' = Y/t^2  gives  Y^2 = t * (t^3 c(r + 1/t)) =: cubic(t)
  const b = pnums(bQ);
  const roots = polyRootsNum(b).filter(z => Math.abs(z[1]) < 1e-9 * (1 + Math.abs(z[0])));
  if (!roots.length) return null;
  roots.sort((u, v) => Math.abs(u[0]) - Math.abs(v[0]));
  let r = roots[0][0];
  for (let i = 0; i < 3; i++) { const f = pevalNumArr(b, r), df = pevalNumArr(pderivNum(b), r); if (df) r -= f / df; }
  // c(x) = b(x)/(x - r): synthetic division
  const c = [0, 0, 0, 0]; let carry = b[4]; c[3] = carry;
  for (let i = 3; i >= 1; i--) { carry = b[i] + r * carry; c[i - 1] = carry; }
  // cubic(t) = t * sum_i c_i (1 + r t)^i t^{3-i}  -- expand into t^1..t^4? No: t^3 c(r + 1/t) = sum_i c_i (1 + r t)^i t^{3 - i} is a cubic in t; times t gives degree 4?
  // Careful: t^4 b(r + 1/t) = t^4 (1/t) c(r + 1/t) = t^3 c(r + 1/t) = sum_i c_i (rt + 1)^i t^{3-i}, a polynomial of degree 3 in t. So Y^2 = cubic(t).
  const cub = [0, 0, 0, 0];
  for (let i = 0; i <= 3; i++) {
    // c_i (r t + 1)^i t^{3-i}
    const binom = [[1], [1, 1], [1, 2, 1], [1, 3, 3, 1]][i];
    for (let j = 0; j <= i; j++) cub[j + 3 - i] += c[i] * binom[j] * r ** j;
  }
  const k = cub[3];
  if (Math.abs(k) < 1e-14) return null;
  const a = [0, cub[2], 0, k * cub[1], k * k * cub[0]];                    // (T, Yc) = (k t, k Y)
  const inv = invariantsNum(a);
  const aP = pnums(aPolyQ), hP = pnums(hPolyQ);
  return { kind: 'quartic', ainvs: a, ainvsQ: null, inv, k, r, aPoly: aP, hPoly: hP, exact: false, singular: inv.disc === 0,
    toOriginal: (T, Yc) => {
      const t = [T[0] / k, T[1] / k], Y = [Yc[0] / k, Yc[1] / k];
      const x = cadd([r, 0], cdiv([1, 0], t));
      const yp = cdiv(Y, cmul(t, t));
      const ax = cpolyEval(aP, x), hx = cpolyEval(hP, x);
      return [x, csub(cmul(ax, yp), [hx[0] / 2, hx[1] / 2])];
    } };
}
function pevalNumArr(p, x) { let r = 0; for (let i = p.length - 1; i >= 0; i--) r = r * x + p[i]; return r; }
function pderivNum(p) { const r = []; for (let i = 1; i < p.length; i++) r.push(i * p[i]); return r; }
// all complex roots of a real polynomial (Durand-Kerner), for small degrees
function polyRootsNum(p) {
  p = p.slice(); while (p.length && Math.abs(p[p.length - 1]) === 0) p.pop();
  const n = p.length - 1; if (n < 1) return [];
  const lc = p[n], q = p.map(c => c / lc);
  let R = 1 + Math.max(...q.slice(0, n).map(Math.abs));
  let z = []; for (let i = 0; i < n; i++) z.push([R * Math.cos(2 * Math.PI * i / n + 0.4), R * Math.sin(2 * Math.PI * i / n + 0.4)]);
  const ev = w => cpolyEval(q, w);
  for (let it = 0; it < 500; it++) {
    let maxd = 0;
    for (let i = 0; i < n; i++) {
      let den = [1, 0]; for (let j = 0; j < n; j++) if (j !== i) den = cmul(den, csub(z[i], z[j]));
      const d = cdiv(ev(z[i]), den); z[i] = csub(z[i], d); maxd = Math.max(maxd, cabs(d));
    }
    if (maxd < 1e-15) break;
  }
  return z;
}

// ---------------------------------------------------------------- analysing an equation
function analyzeEquation(text) {
  const F = parsePolynomial(text);
  const res = { input: text, F, pretty: bformat(F) + ' = 0', degY: bdegY(F), degX: bdegX(F), totalDeg: btotalDeg(F),
                genus: null, genusBound: null, kind: 'plane', model: null, notes: [] };
  const d = res.totalDeg;
  res.genusBound = Math.max(0, (d - 1) * (d - 2) / 2);
  const dy = res.degY;
  if (dy <= 0) { res.kind = 'degenerate'; res.notes.push('no $y$ in the equation: this is a union of vertical lines, not a curve to plot'); return res; }
  if (dy === 1) {
    res.kind = 'rational'; res.genus = 0;
    res.notes.push('linear in $y$: a rational curve (genus 0), $y = -C(x)/A(x)$');
    return res;
  }
  if (dy === 2) {
    const A = bcoeffY(F, 2), B = bcoeffY(F, 1), C = bcoeffY(F, 0);
    if (pdeg(A) === 0) {
      const inv = Q.ONE.div(A[0]);
      const h = pscale(B, inv), f = pscale(C, inv.neg());                    // y^2 + h(x) y = f(x)
      res.kind = 'hyperelliptic-type'; res.h = h; res.f = f;
      // generalised Weierstrass form?  deg h <= 1, deg f = 3
      if (pdeg(h) <= 1 && pdeg(f) === 3) {
        const k = f[3];
        const a1 = h[1] || Q.ZERO, a3 = h[0] || Q.ZERO;
        const aQ = [a1, f[2] || Q.ZERO, k.mul(a3), k.mul(f[1] || Q.ZERO), k.mul(k).mul(f[0] || Q.ZERO)];
        const kn = k.toNumber();
        const m = modelFromAinvs(aQ, { k: kn, toOriginal: (X, Y) => [[X[0] / kn, X[1] / kn], [Y[0] / kn, Y[1] / kn]] });
        m.scaled = !k.isOne();
        if (m.singular) { res.kind = 'singular cubic'; res.genus = 0; res.model = m;
          res.notes.push('Weierstrass cubic with $\\Delta = 0$: a nodal or cuspidal cubic, genus 0'); return res; }
        res.kind = 'weierstrass'; res.genus = 1; res.model = m;
        if (m.scaled) res.notes.push(`leading coefficient $${qtex(k)}$ of the cubic: plotted via $X = ${qtex(k)}\\,x,\\ Y = ${qtex(k)}\\,y$`);
        return res;
      }
      // complete the square: (y + h/2)^2 = g(x) = f + h^2/4
      const g = padd(f, pscale(pmul(h, h), new Q(1n, 4n)));
      if (!g.length) { res.kind = 'degenerate'; res.notes.push('the equation is a perfect square $(y + h(x)/2)^2 = 0$: a double curve'); return res; }
      const { a, b } = squarefreeParts(g);
      const m = pdeg(b);
      res.completedSquare = { g, a, b };
      res.singular = pdeg(a) > 0;
      if (m <= 0) { res.kind = 'reducible'; res.notes.push('$(y + h/2)^2 = c\\,a(x)^2$: two rational components'); res.genus = 0; return res; }
      res.genus = m <= 2 ? 0 : Math.floor((m - 1) / 2);
      if (res.singular) res.notes.push(`singular: $(y + h/2)^2 = a(x)^2\\,b(x)$ with $a = ${ptex(a)}$; the smooth model is $y'^2 = b(x)$`);
      if (m <= 2) { res.kind = 'rational'; res.notes.push('genus 0: rational curve'); return res; }
      if (m === 3) {
        res.kind = 'elliptic (cubic model)';
        res.model = scaledCubicModel(b, a, h);
        res.notes.push(`smooth model $y'^2 = ${ptex(b)}$ with $y' = ` + (h.length ? `\\bigl(y + \\tfrac{${ptex(h)}}{2}\\bigr)` : 'y') + (res.singular ? `/(${ptex(a)})` : '') + '$');
        return res;
      }
      if (m === 4) {
        res.kind = 'elliptic (quartic model)';
        const qm = quarticModel(b, a, h);
        if (qm) { res.model = qm; res.notes.push(`genus 1: $y'^2 = ${ptex(b)}$ is a quartic; plotted through $x = r + 1/t$ with the real root $r = ${qm.r.toPrecision(6)}$, as an elliptic curve over $\\mathbb{Q}(r)$`); }
        else if (b[pdeg(b)].sign() > 0) res.notes.push(`genus 1: $y'^2 = ${ptex(b)}$ is a quartic with no real root but real points everywhere; plotting it through one of them is planned`);
        else res.notes.push(`genus 1: $y'^2 = ${ptex(b)}$ is a quartic that is negative on $\\mathbb{R}$: no real points, and its complex points need a complex Weierstrass model, which this plotter does not have`);
        return res;
      }
      res.kind = 'hyperelliptic';
      res.notes.push(`hyperelliptic of genus ${res.genus}: $y'^2 = ${ptex(b)}$ has degree ${m}`);
      return res;
    }
    res.notes.push('quadratic in $y$ with a non-constant coefficient of $y^2$');
  }
  // general plane curve: only the arithmetic genus bound without singularity analysis
  res.kind = 'plane';
  res.notes.push(`plane curve of degree ${d}: genus $\\le ${res.genusBound}$ (equality iff the projective closure is smooth)`);
  return res;
}

// ---------------------------------------------------------------- input dispatch
const LABEL_RE = /^(\d+)(\.?)([a-z]+)(\d+)$/i;
// What did the user type?  A label comes back unresolved, {type: 'label', label, N, lmfdb}; the app looks it up in
// Cremona's tables (js/cremona.js).  A-invariants and equations are parsed here, exactly.
function parseInput(text) {
  const s = text.trim();
  const lab = s.match(LABEL_RE);
  if (lab) return { type: 'label', label: s.toLowerCase(), N: Number(lab[1]), lmfdb: lab[2] === '.' };
  const av = s.match(/^\[?\s*([+-]?\d+)\s*,\s*([+-]?\d+)\s*,\s*([+-]?\d+)\s*,\s*([+-]?\d+)\s*,\s*([+-]?\d+)\s*\]?$/);
  if (av) return { type: 'ainvs', ainvs: av.slice(1, 6).map(v => Q.of(v)) };
  return { type: 'equation', analysis: analyzeEquation(s) };
}

// ---------------------------------------------------------------- the surface grid
// Returns typed arrays for an n x n grid over the half fundamental domain s in [-1/2, 1/2], t in [0, 1/2]:
//   pos  (3 floats per vertex): (Re x, Im x, Re y) in the ORIGINAL coordinates,
//   uv   (2 floats per vertex): (s, t) for the lattice coloring,
//   ok   (1 byte per vertex):   0 at the pole / non-finite points.
function buildGrid(model, n, lattice, nterms = 10) {
  const lat = lattice || periodLattice(model.ainvs, model.inv.disc instanceof Q ? model.inv.disc.sign() : undefined, model.ainvsQ || undefined);
  const wp = makeWp(lat, nterms);
  const [a1, , a3] = model.ainvs;
  const b2 = lat.inv.b2;
  const pos = new Float32Array(n * n * 3), uv = new Float32Array(n * n * 2), ok = new Uint8Array(n * n);
  const xs = new Float64Array(n * n * 2), ys = new Float64Array(n * n * 2);       // double precision copies (for the real curves)
  let idx = 0;
  for (let i = 0; i < n; i++) {
    const s = -0.5 + i / (n - 1);
    for (let j = 0; j < n; j++, idx++) {
      const t = 0.5 * j / (n - 1);
      uv[2 * idx] = s; uv[2 * idx + 1] = t;
      const [m, nn] = toNormalisedCoords(lat, s, t);
      const w = wp(m, nn);
      if (!w) continue;
      const P = w[0], dP = w[1];
      const X = [P[0] - b2 / 12, P[1]];
      const Y = [(dP[0] - a1 * X[0] - a3) / 2, (dP[1] - a1 * X[1]) / 2];
      const [x, y] = model.toOriginal(X, Y);
      if (!isFinite(x[0]) || !isFinite(x[1]) || !isFinite(y[0]) || !isFinite(y[1])) continue;
      pos[3 * idx] = x[0]; pos[3 * idx + 1] = x[1]; pos[3 * idx + 2] = y[0];
      xs[2 * idx] = x[0]; xs[2 * idx + 1] = x[1]; ys[2 * idx] = y[0]; ys[2 * idx + 1] = y[1];
      ok[idx] = 1;
    }
  }
  return { n, pos, uv, ok, xs, ys, lattice: lat };
}
// The real components: the rows t = 0 and t = 1/2 on which x and y are real.  Returns polylines
// (arrays of [X, Y, Z]) broken at the pole, at jumps and where |P| > radius.
function realComponents(grid, radius) {
  const { n, xs, ys, ok } = grid;
  const comps = [];
  for (const j of [0, n - 1]) {
    // is the row real?
    let imax = 0, cnt = 0;
    for (let i = 0; i < n; i++) { const idx = i * n + j; if (!ok[idx]) continue; if (Math.hypot(xs[2 * idx], ys[2 * idx]) > radius) continue; imax = Math.max(imax, Math.abs(xs[2 * idx + 1]), Math.abs(ys[2 * idx + 1])); cnt++; }
    if (!cnt || imax > 1e-6) continue;
    // walk the row starting just after the pole (which is at s = 0 for row t = 0) so the unbounded component is one polyline
    const order = [];
    const start = j === 0 ? Math.floor(n / 2) + 1 : 0;
    for (let k = 0; k < n; k++) order.push((start + k) % n);
    let cur = [];
    const flush = () => { if (cur.length > 1) comps.push(cur); cur = []; };
    let prev = null;
    const before = comps.length;
    for (const i of order) {
      const idx = i * n + j;
      const p = ok[idx] ? [xs[2 * idx], 0, ys[2 * idx]] : null;
      if (!p || Math.hypot(p[0], p[2]) > radius) { flush(); prev = null; continue; }
      if (prev) {
        const d = Math.hypot(p[0] - prev[0], p[2] - prev[2]);
        if (d > 0.5 * radius) flush();
        // s = 1/2 and s = -1/2 are the same point of the torus (a 2-torsion point, on the x-axis): a zero-length
        // segment there would break the tube's frames, so a repeated point is skipped
        else if (d <= 1e-9 * (1 + Math.hypot(p[0], p[2]))) continue;
      }
      cur.push(p); prev = p;
    }
    flush();
    // close the oval: on the row t = 1/2 every point is real, and the walk is one loop
    if (j === n - 1 && cnt === n && comps.length === before + 1) {
      const c = comps[comps.length - 1], l = c[c.length - 1];
      if (Math.hypot(l[0] - c[0][0], l[2] - c[0][2]) <= 1e-9 * (1 + Math.hypot(l[0], l[2]))) c.pop();   // s = 1/2 repeats s = -1/2
      c.push(c[0]);
    }
  }
  return comps;
}

return { Q, parsePolynomial, analyzeEquation, parseInput, invariantsQ, invariantsNum, formatWeierstrass, formatWeierstrassTeX, bformat, btex, pformat, ptex, qtex,
         periodLattice, normalisePeriods, reduceTau, makeWp, cubicRootsExact, refineRootExact, toNormalisedCoords, modelFromAinvs, buildGrid, realComponents,
         squarefreeParts, polyRootsNum, cubicRoots, agm,
         _poly: { padd, psub, pmul, pdivmod, pgcd, pderiv, ptrim, pdeg } };
});
