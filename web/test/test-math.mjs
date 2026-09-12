// Node tests for ec3d-math.js against Sage/PARI reference values (vectors.json, made by make_vectors.sage).
//   node test/test-math.mjs
import { createRequire } from 'module';
import fs from 'fs';
const require = createRequire(import.meta.url);
const EC3D = require('../js/ec3d-math.js');
const V = JSON.parse(fs.readFileSync(new URL('./vectors.json', import.meta.url)));
const Cremona = require('../js/cremona.js');
Cremona.base = new URL('../data/cremona/', import.meta.url).pathname;
Cremona.fetchJSON = async p => JSON.parse(fs.readFileSync(p, 'utf8'));
const { Q } = EC3D;

let pass = 0, fail = 0; const failures = [];
function check(name, ok, detail = '') { if (ok) pass++; else { fail++; failures.push(name + (detail ? ': ' + detail : '')); } }
const relerr = (a, b) => Math.abs(a - b) / Math.max(1e-300, Math.abs(b), 1);
const crel = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1]) / Math.max(1, Math.hypot(b[0], b[1]));

// ---- the knot boundary ----
{ // on 20.a3 with the clipping radius 3: at depth 0 the boundary is the sphere cut and meets itself; at the automatic
  // depth it is an embedded curve of surface points between the sphere and the inner cutoff, one diameter clear of itself
  const model = EC3D.modelFromAinvs([0, 1, 0, -1, 0].map(v => Q.of(v)));
  const lat = EC3D.periodLattice(model.ainvs, 1, model.ainvsQ), wp = EC3D.makeWp(lat), R = 3, tube = 0.004 * R;
  const flat = EC3D.knotCurve(model, lat, wp, { radius: R, depth: 0, count: 360 });
  const onSphere = [...Array(360).keys()].every(k => Math.abs(Math.hypot(flat.pts[3 * k], flat.pts[3 * k + 1], flat.pts[3 * k + 2]) - R) < 1e-5 * R);
  check('knot: at depth 0 the boundary lies on the clipping sphere', onSphere);
  check('knot: at depth 0 the boundary meets itself', EC3D.knotClearance(flat) < 1e-3 * R, String(EC3D.knotClearance(flat)));
  const d = EC3D.autoKnotDepth(model, lat, wp, R, tube);
  const kc = EC3D.knotCurve(model, lat, wp, { radius: R, depth: d, count: 720 }), c = EC3D.knotClearance(kc);
  check('knot: the automatic depth keeps the tube one diameter clear', c >= 4 * tube && c <= 1.3 * 4 * tube, `d = ${d}, clearance ${c}, wanted ${4 * tube}`);
  const onSurface = [...Array(720).keys()].every(k => { const p = EC3D.surfacePoint(model, lat, wp, kc.st[2 * k], kc.st[2 * k + 1]); return p && Math.abs(p[0][0] - kc.pts[3 * k]) < 1e-9 && Math.abs(p[1][0] - kc.pts[3 * k + 2]) < 1e-9; });
  check('knot: the samples are points of the surface', onSurface);
  const between = [...Array(720).keys()].every(k => { const r = Math.hypot(kc.pts[3 * k], kc.pts[3 * k + 1], kc.pts[3 * k + 2]); return r <= R * (1 + 1e-6) && r >= R * (1 - d) * (1 - 1e-6); });
  check('knot: between the sphere and the inner cutoff R(1 - d)', between);
  const seeded = EC3D.knotCurve(model, lat, wp, { radius: R, depth: d, theta: 0.01, count: 360, prev: kc.rho });
  check('knot: a seeded sample at a nearby θ stays close', Math.max(...[...Array(360).keys()].map(k => Math.abs(seeded.rho[k] - kc.rho[2 * k]) / kc.rho[2 * k])) < 0.05);
}

// ---- labels ----
for (const c of V.curves) {
  const r1 = await Cremona.lookup(c.label), r2 = await Cremona.lookup(c.lmfdb);
  check(`label ${c.label}`, r1 && r1.ainvs.join() === c.ainvs.join() && r1.lmfdb === c.lmfdb, JSON.stringify(r1));
  check(`lmfdb ${c.lmfdb}`, r2 && r2.ainvs.join() === c.ainvs.join() && r2.cremona === c.label);
}
check('no such curve', (await Cremona.lookup('11a4')) === null);
check('beyond the tables', (await Cremona.lookup('999999a1')).missing === true);
{ // the last shard: every curve is found under both labels and the labels round-trip
  const idx = await Cremona.loadIndex(), last = idx.shards[idx.shards.length - 1], shard = await Cremona.ensure(last.lo);
  let n = 0, bad = 0;
  for (const N in shard) for (const [cl, ll, curves] of shard[N]) for (let i = 0; i < curves.length; i++) {
    const a = await Cremona.lookup(`${N}${cl}${i + 1}`), b = await Cremona.lookup(`${N}.${ll}${curves[i][0]}`);
    n++; if (!a || !b || a.lmfdb !== b.lmfdb || a.cremona !== b.cremona || a.ainvs.join() !== curves[i].slice(1, 6).join() || !Number.isInteger(a.rank) || !Array.isArray(a.torsion)) bad++;
    if (n >= 20000) break;
  }
  check(`shard ${last.file}: labels round-trip`, n > 1000 && bad === 0, `${bad} bad of ${n}`);
  const r5077 = await Cremona.lookup('5077a1'), r15 = await Cremona.lookup('15a1'), r108 = await Cremona.lookup('108a1');
  check('Mordell-Weil: 5077a1 has rank 3, no torsion', r5077.rank === 3 && r5077.torsion.length === 0, JSON.stringify(r5077));
  check('Mordell-Weil: 15a1 has torsion Z/2 x Z/4', r15.rank === 0 && r15.torsion.join() === '2,4', JSON.stringify(r15));
  check('Mordell-Weil: 108a1 (y^2 = x^3 + 4) has torsion Z/3', r108.rank === 0 && r108.torsion.join() === '3', JSON.stringify(r108));
  check('index: 50 shards to 499999', idx.shards.length === 50 && idx.max_conductor === 499999 && idx.curves > 3000000, JSON.stringify([idx.shards.length, idx.max_conductor, idx.curves]));
}

// ---- invariants, periods, wp ----
let altBasis = 0;
for (const c of V.curves) {
  const aQ = c.ainvs.map(v => Q.of(v));
  const inv = EC3D.invariantsQ(aQ);
  check(`disc ${c.label}`, inv.disc.n === BigInt(c.disc) && inv.disc.d === 1n, inv.disc.toString());
  check(`j ${c.label}`, inv.j.n === BigInt(c.j[0]) && inv.j.d === BigInt(c.j[1]), inv.j.toString());
  const lat = EC3D.periodLattice(c.ainvs, inv.disc.sign(), aQ);
  check(`w1 ${c.label}`, relerr(lat.w1, c.w1) < 1e-12, `${lat.w1} vs ${c.w1}`);
  { const mine = lat.e.map(r => Array.isArray(r) ? r : [r, 0]).sort((p, q) => (p[0] - q[0]) || (p[1] - q[1]));
    const sage = c.two_division_roots.slice().sort((p, q) => (p[0] - q[0]) || (p[1] - q[1]));
    check(`2-division roots ${c.label}`, mine.every((r, i) => crel(r, sage[i]) < 1e-11), `${JSON.stringify(mine)} vs ${JSON.stringify(sage)}`); }
  check(`w2 ${c.label}`, crel(lat.w2, c.w2) < 1e-12, `${lat.w2} vs ${c.w2}`);
  const { n1, n2, mat, tau } = lat.normalised;
  const same = crel(n1, c.normalised[0]) < 1e-11 && crel(n2, c.normalised[1]) < 1e-11 && mat.join() === c.matrix.join();
  if (!same) {
    altBasis++;
    const det = mat[0] * mat[3] - mat[1] * mat[2];
    check(`normalised basis valid ${c.label}`, Math.abs(det) === 1 && Math.hypot(tau[0], tau[1]) >= 0.999 && Math.abs(tau[0]) <= 0.5 + 1e-9 && tau[1] > 0.86,
          `mat ${mat} tau ${tau} (Sage: ${c.matrix})`);
  } else pass++;
  // wp at the sample points
  const wp = EC3D.makeWp(lat);
  const b2 = inv.b2.toNumber(), a1 = c.ainvs[0], a3 = c.ainvs[2];
  for (const smp of c.samples) {
    // z = m n1 + n n2 : solve the real 2x2 system
    const det = n1[0] * n2[1] - n2[0] * n1[1];
    const m = (smp.z[0] * n2[1] - n2[0] * smp.z[1]) / det, n = (n1[0] * smp.z[1] - smp.z[0] * n1[1]) / det;
    const w = wp(m, n);
    if (!w) { check(`wp ${c.label} pole?`, false, `z=${smp.z}`); continue; }
    check(`wp ${c.label} z=${smp.z}`, crel(w[0], smp.wp) < 1e-9, `${w[0]} vs ${smp.wp}`);
    check(`wp' ${c.label} z=${smp.z}`, crel(w[1], smp.wpp) < 1e-8 || Math.hypot(...w[1]) < 1e-9, `${w[1]} vs ${smp.wpp}`);
    const X = [w[0][0] - b2 / 12, w[0][1]], Y = [(w[1][0] - a1 * X[0] - a3) / 2, (w[1][1] - a1 * X[1]) / 2];
    check(`point ${c.label}`, crel(X, smp.x) < 1e-9 && crel(Y, smp.y) < 1e-8, `${X},${Y} vs ${smp.x},${smp.y}`);
  }
}
console.log(`normalised basis differs from Sage's (but valid) for ${altBasis} of ${V.curves.length} curves`);

// ---- grid + real components ----
async function gridChecks(label, expectComponents, radius = 3) {
  const rec = await Cremona.lookup(label);
  const model = EC3D.modelFromAinvs(rec.ainvs.map(v => Q.of(v)));
  const g = EC3D.buildGrid(model, 101);
  const [a1, a2, a3, a4, a6] = rec.ainvs;
  let maxres = 0, imRow0 = 0, cnt = 0;
  for (let i = 0; i < g.n; i++) for (let j = 0; j < g.n; j++) {
    const idx = i * g.n + j; if (!g.ok[idx]) continue;
    const x = [g.xs[2 * idx], g.xs[2 * idx + 1]], y = [g.ys[2 * idx], g.ys[2 * idx + 1]];
    if (Math.hypot(...x) > 50) continue;
    // y^2 + a1 x y + a3 y - (x^3 + a2 x^2 + a4 x + a6)
    const cm = (p, q) => [p[0] * q[0] - p[1] * q[1], p[0] * q[1] + p[1] * q[0]];
    const x2 = cm(x, x), x3 = cm(x2, x), y2 = cm(y, y), xy = cm(x, y);
    const r = [y2[0] + a1 * xy[0] + a3 * y[0] - (x3[0] + a2 * x2[0] + a4 * x[0] + a6), y2[1] + a1 * xy[1] + a3 * y[1] - (x3[1] + a2 * x2[1] + a4 * x[1])];
    maxres = Math.max(maxres, Math.hypot(...r)); cnt++;
    if (j === 0) imRow0 = Math.max(imRow0, Math.abs(x[1]), Math.abs(y[1]));
  }
  check(`grid ${label} on curve`, cnt > 5000 && maxres < 1e-8, `max residual ${maxres} over ${cnt} points`);
  check(`grid ${label} row t=0 real`, imRow0 < 1e-8, `${imRow0}`);
  check(`grid ${label} pole marked`, g.ok[Math.floor(g.n / 2) * g.n + 0] === 0);
  const comps = EC3D.realComponents(g, radius);
  check(`real components ${label}`, comps.length === expectComponents, `${comps.length} components: ${comps.map(c => c.length)}`);
  return comps;
}
const comps20 = await gridChecks('20.a3', 2);
check('oval closed', comps20.some(c => c.length > 50 && Math.hypot(c[0][0] - c[c.length - 1][0], c[0][2] - c[c.length - 1][2]) < 1e-9));
await gridChecks('11a1', 1, 8);      // its real points start at x = 4.35, outside a radius-3 ball
await gridChecks('11a1', 0, 3);
await gridChecks('37a1', 2);

// ---- equations ----
const eq = (s) => EC3D.analyzeEquation(s);
const ain = (r) => r.model && r.model.ainvsQ && r.model.ainvsQ.map(v => v.toString()).join(',');
check('parse 20.a3 eq', ain(eq('y^2 = x^3 + x^2 - x')) === '0,1,0,-1,0' && eq('y^2 = x^3 + x^2 - x').kind === 'weierstrass');
check('parse 11a1 eq', ain(eq('y^2 + y = x^3 - x^2 - 10x - 20')) === '0,-1,1,-10,-20');
check('parse unicode', ain(eq('y² = x³ − x')) === '0,0,0,-1,0');
check('parse no spaces', ain(eq('y^2=x^3+x')) === '0,0,0,1,0');
check('parse ** and implicit mult', ain(eq('y**2 = x**3 - 3x + 2')) === '0,0,0,-3,2' || eq('y**2 = x**3 - 3x + 2').kind === 'singular cubic');
check('parse product form', ain(eq('(y)^2 = x(x-1)(x+1)')) === '0,0,0,-1,0');
check('parse a1 term', ain(eq('y^2 + x*y = x^3 + 1')) === '1,0,0,0,1');
const sc = eq('2*y^2 = x^3 - x');
check('scaled model', sc.kind === 'weierstrass' && ain(sc) === '0,0,0,-1/4,0' && sc.model.scaled, ain(sc));
check('scaled map back', (() => { const [x, y] = sc.model.toOriginal([1, 0], [2, 0]); return Math.abs(x[0] - 2) < 1e-12 && Math.abs(y[0] - 4) < 1e-12; })());
check('leading 2', ain(eq('y^2 = 2*x^3 + 1')) === '0,0,0,0,4');
for (const [s, gexp] of Object.entries(V.genus)) {
  const r = eq(s);
  if (typeof r.genus === 'number') check(`genus ${s}`, r.genus === gexp, `got ${r.genus} (${r.kind}), Sage ${gexp}`);
  else check(`genus bound ${s}`, r.genusBound >= gexp, `bound ${r.genusBound} < Sage ${gexp}`);
}
check('quartic model exists', eq('y^2 + x^2*y = x^3 + 1').model !== null && eq('y^2 + x^2*y = x^3 + 1').kind === 'elliptic (quartic model)');
check('reducible detected', eq('y^2 + (x^2+1)*y = x^3 - x').kind === 'reducible');
check('genus 2 from completed square', eq('y^2 + x^3*y = x + 2').genus === 2);
check('quartic without real root', eq('y^2 = x^4 + 1').model === null && eq('y^2 = x^4 + 1').genus === 1);
// points of the quartic and cubic-square models satisfy the original equation
function evalF(F, x, y) {
  const cm = (p, q) => [p[0] * q[0] - p[1] * q[1], p[0] * q[1] + p[1] * q[0]];
  let r = [0, 0];
  for (const [k, c] of F) { const [i, j] = k.split(',').map(Number); let t = [c.toNumber(), 0]; for (let a = 0; a < i; a++) t = cm(t, x); for (let b = 0; b < j; b++) t = cm(t, y); r = [r[0] + t[0], r[1] + t[1]]; }
  return Math.hypot(...r);
}
for (const s of ['y^2 + x^2*y = x^3 + 1', 'y^2 = x^4 - 3x^2 + x + 1', 'y^2 + (x^2+1)*y = x^3 + x', 'y^2 = (x-1)^2*(x^3 - x)', '3*y^2 + x*y = 2*x^3 + x - 5']) {
  const r = eq(s);
  if (!r.model) { check(`model for ${s}`, false, r.notes.join('; ')); continue; }
  const g = EC3D.buildGrid(r.model, 61);
  let worst = 0, cnt = 0;
  for (let idx = 0; idx < g.n * g.n; idx++) { if (!g.ok[idx]) continue; const x = [g.xs[2 * idx], g.xs[2 * idx + 1]], y = [g.ys[2 * idx], g.ys[2 * idx + 1]]; if (Math.hypot(...x) > 20 || Math.hypot(...y) > 100) continue; worst = Math.max(worst, evalF(r.F, x, y)); cnt++; }
  check(`model points on original curve: ${s} [${r.kind}]`, cnt > 500 && worst < 1e-6, `worst ${worst} over ${cnt}`);
}
for (const bad of ['y^2 = x^3 + z', 'y^2 = x^(3', 'y^2 / x = 1', 'y^2 = x^3 = 1', 'y^2 = x^-1']) {
  let threw = false; try { eq(bad); } catch (e) { threw = true; }
  check(`rejects ${bad}`, threw);
}
// ---- input dispatch ----
{ const p = EC3D.parseInput('20.a3'); check('input label', p.type === 'label' && p.N === 20 && p.lmfdb === true && (await Cremona.lookup(p.label)).ainvs.join() === '0,1,0,-1,0'); }
{ const p = EC3D.parseInput('20a2'); check('input cremona', p.type === 'label' && p.lmfdb === false && (await Cremona.lookup(p.label)).lmfdb === '20.a3'); }
check('input ainvs', EC3D.parseInput('[0,1,0,-1,0]').ainvs.map(String).join() === '0,1,0,-1,0');
check('input ainvs bare', EC3D.parseInput('0, 1, 0, -1, 0').type === 'ainvs');
check('input label unresolved', EC3D.parseInput('999999a1').type === 'label');
check('input equation', EC3D.parseInput('y^2 = x^3 - x').type === 'equation');
check('findByAinvs over loaded shards', (Cremona.findByAinvs([0, 1, 0, -1, 0]) || {}).lmfdb === '20.a3');
check('format', EC3D.formatWeierstrass([0, -1, 1, -10, -20].map(v => Q.of(v))) === 'y² + y = x³ - x² - 10x - 20', EC3D.formatWeierstrass([0, -1, 1, -10, -20].map(v => Q.of(v))));

{ const seen = {}; for (const f of failures) { const k = f.split(' ')[0]; seen[k] = (seen[k] || 0) + 1; if (seen[k] <= 3) console.log('  FAIL', f); } }
const cats = {}; for (const f of failures) { const k = f.split(' ')[0]; cats[k] = (cats[k] || 0) + 1; }
console.log('failure categories:', JSON.stringify(cats));
console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
