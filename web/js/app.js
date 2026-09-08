/* app.js -- viewer and UI for the EllipticCurve3D web app.  Globals: THREE, katex, EC3D, Cremona, COLORMAPS. */
(function () {
'use strict';
const $ = id => document.getElementById(id);
const isMobile = matchMedia('(max-width: 640px)').matches || navigator.maxTouchPoints > 1;

const state = { grid: isMobile ? 250 : 400, radius: 3, cutoff: 0, opacity: 0.95, lines: 3, soft: 0.33, colormap: 'Greens_r',
                real: true, thick: 0.004, mirror: true, axes: true, sphere: false,
                slice: 're', theta: 0, playing: false, speed: 30 };     // third coordinate: Re y (θ = 0), Im y (θ = π/2), or rotating
let model = null, lattice = null, gridData = null, curveInfo = null, lastText = '';

// ------------------------------------------------------------------ scene
const view = $('view');
const renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
renderer.setSize(view.clientWidth, view.clientHeight);
renderer.setClearColor(0xffffff, 1);
view.appendChild(renderer.domElement);
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(45, view.clientWidth / view.clientHeight, 0.01, 2000);
camera.up.set(0, 0, 1);
const controls = new THREE.OrbitControls(camera, renderer.domElement);
controls.enableDamping = true; controls.dampingFactor = 0.12;
controls.touches = { ONE: THREE.TOUCH.ROTATE, TWO: THREE.TOUCH.DOLLY_PAN };
controls.addEventListener('change', requestRender);
const keyLight = new THREE.DirectionalLight(0xffffff, 0.85); keyLight.position.set(-6, 4, 8); camera.add(keyLight); scene.add(camera);
scene.add(new THREE.AmbientLight(0xffffff, 0.6));
const fill = new THREE.DirectionalLight(0xffffff, 0.3); fill.position.set(4, -5, -7); scene.add(fill);

let dirty = true;
function requestRender() { dirty = true; }
let lastFrame = performance.now();
(function loop(now) {
  requestAnimationFrame(loop);
  const dt = Math.min(0.1, (now - lastFrame) / 1000); lastFrame = now;
  if (state.playing) setTheta(state.theta + state.speed * Math.PI / 180 * dt);
  const moved = controls.update();
  if (moved || dirty) { renderer.render(scene, camera); dirty = false; }
})(lastFrame);
window.addEventListener('resize', () => { renderer.setSize(view.clientWidth, view.clientHeight); camera.aspect = view.clientWidth / view.clientHeight; camera.updateProjectionMatrix(); requestRender(); });

function resetView() {
  const R = state.radius, dist = R / Math.sin(Math.PI / 8) * 0.9;
  const s3 = Math.sin(Math.PI / 3), c3 = Math.cos(Math.PI / 3), s5 = Math.sin(Math.PI / 5), c5 = Math.cos(Math.PI / 5);
  camera.position.set(dist * s3 * c5, -dist * s3 * s5, dist * c3);      // like Sage's default viewpoint, from the -Im x side
  controls.target.set(0, 0, 0); camera.near = 0.01 * R; camera.far = 200 * R; camera.updateProjectionMatrix(); controls.update(); requestRender();
}

// ------------------------------------------------------------------ colormap texture
const cmapTex = new THREE.DataTexture(new Uint8Array(64 * 4), 64, 1, THREE.RGBAFormat);
cmapTex.minFilter = cmapTex.magFilter = THREE.LinearFilter; cmapTex.wrapS = cmapTex.wrapT = THREE.ClampToEdgeWrapping;
function setColormap(name) {
  const t = COLORMAPS[name] || COLORMAPS.Greens_r, d = cmapTex.image.data;
  for (let i = 0; i < 64; i++) { d[4 * i] = Math.round(255 * t[i][0]); d[4 * i + 1] = Math.round(255 * t[i][1]); d[4 * i + 2] = Math.round(255 * t[i][2]); d[4 * i + 3] = 255; }
  cmapTex.needsUpdate = true; requestRender();
}
setColormap(state.colormap);

// ------------------------------------------------------------------ surface material
// The third coordinate is z = Re y cos θ + Im y sin θ, computed in the vertex shader from the attributes position.z = Re y
// and yim = Im y, so Re y (θ = 0), Im y (θ = π/2) and the rotation between them need no rebuild.  The normal is the cross
// product of the surface's tangent vectors dS, dT (finite differences of (Re x, Im x, Re y, Im y) along the grid), mixed
// the same way, so lighting is exact at every θ.  The mirror half (the complex conjugate) has uMirrorY = -1.
// Also here: the lattice coloring and the spherical clipping.
const uniforms = { uRadius: { value: state.radius }, uCutoff: { value: state.cutoff }, uLines: { value: state.lines }, uSoft: { value: state.soft }, uColormap: { value: cmapTex } };
const thetaUniform = { value: 0 };
function makeSurfaceMaterial(mirrorY) {
  const m = new THREE.MeshPhongMaterial({ side: THREE.DoubleSide, transparent: true, opacity: state.opacity, depthWrite: state.opacity >= 1,
                                          shininess: 25, specular: new THREE.Color(0x2a2a2a), color: 0xffffff });
  m.onBeforeCompile = shader => {
    Object.assign(shader.uniforms, uniforms, { uTheta: thetaUniform, uMirrorY: { value: mirrorY } });
    shader.vertexShader = shader.vertexShader
      .replace('#include <common>', '#include <common>\nattribute vec2 st;\nattribute float yim;\nattribute vec4 dS;\nattribute vec4 dT;\nuniform float uTheta;\nuniform float uMirrorY;\nvarying vec2 vST;\nvarying vec3 vPos;')
      .replace('#include <beginnormal_vertex>', 'float cs = cos(uTheta), sn = sin(uTheta) * uMirrorY;\nvec3 tS = vec3(dS.x, dS.y, cs * dS.z + sn * dS.w);\nvec3 tT = vec3(dT.x, dT.y, cs * dT.z + sn * dT.w);\nvec3 objectNormal = normalize(cross(tS, tT));')
      .replace('#include <begin_vertex>', 'vec3 transformed = vec3(position.x, position.y, cs * position.z + sn * yim);\nvST = st;\nvPos = transformed;');
    shader.fragmentShader = shader.fragmentShader
      .replace('#include <common>', '#include <common>\nuniform float uRadius;\nuniform float uCutoff;\nuniform float uLines;\nuniform float uSoft;\nuniform sampler2D uColormap;\nvarying vec2 vST;\nvarying vec3 vPos;')
      .replace('#include <clipping_planes_fragment>', '#include <clipping_planes_fragment>\n{ float rr = length(vPos); if (rr > uRadius || rr < uCutoff) discard; }')
      .replace('#include <color_fragment>', '#include <color_fragment>\n{ float c = pow(abs(sin(3.14159265358979 * vST.x * uLines) * sin(3.14159265358979 * vST.y * uLines)), uSoft);\n  diffuseColor.rgb = texture2D(uColormap, vec2(c, 0.5)).rgb; }');
  };
  m.customProgramCacheKey = () => 'ec3d-surface';
  return m;
}
const surfaceMaterial = makeSurfaceMaterial(1), mirrorMaterial = makeSurfaceMaterial(-1);
const realMaterial = new THREE.MeshPhongMaterial({ color: 0xc8281e, shininess: 40, specular: new THREE.Color(0x442222) });

// ------------------------------------------------------------------ scene objects
const group = new THREE.Group(); scene.add(group);
let surfaceMesh = null, mirrorMesh = null, realGroup = null, axesGroup = null, sphereMesh = null;
let builtBig = 2000;                                          // |P| beyond which grid cells were dropped when the geometry was built

function buildSurfaceGeometry(g) {
  const n = g.n, N = n * n, idx = [];
  const pos = new Float32Array(3 * N), yim = new Float32Array(N);          // (Re x, Im x, Re y) and Im y
  for (let i = 0; i < N; i++) { pos[3 * i] = g.xs[2 * i]; pos[3 * i + 1] = g.xs[2 * i + 1]; pos[3 * i + 2] = g.ys[2 * i]; yim[i] = g.ys[2 * i + 1]; }
  const big = builtBig = Math.max(2000, 50 * state.radius);                 // drop cells at the pole itself (|P| enormous)
  const okv = new Uint8Array(N);
  for (let i = 0; i < N; i++) okv[i] = (g.ok[i] && Math.abs(pos[3 * i]) < big && Math.abs(pos[3 * i + 1]) < big && Math.abs(pos[3 * i + 2]) < big && Math.abs(yim[i]) < big) ? 1 : 0;
  for (let i = 0; i < n - 1; i++) for (let j = 0; j < n - 1; j++) {
    const a = i * n + j, b = (i + 1) * n + j, c = (i + 1) * n + j + 1, d = i * n + j + 1;
    if (okv[a] && okv[b] && okv[c] && okv[d]) idx.push(a, b, c, a, c, d);
  }
  // tangent vectors along the grid, (d/ds, d/dt) of (Re x, Im x, Re y, Im y), by central differences where both neighbours exist
  const dS = new Float32Array(4 * N), dT = new Float32Array(4 * N);
  const comp = (k, c) => c < 2 ? g.xs[2 * k + c] : g.ys[2 * k + c - 2];
  const diff = (out, k, ka, kb) => { for (let c = 0; c < 4; c++) out[4 * k + c] = comp(kb, c) - comp(ka, c); };
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    const k = i * n + j;
    if (!okv[k]) continue;
    const im = i > 0 && okv[k - n] ? k - n : k, ip = i < n - 1 && okv[k + n] ? k + n : k;
    const jm = j > 0 && okv[k - 1] ? k - 1 : k, jp = j < n - 1 && okv[k + 1] ? k + 1 : k;
    diff(dS, k, im, ip); diff(dT, k, jm, jp);
  }
  const geom = new THREE.BufferGeometry();
  geom.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  geom.setAttribute('yim', new THREE.BufferAttribute(yim, 1));
  geom.setAttribute('dS', new THREE.BufferAttribute(dS, 4));
  geom.setAttribute('dT', new THREE.BufferAttribute(dT, 4));
  geom.setAttribute('st', new THREE.BufferAttribute(g.uv, 2));
  geom.setIndex(new THREE.BufferAttribute(new Uint32Array(idx), 1));
  return geom;
}
function buildRealCurves(comps) {
  const grp = new THREE.Group();
  for (const comp of comps) {
    const pts = [];
    for (const p of comp) {                                       // no repeated points: they pinch the tube
      const v = new THREE.Vector3(p[0], p[1], p[2]);
      if (!pts.length || pts[pts.length - 1].distanceTo(v) > 1e-9 * (1 + v.length())) pts.push(v);
    }
    const closed = pts.length > 3 && pts[0].distanceTo(pts[pts.length - 1]) < 1e-9 * (1 + pts[0].length());
    if (closed) pts.pop();
    if (pts.length < 2) continue;
    const curve = new THREE.CatmullRomCurve3(pts, closed, 'centripetal');
    const tube = new THREE.TubeGeometry(curve, Math.min(800, 2 * pts.length), state.thick * state.radius, 8, closed);
    grp.add(new THREE.Mesh(tube, realMaterial));
  }
  return grp;
}
const SERIF = '"STIX Two Text", "STIX Two Math", "Times New Roman", Times, serif';
function makeLabel(sym, variable, x, y, z) {                 // e.g. ℜ x: the Fraktur symbol upright, the variable in italics
  const canvas = document.createElement('canvas'), ctx = canvas.getContext('2d'), pr = 2, fs = 24;
  const f1 = `${fs}px ${SERIF}`, f2 = `italic ${fs}px ${SERIF}`;
  ctx.font = f1; const w1 = ctx.measureText(sym + ' ').width; ctx.font = f2; const w2 = ctx.measureText(variable).width;
  const w = Math.ceil(w1 + w2) + 8, h = fs + 10;
  canvas.width = w * pr; canvas.height = h * pr; ctx.scale(pr, pr);
  ctx.fillStyle = '#222'; ctx.textBaseline = 'middle'; ctx.textAlign = 'left';
  ctx.font = f1; ctx.fillText(sym + ' ', 4, h / 2); ctx.font = f2; ctx.fillText(variable, 4 + w1, h / 2);
  const tex = new THREE.Texture(canvas); tex.needsUpdate = true;
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, sizeAttenuation: false, depthWrite: false, transparent: true }));
  sp.position.set(x, y, z); sp.scale.set(w / 700, h / 700, 1);
  return sp;
}
function buildAxes(R) {
  const grp = new THREE.Group();
  grp.add(new THREE.AxesHelper(1.15 * R));
  grp.add(makeLabel('ℜ', 'x', 1.3 * R, 0, 0)); grp.add(makeLabel('ℑ', 'x', 0, 1.3 * R, 0));
  grp.add(state.slice === 'anim' ? makeLabel('ℜy cos θ + ℑy sin θ', '', 0, 0, 1.3 * R) : makeLabel(state.slice === 'im' ? 'ℑ' : 'ℜ', 'y', 0, 0, 1.3 * R));
  return grp;
}
function rebuildDecorations() {
  const R = state.radius;
  if (axesGroup) group.remove(axesGroup);
  axesGroup = buildAxes(R); axesGroup.visible = state.axes; group.add(axesGroup);
  if (sphereMesh) group.remove(sphereMesh);
  sphereMesh = new THREE.Mesh(new THREE.SphereGeometry(R, 48, 32), new THREE.MeshBasicMaterial({ color: 0x777777, wireframe: true, transparent: true, opacity: 0.12, depthWrite: false }));
  sphereMesh.visible = state.sphere; group.add(sphereMesh);
  if (gridData) {
    if (realGroup) group.remove(realGroup);
    realGroup = buildRealCurves(EC3D.realComponents(gridData, R)); realGroup.visible = state.real && state.slice === 're'; group.add(realGroup);
  }
  requestRender();
}
function computeGrid() {                                       // wp on the grid; returns the time taken
  const t0 = performance.now();
  gridData = EC3D.buildGrid(model, state.grid, lattice);
  return performance.now() - t0;
}
function rebuildSurface() {                                    // the mesh from the current grid
  if (!gridData) return { mesh: 0, faces: 0 };
  const t1 = performance.now();
  const geom = buildSurfaceGeometry(gridData);
  if (surfaceMesh) { group.remove(surfaceMesh); surfaceMesh.geometry.dispose(); group.remove(mirrorMesh); }
  surfaceMesh = new THREE.Mesh(geom, surfaceMaterial); group.add(surfaceMesh);
  // the other half of the torus is the complex conjugate: Im x -> -Im x by the scale, Im y -> -Im y in the shader (uMirrorY)
  mirrorMesh = new THREE.Mesh(geom, mirrorMaterial); mirrorMesh.scale.set(1, -1, 1); mirrorMesh.visible = state.mirror; group.add(mirrorMesh);
  rebuildDecorations();
  return { mesh: performance.now() - t1, faces: geom.index.count / 3 };
}
// The default clipping radius: every real component in view, i.e. 1.5 times the distance from the origin to
// the farthest real 2-torsion point (the turning points of E(R), at s = 0, +-1/2 on the real rows), at least 3,
// rounded up to 2 digits.  11.a1's real points start at x = 103, so it gets 160.
function autoRadius(g) {
  const { n, xs, ys, ok } = g;
  let r = 0;
  const mid = [Math.floor((n - 1) / 2), Math.ceil((n - 1) / 2)];               // s = 0 (the pole on the row t = 0, so only on t = 1/2)
  for (const [j, is] of [[0, [0, n - 1]], [n - 1, [0, ...mid, n - 1]]]) for (const i of is) {
    const idx = i * n + j;
    if (!ok[idx] || Math.abs(xs[2 * idx + 1]) > 1e-6 || Math.abs(ys[2 * idx + 1]) > 1e-6) continue;
    r = Math.max(r, Math.hypot(xs[2 * idx], ys[2 * idx]));
  }
  const R = Math.max(3, 1.5 * r), p = Math.pow(10, Math.floor(Math.log10(R)) - 1);
  return Math.ceil(R / p) * p;
}

// ------------------------------------------------------------------ the picture of the lattice coloring
const latticeCanvas = $('lattice-plot'), lctx = latticeCanvas.getContext('2d');
let latticePlotQueued = false;
function scheduleLatticePlot() { if (!latticePlotQueued) { latticePlotQueued = true; requestAnimationFrame(() => { latticePlotQueued = false; drawLatticePlot(); }); } }
function drawLatticePlot() {
  const W = latticeCanvas.width, H = latticeCanvas.height;
  lctx.clearRect(0, 0, W, H);
  if (!lattice) return;
  const w1 = lattice.w1, w2 = lattice.w2;                      // ω₁ > 0 real, ω₂ = [re, im]
  const xs = [0, w1, w2[0], w1 + w2[0]], x0 = Math.min(...xs), x1 = Math.max(...xs), y1 = w2[1];
  const mx = 0.45 * (x1 - x0), my = 0.45 * y1;                  // show some of the neighbouring translates too
  const scale = Math.min(W / (x1 - x0 + 2 * mx), H / (y1 + 2 * my)), cx = (x0 + x1) / 2, cy = y1 / 2;
  const toPx = (x, y) => [W / 2 + (x - cx) * scale, H / 2 - (y - cy) * scale];
  const dark = matchMedia('(prefers-color-scheme: dark)').matches && document.documentElement.dataset.theme !== 'light';
  const bg = dark ? 24 : 255, cm = COLORMAPS[state.colormap] || COLORMAPS.Greens_r, n = state.lines, soft = state.soft;
  const img = lctx.createImageData(W, H), d = img.data;
  for (let py = 0; py < H; py++) for (let px = 0; px < W; px++) {
    const x = cx + (px + 0.5 - W / 2) / scale, y = cy - (py + 0.5 - H / 2) / scale;
    const t = y / w2[1], s = (x - t * w2[0]) / w1;               // z = s ω₁ + t ω₂
    const c = Math.pow(Math.abs(Math.sin(Math.PI * s * n) * Math.sin(Math.PI * t * n)), soft);
    const k = Math.min(63, Math.max(0, Math.round(c * 63))), o = 4 * (py * W + px);
    const f = (s >= 0 && s <= 1 && t >= 0 && t <= 1) ? 1 : 0.45;                 // dim outside the period parallelogram
    d[o] = 255 * cm[k][0] * f + bg * (1 - f); d[o + 1] = 255 * cm[k][1] * f + bg * (1 - f); d[o + 2] = 255 * cm[k][2] * f + bg * (1 - f); d[o + 3] = 255;
  }
  lctx.putImageData(img, 0, 0);
  const accent = getComputedStyle(document.documentElement).getPropertyValue('--accent-2').trim() || '#b3261e';
  const text = getComputedStyle(document.documentElement).getPropertyValue('--text').trim() || '#222';
  const corners = [[0, 0], [w1, 0], [w1 + w2[0], w2[1]], [w2[0], w2[1]]].map(p => toPx(p[0], p[1]));
  lctx.lineWidth = 4; lctx.strokeStyle = accent; lctx.lineJoin = 'round';
  lctx.beginPath(); corners.forEach((p, i) => (i ? lctx.lineTo : lctx.moveTo).call(lctx, p[0], p[1])); lctx.closePath(); lctx.stroke();
  lctx.setLineDash([8, 8]); lctx.lineWidth = 2.5;                 // the half 0 <= t <= 1/2 actually drawn
  const half = [[0, 0], [w1, 0], [w1 + w2[0] / 2, w2[1] / 2], [w2[0] / 2, w2[1] / 2]].map(p => toPx(p[0], p[1]));
  lctx.beginPath(); half.forEach((p, i) => (i ? lctx.lineTo : lctx.moveTo).call(lctx, p[0], p[1])); lctx.closePath(); lctx.stroke();
  lctx.setLineDash([]);
  lctx.fillStyle = text; lctx.font = `italic 24px ${SERIF}`; lctx.textBaseline = 'middle';
  const label = (txt, x, y, dx, dy) => { const p = toPx(x, y); lctx.beginPath(); lctx.arc(p[0], p[1], 5, 0, 2 * Math.PI); lctx.fill(); lctx.fillText(txt, p[0] + dx, p[1] + dy); };
  label('0', 0, 0, 10, 16); label('ω₁', w1, 0, 10, 16); label('ω₂', w2[0], w2[1], -34, -14); label('ω₁+ω₂', w1 + w2[0], w2[1], 10, -14);
}
function clearSurface() {
  if (surfaceMesh) { group.remove(surfaceMesh); group.remove(mirrorMesh); surfaceMesh = mirrorMesh = null; }
  if (realGroup) { group.remove(realGroup); realGroup = null; }
  gridData = null; requestRender();
}

// ------------------------------------------------------------------ input handling
const fmt = (v, d = 5) => (Math.abs(v) < 1e-12 ? '0' : Number(v.toPrecision(d)).toString());
const cfmt = z => { const re = fmt(z[0]), im = fmt(Math.abs(z[1])); if (Math.abs(z[1]) < 1e-12) return re; if (Math.abs(z[0]) < 1e-12) return (z[1] < 0 ? '-' : '') + im + 'i'; return `${re} ${z[1] < 0 ? '-' : '+'} ${im}i`; };
// KaTeX: T renders TeX, mixed renders "text $tex$ text", qtex/ctex format numbers for TeX
const T = tex => katex.renderToString(tex, { throwOnError: false, output: 'html' });
const esc = t => String(t).replace(/[&<>]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));
const mixed = str => String(str).split('$').map((part, i) => i % 2 ? T(part) : esc(part)).join('');
const qtex = v => (v instanceof EC3D.Q) ? EC3D.qtex(v) : fmt(v);
const ctex = z => { const re = fmt(z[0]), im = fmt(Math.abs(z[1])); if (Math.abs(z[1]) < 1e-12) return re; if (Math.abs(z[0]) < 1e-12) return (z[1] < 0 ? '-' : '') + im + '\\,i'; return `${re} ${z[1] < 0 ? '-' : '+'} ${im}\\,i`; };
for (const el of document.querySelectorAll('.tex')) katex.render(el.textContent, el, { throwOnError: false, output: 'html' });
function setInfo(html, cls) { const el = $('info'); el.innerHTML = html; el.className = cls || ''; }
let infoParts = null;                                          // the info line, re-rendered when the slice changes
const sliceText = () => mixed(state.slice === 'im' ? 'imaginary slice $(\\Re x, \\Im x, \\Im y)$'
                            : state.slice === 'anim' ? 'rotating $(\\Re x, \\Im x, \\Re y\\cos\\theta + \\Im y\\sin\\theta)$'
                            : 'surface $(\\Re x, \\Im x, \\Re y)$');
function renderInfo() { if (infoParts) setInfo([...infoParts.base, sliceText()].join(' · ')); }
function updateHash(text) { try { history.replaceState(null, '', '#' + (state.slice === 'im' ? 'im:' : state.slice === 'anim' ? 'anim:' : '') + encodeURIComponent(text)); } catch (e) {} }
// A typed curve is identified up to isomorphism over Q: the tables list minimal models, so the a-invariants are
// first made integral by (x, y) -> (u^2 x, u^3 y), a_i -> u^i a_i, and then reduced by the same scaling where possible.
function integralModel(aQ) {
  const w = [1, 2, 3, 4, 6];
  let a = null;
  for (let u = 1; u <= 100 && !a; u++) { const b = aQ.map((v, i) => v.mul(EC3D.Q.of(u).pow(w[i]))); if (b.every(v => v.isInteger())) a = b; }
  if (!a) return null;
  for (let u = 2; u <= 30; ) { const b = a.map((v, i) => v.div(EC3D.Q.of(u).pow(w[i]))); if (b.every(v => v.isInteger())) a = b; else u++; }
  return a.map(v => Number(v.toString()));
}
function identify(aQ) {
  const nums = aQ.map(v => Number(v.toString()));
  let rec = nums.every(Number.isInteger) ? Cremona.findByAinvs(nums) : null;
  if (!rec) { const m = integralModel(aQ); if (m) rec = Cremona.findByAinvs(m); }
  return rec;
}
// the Mordell-Weil group from Cremona's tables: E(Q) = Z^r (+) torsion
function mwText(rec) {
  if (!rec || rec.rank === null || rec.rank === undefined) return null;
  const parts = [];
  if (rec.rank > 0) parts.push(rec.rank === 1 ? '\\mathbb{Z}' : `\\mathbb{Z}^{${rec.rank}}`);
  for (const n of rec.torsion || []) parts.push(`\\mathbb{Z}/${n}\\mathbb{Z}`);
  return T(`E(\\mathbb{Q}) \\cong ` + (parts.length ? parts.join(' \\oplus ') : '0'));
}
const mwUnknown = () => mixed('$E(\\mathbb{Q})$: not in Cremona’s tables');
function nextFrame() {                                         // let the busy marker paint first; a hidden tab gets no frames, so fall back to a timer
  return new Promise(r => { let done = false; const go = () => { if (!done) { done = true; setTimeout(r, 0); } }; requestAnimationFrame(go); setTimeout(go, 100); });
}

// a label is looked up in Cremona's tables, fetched shard by shard (js/cremona.js)
async function resolveLabel(p) {
  let r;
  try { r = await Cremona.lookup(p.label); }
  catch (e) { throw new Error(`the curve tables could not be loaded (${e.message}); serve the app over http, e.g. python3 -m http.server`); }
  if (r && r.missing) throw new Error(`conductor ${p.N} is beyond Cremona's tables (conductor ≤ ${Cremona.maxConductor()}); paste the a-invariants from https://www.lmfdb.org/EllipticCurve/Q/${p.label}/`);
  if (!r) throw new Error(`there is no curve ${p.label} in Cremona's tables`);
  return r;
}
async function plot(text) {
  text = text.trim(); if (!text) return;
  lastText = text; $('input').value = text;
  updateHash(text);
  $('busy').hidden = false; await nextFrame();
  try {
    const parsed = EC3D.parseInput(text);
    let desc = [];
    model = null; curveInfo = null;
    if (parsed.error) throw new Error(parsed.error);
    if (parsed.type === 'label' || parsed.type === 'ainvs') {
      let ainvs = parsed.ainvs, rec = null;
      if (parsed.type === 'label') { rec = await resolveLabel(parsed); ainvs = rec.ainvs.map(v => EC3D.Q.of(v)); }
      else rec = identify(ainvs);
      const m = EC3D.modelFromAinvs(ainvs);
      if (m.singular) throw new Error(`singular Weierstrass cubic ($\\Delta = 0$): $${EC3D.formatWeierstrassTeX(ainvs)}$`);
      model = m;
      desc.push(`<span class="eq">${T(EC3D.formatWeierstrassTeX(ainvs))}</span>`);
      if (rec) desc.push(esc(`${rec.cremona} = ${rec.lmfdb}, conductor ${rec.conductor}`));
      desc.push(mwText(rec) || mwUnknown());
    } else {
      const an = parsed.analysis;
      if (!an.model) {
        const g = an.genus === null ? mixed(`genus $\\le ${an.genusBound}$ (not computed)`) : `genus ${an.genus}`;
        clearSurface(); infoParts = null;
        setInfo(`<span class="eq">${T(EC3D.btex(an.F) + ' = 0')}</span> — <span class="err">not an elliptic curve I can plot</span>: ${esc(an.kind)}, ${g}. ${an.notes.map(mixed).join(' · ')}`, '');
        return;
      }
      model = an.model;
      desc.push(`<span class="eq">${T(EC3D.btex(an.F) + ' = 0')}</span>`);
      if (an.kind === 'weierstrass' && model.ainvsQ) {
        const rec = identify(model.ainvsQ);
        if (rec) desc.push(esc(`${rec.cremona} = ${rec.lmfdb}, conductor ${rec.conductor}`));
        desc.push(mwText(rec) || mwUnknown());
      } else desc.push(esc(an.kind));
      desc.push(...an.notes.map(mixed));
    }
    lattice = EC3D.periodLattice(model.ainvs, model.inv.disc instanceof EC3D.Q ? model.inv.disc.sign() : undefined, model.ainvsQ || undefined);
    computeGrid();
    setRadius(autoRadius(gridData));                          // per curve: the real points must be inside the ball
    rebuildSurface();
    drawLatticePlot();
    resetView();
    const inv = model.inv;
    desc.push(T(`\\Delta = ${qtex(inv.disc)},\\; j = ${inv.j === null ? '\\infty' : qtex(inv.j)}`));
    desc.push(T(`\\omega_1 = ${fmt(lattice.w1)},\\ \\omega_2 = ${ctex(lattice.w2)},\\ \\tau = ${ctex(lattice.normalised.tau)}`));
    desc.push(mixed(lattice.discSign > 0 ? '$\\Delta > 0$: two real components (rows $t = 0$ and $t = \\tfrac12$)' : '$\\Delta < 0$: one real component (row $t = 0$)'));
    desc.push(mixed(`clipped to $|(x, y)| < ${fmt(state.radius, 3)}$`));
    infoParts = { base: desc };
    renderInfo();
  } catch (e) {
    clearSurface(); infoParts = null; setInfo(`<span class="err">${mixed(e.message)}</span>`, '');
  } finally { $('busy').hidden = true; requestRender(); }
}

// ------------------------------------------------------------------ UI wiring
// The last three examples are not elliptic curves; they exercise the error paths and are only offered with ?dev in the URL.
const DEV = new URLSearchParams(location.search).has('dev') || (typeof window !== 'undefined' && !!window.EC3D_DEV);
const EXAMPLES = [
  ['20.a3  (two real components)', '20.a3'], ['11a1  (one real component)', '11a1'], ['37a1  (rank 1)', '37a1'], ['389a1', '389a1'], ['5077a1', '5077a1'],
  ['y² = x³ − x', 'y^2 = x^3 - x'], ['y² = x³ + x', 'y^2 = x^3 + x'], ['y² = x³ − x + 1', 'y^2 = x^3 - x + 1'],
  ['y² = x³ + 4  (E(ℚ) ≅ ℤ/3)', 'y^2 = x^3 + 4'], ['y² = x³ − 7x + 10', 'y^2 = x^3 - 7x + 10'], ['y² = x³ − x/4  (Arapura; try ∞)', 'y^2 = x^3 - x/4'],
  ['y² + y = x³ − x²', 'y^2 + y = x^3 - x^2'], ['2y² = x³ − x  (scaled)', '2*y^2 = x^3 - x'],
  ['y² + x²y = x³ + 1  (quartic model)', 'y^2 + x^2*y = x^3 + 1'], ['y² = x⁴ − 3x² + x + 1  (quartic)', 'y^2 = x^4 - 3x^2 + x + 1'],
  ['y² = x⁵ + 1  (genus 2)', 'y^2 = x^5 + 1', 'dev'], ['y² = x²(x + 1)  (nodal cubic)', 'y^2 = x^2*(x+1)', 'dev'], ['x³ + y³ = 1  (not Weierstrass)', 'x^3 + y^3 = 1', 'dev'],
];
for (const [label, value, dev] of EXAMPLES) {
  if (dev && !DEV) continue;
  const o = document.createElement('option'); o.value = value; o.textContent = label; $('examples').appendChild(o);
}
$('examples').addEventListener('change', e => { if (e.target.value) plot(e.target.value); e.target.value = ''; });
$('plot').addEventListener('click', () => plot($('input').value));
$('input').addEventListener('keydown', e => { if (e.key === 'Enter') plot($('input').value); });
// the third coordinate: Re y (θ = 0, a neighbourhood of the real points), Im y (θ = π/2, the imaginary slice), or the
// rotation z = Re y cos θ + Im y sin θ animated (Donu Arapura's picture); the real points are drawn only at θ = 0
const sliceButtons = [...document.querySelectorAll('#topbar .seg button')];
const thetaSlider = $('theta'), speedSlider = $('speed'), playButton = $('anim-play');
function setTheta(v) {
  v = ((v % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
  state.theta = v; thetaUniform.value = v;
  const deg = v * 180 / Math.PI;
  thetaSlider.value = deg.toFixed(1); $('v-theta').textContent = deg.toFixed(0) + '°';
  dirty = true;
}
function setPlaying(on) { state.playing = on; playButton.textContent = on ? '❚❚' : '▶'; playButton.setAttribute('aria-label', on ? 'pause' : 'play'); }
function setSlice(v, silent) {
  state.slice = v;
  for (const b of sliceButtons) { const on = b.dataset.slice === v; b.classList.toggle('active', on); b.setAttribute('aria-pressed', String(on)); }
  $('animbar').hidden = v !== 'anim';
  if (v === 're') { setTheta(0); setPlaying(false); }
  else if (v === 'im') { setTheta(Math.PI / 2); setPlaying(false); }
  else setPlaying(true);
  if (silent) return;
  if (gridData) { rebuildDecorations(); renderInfo(); }
  if (lastText) updateHash(lastText);
  requestRender();
}
for (const b of sliceButtons) b.addEventListener('click', () => { if (state.slice !== b.dataset.slice) setSlice(b.dataset.slice); });
playButton.addEventListener('click', () => setPlaying(!state.playing));
thetaSlider.addEventListener('input', () => { setPlaying(false); setTheta(+thetaSlider.value * Math.PI / 180); });
speedSlider.value = state.speed; $('v-speed').textContent = state.speed + '°/s';
speedSlider.addEventListener('input', () => { state.speed = +speedSlider.value; $('v-speed').textContent = state.speed + '°/s'; });
// the options drawer slides in from the left, below the top bar; its two tabs ride on its right edge and pick the section
const drawer = $('drawer'), tabs = [...drawer.querySelectorAll('.tab')];
let openSection = null;
function showSection(name) {
  openSection = name;
  drawer.classList.toggle('open', !!name);
  for (const t of tabs) { const on = t.dataset.section === name; t.classList.toggle('active', on); t.setAttribute('aria-expanded', String(on)); }
  if (name) for (const sec of drawer.querySelectorAll('#panel > section')) sec.hidden = sec.id !== 'section-' + name;
  if (name === 'color') scheduleLatticePlot();
}
for (const t of tabs) t.addEventListener('click', () => showSection(openSection === t.dataset.section ? null : t.dataset.section));
const topbar = $('topbar');                                    // the drawer sits below the top bar, whatever its height
function measureTopbar() { document.documentElement.style.setProperty('--topbar-h', topbar.offsetHeight + 'px'); }
new ResizeObserver(measureTopbar).observe(topbar); measureTopbar();
for (const name of Object.keys(COLORMAPS)) { const o = document.createElement('option'); o.value = o.textContent = name; $('colormap').appendChild(o); }
$('colormap').value = state.colormap;
$('colormap').addEventListener('change', e => { state.colormap = e.target.value; setColormap(state.colormap); scheduleLatticePlot(); });

function bindRange(id, key, show, onChange) {
  const el = $(id); el.value = state[key]; $('v-' + id).textContent = show(state[key]);
  el.addEventListener('input', () => { state[key] = Number(el.value); $('v-' + id).textContent = show(state[key]); onChange('input'); });
  el.addEventListener('change', () => { onChange('change'); });
}
bindRange('grid', 'grid', v => v, kind => { if (kind === 'change' && model) { $('busy').hidden = false; setTimeout(() => { computeGrid(); rebuildSurface(); $('busy').hidden = true; requestRender(); }, 20); } });
bindRange('cutoff', 'cutoff', v => v.toFixed(2), () => { uniforms.uCutoff.value = state.cutoff; requestRender(); });
bindRange('opacity', 'opacity', v => v.toFixed(2), () => { for (const m of [surfaceMaterial, mirrorMaterial]) { m.opacity = state.opacity; m.depthWrite = state.opacity >= 1; } requestRender(); });
bindRange('thick', 'thick', v => v.toFixed(3), kind => { if (kind === 'change') rebuildDecorations(); });
bindRange('lines', 'lines', v => v, () => { uniforms.uLines.value = state.lines; scheduleLatticePlot(); requestRender(); });
bindRange('soft', 'soft', v => v.toFixed(2), () => { uniforms.uSoft.value = state.soft; scheduleLatticePlot(); requestRender(); });
// the clipping radius has no upper bound: a logarithmic slider plus a free number field
const radiusSlider = $('radius'), radiusNum = $('radius-num');
function setRadius(R, from) {
  R = Math.max(0.1, +R || state.radius); state.radius = R;
  if (from !== 'slider') radiusSlider.value = Math.max(+radiusSlider.min, Math.min(+radiusSlider.max, Math.log10(R))).toFixed(2);
  if (from !== 'num') radiusNum.value = R >= 100 ? Math.round(R) : Number(R.toPrecision(3));
  uniforms.uRadius.value = R; $('cutoff').max = R; if (state.cutoff > R) { state.cutoff = 0; $('cutoff').value = 0; uniforms.uCutoff.value = 0; $('v-cutoff').textContent = '0.00'; }
  requestRender();
}
function radiusChanged() { if (!gridData) return; if (50 * state.radius > builtBig) rebuildSurface(); else rebuildDecorations(); }
radiusSlider.addEventListener('input', () => setRadius(Math.pow(10, +radiusSlider.value), 'slider'));
radiusSlider.addEventListener('change', radiusChanged);
radiusNum.addEventListener('change', () => { setRadius(radiusNum.value, 'num'); radiusChanged(); });
setRadius(state.radius);
function bindCheck(id, key, onChange) { const el = $(id); el.checked = state[key]; el.addEventListener('change', () => { state[key] = el.checked; onChange(); requestRender(); }); }
bindCheck('real', 'real', () => { if (realGroup) realGroup.visible = state.real && state.slice === 're'; });
bindCheck('mirror', 'mirror', () => { if (mirrorMesh) mirrorMesh.visible = state.mirror; });
bindCheck('axes', 'axes', () => { if (axesGroup) axesGroup.visible = state.axes; });
bindCheck('sphere', 'sphere', () => { if (sphereMesh) sphereMesh.visible = state.sphere; });
$('reset').addEventListener('click', resetView);
$('snapshot').addEventListener('click', () => { renderer.render(scene, camera); const a = document.createElement('a'); a.href = renderer.domElement.toDataURL('image/png'); a.download = 'elliptic-curve-3d.png'; document.body.appendChild(a); a.click(); a.remove(); });
$('share').addEventListener('click', async () => { try { await navigator.clipboard.writeText(location.href); $('share').textContent = 'Copied'; setTimeout(() => $('share').textContent = 'Copy link', 1200); } catch (e) { prompt('Link:', location.href); } });
setTimeout(() => { $('hint').hidden = true; }, 9000);

Cremona.ensure(11).catch(() => {});                              // warm the first shard: examples and curve identification
if (DEV) window.EC3D_DEBUG = { state, render: () => renderer.render(scene, camera), setTheta, canvas: renderer.domElement };   // for tests
let initial = decodeURIComponent((location.hash || '').slice(1));
if (initial.startsWith('im:')) { initial = initial.slice(3); setSlice('im', true); }
else if (initial.startsWith('anim:')) { initial = initial.slice(5); setSlice('anim', true); }
plot(initial || '20.a3');
})();
