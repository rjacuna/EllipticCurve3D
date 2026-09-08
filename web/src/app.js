/* app.js -- viewer and UI for the EllipticCurve3D web app.  Globals: THREE, EC3D, CURVE_TABLE, COLORMAPS. */
(function () {
'use strict';
const $ = id => document.getElementById(id);
const isMobile = matchMedia('(max-width: 640px)').matches || navigator.maxTouchPoints > 1;

const state = { grid: isMobile ? 250 : 400, radius: 3, cutoff: 0, opacity: 0.5, lines: 8, soft: 0.12, colormap: 'Greens_r',
                real: true, mirror: false, axes: true, sphere: false };
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
(function loop() { requestAnimationFrame(loop); const moved = controls.update(); if (moved || dirty) { renderer.render(scene, camera); dirty = false; } })();
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

// ------------------------------------------------------------------ surface material (lattice colouring + spherical clipping in the shader)
const uniforms = { uRadius: { value: state.radius }, uCutoff: { value: state.cutoff }, uLines: { value: state.lines }, uSoft: { value: state.soft }, uColormap: { value: cmapTex } };
function makeSurfaceMaterial() {
  const m = new THREE.MeshPhongMaterial({ side: THREE.DoubleSide, transparent: true, opacity: state.opacity, depthWrite: state.opacity >= 1,
                                          shininess: 25, specular: new THREE.Color(0x2a2a2a), color: 0xffffff });
  m.onBeforeCompile = shader => {
    Object.assign(shader.uniforms, uniforms);
    shader.vertexShader = shader.vertexShader
      .replace('#include <common>', '#include <common>\nattribute vec2 st;\nvarying vec2 vST;\nvarying vec3 vPos;')
      .replace('#include <begin_vertex>', '#include <begin_vertex>\nvST = st;\nvPos = position;');
    shader.fragmentShader = shader.fragmentShader
      .replace('#include <common>', '#include <common>\nuniform float uRadius;\nuniform float uCutoff;\nuniform float uLines;\nuniform float uSoft;\nuniform sampler2D uColormap;\nvarying vec2 vST;\nvarying vec3 vPos;')
      .replace('#include <clipping_planes_fragment>', '#include <clipping_planes_fragment>\n{ float rr = length(vPos); if (rr > uRadius || rr < uCutoff) discard; }')
      .replace('#include <color_fragment>', '#include <color_fragment>\n{ float c = pow(abs(sin(3.14159265358979 * vST.x * uLines) * sin(3.14159265358979 * vST.y * uLines)), uSoft);\n  diffuseColor.rgb = texture2D(uColormap, vec2(c, 0.5)).rgb; }');
  };
  m.customProgramCacheKey = () => 'ec3d-surface';
  return m;
}
const surfaceMaterial = makeSurfaceMaterial();
const realMaterial = new THREE.MeshPhongMaterial({ color: 0xc8281e, shininess: 40, specular: new THREE.Color(0x442222) });

// ------------------------------------------------------------------ scene objects
const group = new THREE.Group(); scene.add(group);
let surfaceMesh = null, mirrorMesh = null, realGroup = null, axesGroup = null, sphereMesh = null;

function buildSurfaceGeometry(g) {
  const n = g.n, idx = [];
  const big = 2000;                                                         // drop cells at the pole itself (|P| enormous)
  const ok = i => g.ok[i] && Math.abs(g.pos[3 * i]) < big && Math.abs(g.pos[3 * i + 1]) < big && Math.abs(g.pos[3 * i + 2]) < big;
  for (let i = 0; i < n - 1; i++) for (let j = 0; j < n - 1; j++) {
    const a = i * n + j, b = (i + 1) * n + j, c = (i + 1) * n + j + 1, d = i * n + j + 1;
    if (ok(a) && ok(b) && ok(c) && ok(d)) idx.push(a, b, c, a, c, d);
  }
  const geom = new THREE.BufferGeometry();
  geom.setAttribute('position', new THREE.BufferAttribute(g.pos, 3));
  geom.setAttribute('st', new THREE.BufferAttribute(g.uv, 2));
  geom.setIndex(new THREE.BufferAttribute(new Uint32Array(idx), 1));
  geom.computeVertexNormals();
  return geom;
}
function buildRealCurves(comps) {
  const grp = new THREE.Group();
  for (const comp of comps) {
    const pts = comp.map(p => new THREE.Vector3(p[0], p[1], p[2]));
    const closed = pts.length > 3 && pts[0].distanceTo(pts[pts.length - 1]) < 1e-9;
    if (closed) pts.pop();
    if (pts.length < 2) continue;
    const curve = new THREE.CatmullRomCurve3(pts, closed, 'centripetal');
    const tube = new THREE.TubeGeometry(curve, Math.min(800, 2 * pts.length), 0.011 * state.radius, 8, closed);
    grp.add(new THREE.Mesh(tube, realMaterial));
  }
  return grp;
}
function makeLabel(text, x, y, z) {
  const canvas = document.createElement('canvas'), ctx = canvas.getContext('2d'), pr = 2, fs = 22;
  ctx.font = `${fs}px system-ui, sans-serif`;
  const w = Math.ceil(ctx.measureText(text).width) + 8, h = fs + 8;
  canvas.width = w * pr; canvas.height = h * pr; ctx.scale(pr, pr);
  ctx.font = `${fs}px system-ui, sans-serif`; ctx.fillStyle = '#222'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText(text, w / 2, h / 2);
  const tex = new THREE.Texture(canvas); tex.needsUpdate = true;
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, sizeAttenuation: false, depthWrite: false, transparent: true }));
  sp.position.set(x, y, z); sp.scale.set(w / 700, h / 700, 1);
  return sp;
}
function buildAxes(R) {
  const grp = new THREE.Group();
  grp.add(new THREE.AxesHelper(1.15 * R));
  grp.add(makeLabel('Re x', 1.3 * R, 0, 0)); grp.add(makeLabel('Im x', 0, 1.3 * R, 0)); grp.add(makeLabel('Re y', 0, 0, 1.3 * R));
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
    realGroup = buildRealCurves(EC3D.realComponents(gridData, R)); realGroup.visible = state.real; group.add(realGroup);
  }
  requestRender();
}
function rebuildSurface() {
  if (!model) return;
  const t0 = performance.now();
  gridData = EC3D.buildGrid(model, state.grid, lattice);
  const t1 = performance.now();
  const geom = buildSurfaceGeometry(gridData);
  if (surfaceMesh) { group.remove(surfaceMesh); surfaceMesh.geometry.dispose(); group.remove(mirrorMesh); }
  surfaceMesh = new THREE.Mesh(geom, surfaceMaterial); group.add(surfaceMesh);
  mirrorMesh = new THREE.Mesh(geom, surfaceMaterial); mirrorMesh.scale.set(1, -1, 1); mirrorMesh.visible = state.mirror; group.add(mirrorMesh);
  rebuildDecorations();
  const t2 = performance.now();
  return { wp: t1 - t0, mesh: t2 - t1, faces: geom.index.count / 3 };
}
function clearSurface() {
  if (surfaceMesh) { group.remove(surfaceMesh); group.remove(mirrorMesh); surfaceMesh = mirrorMesh = null; }
  if (realGroup) { group.remove(realGroup); realGroup = null; }
  gridData = null; requestRender();
}

// ------------------------------------------------------------------ input handling
const fmt = (v, d = 5) => (Math.abs(v) < 1e-12 ? '0' : Number(v.toPrecision(d)).toString());
const cfmt = z => { const re = fmt(z[0]), im = fmt(Math.abs(z[1])); if (Math.abs(z[1]) < 1e-12) return re; if (Math.abs(z[0]) < 1e-12) return (z[1] < 0 ? '-' : '') + im + 'i'; return `${re} ${z[1] < 0 ? '-' : '+'} ${im}i`; };
function setInfo(html, cls) { const el = $('info'); el.innerHTML = html; el.className = cls || ''; }
function nextFrame() { return new Promise(r => requestAnimationFrame(() => setTimeout(r, 0))); }

async function plot(text) {
  text = text.trim(); if (!text) return;
  lastText = text; $('input').value = text;
  try { history.replaceState(null, '', '#' + encodeURIComponent(text)); } catch (e) {}
  $('busy').hidden = false; await nextFrame();
  try {
    const parsed = EC3D.parseInput(text, CURVE_TABLE);
    let desc = [];
    model = null; curveInfo = null;
    if (parsed.error) throw new Error(parsed.error);
    if (parsed.type === 'label' || parsed.type === 'ainvs') {
      const m = EC3D.modelFromAinvs(parsed.ainvs);
      if (m.singular) throw new Error(`singular Weierstrass cubic (Δ = 0): ${EC3D.formatWeierstrass(parsed.ainvs)}`);
      model = m;
      const rec = parsed.rec || EC3D.findByAinvs(parsed.ainvs.map(v => Number(v.toString())), CURVE_TABLE);
      desc.push(`<b>${EC3D.formatWeierstrass(parsed.ainvs)}</b>`);
      if (rec) desc.push(`${rec.cremona} = ${rec.lmfdb}, conductor ${rec.conductor}`);
    } else {
      const an = parsed.analysis;
      if (!an.model) {
        const g = an.genus === null ? `genus ≤ ${an.genusBound} (not computed)` : `genus ${an.genus}`;
        clearSurface();
        setInfo(`<b>${an.pretty}</b> — <span class="err">not an elliptic curve I can plot</span>: ${an.kind}, ${g}. ${an.notes.join(' · ')}`, '');
        return;
      }
      model = an.model;
      desc.push(`<b>${an.pretty.replace(/ = 0$/, '')} = 0</b>`);
      if (an.kind === 'weierstrass' && model.ainvsQ) {
        const rec = EC3D.findByAinvs(model.ainvsQ.map(v => Number(v.toString())), CURVE_TABLE);
        if (rec) desc.push(`${rec.cremona} = ${rec.lmfdb}, conductor ${rec.conductor}`);
      } else desc.push(an.kind);
      desc.push(...an.notes);
    }
    lattice = EC3D.periodLattice(model.ainvs, model.inv.disc instanceof EC3D.Q ? model.inv.disc.sign() : undefined, model.ainvsQ || undefined);
    const timing = rebuildSurface();
    resetView();
    const inv = model.inv;
    const discStr = inv.disc instanceof EC3D.Q ? inv.disc.toString() : fmt(inv.disc);
    const jStr = inv.j === null ? '∞' : (inv.j instanceof EC3D.Q ? inv.j.toString() : fmt(inv.j));
    desc.push(`Δ = ${discStr}, j = ${jStr}`);
    desc.push(`ω₁ = ${fmt(lattice.w1)}, ω₂ = ${cfmt(lattice.w2)}, τ = ${cfmt(lattice.normalised.tau)}`);
    desc.push(lattice.discSign > 0 ? 'Δ > 0: two real components (rows t = 0 and t = ½)' : 'Δ < 0: one real component (row t = 0)');
    desc.push(`<span class="ok">${timing.faces.toLocaleString()} triangles in ${(timing.wp + timing.mesh).toFixed(0)} ms</span>`);
    setInfo(desc.join(' · '));
  } catch (e) {
    clearSurface(); setInfo(`<span class="err">${e.message}</span>`, '');
  } finally { $('busy').hidden = true; requestRender(); }
}

// ------------------------------------------------------------------ UI wiring
// The last three examples are not elliptic curves; they exercise the error paths and are only offered when the
// page was built with  build.py --dev  (which sets window.EC3D_DEV).
const DEV = typeof window !== 'undefined' && !!window.EC3D_DEV;
const EXAMPLES = [
  ['20.a3  (two real components)', '20.a3'], ['11a1  (one real component)', '11a1'], ['37a1  (rank 1)', '37a1'], ['389a1', '389a1'], ['5077a1', '5077a1'],
  ['y² = x³ − x', 'y^2 = x^3 - x'], ['y² + y = x³ − x²', 'y^2 + y = x^3 - x^2'], ['2y² = x³ − x  (scaled)', '2*y^2 = x^3 - x'],
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
// the options drawer slides in from the left; its tab rides on its right edge
const drawer = $('drawer'), optionsTab = $('options-tab');
function setOptionsOpen(open) { drawer.classList.toggle('open', open); optionsTab.setAttribute('aria-expanded', String(open)); }
optionsTab.addEventListener('click', () => setOptionsOpen(!drawer.classList.contains('open')));
for (const name of Object.keys(COLORMAPS)) { const o = document.createElement('option'); o.value = o.textContent = name; $('colormap').appendChild(o); }
$('colormap').value = state.colormap;
$('colormap').addEventListener('change', e => { state.colormap = e.target.value; setColormap(state.colormap); });

function bindRange(id, key, show, onChange) {
  const el = $(id); el.value = state[key]; $('v-' + id).textContent = show(state[key]);
  el.addEventListener('input', () => { state[key] = Number(el.value); $('v-' + id).textContent = show(state[key]); onChange('input'); });
  el.addEventListener('change', () => { onChange('change'); });
}
bindRange('grid', 'grid', v => v, kind => { if (kind === 'change' && model) { $('busy').hidden = false; setTimeout(() => { rebuildSurface(); $('busy').hidden = true; requestRender(); }, 20); } });
bindRange('radius', 'radius', v => v.toFixed(1), kind => { uniforms.uRadius.value = state.radius; $('cutoff').max = state.radius; if (kind === 'change') rebuildDecorations(); requestRender(); });
bindRange('cutoff', 'cutoff', v => v.toFixed(2), () => { uniforms.uCutoff.value = state.cutoff; requestRender(); });
bindRange('opacity', 'opacity', v => v.toFixed(2), () => { surfaceMaterial.opacity = state.opacity; surfaceMaterial.depthWrite = state.opacity >= 1; requestRender(); });
bindRange('lines', 'lines', v => v, () => { uniforms.uLines.value = state.lines; requestRender(); });
bindRange('soft', 'soft', v => v.toFixed(2), () => { uniforms.uSoft.value = state.soft; requestRender(); });
function bindCheck(id, key, onChange) { const el = $(id); el.checked = state[key]; el.addEventListener('change', () => { state[key] = el.checked; onChange(); requestRender(); }); }
bindCheck('real', 'real', () => { if (realGroup) realGroup.visible = state.real; });
bindCheck('mirror', 'mirror', () => { if (mirrorMesh) mirrorMesh.visible = state.mirror; });
bindCheck('axes', 'axes', () => { if (axesGroup) axesGroup.visible = state.axes; });
bindCheck('sphere', 'sphere', () => { if (sphereMesh) sphereMesh.visible = state.sphere; });
$('reset').addEventListener('click', resetView);
$('snapshot').addEventListener('click', () => { renderer.render(scene, camera); const a = document.createElement('a'); a.href = renderer.domElement.toDataURL('image/png'); a.download = 'elliptic-curve-3d.png'; document.body.appendChild(a); a.click(); a.remove(); });
$('share').addEventListener('click', async () => { try { await navigator.clipboard.writeText(location.href); $('share').textContent = 'Copied'; setTimeout(() => $('share').textContent = 'Copy link', 1200); } catch (e) { prompt('Link:', location.href); } });
setTimeout(() => { $('hint').hidden = true; }, 9000);

const initial = decodeURIComponent((location.hash || '').slice(1)) || '20.a3';
plot(initial);
})();
