# EllipticCurve3D

The complex points of an elliptic curve `E` over ℝ form a torus `E(ℂ) ⊂ ℂ² = ℝ⁴`. This project draws its
projection `(Re x, Im x, Re y)` to ℝ³, colored by the period lattice, so that the real points `E(ℝ)` are the
slice `Im x = 0` and the colored lines are the cycles `s·ω₁ + t·ω₂` with `s` or `t` constant. It exists in
two forms that share the same mathematics: a SageMath module with a notebook, and a single-page web app.

## Layout

| path | what it is |
|---|---|
| `elliptic_curve_3d.sage` | the Sage module: `EllipticCurve3D`, with doctests. `load("elliptic_curve_3d.sage")` |
| `Elliptic Curve 3D.ipynb` | the notebook (SageMath 10.9 kernel): builds the surface with the module and runs the two checks below |
| `run_doctests.py` | `sage -python run_doctests.py` runs the module's doctests |
| `verify.sage` | `sage verify.sage`: algebraic checks of the parametrisation on several curves |
| `slice_check.sage` | the planar diagnosis of the original bug, writes `figures/slice_*.png` |
| `render_figures.sage` | tachyon renders of the fixed and the original surface, writes `figures/render_*.png` |
| `original_curve.sage` | cell 0 of the original notebook, verbatim, for reference |
| `web/` | the web app, see below |

## The mathematics, and the bug that started this

Sage's `E.weierstrass_p()` is the Weierstrass function of the period lattice of `E`, and `(℘, ℘′)` satisfy
`℘′² = 4℘³ − g₂℘ − g₃` with `g₂ = c₄/12`, `g₃ = c₆/216`. That is the c₄/c₆ model, not the Weierstrass equation
`y² + a₁xy + a₃y = x³ + a₂x² + a₄x + a₆` one starts from. The two are related by

    ℘ = x + b₂/12,      ℘′ = 2y + a₁x + a₃,      b₂ = a₁² + 4a₂,

so plotting `(Re ℘, Im ℘, Re ℘′)` gives `E` shifted by `b₂/12` in `x` and stretched by 2 in `y`, and the shift
changes from curve to curve. For `20.a3`, `y² = x³ + x² − x`, the shift is `1/3`. The uniformisation used now is

    x = ℘(z) − b₂/12,      y = (℘′(z) − a₁x − a₃)/2,

the inverse of the Abel–Jacobi map of the given model. A second, independent problem was numerical: the Laurent
series of `℘` at `z = 0`, truncated at order 10, is badly wrong at the half-periods, which are exactly the real
points. `℘` and `℘′` are now evaluated on the whole grid by the q-series

    ℘(z) = (2πi)² [ 1/12 + Σₙ qⁿu/(1 − qⁿu)² − 2 Σₙ≥₁ qⁿ/(1 − qⁿ)² ],   u = e^{2πiz}, q = e^{2πiτ},

in a reduced basis of the lattice. Reducing the basis matters: Sage's basis has `ω₁` real but not short, and
for `1001a1` its `τ` has imaginary part 0.13, so `|q| = 0.45` and ten terms would give four digits. In the
reduced basis `|q| ≤ e^{−π√3} ≈ 0.0043` for every lattice and ten terms give double precision.

The surface is drawn over the half fundamental domain `z = s·ω₁ + t·ω₂`, `−½ ≤ s ≤ ½`, `0 ≤ t ≤ ½`, a
fundamental domain for `z ↦ −z`, so every `x` occurs once. It is an annulus whose boundary circles are the rows
`t = 0` and `t = ½`, the two real components when `Δ > 0`. The lines `t = k/n` are complete cycles in the class
of `ω₁`; the lines `s = k/n` are half-cycles from one real component to the other, and two of them, `s = 0` and
`s = ½`, are squashed onto the `x`-axis because there `y` is purely imaginary and the projection drops `Im y`.
The other half of the torus is the mirror image under `Im x ↦ −Im x`.

`figures/slice_bug.png` shows the original, shifted slice against the real points, `figures/slice_truncation.png`
the truncated series against exact `℘`, `figures/slice_fixed.png` the fixed slice, and `figures/render_*.png`
the 3D pictures with the real points drawn in red from the equation alone.

## The Sage module

```python
load("elliptic_curve_3d.sage")
S = EllipticCurve3D("20.a3")            # a curve, a label, a-invariants, a Weierstrass polynomial, or a cubic f(x)
S.plot(colormaps.Greens_r, radius=3, opacity=1/2, plot_points=600)      # surface + real points
S.point(z), S.grid(n), S.real_components(radius), S.basis(), S.normalised_basis()
```

`point(z)` returns the point of the given model at `z`; `grid(n)` the parametrisation on an `n × n` grid;
`surface(...)` the colored, clipped surface with `cutoff ≤ |(X, Y, Z)| < radius`, its opacity restored after
the clipping (Sage's `add_condition` drops it for per-face colors). Curves over a number field are taken with
a real embedding, the first one by default or the `embedding` argument.

Checks: `sage -python run_doctests.py` (81 doctests) and `sage verify.sage`.

## The web app

`web/dist/index.html` is the whole app in one file: open it directly or host it anywhere. It takes an LMFDB
label (`20.a3`), a Cremona label (`11a1`), a-invariants (`[0,1,0,-1,0]`), or an equation, draws the surface with
three.js, and reports the equation, labels, conductor, `Δ`, `j`, `ω₁`, `ω₂`, `τ` and the number of real components.

* Labels resolve from an offline table of every curve of conductor below 10000 from Sage's copy of Cremona's
  tables, with both label systems (`web/src/curves-data.js`, 64,687 curves). lmfdb.org's API sends no CORS
  header, so a browser page cannot query it; for a larger conductor paste the a-invariants from the LMFDB page.
* Equations are parsed with exact rational arithmetic. Weierstrass equations are used as given; a scaled one
  such as `2y² = x³ − x` is rescaled to a monic model and the surface drawn in the original coordinates; a
  quartic `y² = f₄(x)` with a real root, or a cubic with an `x²y` term, is brought to Weierstrass form and drawn
  through the corresponding model; anything else is reported with its genus.
* Next to Plot, **Re y / Im y** chooses the third coordinate: `(Re x, Im x, Re y)`, a neighbourhood of the real
  points, or the imaginary slice `(Re x, Im x, Im y)`, in which the real points are not drawn (there they lie in the
  plane `Im y = 0`, not on a curve of the picture). The mirror half is the complex conjugate, so in the imaginary
  slice it is reflected in both `Im x` and `Im y`. A link to the imaginary slice reads `#im:20.a3`.
* Controls live in a drawer that slides in from the left, below the top bar, with two tabs on its edge.
  **Surface**: grid size, clipping radius and inner cutoff (done in the shader), opacity, the real points and their
  thickness, the mirror half (on by default, so the whole torus is shown), axes, the clipping sphere, PNG export
  and a shareable link. **Coloring**: a picture of the z-plane colored the same way as the surface, with the period
  parallelogram outlined and the drawn half dashed, the subdivision count n (dark lines where s or t is a multiple
  of 1/n), the line softness and the colormap.
* The examples menu lists elliptic curves only. `python3 web/build.py --dev` also adds three inputs that are not
  elliptic curves (genus 2, a nodal cubic, a cubic not in Weierstrass form), which exercise the error messages.

Development layout, assembled by `python3 web/build.py` (add `--dev` for the extra examples):

| path | |
|---|---|
| `web/src/ec3d-math.js` | parsing, invariants, genus, period lattice by the AGM (a port of Sage's `_compute_periods_real` and `normalise_periods`), `℘` by the q-series, the grid, real components |
| `web/src/app.js`, `app.css`, `body.html`, `index.template.html` | viewer and UI |
| `web/src/curves-data.js`, `colormaps.js` | the label table and 64-sample matplotlib colormaps |
| `web/vendor/` | three.js r128 and OrbitControls, vendored so the page works offline |
| `web/index.html` | development page loading the files above separately; `python3 -m http.server -d web` |
| `web/dist/index.html`, `dist/artifact.html` | the single-file build, and the same without the document wrapper |
| `web/test/test-math.mjs` | `node web/test/test-math.mjs`: 8,238 checks of the JavaScript against Sage/PARI values in `vectors.json` (made by `make_vectors.sage`) |

Rebuild after editing anything in `web/src`, and rerun the test after touching `ec3d-math.js`.
