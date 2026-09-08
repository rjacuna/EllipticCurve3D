# render_figures.sage -- static (tachyon) renders of the surface together with the real points of E,
# drawn independently of ℘ in the plane Im x = 0.  Run from the project directory:  sage render_figures.sage
#
# Note: tachyon builds per-face textures from the colour only, so lattice-coloured surfaces render
# opaque here; the three.js viewer in the notebook honours the opacity.
import os, json
FIG = "figures"; os.makedirs(FIG, exist_ok=True)

E = EllipticCurve("20.a3")
_, a2, _, a4, a6 = E.ainvs()
x = var('x'); f = x^3 + a2*x^2 + a4*x + a6
e3, e2, e1 = sorted(E.two_division_polynomial().roots(RR, multiplicities=False))

def real_points(xmax=1.8):
    """the real locus of E as (x, 0, ±sqrt f(x)) in (Re x, Im x, Re y)-space -- no ℘ involved"""
    P = None
    for lo, hi in [(e1, xmax), (e3, e2)]:
        for s in (1, -1):
            c = parametric_plot3d((x, 0, s*sqrt(abs(f))), (x, lo, hi), color='red', thickness=7)
            P = c if P is None else P + c
    return P

VIEW = dict(viewer='tachyon', camera_position=(1.0, -2.4, 1.9), zoom=1.2, aspect_ratio=1, figsize=9)

# 1) the fixed surface: cell 0 of the notebook
load("elliptic_curve_3d.sage")
G = curve("20.a3", colormaps.Greens_r, spacing=1/8, cutoff=0, radius=3, opacity=1/2, plot_points=300)
(G + real_points()).save(os.path.join(FIG, "render_fixed.png"), **VIEW)
print("wrote figures/render_fixed.png")

# 2) the original construction (Re ℘, Im ℘, Re ℘') with the truncated Laurent series, for comparison
ns = globals().copy()
exec(preparse(open("original_curve.sage").read()), ns)
G0 = ns['curve']("20.a3", colormaps.Greens_r, spacing=1/8, cutoff=0, radius=3, opacity=1/2, plot_points=250)
(G0 + real_points()).save(os.path.join(FIG, "render_original.png"), **VIEW)
print("wrote figures/render_original.png")
