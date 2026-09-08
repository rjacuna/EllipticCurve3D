# verify.sage -- checks the notebook's own code.  Run from the project directory:   sage verify.sage
#
# It loads elliptic_curve_3d.sage (the notebook uses the same file) and checks that
#   * every point of the parametrisation satisfies the Weierstrass equation of the given model,
#   * the rows t = 0 (and t = 1/2 when disc > 0) are the real components of E,
#   * the corners of the grid are the 2-torsion points (roots of the 2-division polynomial),
#   * the surface builds, keeps its opacity through the clipping, and has no non-finite vertices
#     even when the pole z = 0 is a grid point.
import json, sys, numpy as np
load("elliptic_curve_3d.sage")

fails = 0
def check(name, ok, detail=""):
    global fails
    fails += (not ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")

for label in ["20.a3", "11a1", "37a1", "389a1", "5077a1"]:
    E = EllipticCurve(label)
    a1, a2, a3, a4, a6 = (float(a) for a in E.ainvs())
    S, T, x, y = parametrize(E, plot_points=301)             # odd: the pole is a grid point
    res = y**2 + a1*x*y + a3*y - (x**3 + a2*x**2 + a4*x + a6)
    m = np.isfinite(res) & (np.abs(x) < 50)
    check(f"{label}: grid points lie on E", np.abs(res[m]).max() < 1e-8, f"(max residual {np.abs(res[m]).max():.1e})")
    mask = np.abs(x[:, 0]) < 50
    check(f"{label}: row t=0 is real", np.abs(x[:, 0].imag[mask]).max() < 1e-8 and np.abs(y[:, 0].imag[mask]).max() < 1e-7)
    roots = sorted(float(r) for r in E.two_division_polynomial().roots(RR, multiplicities=False))
    if E.discriminant() > 0:
        check(f"{label}: row t=1/2 is real (second component)", np.abs(x[:, -1].imag).max() < 1e-8 and np.abs(y[:, -1].imag).max() < 1e-8)
        mid = len(S)//2                                             # s = 0 (plot_points is odd)
        corners = sorted([x[mid, -1].real, x[-1, -1].real, x[-1, 0].real])   # z = w2/2, (w1+w2)/2, w1/2
    else:
        corners = [x[-1, 0].real]
    check(f"{label}: half-periods map to the 2-torsion points", np.allclose(corners, roots, atol=1e-8), f"{[round(c, 6) for c in corners]}")

G = curve("20.a3", colormaps.Greens_r, spacing=1/8, cutoff=1/2, radius=3, opacity=1/2, plot_points=151)
V = np.array(G.vertex_list()); r = np.sqrt((V**2).sum(axis=1))
check("surface builds with the pole on the grid, all vertices finite", np.isfinite(V).all())
check("clipping: cutoff <= |P| <= radius", r.min() > 0.5 - 1e-6 and r.max() < 3 + 1e-6, f"(|P| in [{r.min():.4f}, {r.max():.4f}])")
surf = G.threejs_repr(G.default_render_params())[0][1]
check("opacity survives add_condition (three.js)", surf["opacity"] == 0.5, f"(opacity {surf['opacity']})")
check("per-face lattice colors present", len(surf.get("faceColors", [])) == len(G.face_list()))

print("\n" + ("ALL CHECKS PASSED" if not fails else f"{fails} CHECK(S) FAILED"))
sys.exit(fails)
