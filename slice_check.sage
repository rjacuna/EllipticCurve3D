# slice_check.sage -- planar diagnosis of the "surface is off from the real points" problem.
#
# Run from the project directory:   sage slice_check.sage      (writes figures/slice_*.png)
#
# The surface drawn by the notebook is (Re X, Im X, Re Z) for a parametrisation
# z -> (X(z), Z(z)) of the curve.  Its slice by the plane Im X = 0 is therefore a
# 1-parameter curve u -> (X(u), Z(u)) with u running over the lines where X is real
# (the real axis and, for a rectangular lattice, the line Im z = Im(w2)/2).  We take
# that slice algebraically and compare it with the implicit plot of the real points.
import os, numpy as np
FIG = "figures"; os.makedirs(FIG, exist_ok=True)

E = EllipticCurve("20.a3")                       # y^2 = x^3 + x^2 - x
a1, a2, a3, a4, a6 = E.ainvs()
b2, g2, g3 = E.b2(), E.c4()/12, E.c6()/216
L = E.period_lattice()
w1, w2 = L.basis()                               # w1 real, w2 imaginary (disc > 0)
print("E =", E, "\n b2/12 =", b2/12, "  g2 =", g2, " g3 =", g3, "\n w1 =", w1, " w2 =", w2)

R = PolynomialRing(RR, 'X'); X_ = R.gen()
f_roots = sorted((X_**3 + a2*X_**2 + a4*X_ + a6).roots(multiplicities=False))       # roots of f: 2-torsion of E
e_roots = sorted((4*X_**3 - g2*X_ - g3).roots(multiplicities=False))               # e_i = roots of 4X^3 - g2 X - g3
print(" roots of f      :", [round(r, 5) for r in f_roots])
print(" e_i             :", [round(r, 5) for r in e_roots], "  = roots of f + b2/12")

def wp(z):
    """exact (wp(z), wp'(z)) for the lattice of E, via PARI"""
    P, dP = L.elliptic_exponential(CC(z), to_curve=False)
    return complex(P), complex(dP)

def algebraic_curve(g, lo_hi_list, n=600):
    """real locus of Z^2 = g(X) as a list of polylines, sampling X on the given intervals"""
    polys = []
    for lo, hi in lo_hi_list:
        xs = np.linspace(float(lo), float(hi), n)
        ys = np.sqrt(np.maximum(np.array([float(g(v)) for v in xs]), 0))
        polys.append(list(zip(xs, ys)) + list(zip(xs[::-1], -ys[::-1])))
    return polys

def slice_points(zs, transform):
    """the 1-parameter curve u -> (Re X, Re Z) traced by the surface on the plane Im X = 0"""
    pts = []
    for z in zs:
        P, dP = wp(z)
        X, Z = transform(P, dP)
        if abs(X) < 12 and abs(Z) < 30:
            pts.append((X.real, Z.real))
    return pts

W1, IW2 = float(w1), float(w2.imag())
us = np.linspace(0.05, W1 - 0.05, 900)                        # real axis (avoiding the pole)
vs = np.linspace(0, W1, 900)                                  # the line Im z = Im(w2)/2
real_axis  = [complex(u, 0) for u in us]
oval_line  = [complex(v, IW2/2) for v in vs]
imag_axis  = [complex(0, t) for t in np.linspace(0.05, IW2/2, 300)]        # X real, Z imaginary: lands on the X-axis
half_line  = [complex(W1/2, t) for t in np.linspace(0, IW2/2, 300)]

xr, yr, XMAX = (-2.7, 3.7), (-6.5, 6.5), 3.7
fX = lambda v: v**3 + a2*v**2 + a4*v + a6
gX = lambda v: 4*v**3 - g2*v - g3
E_real  = algebraic_curve(fX, [(f_roots[2], XMAX), (f_roots[0], f_roots[1])])
P_curve = algebraic_curve(gX, [(e_roots[2], XMAX), (e_roots[0], e_roots[1])])

def save(g, name, title):
    g.set_legend_options(loc='lower right', font_size=9)
    g.axes_labels(['Re x', 'Re y'])
    g.save(os.path.join(FIG, name), figsize=[9, 6.5], dpi=110, title=title,
           xmin=xr[0], xmax=xr[1], ymin=yr[0], ymax=yr[1], title_pos=(0.5, 1.05))
    print("wrote", os.path.join(FIG, name))

# ---- Figure 1: the coordinate mismatch (exact wp, so this is purely algebraic) ----
g = Graphics()
for i, poly in enumerate(E_real):
    g += line2d(poly, color='blue', thickness=2.5, legend_label=None if i else 'real points of E:  y² = x³ + x² − x')
for i, poly in enumerate(P_curve):
    g += line2d(poly, color='darkorange', thickness=2, linestyle='--', legend_label=None if i else 'η² = 4X³ − g₂X − g₃   (what (℘, ℘′) satisfies)')
identity = lambda P, dP: (P, dP)
g += point2d(slice_points(real_axis, identity)[::6] + slice_points(oval_line, identity)[::6], color='red', size=9,
             legend_label='Im x = 0 slice of the original surface:  u ↦ (℘(u), ℘′(u))')
g += line2d([(p[0], 0) for p in slice_points(imag_axis, identity)], color='gray', thickness=1.5)
g += line2d([(p[0], 0) for p in slice_points(half_line, identity)], color='gray', thickness=1.5,
            legend_label='slice where ℘′ is imaginary (surface meets the x-axis)')
g += point2d([(r, 0) for r in f_roots], color='blue', size=45, zorder=5, legend_label='2-torsion of E:  roots of f')
g += point2d([(r, 0) for r in e_roots], color='red', size=45, zorder=5, legend_label='roots of 4X³ − g₂X − g₃  =  roots of f + b₂/12')
for r, e in zip(f_roots, e_roots):
    g += arrow2d((r, 0.9), (e, 0.9), color='black', width=1, arrowsize=2)
g += text('shift b₂/12 = 1/3;  y-scale ×2', (1.6, 1.4), color='black', fontsize=10)
save(g, 'slice_bug.png', 'Original surface (Re ℘, Im ℘, Re ℘′): the Im x = 0 slice is E shifted by b₂/12 and stretched ×2')

# ---- Figure 2: the truncated Laurent series the notebook actually evaluated ----
p10 = E.weierstrass_p(prec=10).truncate(10); pp10 = p10.derivative()
def series_slice(zs):
    pts = []
    for z in zs:
        P, dP = complex(p10(CC(z))), complex(pp10(CC(z)))
        if abs(P) < 12 and abs(dP) < 30:
            pts.append((P.real, dP.real))
    return pts
g = Graphics()
for i, poly in enumerate(P_curve):
    g += line2d(poly, color='darkorange', thickness=2, legend_label=None if i else 'η² = 4X³ − g₂X − g₃   (exact ℘: the slice lies on this)')
g += line2d(series_slice(real_axis), color='red', thickness=2, linestyle='--', legend_label='same slice with ℘ ≈ Laurent series truncated at z¹⁰ (the notebook)')
g += line2d(series_slice(oval_line), color='red', thickness=2, linestyle='--')
g += point2d([(r, 0) for r in e_roots], color='darkorange', size=45, zorder=5)
zc = complex(W1/2, IW2/2)
g += point2d([(wp(zc)[0].real, 0)], color='black', size=60, zorder=6, legend_label='℘(ω₁/2 + ω₂/2) = e₂ = 1/3 exactly')
g += point2d([(complex(p10(CC(zc))).real, complex(pp10(CC(zc))).real)], color='red', marker='x', size=80, zorder=6,
             legend_label='series value there: %.3f (error %.2f)' % (complex(p10(CC(zc))).real, abs(complex(p10(CC(zc))) - wp(zc)[0])))
save(g, 'slice_truncation.png', 'Second problem: the order-10 Laurent series about z = 0 is wrong near the half-periods, i.e. at the real points')

# ---- Figure 3: the fix ----
fix = lambda P, dP: (P - float(b2)/12, (dP - float(a1)*(P - float(b2)/12) - float(a3))/2)
g = Graphics()
for i, poly in enumerate(E_real):
    g += line2d(poly, color='blue', thickness=2.5, legend_label=None if i else 'real points of E:  y² = x³ + x² − x')
g += point2d(slice_points(real_axis, fix)[::6] + slice_points(oval_line, fix)[::6], color='red', size=9,
             legend_label='Im x = 0 slice of the fixed surface:  u ↦ (℘(u) − b₂/12, ℘′(u)/2)')
# and the rows t = 0, t = 1/2 of the notebook's own vectorised parametrisation (q-series)
import json
load("elliptic_curve_3d.sage")
S, T, xg, yg = parametrize(E, plot_points=300)
rows = np.concatenate([np.stack([xg[:, j].real, yg[:, j].real], axis=1) for j in (0, -1)])
rows = rows[np.abs(rows[:, 0]) < 12]
g += point2d([tuple(p) for p in rows[::3]], color='limegreen', size=4, zorder=4, legend_label='rows t = 0, ½ of parametrize() in the notebook (q-series)')
g += point2d([(r, 0) for r in f_roots], color='blue', size=45, zorder=5)
save(g, 'slice_fixed.png', 'Fixed surface: (x, y) = (℘ − b₂/12, (℘′ − a₁x − a₃)/2) lies on E, and the q-series agrees')
