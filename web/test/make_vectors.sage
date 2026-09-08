# Reference values from Sage/PARI for the JavaScript tests.  Run:  sage make_vectors.sage
import json
from sage.schemes.elliptic_curves.period_lattice import normalise_periods
from sage.databases.cremona import CremonaDatabase, cremona_to_lmfdb
db = CremonaDatabase()
labels = []
for N in range(11, 101):
    labels += [f"{N}{lab}" for lab in db.allcurves(N)]
for N in [389, 5077, 1000, 2310, 9999, 9998, 6000, 3000, 128, 1728, 4000, 7777]:
    labs = sorted(db.allcurves(N))
    if labs:
        labels.append(f"{N}{labs[0]}")
curves = []
for lab in labels:
    E = EllipticCurve(lab)
    L = E.period_lattice()
    w1, w2 = L.basis(prec=53)
    (n1, n2), mat = normalise_periods(CC(w1), CC(w2))
    zs = [CC(0.37, 0.11), CC(w1/2), CC(w2/2), CC((w1 + w2)/2), CC(0.13*w1 + 0.41*w2), CC(-0.45*w1 + 0.05*w2)]
    samples = []
    for z in zs:
        P, dP = L.elliptic_exponential(z, to_curve=False)
        pt = L.elliptic_exponential(z, to_curve=True)
        samples.append({"z": [float(z.real()), float(z.imag())],
                        "wp": [float(P.real()), float(P.imag())], "wpp": [float(dP.real()), float(dP.imag())],
                        "x": [float(pt[0].real()), float(pt[0].imag())], "y": [float(pt[1].real()), float(pt[1].imag())]})
    ei = sorted(E.two_division_polynomial().roots(CC, multiplicities=False), key=lambda r: (r.real(), r.imag()))
    curves.append({"label": lab, "lmfdb": cremona_to_lmfdb(lab, db), "ainvs": [int(a) for a in E.ainvs()], "disc": str(E.discriminant()),
                   "j": [str(E.j_invariant().numerator()), str(E.j_invariant().denominator())],
                   "w1": float(w1), "w2": [float(w2.real()), float(w2.imag())],
                   "normalised": [[float(n1.real()), float(n1.imag())], [float(n2.real()), float(n2.imag())]],
                   "matrix": [int(m) for m in mat],
                   "two_division_roots": [[float(r.real()), float(r.imag())] for r in ei],
                   "samples": samples})
print("curves:", len(curves))

# genus of plane curves (geometric genus via Singular)
R = PolynomialRing(QQ, 'x,y'); x, y = R.gens()
plane = {
    "y^2 = x^3 + x^2 - x": y**2 - (x**3 + x**2 - x),
    "y^2 + y = x^3 - x^2": y**2 + y - x**3 + x**2,
    "y^2 = x^3": y**2 - x**3,
    "y^2 = x^2*(x+1)": y**2 - x**2*(x+1),
    "y^2 = x^5 + 1": y**2 - x**5 - 1,
    "y^2 = x^6 + 1": y**2 - x**6 - 1,
    "y^2 = x^4 + 1": y**2 - x**4 - 1,
    "y^2 = (x^2-1)^2*(x+2)": y**2 - (x**2-1)**2*(x+2),
    "x^3 + y^3 = 1": x**3 + y**3 - 1,
    "x^4 + y^4 = 1": x**4 + y**4 - 1,
    "y^3 = x^4 - 1": y**3 - x**4 + 1,
    "(x^2+y^2)^2 = x^2 - y^2": (x**2+y**2)**2 - x**2 + y**2,
    "x*y = 1": x*y - 1,
    "y = x^2": y - x**2,
    "x^2 + y^2 = 1": x**2 + y**2 - 1,
    "y^2 = x^7 - x": y**2 - x**7 + x,
    "x^3*y + y^3 + x = 0": x**3*y + y**3 + x,
    "y^2*x + y = x^3 + 1": y**2*x + y - x**3 - 1,
    "x^2*y^2 = x^3 + y^3": x**2*y**2 - x**3 - y**3,
    "y^4 = x^3 - x": y**4 - x**3 + x,
    "2*y^2 = x^3 - x": 2*y**2 - x**3 + x,
    "y^2 = 2*x^3 + 1": y**2 - 2*x**3 - 1,
    "y^2 + x^2*y = x^3 + 1": y**2 + x**2*y - x**3 - 1,
    "x^5 + y^5 = 1": x**5 + y**5 - 1,
    "y^2 = x^3 - x + (y - x)^2*0": y**2 - x**3 + x,
}
genus = {}
for k, F in plane.items():
    try:
        genus[k] = int(Curve(F).genus())
    except Exception as e:
        genus[k] = f"error: {type(e).__name__}: {e}"
print(genus)
json.dump({"curves": curves, "genus": genus}, open("vectors.json", "w"))
print("wrote vectors.json")
