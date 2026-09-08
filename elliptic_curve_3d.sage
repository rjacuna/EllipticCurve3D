r"""
Elliptic curves as surfaces in `\RR^3`

The complex points of an elliptic curve `E` defined over `\RR` form a torus
`E(\CC) \subset \CC^2 = \RR^4`.  This module draws its projection
`(\operatorname{Re} x, \operatorname{Im} x, \operatorname{Re} y)` to `\RR^3`,
colored by the period lattice, so that the real points `E(\RR)` are the slice
`\operatorname{Im} x = 0`, and the colored lines are the images of the cycles
`s\,\omega_1 + t\,\omega_2`, `s` or `t` constant.

The uniformisation is the inverse of the Abel--Jacobi map of the *given* model
of `E`: if `\wp` is the Weierstrass function of the period lattice then

.. MATH::

    x = \wp(z) - b_2/12, \qquad y = (\wp'(z) - a_1 x - a_3)/2 .

`\wp` and `\wp'` are evaluated on a whole grid at once by their `q`-series in a
reduced basis of the lattice, which gives double precision on the whole
fundamental domain (a Laurent series about `z = 0` does not).

EXAMPLES::

    sage: S = EllipticCurve3D("20.a3"); S
    Surface in R^3 of Elliptic Curve defined by y^2 = x^3 + x^2 - x over Rational Field
    sage: w1, w2 = S.basis()
    sage: x, y = S.point(w1/2)                      # the 2-torsion point (\sqrt5 - 1)/2
    sage: abs(x - (sqrt(5) - 1)/2) < 1e-12 and abs(y) < 1e-12
    True
    sage: S.plot(radius=3, plot_points=50)
    Graphics3d Object

Load it in a notebook or script with ``load("elliptic_curve_3d.sage")``.

AUTHORS:

- drakezhard (2026-06): original notebook
- drakezhard (2026-09): class, coordinates of the given model, `q`-series
"""
# ****************************************************************************
#       Copyright (C) 2026 drakezhard <drakezhard@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import math
import numpy as np
from sage.misc.cachefunc import cached_method
from sage.plot.colors import colormaps
from sage.plot.plot3d.parametric_surface import ParametricSurface
from sage.plot.plot3d.shapes2 import line3d
from sage.rings.cc import CC
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.schemes.elliptic_curves.period_lattice import normalise_periods
from sage.structure.sage_object import SageObject


def elliptic_curve(f):
    r"""
    Return the elliptic curve described by ``f``.

    INPUT:

    - ``f`` -- an elliptic curve; a Cremona or LMFDB label; a list of
      `a`-invariants; a Weierstrass polynomial in `x` and `y`; or a monic cubic
      `f(x)` in one variable, meaning the curve `y^2 = f(x)`

    OUTPUT: an elliptic curve

    EXAMPLES::

        sage: elliptic_curve("20.a3")
        Elliptic Curve defined by y^2 = x^3 + x^2 - x over Rational Field
        sage: x = polygen(QQ)
        sage: elliptic_curve(x^3 + x^2 - x)
        Elliptic Curve defined by y^2 = x^3 + x^2 - x over Rational Field
        sage: elliptic_curve([0, 1, 0, -1, 0]) == elliptic_curve("20a2")
        True
        sage: R.<u, v> = QQ[]
        sage: elliptic_curve(v^2 + v - u^3 + u^2)
        Elliptic Curve defined by y^2 + y = x^3 - x^2 over Rational Field

    TESTS::

        sage: elliptic_curve(2*x^3 + 1)
        Traceback (most recent call last):
        ...
        ValueError: y^2 = f(x) needs a monic cubic f
    """
    if isinstance(f, EllipticCurve_generic):
        return f
    parent = getattr(f, 'parent', None)
    if parent is not None and hasattr(parent(), 'ngens') and parent().ngens() == 1:
        if f.degree() != 3 or not f.is_monic():
            raise ValueError("y^2 = f(x) needs a monic cubic f")
        a6, a4, a2, _ = f.list()
        return EllipticCurve([0, a2, 0, a4, a6])
    return EllipticCurve(f)


def weierstrass_p_series(zeta, tau, nterms=10):
    r"""
    The Weierstrass function of the lattice `\ZZ + \tau\ZZ` and its derivative,
    evaluated by the `q`-series.

    INPUT:

    - ``zeta`` -- a complex number, or a numpy array of complex numbers, with
      `|\operatorname{Im} \zeta| < \operatorname{Im} \tau`
    - ``tau`` -- a complex number with `\operatorname{Im} \tau > 0`
    - ``nterms`` -- (default: 10) number of terms; for `\tau` in the standard
      fundamental domain `|q| < 0.0044`, and 10 terms give double precision

    OUTPUT: the pair `(\wp(\zeta), \wp'(\zeta))`, numpy arrays if ``zeta`` is one

    ALGORITHM:

    With `u = e^{2\pi i\zeta}` and `q = e^{2\pi i\tau}`,

    .. MATH::

        \wp(\zeta) = (2\pi i)^2 \Bigl[\frac{1}{12} + \sum_{n\in\ZZ} \frac{q^n u}{(1-q^n u)^2}
                     - 2\sum_{n\ge 1} \frac{q^n}{(1-q^n)^2}\Bigr],

    the terms with `n < 0` being rewritten with `w = q^{|n|}/u` so that nothing
    overflows.  At a lattice point `u = 1` and the result is not finite.

    EXAMPLES:

    Against PARI, for the square lattice of `y^2 = x^3 - x`::

        sage: E = EllipticCurve([0, 0, 0, -1, 0])
        sage: L = E.period_lattice(); w1, w2 = L.basis()
        sage: zeta = complex(0.3, 0.2)
        sage: P, dP = weierstrass_p_series(zeta, complex(w2/w1))
        sage: P0, dP0 = L.elliptic_exponential(CC(zeta)*w1, to_curve=False)
        sage: abs(P/w1^2 - P0) < 1e-12 and abs(dP/w1^3 - dP0) < 1e-12
        True

    On an array::

        sage: import numpy as np
        sage: P, dP = weierstrass_p_series(np.array([0.1 + 0.1j, 0.5]), 1j)
        sage: P.shape, np.isfinite(P).all()
        ((2,), True)
    """
    two_pi_i = complex(0, 2*np.pi)
    u = np.exp(two_pi_i*zeta)
    q = np.exp(two_pi_i*complex(tau))
    P = 1/12 + u/(1 - u)**2
    dP = u*(1 + u)/(1 - u)**3
    qn = complex(1)
    for n in range(1, nterms + 1):
        qn = qn*q
        v, w = qn*u, qn/u
        P += v/(1 - v)**2 + w/(1 - w)**2 - 2*qn/(1 - qn)**2
        dP += v*(1 + v)/(1 - v)**3 - w*(1 + w)/(1 - w)**3
    return two_pi_i*two_pi_i*P, two_pi_i*two_pi_i*two_pi_i*dP


class EllipticCurve3D(SageObject):
    r"""
    The complex points of an elliptic curve over `\RR` as a surface in `\RR^3`.

    The surface is `(\operatorname{Re} x, \operatorname{Im} x, \operatorname{Re} y)`
    for the points `(x, y)` of the given model of `E`, parametrised over the half
    fundamental domain `z = s\,\omega_1 + t\,\omega_2`, `-\tfrac12 \le s \le \tfrac12`,
    `0 \le t \le \tfrac12` of the period lattice (a fundamental domain for
    `z \mapsto -z`, so every `x` occurs once).  The rows `t = 0` and, when
    `\Delta > 0`, `t = \tfrac12` are the real components of `E`.

    INPUT:

    - ``E`` -- anything accepted by :func:`elliptic_curve`; the curve must be
      defined over `\QQ` or over a number field with a real embedding
    - ``nterms`` -- (default: 10) number of terms of the `q`-series
    - ``embedding`` -- (default: ``None``) for a curve over a number field, the
      real embedding to use; by default the first one of ``K.embeddings(RR)``

    EXAMPLES::

        sage: S = EllipticCurve3D("11a1"); S
        Surface in R^3 of Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
        sage: S.curve().discriminant() < 0                 # one real component
        True
        sage: len(S.real_components(radius=8))
        1
        sage: len(EllipticCurve3D("20.a3").real_components(radius=3))
        2

    A curve over a number field is taken with a real embedding (by default the
    first one, which for `\QQ(\sqrt2)` sends `\sqrt2` to `-1.414`)::

        sage: K.<r> = QuadraticField(2)
        sage: S = EllipticCurve3D(EllipticCurve(K, [0, 0, 0, 0, r]))   # y^2 = x^3 + r
        sage: x, y = S.point(S.basis()[0]/2)               # the real 2-torsion point: x^3 = -r
        sage: abs(x^6 - 2) < 1e-9 and abs(y) < 1e-9
        True
        sage: S = EllipticCurve3D(EllipticCurve(K, [0, 0, 0, 0, r]), embedding=K.embeddings(RR)[1])
        sage: x, y = S.point(S.basis()[0]/2); abs(x + float(2^(1/6))) < 1e-9
        True

    TESTS::

        sage: EllipticCurve3D(EllipticCurve(QuadraticField(-1, 'i'), [0, 0, 0, 1, 0]))
        Traceback (most recent call last):
        ...
        ValueError: the curve must be defined over a subfield of the reals
    """
    def __init__(self, E, nterms=10, embedding=None):
        r"""
        Initialise ``self``.

        EXAMPLES::

            sage: S = EllipticCurve3D([0, 0, 0, -1, 0])
            sage: TestSuite(S).run(skip=["_test_pickling"])
        """
        E = elliptic_curve(E)
        K = E.base_ring()
        if K is QQ:
            L = E.period_lattice()
            conv = float
        else:                                       # a number field: use a real embedding
            if embedding is None:
                embeddings = K.embeddings(RR) if hasattr(K, 'embeddings') else []
                if not embeddings:
                    raise ValueError("the curve must be defined over a subfield of the reals")
                embedding = embeddings[0]
            L = E.period_lattice(embedding)
            conv = lambda a: float(embedding(a))
        if not L.is_real():
            raise ValueError("the curve must be defined over a subfield of the reals")
        w1, w2 = L.basis()
        (n1, n2), mat = normalise_periods(CC(w1), CC(w2))
        self._E = E
        self._L = L
        self._ainvs = tuple(conv(a) for a in E.ainvs())
        self._b2 = conv(E.b2())
        self._w1, self._w2 = w1, w2
        self._n1, self._n2 = complex(n1), complex(n2)
        self._mat = tuple(int(m) for m in mat)
        self._tau = self._n1/self._n2
        self._nterms = int(nterms)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: EllipticCurve3D("37a1")
            Surface in R^3 of Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        """
        return "Surface in R^3 of %s" % self._E

    def __eq__(self, other):
        r"""
        Two surfaces are equal when their curves are.

        EXAMPLES::

            sage: EllipticCurve3D("20.a3") == EllipticCurve3D("20a2")
            True
            sage: EllipticCurve3D("20.a3") == EllipticCurve3D("11a1")
            False
        """
        return isinstance(other, EllipticCurve3D) and self._E == other._E

    def __ne__(self, other):
        r"""
        EXAMPLES::

            sage: EllipticCurve3D("20.a3") != EllipticCurve3D("11a1")
            True
        """
        return not self == other

    def __hash__(self):
        r"""
        EXAMPLES::

            sage: hash(EllipticCurve3D("20.a3")) == hash(EllipticCurve3D("20a2"))
            True
        """
        return hash(self._E)

    def curve(self):
        r"""
        Return the elliptic curve.

        EXAMPLES::

            sage: EllipticCurve3D("20.a3").curve()
            Elliptic Curve defined by y^2 = x^3 + x^2 - x over Rational Field
        """
        return self._E

    def period_lattice(self):
        r"""
        Return the period lattice of the curve.

        EXAMPLES::

            sage: EllipticCurve3D("20.a3").period_lattice()
            Period lattice associated to Elliptic Curve defined by y^2 = x^3 + x^2 - x over Rational Field
        """
        return self._L

    def basis(self):
        r"""
        Return the basis `(\omega_1, \omega_2)` of the period lattice used for
        the parametrisation: `\omega_1 > 0` is real and `\operatorname{Im}\omega_2 > 0`,
        with `\operatorname{Re}\omega_2 \in \{0, \omega_1/2\}`.

        EXAMPLES::

            sage: w1, w2 = EllipticCurve3D("20.a3").basis(); w1, w2
            (2.82437514195911, 2.27416519904108*I)
        """
        return (self._w1, self._w2)

    def normalised_basis(self):
        r"""
        Return the reduced basis `(n_1, n_2)` in which the `q`-series is
        evaluated, and the matrix `(a, b, c, d)` with
        `(n_1, n_2) = (a\omega_1 + b\omega_2, c\omega_1 + d\omega_2)`.
        Its ratio `\tau = n_1/n_2` lies in the standard fundamental domain.

        EXAMPLES::

            sage: (n1, n2), mat = EllipticCurve3D("20.a3").normalised_basis()
            sage: tau = n1/n2; tau.imag() > 0.86 and abs(tau.real()) <= 0.5
            True
            sage: mat
            (1, 0, 0, -1)
        """
        return (CC(self._n1), CC(self._n2)), self._mat

    def weierstrass_p(self, z):
        r"""
        Return `(\wp(z), \wp'(z))` for the period lattice of the curve.

        INPUT:

        - ``z`` -- a complex number or a numpy array of complex numbers

        OUTPUT: two complex numbers, or two numpy arrays; not finite at lattice points

        EXAMPLES::

            sage: S = EllipticCurve3D("20.a3"); L = S.period_lattice()
            sage: z = CC(0.37, 0.11)
            sage: P, dP = S.weierstrass_p(complex(z))
            sage: P0, dP0 = L.elliptic_exponential(z, to_curve=False)
            sage: abs(P - P0) < 1e-12 and abs(dP - dP0) < 1e-12
            True

        The pole::

            sage: import numpy as np
            sage: with np.errstate(all='ignore'):
            ....:     P, dP = S.weierstrass_p(np.array([0, 0.5]))
            sage: np.isfinite(P)
            array([False,  True])
        """
        z = np.asarray(z, dtype=complex)
        n1, n2 = self._n1, self._n2
        det = n1.real*n2.imag - n2.real*n1.imag
        m = (z.real*n2.imag - n2.real*z.imag)/det
        n = (n1.real*z.imag - z.real*n1.imag)/det
        zeta = (m - np.round(m))*self._tau + n                    # translate so that |Im zeta| <= Im(tau)/2
        with np.errstate(all='ignore'):
            P, dP = weierstrass_p_series(zeta, self._tau, self._nterms)
        P, dP = P/(n2*n2), dP/(n2*n2*n2)
        if P.ndim == 0:
            return complex(P), complex(dP)
        return P, dP

    def point(self, z):
        r"""
        Return the point `(x(z), y(z))` of the given model of the curve,

        .. MATH::

            x = \wp(z) - b_2/12, \qquad y = (\wp'(z) - a_1 x - a_3)/2 .

        INPUT:

        - ``z`` -- a complex number or a numpy array of complex numbers

        OUTPUT: `(x, y)` as complex numbers or numpy arrays

        EXAMPLES::

            sage: S = EllipticCurve3D("11a1"); E = S.curve()
            sage: z = CC(0.37, 0.11)
            sage: x, y = S.point(complex(z))
            sage: P = E.period_lattice().elliptic_exponential(z)
            sage: abs(x - P[0]) < 1e-11 and abs(y - P[1]) < 1e-11
            True
            sage: a1, a2, a3, a4, a6 = E.ainvs()
            sage: abs(y^2 + a1*x*y + a3*y - (x^3 + a2*x^2 + a4*x + a6)) < 1e-9
            True
        """
        P, dP = self.weierstrass_p(z)
        a1, _, a3, _, _ = self._ainvs
        b2 = self._b2
        x = P - b2/12
        y = (dP - a1*x - a3)/2
        return x, y

    def grid(self, plot_points=300):
        r"""
        Return the parametrisation on a grid over the half fundamental domain.

        INPUT:

        - ``plot_points`` -- (default: 300) number of points in each direction

        OUTPUT: ``(S, T, x, y)`` where ``S`` and ``T`` are the 1-dimensional
        arrays of lattice coordinates `s \in [-\tfrac12, \tfrac12]`,
        `t \in [0, \tfrac12]`, and ``x[i, j]``, ``y[i, j]`` are the complex
        coordinates of the point `z = S[i]\,\omega_1 + T[j]\,\omega_2`

        EXAMPLES::

            sage: import numpy as np
            sage: S, T, x, y = EllipticCurve3D("20.a3").grid(101)
            sage: x.shape
            (101, 101)
            sage: finite = np.isfinite(x[:, 0])                  # the pole z = 0 is a grid point
            sage: np.abs(x[finite, 0].imag).max() < 1e-9 and np.abs(x[:, -1].imag).max() < 1e-9    # real rows
            True
            sage: sorted(round(float(v), 6) for v in (x[50, -1].real, x[-1, -1].real, x[-1, 0].real))
            [-1.618034, 0.0, 0.618034]
        """
        S = np.linspace(-0.5, 0.5, int(plot_points))
        T = np.linspace(0, 0.5, int(plot_points))
        SS, TT = np.meshgrid(S, T, indexing='ij')
        x, y = self.point(SS*float(self._w1) + TT*complex(self._w2))
        return S, T, x, y

    def real_components(self, radius=None, plot_points=600):
        r"""
        Return the real points of the curve as polylines.

        These are the rows `t = 0` and `t = \tfrac12` of :meth:`grid` on which
        `x` and `y` are real, restricted to `|(x, 0, y)| \le` ``radius`` and
        ordered along the curve; the unbounded component starts and ends at the
        pole, and a compact component is closed.

        INPUT:

        - ``radius`` -- (default: ``None``, meaning no restriction) clip to the ball of this radius
        - ``plot_points`` -- (default: 600) number of sample points along each component

        OUTPUT: a list of numpy arrays of shape ``(k, 2)`` with columns `x, y`

        EXAMPLES::

            sage: comps = EllipticCurve3D("20.a3").real_components(radius=3)
            sage: [len(c) > 100 for c in comps]
            [True, True]
            sage: oval = min(comps, key=lambda c: c[:, 0].max())
            sage: (oval[0] == oval[-1]).all()                    # closed
            True
            sage: round(float(oval[:, 0].min()), 4), abs(float(oval[:, 0].max())) < 1e-9
            (-1.618, True)
        """
        S, T, x, y = self.grid(plot_points)
        n = len(S)
        comps = []
        rad = float(radius) if radius is not None else float('inf')
        for j in (0, n - 1):
            xr, yr = x[:, j], y[:, j]
            finite = np.isfinite(xr) & np.isfinite(yr)
            inside = finite & (np.hypot(xr.real, yr.real) <= rad)
            real = finite & (np.abs(xr.imag) < 1e-8) & (np.abs(yr.imag) < 1e-8)
            if not np.any(inside) or np.any(inside & ~real):
                continue
            order = np.roll(np.arange(n), -(n//2 + 1)) if j == 0 else np.arange(n)
            cur, prev = [], None
            for i in order:
                if not inside[i]:
                    if len(cur) > 1:
                        comps.append(cur)
                    cur, prev = [], None
                    continue
                p = (float(xr[i].real), float(yr[i].real))
                if prev is not None and math.hypot(p[0] - prev[0], p[1] - prev[1]) > 0.5*min(rad, 1e300):
                    if len(cur) > 1:
                        comps.append(cur)
                    cur = []
                cur.append(p)
                prev = p
            if len(cur) > 1:
                comps.append(cur)
            if j == n - 1 and comps and len(comps[-1]) == n:
                comps[-1].append(comps[-1][0])
        return [np.array(c) for c in comps]

    def surface(self, colormap=None, spacing=1/8, cutoff=0, radius=1, opacity=1, plot_points=300):
        r"""
        Return the surface `(\operatorname{Re} x, \operatorname{Im} x, \operatorname{Re} y)`
        colored by the period lattice.

        INPUT:

        - ``colormap`` -- (default: ``colormaps.Greens_r``) a matplotlib colormap
        - ``spacing`` -- (default: `1/8`) period of the lattice coloring in the
          lattice coordinates `(s, t)`; the dark lines are the cycles
          `s \in` ``spacing`` `\ZZ` and `t \in` ``spacing`` `\ZZ`
        - ``cutoff``, ``radius`` -- (default: 0 and 1) keep the part of the
          surface with ``cutoff`` `\le |(X, Y, Z)| <` ``radius``
        - ``opacity`` -- (default: 1) opacity of the surface, kept through the clipping
        - ``plot_points`` -- (default: 300) grid points in each lattice direction

        OUTPUT: a Graphics3d object

        EXAMPLES::

            sage: S = EllipticCurve3D("20.a3")
            sage: G = S.surface(colormaps.Greens_r, radius=3, opacity=1/2, plot_points=60); G
            Graphics3d Object
            sage: import numpy as np
            sage: V = np.array(G.vertex_list()); np.isfinite(V).all() and np.abs(V).max() <= 3 + 1e-6
            True
            sage: G.threejs_repr(G.default_render_params())[0][1]['opacity']
            0.5
        """
        if colormap is None:
            colormap = colormaps.Greens_r
        S, T, x, y = self.grid(plot_points)
        X, Y, Z = (np.nan_to_num(A, nan=1e150, posinf=1e150, neginf=-1e150)
                   for A in (x.real, x.imag, y.real))            # the pole is clipped away below
        n, half, sp = int(plot_points) - 1, 0.5, float(spacing)

        def point(s, t):                                         # look the precomputed grid up by index
            i, j = int(round((s + half)*n)), int(round(2*t*n))
            return X[i, j], Y[i, j], Z[i, j]

        def lattice_cf(s, t):
            sf, tf = s - math.floor(s), t - math.floor(t)
            return abs(math.sin(math.pi*sf/sp)*math.sin(math.pi*tf/sp))**0.12

        G = ParametricSurface(point, (S.tolist(), T.tolist()), color=(lattice_cf, colormap), opacity=opacity)
        r2, c2 = float(radius)**2, float(cutoff)**2
        G = G.add_condition(lambda X, Y, Z: c2 <= X*X + Y*Y + Z*Z < r2)
        # add_condition rebuilds the mesh and forgets the viewing options; restore them
        G._extra_kwds.update(opacity=float(opacity), aspect_ratio=1)
        return G

    def real_points_plot(self, radius=1, plot_points=600, color='red', thickness=6):
        r"""
        Return the real points of the curve drawn in the plane `\operatorname{Im} x = 0`.

        INPUT:

        - ``radius`` -- (default: 1) clip to the ball of this radius
        - ``plot_points`` -- (default: 600) sample points along each component
        - ``color``, ``thickness`` -- passed to :func:`line3d`

        OUTPUT: a Graphics3d object, empty if no real point lies in the ball

        EXAMPLES::

            sage: EllipticCurve3D("20.a3").real_points_plot(radius=3)
            Graphics3d Object
            sage: EllipticCurve3D("11a1").real_points_plot(radius=3)     # nothing inside
            Graphics3d Object
        """
        from sage.plot.plot3d.base import Graphics3dGroup
        parts = [line3d([(p[0], 0, p[1]) for p in comp], color=color, thickness=thickness)
                 for comp in self.real_components(radius, plot_points)]
        return Graphics3dGroup(parts)

    def plot(self, colormap=None, spacing=1/8, cutoff=0, radius=1, opacity=1, plot_points=300,
             real_points=True, color='red', thickness=6):
        r"""
        Return the surface together with the real points of the curve.

        INPUT: as for :meth:`surface` and :meth:`real_points_plot`; ``real_points``
        (default: ``True``) says whether to draw the latter.

        EXAMPLES::

            sage: EllipticCurve3D("37a1").plot(radius=3, plot_points=50)
            Graphics3d Object
        """
        G = self.surface(colormap, spacing, cutoff, radius, opacity, plot_points)
        if real_points:
            G = G + self.real_points_plot(radius, max(plot_points, 400), color, thickness)
            G._extra_kwds.update(opacity=float(opacity), aspect_ratio=1)
        return G


def parametrize(E, plot_points=300, nterms=10):
    r"""
    Return ``EllipticCurve3D(E, nterms).grid(plot_points)``; kept for the scripts
    written against the notebook version.

    EXAMPLES::

        sage: S, T, x, y = parametrize("20.a3", 51); x.shape
        (51, 51)
    """
    return EllipticCurve3D(E, nterms).grid(plot_points)


def curve(f, colormap=None, spacing=1/8, cutoff=0, radius=1, opacity=1, plot_points=300, nterms=10):
    r"""
    Return ``EllipticCurve3D(f, nterms).surface(...)``; the interface of the
    original notebook function.

    EXAMPLES::

        sage: curve("20.a3", colormaps.Greens_r, radius=3, opacity=1/2, plot_points=40)
        Graphics3d Object
    """
    return EllipticCurve3D(f, nterms).surface(colormap, spacing, cutoff, radius, opacity, plot_points)
