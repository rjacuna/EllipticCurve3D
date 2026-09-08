# Cell 0 of Plotting/Elliptic Curve 3D.ipynb, verbatim, kept for reference and for the
# 'before' picture in render_figures.sage.  See README.md for what is wrong with it.

from sage.schemes.riemann_surfaces.riemann_surface import RiemannSurface

R = PolynomialRing(QQ, 2, names=('x','y'))
x, y = R.gens()

def curve(f, colormap, 
          spacing=1/8, 
          precision=10,
          cutoff = 0,
          radius = 1,
          opacity =1, plot_points= 300 ):
    """
    Given a defining polynomial f for an elliptic curve y^2 = f(x),
    and a Sage matplotlib-style colormap object, this routine
    plots the ℘-surface patch over one fundamental parallelogram
    together with the colored plane beneath it.
    """
    # 4) Build ℘ and its half‐derivative:
    E = EllipticCurve(f)  # or EllipticCurve([…])
    p  = E.weierstrass_p(prec=precision).truncate(precision)
    pp = p.derivative()

        # 1) Build Riemann surface and get periods:
    S     = E.riemann_surface()
    Ω     = S.period_matrix()
    alpha, beta = Ω[0,0], Ω[0,1]
    # 2) Compute basis‐change for lattice coloring:
    M  = matrix(RR, [[alpha.real(), beta.real()],
                     [alpha.imag(), beta.imag()]])
    Mi = M.inverse()
    
    def lattice_cf(x, y, spacing=spacing):    
        u, v    = Mi * vector([x, y])
        uf, vf  = u - math.floor(u), v - math.floor(v)
        return abs(math.sin(math.pi*uf/spacing)* math.sin(math.pi*vf/spacing))**0.12

    # 3) Set up radii a,b from periods:
    a = max(abs(alpha.real()), abs(beta.real()))
    b = max(abs(alpha.imag()), abs(beta.imag()))

    # 5) Parametric surface patch:
    X = lambda u,v: float(p(u + 1j*v).real()) 
    Y = lambda u,v: float(p(u + 1j*v).imag()) 
    Z = lambda u,v: float(pp(u + 1j*v).real())
    E_C = parametric_plot3d(
         [X,Y,Z],
        (-(a/2), (a/2)), (0,b/2 ),
        aspect_ratio=1,
        color=(lattice_cf, colormap),
        plot_points=plot_points,
        opacity = 0.1
    ).add_condition(lambda X,Y,Z: X**2 + Y**2 + Z**2 < (radius)**2)
    
    return E_C

