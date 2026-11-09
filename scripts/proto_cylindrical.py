#!/usr/bin/env python3
"""Prototyp zu Plan v2.14, Abschnitte M20/M21 (Zylinderkoordinaten).

Prueft die Fourier-Bessel-Formulierung an exakten Fasermoden: in einer y-invarianten Scheibe mit
Kern (eps1, Radius a) und Mantel eps2 muessen die gefuehrten Moden Eigenwerte +-i beta der
Propagationsmatrix (Plan eq:cyl_M) sein.

  python3 scripts/proto_cylindrical.py fiber     # V21b: Laurent vs. Li in r (Plan eq:cyl_li)
  python3 scripts/proto_cylindrical.py offaxis   # V22b: exzentrischer Kern, P = 1 (Plan eq:cylP_Q)

Normierte Einheiten (k0 = 1). Nur numpy/scipy. Kein Bestandteil des Loesers.
"""
import sys
import numpy as np
from scipy.special import jv, kv, jvp, kvp
from scipy.optimize import brentq

def J(n, x):
    return jv(n, x)

def gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

def knodes(kmax, nk, panels):
    # zusammengesetzte Gauss-Legendre-Knoten, Panelgrenzen bei kritischen k
    edges = np.unique(np.concatenate(([0.0], panels, [kmax])))
    edges = edges[edges <= kmax]
    lens = np.diff(edges)
    counts = np.maximum(4, np.round(nk * lens / kmax).astype(int))
    ks, ws = [], []
    for a, b, c in zip(edges[:-1], edges[1:], counts):
        x, w = gl(a, b, c)
        ks.append(x); ws.append(w)
    return np.concatenate(ks), np.concatenate(ws)

class Grid:
    def __init__(self, k, w, rbreaks, nr=60):
        self.k, self.w = k, w
        # r-Quadratur fuer kompakte Traeger: Panels zwischen den Bruchstellen
        rs, rw = [], []
        kmax = k.max()
        for a, b in zip(rbreaks[:-1], rbreaks[1:]):
            npan = max(1, int(np.ceil((b - a) * kmax / 3.0)))
            for p in range(npan):
                x, ww = gl(a + (b - a) * p / npan, a + (b - a) * (p + 1) / npan, nr)
                rs.append(x); rw.append(ww)
        self.r, self.rw = np.concatenate(rs), np.concatenate(rw)

    def gal(self, g, nu, nup):
        """[[g]]_{nu<-nup}(j,l) = w_l k_l int g(r) J_nu(k_j r) J_nup(k_l r) r dr (g kompakt)."""
        k, r = self.k, self.r
        A = J(nu, np.outer(k, r)) * (g(r) * r * self.rw)[None, :]
        B = J(nup, np.outer(r, k))
        return A @ B * (self.w * self.k)[None, :]

def Mmatrix(G, m, Qpp, Qpm, Qmp, Qmm, Ey_inv):
    k = G.k; n = len(k); I = np.eye(n); K = np.diag(k)
    KEK = K @ Ey_inv @ K
    Z = np.zeros((n, n))
    top = np.block([[Z, Z, 0.5 * KEK - I, 0.5 * KEK],
                    [Z, Z, -0.5 * KEK, -0.5 * KEK + I]])
    bot = np.block([[Qpp - 0.5 * K @ K, Qpm - 0.5 * K @ K, Z, Z],
                    [-Qmp + 0.5 * K @ K, -Qmm + 0.5 * K @ K, Z, Z]])
    return np.vstack([top, bot])

def build(m, eps1, eps2, a, kmax, nk, method, ref_split=0.5):
    k, w = knodes(kmax, nk, [np.sqrt(eps2), np.sqrt(eps1)])
    G = Grid(k, w, [0.0, ref_split * a, a])
    n = len(k); I = np.eye(n)
    d = eps1 - eps2
    chi = lambda r: (r <= a).astype(float)
    p, q = m + 1, m - 1
    Ey = eps2 * I + G.gal(lambda r: d * chi(r), m, m)          # Laurent fuer D_y
    if method == "laurent":
        Qpp = eps2 * I + G.gal(lambda r: d * chi(r), p, p)
        Qmm = eps2 * I + G.gal(lambda r: d * chi(r), q, q)
        Qpm = np.zeros((n, n)); Qmp = np.zeros((n, n))
    else:
        # Li in r: G+- = D_r/eps_ref +- i E_phi, eps_ref glatt: eps1 auf der Achse, eps2 aussen
        r1, r2 = ref_split * a, a
        def s(r):
            t = np.clip((r - r1) / (r2 - r1), 0, 1)
            return t**3 * (10 - 15 * t + 6 * t**2)
        if method == "li_r_axis_ignored":
            eref = lambda r: eps2 + 0 * r
        else:
            eref = lambda r: np.exp((1 - s(r)) * np.log(eps1 + 0j) + s(r) * np.log(eps2 + 0j))
        eps = lambda r: np.where(r <= a, eps1, eps2)
        rho = lambda r: eref(r) / eps(r)
        # Bloecke: Hintergrund analytisch (1 bzw. eps2), Rest kompakt auf [0, a]
        A11 = I + G.gal(lambda r: 0.5 * (rho(r) - 1) * chi(r), p, p)
        A22 = I + G.gal(lambda r: 0.5 * (rho(r) - 1) * chi(r), q, q)
        A12 = G.gal(lambda r: 0.5 * (rho(r) - 1) * chi(r), p, q)
        A21 = G.gal(lambda r: 0.5 * (rho(r) - 1) * chi(r), q, p)
        B11 = eps2 * I + G.gal(lambda r: 0.5 * ((eref(r) - eps2) + (eps(r) - eps2)) * chi(r), p, p)
        B22 = eps2 * I + G.gal(lambda r: 0.5 * ((eref(r) - eps2) + (eps(r) - eps2)) * chi(r), q, q)
        B12 = G.gal(lambda r: 0.5 * (eref(r) - eps(r)) * chi(r), p, q)
        B21 = G.gal(lambda r: 0.5 * (eref(r) - eps(r)) * chi(r), q, p)
        Aop = np.block([[A11, A12], [A21, A22]])
        Bop = np.block([[B11, B12], [B21, B22]])
        Q = Bop @ np.linalg.inv(Aop)
        Qpp, Qpm, Qmp, Qmm = Q[:n, :n], Q[:n, n:], Q[n:, :n], Q[n:, n:]
    M = Mmatrix(G, m, Qpp, Qpm, Qmp, Qmm, np.linalg.inv(Ey))
    lam = np.linalg.eigvals(M)
    return lam

# exakte Fasermoden
def fiber_beta(m, eps1, eps2, a, kind):
    n1, n2 = np.sqrt(eps1), np.sqrt(eps2)
    def UW(b):
        return a * np.sqrt(eps1 - b * b), a * np.sqrt(b * b - eps2)
    def f(b):
        U, W = UW(b)
        jj = jvp(m, U) / (U * jv(m, U)); kk = kvp(m, W) / (W * kv(m, W))
        if m == 0 and kind == "TE":
            return jj + kk
        if m == 0 and kind == "TM":
            return eps1 * jj + eps2 * kk
        return (jj + kk) * (eps1 * jj + eps2 * kk) - (m * b) ** 2 * (1 / U**2 + 1 / W**2) ** 2
    bs = np.linspace(n2 + 1e-9, n1 - 1e-9, 20001)
    vals = np.array([f(b) for b in bs])
    roots = []
    for i in range(len(bs) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] * vals[i + 1] < 0:
            r = brentq(f, bs[i], bs[i + 1])
            if abs(f(r)) < 1e-6:
                roots.append(r)
    return sorted(roots, reverse=True)

def build2(m, eps1, eps2, a, kmax, nk, method, r1f=0.5, r2f=1.0):
    # Variante mit waehlbarer Uebergangszone [r1f a, r2f a] fuer eps_ref
    k, w = knodes(kmax, nk, [np.sqrt(eps2), np.sqrt(eps1)])
    rb = sorted(set([0.0, r1f*a, a, r2f*a]))
    G = Grid(k, w, rb)
    n = len(k); I = np.eye(n); d = eps1-eps2
    R = max(a, r2f*a)
    chi = lambda r: (r <= R+1e-12).astype(float)
    p, q = m+1, m-1
    Ey = eps2*I + G.gal(lambda r: d*(r<=a), m, m)
    if method == "laurent":
        Qpp = eps2*I + G.gal(lambda r: d*(r<=a), p, p); Qmm = eps2*I + G.gal(lambda r: d*(r<=a), q, q)
        Qpm = Qmp = np.zeros((n,n))
    else:
        def s(r):
            t = np.clip((r - r1f*a)/((r2f-r1f)*a), 0, 1); return t**3*(10-15*t+6*t**2)
        if method == "ref_inf": eref = lambda r: eps2+0j*r
        else: eref = lambda r: np.exp((1-s(r))*np.log(eps1+0j)+s(r)*np.log(eps2+0j))
        eps = lambda r: np.where(r <= a, eps1, eps2)
        rho = lambda r: eref(r)/eps(r)
        g1 = lambda r: 0.5*(rho(r)-1)*chi(r)
        gB = lambda r: 0.5*((eref(r)-eps2)+(eps(r)-eps2))*chi(r)
        gC = lambda r: 0.5*(eref(r)-eps(r))*chi(r)
        A = np.block([[I+G.gal(g1,p,p), G.gal(g1,p,q)],[G.gal(g1,q,p), I+G.gal(g1,q,q)]])
        B = np.block([[eps2*I+G.gal(gB,p,p), G.gal(gC,p,q)],[G.gal(gC,q,p), eps2*I+G.gal(gB,q,q)]])
        Q = B@np.linalg.inv(A)
        Qpp,Qpm,Qmp,Qmm = Q[:n,:n],Q[:n,n:],Q[n:,:n],Q[n:,n:]
    M = Mmatrix(G, m, Qpp,Qpm,Qmp,Qmm, np.linalg.inv(Ey))
    return np.linalg.eigvals(M)

def coef(q, r, a, x0):
    # Winkelkoeffizient der Indikatorfunktion einer Kreisscheibe (Radius a, Mitte (x0,0)) auf dem Kreis r
    c = (r**2 + x0**2 - a**2) / (2*r*x0)
    phi0 = np.where(c <= -1, np.pi, np.where(c >= 1, 0.0, np.arccos(np.clip(c, -1, 1))))
    return phi0/np.pi if q == 0 else np.sin(q*phi0)/(np.pi*q)

def run(eps1, eps2, a, x0, m0, L, kmax, nk):
    k, w = knodes(kmax, nk, [np.sqrt(eps2), np.sqrt(eps1)])
    rb = [0.0] + sorted(set([abs(a-x0), a+x0]))
    G = Grid(k, w, rb, nr=40)
    n = len(k); I = np.eye(n); K = np.diag(k); d = eps1-eps2
    ms = [m0 + l for l in range(-L, L+1)]; nh = len(ms)
    def Qblock(shift):
        Qm = np.zeros((nh*n, nh*n), complex)
        for i, mi in enumerate(ms):
            for j, mj in enumerate(ms):
                q = mi - mj
                g = lambda r, q=q: d*coef(q, r, a, x0)
                blk = G.gal(g, mi+shift, mj+shift)
                if i == j: blk = blk + eps2*I
                Qm[i*n:(i+1)*n, j*n:(j+1)*n] = blk
        return Qm
    Qp, Qm_, Qy = Qblock(+1), Qblock(-1), Qblock(0)
    Ny = nh*n; Id = np.eye(Ny); KK = np.kron(np.eye(nh), K)
    KEK = KK @ np.linalg.inv(Qy) @ KK
    Z = np.zeros((Ny, Ny))
    top = np.block([[Z, Z, 0.5*KEK-Id, 0.5*KEK],[Z, Z, -0.5*KEK, -0.5*KEK+Id]])
    bot = np.block([[Qp-0.5*KK@KK, -0.5*KK@KK, Z, Z],[0.5*KK@KK, -Qm_+0.5*KK@KK, Z, Z]])
    lam = np.linalg.eigvals(np.vstack([top, bot]))
    b = np.imag(lam[np.abs(np.real(lam)) < 1e-6*np.abs(lam)]); return np.sort(b[b > np.sqrt(eps2)+1e-6])[::-1]


def study_fiber():
    eps1, eps2, a = 4.0, 1.0, 2.0
    for m, kind in [(0, "TM"), (1, "HE"), (2, "HE")]:
        b0 = fiber_beta(m, eps1, eps2, a, kind)[0]
        print(f"m={m} {kind}: beta_exakt={b0:.12f}")
        for meth, label, r1, r2 in [("laurent", "Laurent", .5, 1), ("ref_inf", "Li r, eps_inf", .5, 1),
                                     ("ref", "Li r, glatt", .5, 1)]:
            errs = []
            for kmax in [4, 8, 16]:
                lam = build2(m, eps1, eps2, a, kmax, int(20 * kmax), meth, r1, r2)
                bb = np.imag(lam[np.abs(np.real(lam)) < 1e-6 * np.abs(lam)]); bb = bb[bb > np.sqrt(eps2) + 1e-6]
                errs.append(np.min(np.abs(bb - b0)) if len(bb) else np.nan)
            print(f"  {label:14s} k_max=4/8/16: " + "  ".join(f"{e:.1e}" for e in errs))

def study_offaxis():
    eps1, eps2, a = 4.0, 1.0, 2.0
    bc = run(eps1, eps2, a, 1e-9, 1, 0, 3.0, 60)[0]
    print("zentriert, Laurent m=1:", bc)
    for x0, L in [(0.3, 3), (0.3, 5), (0.6, 5), (0.6, 7)]:
        b = run(eps1, eps2, a, x0, 1, L, 3.0, 60)
        print(f"x0={x0} L={L}: max. Abweichung HE11-Paar = {np.max(np.abs(b[:2] - bc)):.1e}")

if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "fiber"
    {"fiber": study_fiber, "offaxis": study_offaxis}[what]()
