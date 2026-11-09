"""
Background model and likelihoods used for the published next-release prediction
(Deshmukh, Annalen der Physik 538, e70321 (2026), Section 11 and Table 15).

Same numerics as the paper pipeline:
  flat, radiation included, r_d fixed at 147.09 Mpc,
  lcdm : E^2 = Om(1+z)^3 + Or(1+z)^4 + OL
  frac : E^2 = lcdm + Om(1+z)^3 f(z)/2,  f(z) = exp(-z/z_char)
Distances are integrated on a fixed grid z in [0, 2.5] with dz = 1e-3.
All functions are vectorised over parameter samples.
"""
import os
import numpy as np
from scipy.linalg import cholesky, solve_triangular

C_KMS = 299792.458
RD = 147.09                                      # Planck 2018 r_d [Mpc], as in the paper
OMEGA_R_H2 = 2.469e-5 * (1 + 0.2271 * 3.044)     # photons + massless neutrinos

ZMAX = 2.5
ZG = np.linspace(0.0, ZMAX, 2501)
DZ = ZG[1] - ZG[0]

QUANTITIES = ('DM_over_rs', 'DH_over_rs', 'DV_over_rs')


def E2(z, model, H, Om, zc):
    """z: (nz,), H/Om/zc: (n,) -> (n, nz)"""
    H = np.atleast_1d(H)[:, None]
    Om = np.atleast_1d(Om)[:, None]
    Or = OMEGA_R_H2 / (H / 100.0) ** 2
    OL = 1.0 - Om - Or
    zz = np.asarray(z, float)[None, :]
    e2 = Om * (1 + zz) ** 3 + Or * (1 + zz) ** 4 + OL
    if model == 'frac':
        zc = np.atleast_1d(zc)[:, None]
        e2 = e2 + Om * (1 + zz) ** 3 * np.exp(-zz / zc) / 2.0
    elif model != 'lcdm':
        raise ValueError(model)
    return e2


def comoving(zt, model, H, Om, zc):
    """D_M(zt) in Mpc, shape (n, len(zt))"""
    zt = np.asarray(zt, float)
    if zt.max() > ZMAX:
        raise ValueError(f'redshift {zt.max()} beyond the integration grid (z <= {ZMAX})')
    ig = 1.0 / np.sqrt(E2(ZG, model, H, Om, zc))
    cum = np.concatenate([np.zeros((ig.shape[0], 1)),
                          np.cumsum(0.5 * (ig[:, 1:] + ig[:, :-1]) * DZ, axis=1)], axis=1)
    out = np.empty((ig.shape[0], len(zt)))
    for i in range(ig.shape[0]):
        out[i] = np.interp(zt, ZG, cum[i])
    return C_KMS / np.atleast_1d(H)[:, None] * out


def Hz(zt, model, H, Om, zc):
    return np.atleast_1d(H)[:, None] * np.sqrt(E2(zt, model, H, Om, zc))


def bao_observables(z, qty, model, H, Om, zc):
    """BAO observables D/r_d for arrays z, qty -> (n, len(z))"""
    z = np.asarray(z, float)
    qty = np.asarray(qty)
    bad = [q for q in qty if q not in QUANTITIES]
    if bad:
        raise ValueError(f'unsupported BAO quantities: {sorted(set(bad))}')
    DM = comoving(z, model, H, Om, zc)
    DH = C_KMS / Hz(z, model, H, Om, zc)
    DV = (z[None, :] * DM ** 2 * DH) ** (1 / 3)
    return np.where(qty == 'DM_over_rs', DM, np.where(qty == 'DH_over_rs', DH, DV)) / RD


def H_local(H, Om):
    """local expansion rate of the fractal model, H_E sqrt(1 + Om/2)"""
    return np.asarray(H) * np.sqrt(1 + np.asarray(Om) / 2)


# ---------------------------------------------------------------- data
class BAO:
    """BAO data in Cobaya bao_data format: mean file '# z value quantity', covariance matrix file."""

    def __init__(self, mean_file, cov_file):
        rows = [l.split() for l in open(mean_file) if l.strip() and not l.startswith('#')]
        self.z = np.array([float(r[0]) for r in rows])
        self.obs = np.array([float(r[1]) for r in rows])
        self.qty = np.array([r[2] for r in rows])
        self.cov = np.atleast_2d(np.loadtxt(cov_file))
        if self.cov.shape != (len(self.z), len(self.z)):
            raise ValueError(f'covariance shape {self.cov.shape} does not match {len(self.z)} data points')
        L = cholesky(self.cov, lower=True)
        self.W = solve_triangular(L, np.eye(len(self.z)), lower=True)
        self.err = np.sqrt(np.diag(self.cov))
        self.n = len(self.z)

    def chi2(self, model, H, Om, zc):
        r = self.obs[None, :] - bao_observables(self.z, self.qty, model, H, Om, zc)
        u = r @ self.W.T
        return np.sum(u * u, axis=1)


class SN:
    """Pantheon+ (non-calibrator SNe, full STAT+SYS covariance, M_B marginalised analytically)."""

    def __init__(self, dat_file, cov_file):
        import pandas as pd
        d = pd.read_csv(dat_file, sep=r'\s+')
        m = (d['IS_CALIBRATOR'] == 0).values
        self.zhd = d['zHD'].values[m]
        self.zhel = d['zHEL'].values[m]
        self.mb = d['m_b_corr'].values[m]
        with open(cov_file) as f:
            N = int(f.readline())
            C = np.array(f.read().split(), dtype=float).reshape(N, N)
        C = C[np.ix_(m, m)]
        L = cholesky(C, lower=True)
        self.W = solve_triangular(L, np.eye(len(self.zhd)), lower=True)
        self.e = self.W @ np.ones(len(self.zhd))
        self.ee = self.e @ self.e
        self.n = len(self.zhd)

    def chi2(self, model, H, Om, zc):
        DM = comoving(self.zhd, model, H, Om, zc)
        mu = 5 * np.log10((1 + self.zhel)[None, :] * DM) + 25
        u = (self.mb[None, :] - mu) @ self.W.T
        eu = u @ self.e
        return np.sum(u * u, axis=1) - eu ** 2 / self.ee


# ---------------------------------------------------------------- grid posterior
def grid_chi2(bao, sn, model, nH=61, nO=61, nX=70, batch=256):
    """chi^2 of BAO + SN on the grid used for the published prediction (flat priors)."""
    if model == 'frac':
        Hs = np.linspace(52, 78, nH); Os = np.linspace(0.20, 0.50, nO)
        Xs = np.linspace(np.log(0.02), np.log(3.0), nX)
    else:
        Hs = np.linspace(52, 78, 121); Os = np.linspace(0.20, 0.50, 121)
        Xs = np.array([0.0])
    HH, OO = np.meshgrid(Hs, Os, indexing='ij')
    HH = HH.ravel(); OO = OO.ravel()
    chi = np.empty((len(Xs), HH.size))
    for i, x in enumerate(Xs):
        zc = np.full(HH.size, np.exp(x)) if model == 'frac' else None
        c = np.zeros(HH.size)
        for s in range(0, HH.size, batch):
            sl = slice(s, s + batch)
            zz = None if zc is None else zc[sl]
            c[sl] = bao.chi2(model, HH[sl], OO[sl], zz) + sn.chi2(model, HH[sl], OO[sl], zz)
        chi[i] = c
    return Hs, Os, Xs, chi.reshape(len(Xs), len(Hs), len(Os))


def sample_grid(Hs, Os, Xs, chi, rng, n=6000):
    """draw posterior samples from the gridded likelihood, jittered within cells"""
    L = np.exp(-0.5 * (chi - chi.min()))
    w = L.ravel() / L.sum()
    idx = rng.choice(w.size, size=n, p=w)
    ix, ih, io = np.unravel_index(idx, chi.shape)
    dH = Hs[1] - Hs[0]; dO = Os[1] - Os[0]; dX = (Xs[1] - Xs[0]) if len(Xs) > 1 else 0
    H = Hs[ih] + (rng.random(n) - 0.5) * dH
    Om = Os[io] + (rng.random(n) - 0.5) * dO
    X = Xs[ix] + (rng.random(n) - 0.5) * dX
    return H, Om, np.exp(X), float(chi.min())
