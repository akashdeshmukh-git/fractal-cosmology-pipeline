"""Fractal-correction H(z) fit to DESI DR2 BAO + Pantheon+ + compressed CMB.
Same equations and data as the paper / Cobaya reproduction, written in plain numpy
so r_d can be treated as fixed or free."""
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid

C = 299792.458
NEFF = 3.044
ORH2 = 2.469e-5 * (1 + 0.2271 * NEFF)
RD_FID = 147.09
ZSTAR = 1089.92
DMSTAR_OBS, DMSTAR_ERR = 13873.0, 25.0
WM_OBS, WM_ERR = 0.1430, 0.0011

# ---------- grid for distance integrals ----------
ZG = np.concatenate([np.linspace(0, 3, 30001)[:-1], np.geomspace(3, 1200, 20001)])


def E_of_z(z, H, Om, shape):
    h = H / 100.0
    Or = ORH2 / h**2
    OL = 1 - Om - Or
    kind, p = shape
    if kind == 'lcdm':
        f = 0.0
    elif kind == 'exp':
        f = np.exp(-z / p)
    elif kind == 'pl':
        f = (1 + z) ** (-p)
    a3 = (1 + z) ** 3
    return np.sqrt(Om * a3 + Or * (1 + z) ** 4 + OL + Om * a3 * f / 2)


def comoving(H, Om, shape):
    """Return grid z and D_M(z) in Mpc (flat)."""
    invE = 1.0 / E_of_z(ZG, H, Om, shape)
    DC = C / H * cumulative_trapezoid(invE, ZG, initial=0)
    return DC


# ---------- BAO ----------
_b = np.loadtxt('bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt', dtype=str)
BZ = _b[:, 0].astype(float)
BV = _b[:, 1].astype(float)
BQ = _b[:, 2]
BICOV = np.linalg.inv(np.loadtxt('bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt'))


def bao_pred(H, Om, shape, rd, DC=None):
    if DC is None:
        DC = comoving(H, Om, shape)
    DM = np.interp(BZ, ZG, DC)
    DH = C / (H * E_of_z(BZ, H, Om, shape))
    DV = (BZ * DM**2 * DH) ** (1 / 3)
    out = np.where(BQ == 'DM_over_rs', DM, np.where(BQ == 'DH_over_rs', DH, DV))
    return out / rd


def chi2_bao(H, Om, shape, rd, DC=None):
    r = bao_pred(H, Om, shape, rd, DC) - BV
    return r @ BICOV @ r


# ---------- Pantheon+ ----------
_sn = pd.read_csv('sn_data/PantheonPlus/Pantheon+SH0ES.dat', sep=r'\s+')
_cov_full = np.loadtxt('sn_data/PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov', skiprows=1)
_n = len(_sn)
_cov_full = _cov_full.reshape(_n, _n)


def make_sn(cut):
    if cut == 'paper':
        m = (_sn.IS_CALIBRATOR == 0).values
    else:  # standard Pantheon+ cosmology cut
        m = (_sn.zHD > 0.01).values
    d = _sn[m]
    icov = np.linalg.inv(_cov_full[np.ix_(m, m)])
    return dict(zhd=d.zHD.values, zhel=d.zHEL.values, mb=d.m_b_corr.values,
                icov=icov, s=icov.sum(), N=m.sum())


SN = {'paper': make_sn('paper'), 'standard': make_sn('standard')}


def chi2_sn(H, Om, shape, cut, DC=None):
    S = SN[cut]
    if DC is None:
        DC = comoving(H, Om, shape)
    dL = (1 + S['zhel']) * np.interp(S['zhd'], ZG, DC)
    r = S['mb'] - 5 * np.log10(dL)
    a = r @ S['icov'] @ r
    b = (S['icov'] @ r).sum()
    return a - b * b / S['s']  # analytic marginalisation over M_B


# ---------- compressed CMB ----------
def _dmstar(H, Om, shape, DC=None):
    if DC is None:
        DC = comoving(H, Om, shape)
    return np.interp(ZSTAR, ZG, DC)


CAL = DMSTAR_OBS / _dmstar(67.36, 0.3153, ('lcdm', None))


def chi2_cmb(H, Om, shape, rd=None, DC=None, rd_consistent=False):
    """D_M(z*) = r*/theta*. If rd_consistent, r* scales with r_d (r* = r_d * 144.43/147.09)
    so the CMB target D_M(z*) scales by r_d/147.09."""
    dm = CAL * _dmstar(H, Om, shape, DC)
    target = DMSTAR_OBS
    if rd_consistent and rd is not None:
        target = DMSTAR_OBS * rd / RD_FID
    h = H / 100
    return ((dm - target) / DMSTAR_ERR) ** 2 + ((Om * h * h - WM_OBS) / WM_ERR) ** 2


def total(H, Om, shape, rd, data, cut='paper', rd_consistent=False):
    DC = comoving(H, Om, shape)
    c = chi2_bao(H, Om, shape, rd, DC)
    if 'SN' in data:
        c += chi2_sn(H, Om, shape, cut, DC)
    if 'CMB' in data:
        c += chi2_cmb(H, Om, shape, rd, DC, rd_consistent)
    return c
