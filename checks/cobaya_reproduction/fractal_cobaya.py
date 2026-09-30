"""
Independent Cobaya implementation of the fractal-cosmic-web background model
(Deshmukh, Annalen der Physik ms 9212113, revised).

Written from the equations in the paper, NOT imported from the paper's own code
(fraccosmo.py), so that this is a genuine cross-check of the numerics.

Model (flat, radiation included):
    E^2(z) = Om(1+z)^3 + Or(1+z)^4 + OL + X(z),   OL = 1 - Om - Or   (flatness imposed without X,
                                                                        as in the paper)
    X(z)   = Om (1+z)^3 f(z) / 2
    f(z)   = exp(-z/z_char)      (kind='exp')      or   (1+z)^(-n)   (kind='pl')
    H(z)   = H_E E(z),  H_E is the sampled parameter 'H0'
    Or     = 2.469e-5 (1 + 0.2271 * 3.044) / h_E^2
r_d is held fixed at 147.09 Mpc (Planck 2018), exactly as in the paper.

Components
    FractalBackground : cobaya Theory providing Hubble, angular_diameter_distance,
                        comoving_radial_distance, and derived rdrag, omegamh2, H0_local
    CompressedCMB     : D_M(z*) = r*/theta* and omega_m, with the paper's calibration
    PantheonPlusPaperCut : Cobaya Pantheon+ with the paper's sample selection
                           (all non-calibrator SNe, 1624) instead of Cobaya's z_HD > 0.01 (1590)
"""
import os
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_simpson, simpson
from scipy.interpolate import CubicSpline
from cobaya.theory import Theory
from cobaya.likelihood import Likelihood
from cobaya.likelihoods.sn.pantheonplus import PantheonPlus

C_KMS = 299792.458
RD_FIXED = 147.09
OMEGA_R_H2 = 2.469e-5 * (1 + 0.2271 * 3.044)
ZSTAR = 1089.92


def e2(z, H0, Om, family, kind, shape):
    h = H0 / 100.0
    Or = OMEGA_R_H2 / h ** 2
    OL = 1.0 - Om - Or
    out = Om * (1 + z) ** 3 + Or * (1 + z) ** 4 + OL
    if family == 'lcdm':
        return out
    f = np.exp(-z / shape) if kind == 'exp' else (1 + z) ** (-shape)
    return out + 0.5 * Om * (1 + z) ** 3 * f


def comoving_distance_grid(H0, Om, family, kind, shape, zmax=4.0, n=16001):
    """D_C(z) on a uniform grid 0..zmax with Simpson cumulative integration."""
    z = np.linspace(0.0, zmax, n)
    ig = 1.0 / np.sqrt(e2(z, H0, Om, family, kind, shape))
    dc = np.concatenate([[0.0], cumulative_simpson(ig, x=z)]) * C_KMS / H0
    return z, dc


def comoving_distance_star(H0, Om, family, kind, shape, zstar=ZSTAR, n=40001):
    """D_C(z*) integrating in x = ln(1+z)."""
    x = np.linspace(0.0, np.log1p(zstar), n)
    z = np.expm1(x)
    ig = (1 + z) / np.sqrt(e2(z, H0, Om, family, kind, shape))
    return simpson(ig, x=x) * C_KMS / H0


class FractalBackground(Theory):
    family: str = 'lcdm'       # 'lcdm' or 'A'
    kind: str = 'exp'          # 'exp' or 'pl'
    shape_fixed: float = None  # fix z_char (exp) or n (pl); None -> sampled
    # sampled shape parameter names: 'lnzc' (exp) or 'n' (pl)

    def initialize(self):
        self._z_req = set()

    def get_requirements(self):
        req = {'H0': None, 'omegam': None}
        if self.family != 'lcdm' and self.shape_fixed is None:
            req['lnzc' if self.kind == 'exp' else 'n'] = None
        return req

    def must_provide(self, **requirements):
        for k, v in requirements.items():
            if isinstance(v, dict) and 'z' in v:
                self._z_req.update(np.atleast_1d(v['z']).tolist())
        return {}

    def get_can_provide(self):
        return ['Hubble', 'angular_diameter_distance', 'comoving_radial_distance']

    def get_can_provide_params(self):
        return ['rdrag', 'omegamh2', 'H0_local', 'DMstar']

    def _shape(self, params):
        if self.family == 'lcdm':
            return None
        if self.shape_fixed is not None:
            return float(self.shape_fixed)
        return float(np.exp(params['lnzc'])) if self.kind == 'exp' else float(params['n'])

    def calculate(self, state, want_derived=True, **params):
        H0, Om = params['H0'], params['omegam']
        shape = self._shape(params)
        z, dc = comoving_distance_grid(H0, Om, self.family, self.kind, shape)
        state['spline_dc'] = CubicSpline(z, dc)
        state['H0'], state['Om'], state['shape'] = H0, Om, shape
        dstar = comoving_distance_star(H0, Om, self.family, self.kind, shape)
        state['dc_star'] = dstar
        if want_derived:
            e0 = e2(0.0, H0, Om, self.family, self.kind, shape)
            state['derived'] = {'rdrag': RD_FIXED, 'omegamh2': Om * (H0 / 100) ** 2,
                                'H0_local': H0 * np.sqrt(e0), 'DMstar': dstar}

    def _dc(self, z):
        z = np.atleast_1d(np.asarray(z, float))
        s = self.current_state
        out = np.empty_like(z)
        lo = z <= 4.0
        out[lo] = s['spline_dc'](z[lo])
        for i in np.where(~lo)[0]:
            out[i] = comoving_distance_star(s['H0'], s['Om'], self.family, self.kind, s['shape'], zstar=z[i])
        return out

    def get_Hubble(self, z, units='km/s/Mpc'):
        s = self.current_state
        z = np.atleast_1d(np.asarray(z, float))
        H = s['H0'] * np.sqrt(e2(z, s['H0'], s['Om'], self.family, self.kind, s['shape']))
        if units == '1/Mpc':
            return H / C_KMS
        return H

    def get_comoving_radial_distance(self, z):
        return self._dc(z)

    def get_angular_diameter_distance(self, z):
        z = np.atleast_1d(np.asarray(z, float))
        return self._dc(z) / (1 + z)


class CompressedCMB(Likelihood):
    """Planck 2018 compressed: D_M(z*) = r*/theta* and omega_m = Om h^2 (h = H_E/100).
    Same numbers and same LCDM calibration as the paper: the model D_M(z*) is rescaled so that
    the Planck CMB-only LCDM best fit (H0=67.36, Om=0.3153) reproduces r*/theta* exactly."""
    theta_star: float = 0.0104110
    r_star: float = 144.43
    sig_rstar: float = 0.26
    omh2: float = 0.1430
    sig_omh2: float = 0.0011

    def initialize(self):
        self.dm_obs = self.r_star / self.theta_star
        self.dm_sig = self.dm_obs * np.hypot(0.00031 / 1.04110, self.sig_rstar / self.r_star)
        self.cal = self.dm_obs / comoving_distance_star(67.36, 0.3153, 'lcdm', None, None)

    def get_requirements(self):
        return {'DMstar': None, 'omegamh2': None}

    def logp(self, **params_values):
        dm = self.cal * self.provider.get_param('DMstar')
        om = self.provider.get_param('omegamh2')
        chi2 = ((dm - self.dm_obs) / self.dm_sig) ** 2 + ((om - self.omh2) / self.sig_omh2) ** 2
        return -0.5 * chi2


class PantheonPlusPaperCut(PantheonPlus):
    """Cobaya's Pantheon+ likelihood, but with the paper's selection: every SN with
    IS_CALIBRATOR == 0 (1624 SNe), instead of Cobaya's default z_HD > 0.01 (1590 SNe)."""

    def configure(self):
        data_file = os.path.join(self.path, 'PantheonPlus', 'Pantheon+SH0ES.dat')
        d = pd.read_csv(data_file, sep=r'\s+')
        mask = (d['IS_CALIBRATOR'] == 0).values
        assert len(mask) == len(self.mag)
        self._apply_mask(zmask=mask)
        self.pre_vars = 0.0
