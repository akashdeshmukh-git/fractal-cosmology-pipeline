"""r_d computed from omega_m (early universe unchanged in the fractal model), omega_b from BBN.
Aubourg et al. 2015 fitting formula, normalised so rd(0.1430) reproduces 147.09.
CMB D_M(z*) target scales with r* = r_d/1.0184 (theta* is what the CMB measures)."""
import json, sys, warnings
import numpy as np
from scipy.optimize import minimize
warnings.filterwarnings('ignore')
from model import *

WB, WNU = 0.02237, 0.0006


def rd_raw(wm):
    return 55.154 * np.exp(-72.3 * (WNU + 0.0006) ** 2) / ((wm - WNU) ** 0.25351 * WB ** 0.12807)


NORM = RD_FID / rd_raw(0.1430)


def rd_of(H, Om):
    return NORM * rd_raw(Om * (H / 100) ** 2)


cut = sys.argv[1]
models = [('lcdm', None, 'LCDM'), ('exp', 1 / 2.262, 'exp_fix'), ('pl', 2.262, 'pl_fix'),
          ('exp', None, 'exp_free'), ('pl', None, 'pl_free')]
out = {}
for kind, fs, name in models:
    free_shape = kind != 'lcdm' and fs is None

    def unpack(x):
        H, Om = x[0], x[1]
        if kind == 'lcdm':
            sh = ('lcdm', None)
        elif fs is not None:
            sh = (kind, fs)
        else:
            sh = (kind, np.exp(x[2]) if kind == 'exp' else x[2])
        return H, Om, sh

    def f(x):
        H, Om, sh = unpack(x)
        if not (40 < H < 100 and 0.05 < Om < 0.7):
            return 1e10
        if free_shape:
            v = x[2]
            if kind == 'exp' and not (np.log(0.02) <= v <= np.log(3)):
                return 1e10
            if kind == 'pl' and not (0.1 <= v <= 10):
                return 1e10
        return total(H, Om, sh, rd_of(H, Om), 'BAO+SN+CMB', cut, rd_consistent=True)

    base = [[68.5, 0.30], [63.0, 0.36]]
    if free_shape:
        grid = np.linspace(np.log(0.02), np.log(3), 9) if kind == 'exp' else np.geomspace(0.15, 9, 9)
        base = [b + [g] for b in base for g in grid]
    best = None
    for x0 in base:
        r = minimize(f, x0, method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-8, maxiter=20000))
        r = minimize(f, r.x, method='Nelder-Mead', options=dict(xatol=1e-8, fatol=1e-9, maxiter=20000))
        if best is None or r.fun < best.fun:
            best = r
    H, Om, sh = unpack(best.x)
    res = dict(chi2=float(best.fun), H_E=float(H), Om=float(Om), rd=float(rd_of(H, Om)),
               wm=float(Om * (H / 100) ** 2), shape=sh[1],
               H_local=float(H * np.sqrt(1 + (Om / 2 if kind != 'lcdm' else 0))))
    out[name] = res
    print(f'{cut:8s} computed {name:9s} chi2={res["chi2"]:.3f} H_E={H:.2f} Om={Om:.4f} '
          f'wm={res["wm"]:.4f} rd={res["rd"]:.2f} shape={sh[1]} Hloc={res["H_local"]:.2f}', flush=True)
json.dump(out, open(f'computed_{cut}.json', 'w'), indent=1)
