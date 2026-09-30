import json, warnings
import numpy as np
from scipy.optimize import minimize
warnings.filterwarnings('ignore')
from model import *

LZ = (np.log(0.02), np.log(3.0))
NB = (0.1, 10.0)


def fit(kind, data, rdmode, cut='paper', fixed_shape=None, starts=None):
    """rdmode: 'fixed' | 'free' | 'free_consistent'.  Returns best chi2 and params."""
    free_shape = kind in ('exp', 'pl') and fixed_shape is None

    def unpack(x):
        i = 0
        H, Om = x[0], x[1]; i = 2
        rd = RD_FID
        if rdmode != 'fixed':
            rd = x[i]; i += 1
        if kind == 'lcdm':
            shape = ('lcdm', None)
        elif fixed_shape is not None:
            shape = (kind, fixed_shape)
        else:
            shape = (kind, np.exp(x[i]) if kind == 'exp' else x[i])
        return H, Om, rd, shape

    def f(x):
        H, Om, rd, shape = unpack(x)
        if not (40 < H < 100 and 0.05 < Om < 0.7 and 100 < rd < 200):
            return 1e10
        if free_shape:
            v = x[-1]
            if kind == 'exp' and not (LZ[0] <= v <= LZ[1]):
                return 1e10
            if kind == 'pl' and not (NB[0] <= v <= NB[1]):
                return 1e10
        return total(H, Om, shape, rd, data, cut, rd_consistent=(rdmode == 'free_consistent'))

    base = [[68.5, 0.30], [63.0, 0.36]]
    if rdmode != 'fixed':
        base = [b + [RD_FID] for b in base] + [b + [140.0] for b in base] + [b + [155.0] for b in base]
    if free_shape:
        grid = np.linspace(*LZ, 9) if kind == 'exp' else np.geomspace(0.15, 9, 9)
        base = [b + [g] for b in base for g in grid]
    if starts:
        base = starts + base
    best = None
    for x0 in base:
        r = minimize(f, x0, method='Nelder-Mead',
                     options=dict(xatol=1e-6, fatol=1e-7, maxiter=20000, maxfev=20000))
        r = minimize(f, r.x, method='Nelder-Mead',
                     options=dict(xatol=1e-7, fatol=1e-8, maxiter=20000, maxfev=20000))
        if best is None or r.fun < best.fun:
            best = r
    H, Om, rd, shape = unpack(best.x)
    return dict(chi2=float(best.fun), H_E=float(H), Om=float(Om), rd=float(rd),
                shape=None if shape[1] is None else float(shape[1]),
                H_local=float(H * np.sqrt(1 + (Om / 2 if kind != 'lcdm' else 0))),
                hrd=float(H / 100 * rd))


if __name__ == '__main__':
    import sys
    cut = sys.argv[1] if len(sys.argv) > 1 else 'paper'
    models = [('lcdm', None, 'LCDM'), ('exp', 1 / 2.262, 'exp_fix'), ('pl', 2.262, 'pl_fix'),
              ('exp', None, 'exp_free'), ('pl', None, 'pl_free')]
    out = {}
    for data in ['BAO', 'BAO+SN', 'BAO+SN+CMB']:
        modes = ['fixed', 'free'] + (['free_consistent'] if 'CMB' in data else [])
        for rdmode in modes:
            for kind, fs, name in models:
                r = fit(kind, data, rdmode, cut, fixed_shape=fs)
                out[f'{name}|{data}|{rdmode}'] = r
                print(f'{cut:8s} {data:11s} {rdmode:16s} {name:9s} chi2={r["chi2"]:.3f} '
                      f'H_E={r["H_E"]:.2f} Om={r["Om"]:.4f} rd={r["rd"]:.2f} hrd={r["hrd"]:.2f} '
                      f'shape={r["shape"]} Hloc={r["H_local"]:.2f}', flush=True)
    json.dump(out, open(f'results_{cut}.json', 'w'), indent=1)
