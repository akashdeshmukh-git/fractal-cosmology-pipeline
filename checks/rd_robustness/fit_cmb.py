import json, sys
from fit import *
cut = sys.argv[1]
models = [('lcdm', None, 'LCDM'), ('exp', 1 / 2.262, 'exp_fix'), ('pl', 2.262, 'pl_fix'),
          ('exp', None, 'exp_free'), ('pl', None, 'pl_free')]
out = {}
for rdmode in ['fixed', 'free', 'free_consistent']:
    for kind, fs, name in models:
        r = fit(kind, 'BAO+SN+CMB', rdmode, cut, fixed_shape=fs)
        out[f'{name}|{rdmode}'] = r
        print(f'{cut:8s} {rdmode:16s} {name:9s} chi2={r["chi2"]:.3f} H_E={r["H_E"]:.2f} Om={r["Om"]:.4f} '
              f'rd={r["rd"]:.2f} shape={r["shape"]} Hloc={r["H_local"]:.2f}', flush=True)
json.dump(out, open(f'cmb_{cut}.json', 'w'), indent=1)
