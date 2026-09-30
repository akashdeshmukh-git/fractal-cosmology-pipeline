"""Check that this likelihood reproduces the paper / Cobaya chi2 values at their best-fit points."""
from model import *

checks = [
    ('LCDM, BAO', 69.04, 0.2972, ('lcdm', None), 'BAO', 10.283),
    ('LCDM, BAO+SN', 68.65, 0.3045, ('lcdm', None), 'BAO+SN', 1470.098),
    ('exp free, BAO+SN', 62.28, 0.3649, ('exp', 0.566), 'BAO+SN', 1463.958),
    ('LCDM, BAO+SN+CMB', 68.72, 0.3011, ('lcdm', None), 'BAO+SN+CMB', 1471.774),
    ('exp free, BAO+SN+CMB', 62.25, 0.3669, ('exp', 0.542), 'BAO+SN+CMB', 1464.889),
]
print(f'SNe: paper cut {SN["paper"]["N"]}, standard cut {SN["standard"]["N"]}')
for name, H, Om, sh, data, ref in checks:
    c = total(H, Om, sh, RD_FID, data)
    print(f'{name:24s} this code {c:9.3f}   Cobaya {ref:9.3f}   diff {c - ref:+.3f}')
print('(CMB rows differ by ~0.02 only because the reference best-fit parameters are rounded; '
      'the minimiser in fit_cmb.py gives 1471.788 and 1464.895.)')
