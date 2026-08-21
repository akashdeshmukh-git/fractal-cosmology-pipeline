"""
Step 1: LCDM Benchmark
Verifies our likelihood implementation against DESI's published results.
If H0 and Omega_m match DESI's published values, the implementation is correct.
"""

import numpy as np
from scipy.linalg import cho_factor, cho_solve
import emcee, json

c_light = 299792.458
rs_fixed = 147.09
z_grid = np.linspace(0, 3.5, 3000)
dz = np.diff(z_grid)

# Load DESI data
rows = []
with open('data/desi_2024_gaussian_bao_ALL_GCcomb_mean.txt') as f:
    for line in f:
        if line.startswith('#') or not line.strip(): continue
        z, val, qty = line.split()
        rows.append((float(z), float(val), qty.strip()))

z_bao = np.array([r[0] for r in rows])
obs_bao = np.array([r[1] for r in rows])
qty_bao = np.array([r[2] for r in rows])
cov_bao = np.loadtxt('data/desi_2024_gaussian_bao_ALL_GCcomb_cov.txt')
chol_bao = cho_factor(cov_bao)

def chi_fast(z_t, H0, Om, zc=None):
    OL = 1 - Om
    if zc is not None:
        f = np.exp(-z_grid / zc)
        E = np.sqrt(Om*(1+z_grid)**3*(1+f/2) + OL)
    else:
        E = np.sqrt(Om*(1+z_grid)**3 + OL)
    ig = c_light / (H0 * E)
    cum = np.concatenate([[0], np.cumsum(0.5*(ig[:-1]+ig[1:])*dz)])
    return np.interp(z_t, z_grid, cum)

def Hz_arr(z_t, H0, Om, zc=None):
    OL = 1 - Om
    if zc is not None:
        f = np.exp(-z_t / zc)
        return H0 * np.sqrt(Om*(1+z_t)**3*(1+f/2) + OL)
    return H0 * np.sqrt(Om*(1+z_t)**3 + OL)

def bao_pred(H0, Om, zc=None):
    DM = chi_fast(z_bao, H0, Om, zc)
    Hz = Hz_arr(z_bao, H0, Om, zc)
    DH = c_light / Hz
    DV = (z_bao * DM**2 * DH)**(1/3)
    return np.where(qty_bao=='DM_over_rs', DM/rs_fixed,
           np.where(qty_bao=='DH_over_rs', DH/rs_fixed, DV/rs_fixed))

def c2_bao(H0, Om, zc=None):
    res = obs_bao - bao_pred(H0, Om, zc)
    return float(res @ cho_solve(chol_bao, res))

def log_posterior_lcdm(theta):
    H0, Om = theta
    if not (55 < H0 < 85) or not (0.1 < Om < 0.6): return -np.inf
    return -0.5 * c2_bao(H0, Om)

print("Running LCDM MCMC (32 walkers, 6000 steps)...")
nw = 32
p0 = np.array([68.5, 0.295]) + np.array([0.5, 0.01]) * np.random.randn(nw, 2)
sampler = emcee.EnsembleSampler(nw, 2, log_posterior_lcdm)
sampler.run_mcmc(p0, 6000, progress=True)
chain = sampler.get_chain(discard=2000, flat=True)

H0_med = np.median(chain[:,0]); Om_med = np.median(chain[:,1])
H0_err = np.std(chain[:,0]);    Om_err = np.std(chain[:,1])

print(f"\nOur LCDM:  H0 = {H0_med:.2f} ± {H0_err:.2f}")
print(f"DESI pub:  H0 = 68.52 ± 0.62")
print(f"Our LCDM:  Om = {Om_med:.4f} ± {Om_err:.4f}")
print(f"DESI pub:  Om = 0.2941 ± 0.0095")
H0_match = abs(H0_med - 68.52) < 1.0
Om_match = abs(Om_med - 0.2941) < 0.015
print(f"H0 match:  {'PASS ✓' if H0_match else 'FAIL ✗'}")
print(f"Om match:  {'PASS ✓' if Om_match else 'FAIL ✗'}")

results = {
    'our_H0': float(H0_med), 'our_H0_err': float(H0_err),
    'our_Om': float(Om_med), 'our_Om_err': float(Om_err),
    'desi_H0': 68.52, 'desi_Om': 0.2941,
    'H0_match': bool(H0_match), 'Om_match': bool(Om_match),
    'chi2_at_median': float(c2_bao(H0_med, Om_med))
}
# Save for use by step 2
np.save('results/lcdm_best_fit.npy', np.array([H0_med, Om_med]))
with open('results/step1_lcdm_benchmark.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Saved: results/step1_lcdm_benchmark.json")
