"""
Step 2: Full 3-parameter Joint MCMC
Fractal model (H0, Omega_m, z_char) against DESI DR1 + Pantheon+
Full covariance matrices for both datasets.
Estimated runtime: 20-40 minutes.
"""

import numpy as np
from scipy.linalg import cho_factor, cho_solve
import emcee, corner, json, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

c_light = 299792.458
rs_fixed = 147.09
z_grid = np.linspace(0, 3.5, 3000)
dz = np.diff(z_grid)

# ── Model functions ───────────────────────────────────────────
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

# ── DESI ─────────────────────────────────────────────────────
print("Loading DESI DR1...")
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

def bao_pred(H0, Om, zc):
    DM = chi_fast(z_bao, H0, Om, zc)
    Hz = Hz_arr(z_bao, H0, Om, zc)
    DH = c_light / Hz
    DV = (z_bao * DM**2 * DH)**(1/3)
    return np.where(qty_bao=='DM_over_rs', DM/rs_fixed,
           np.where(qty_bao=='DH_over_rs', DH/rs_fixed, DV/rs_fixed))

def c2_bao(H0, Om, zc):
    res = obs_bao - bao_pred(H0, Om, zc)
    return float(res @ cho_solve(chol_bao, res))

# ── Pantheon+ ─────────────────────────────────────────────────
print("Loading Pantheon+ (this takes ~30 seconds)...")
pdata = pd.read_csv('data/Pantheon+SH0ES.dat', sep=' ')
mask = (pdata['IS_CALIBRATOR'] == 0).values
z_sn = pdata[mask]['zHD'].values
mu_obs = pdata[mask]['MU_SH0ES'].values
n_sn = len(z_sn)
print(f"  SNe (non-calibrator): {n_sn}")

print("  Loading covariance matrix...")
with open('data/Pantheon+SH0ES_STAT+SYS.cov') as f:
    N = int(f.readline())
    cov_full = np.array(f.read().split(), dtype=float).reshape(N, N)
all_mask = (pd.read_csv('data/Pantheon+SH0ES.dat', sep=' ')['IS_CALIBRATOR'] == 0).values
cov_sn = cov_full[np.ix_(all_mask, all_mask)]
print("  Pre-factoring covariance (done once)...")
chol_sn = cho_factor(cov_sn)
Ci_ones = cho_solve(chol_sn, np.ones(n_sn))
A_sn = float(np.ones(n_sn) @ Ci_ones)
print(f"  Ready. {n_sn}x{n_sn} covariance loaded.")

def c2_sn(H0, Om, zc):
    DM = chi_fast(z_sn, H0, Om, zc)
    mu_th = 5 * np.log10(DM * (1+z_sn)) + 25
    delta = mu_obs - mu_th
    Ci_d = cho_solve(chol_sn, delta)
    B = float(np.ones(n_sn) @ Ci_d)
    res = delta - B / A_sn
    return float(res @ cho_solve(chol_sn, res))

# ── LCDM baseline ─────────────────────────────────────────────
try:
    lcdm_bf = np.load('results/lcdm_best_fit.npy')
    H0_lcdm, Om_lcdm = lcdm_bf[0], lcdm_bf[1]
    print(f"\nUsing LCDM best fit from Step 1: H0={H0_lcdm:.2f}, Om={Om_lcdm:.4f}")
except:
    H0_lcdm, Om_lcdm = 69.27, 0.2947
    print(f"\nUsing default LCDM best fit: H0={H0_lcdm:.2f}, Om={Om_lcdm:.4f}")

c2_bao_lcdm = c2_bao(H0_lcdm, Om_lcdm, 0.001)  # effectively lcdm
c2_sn_lcdm  = c2_sn(H0_lcdm, Om_lcdm, 0.001)
# Actually recompute lcdm properly
def c2_bao_l(H0, Om):
    OL = 1-Om
    E_arr = np.sqrt(Om*(1+z_grid)**3 + OL)
    ig = c_light/(H0*E_arr)
    cum = np.concatenate([[0], np.cumsum(0.5*(ig[:-1]+ig[1:])*dz)])
    DM = np.interp(z_bao, z_grid, cum)
    Hz = H0*np.sqrt(Om*(1+z_bao)**3+OL)
    DH = c_light/Hz; DV = (z_bao*DM**2*DH)**(1/3)
    pred = np.where(qty_bao=='DM_over_rs',DM/rs_fixed,np.where(qty_bao=='DH_over_rs',DH/rs_fixed,DV/rs_fixed))
    res = obs_bao - pred
    return float(res@cho_solve(chol_bao,res))

def c2_sn_l(H0, Om):
    OL=1-Om
    E_arr=np.sqrt(Om*(1+z_grid)**3+OL)
    ig=c_light/(H0*E_arr)
    cum=np.concatenate([[0],np.cumsum(0.5*(ig[:-1]+ig[1:])*dz)])
    DM=np.interp(z_sn,z_grid,cum)
    mu_th=5*np.log10(DM*(1+z_sn))+25
    d=mu_obs-mu_th; Cd=cho_solve(chol_sn,d); B=float(np.ones(n_sn)@Cd)
    res=d-B/A_sn; return float(res@cho_solve(chol_sn,res))

cb_l = c2_bao_l(H0_lcdm, Om_lcdm)
cs_l = c2_sn_l(H0_lcdm, Om_lcdm)
print(f"LCDM chi2: BAO={cb_l:.3f}, SN={cs_l:.3f}, total={cb_l+cs_l:.3f}")

# ── Joint fractal MCMC ────────────────────────────────────────
H0_prior_mu, H0_prior_sig = 67.66, 0.42

def log_posterior(theta):
    H0, Om, zc = theta
    if not (55 < H0 < 85): return -np.inf
    if not (0.1 < Om < 0.6): return -np.inf
    if not (0.02 < zc < 2.0): return -np.inf
    lp = -0.5 * ((H0 - H0_prior_mu) / H0_prior_sig)**2
    return lp - 0.5 * (c2_bao(H0, Om, zc) + c2_sn(H0, Om, zc))

nwalkers = 32
ndim = 3
nsteps = 6000
burnin = 2000

print(f"\nRunning joint fractal MCMC ({nwalkers} walkers, {nsteps} steps)...")
print("Expected time: 20-40 minutes. Progress bar below:")

p0 = np.array([67.66, 0.295, 0.40]) + \
     np.array([0.3, 0.01, 0.05]) * np.random.randn(nwalkers, ndim)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior)
sampler.run_mcmc(p0, nsteps, progress=True)

chain = sampler.get_chain(discard=burnin, flat=True)
labels = ['H₀', 'Ω_m', 'z_char']

print("\n=== RESULTS ===")
medians = {}
for i, label in enumerate(labels):
    med = np.median(chain[:,i])
    lo, hi = np.percentile(chain[:,i], [16, 84])
    medians[label] = med
    print(f"{label}: {med:.4f}  +{hi-med:.4f} / -{med-lo:.4f}  (68% CI: [{lo:.4f}, {hi:.4f}])")

zc_med = medians['z_char']
zc_lo, zc_hi = np.percentile(chain[:,2], [16, 84])
within = np.mean((chain[:,2] > 0.3) & (chain[:,2] < 0.7))
print(f"\nz_char within theoretical bounds [0.3, 0.7]: {within*100:.1f}%")

cb_f = c2_bao(medians['H₀'], medians['Ω_m'], zc_med)
cs_f = c2_sn(medians['H₀'], medians['Ω_m'], zc_med)
print(f"\nchi2 BAO:   fractal={cb_f:.3f},  LCDM={cb_l:.3f},  Δ={cb_l-cb_f:.3f}")
print(f"chi2 SNe:   fractal={cs_f:.3f}, LCDM={cs_l:.3f}, Δ={cs_l-cs_f:.3f}")
print(f"chi2 joint: fractal={cb_f+cs_f:.3f}, LCDM={cb_l+cs_l:.3f}, Δ={cb_l+cs_l-cb_f-cs_f:.3f}")

n_tot = len(obs_bao) + n_sn
k_f, k_l = 3, 2
AIC_f = cb_f+cs_f+2*k_f; AIC_l = cb_l+cs_l+2*k_l
BIC_f = cb_f+cs_f+k_f*np.log(n_tot); BIC_l = cb_l+cs_l+k_l*np.log(n_tot)
print(f"AIC: fractal={AIC_f:.2f}, LCDM={AIC_l:.2f}, ΔAIC={AIC_l-AIC_f:.2f}")
print(f"BIC: fractal={BIC_f:.2f}, LCDM={BIC_l:.2f}, ΔBIC={BIC_l-BIC_f:.2f}")

# Corner plot
fig = corner.corner(chain, labels=labels,
                    quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_fmt='.4f')
fig.suptitle('Fractal model — Joint DESI DR1 + Pantheon+ (full MCMC)', y=1.02)
fig.savefig('results/corner_joint_full.png', dpi=150, bbox_inches='tight')
print("Saved: results/corner_joint_full.png")

# z_char posterior
fig2, ax = plt.subplots(figsize=(7,5))
ax.hist(chain[:,2], bins=80, density=True, color='steelblue', alpha=0.75)
ax.axvline(zc_med, color='black', lw=1.5, label=f'median = {zc_med:.3f}')
ax.axvline(zc_lo, color='black', lw=1, ls='--', alpha=0.6)
ax.axvline(zc_hi, color='black', lw=1, ls='--', alpha=0.6)
ax.axvspan(0.3, 0.7, alpha=0.15, color='green', label='Theoretical bounds [0.3, 0.7]')
ax.set_xlabel('z_char', fontsize=12)
ax.set_ylabel('Posterior density', fontsize=12)
ax.set_title('z_char posterior — Joint DESI + Pantheon+\n(Full 3-parameter MCMC)', fontsize=11)
ax.legend()
fig2.tight_layout()
fig2.savefig('results/zchar_posterior_full.png', dpi=150)
print("Saved: results/zchar_posterior_full.png")

results = {
    'datasets': 'DESI DR1 BAO (12 pts) + Pantheon+ ({} SNe)'.format(n_sn),
    'H0':    {'median': float(medians['H₀']),   '16th': float(np.percentile(chain[:,0],16)), '84th': float(np.percentile(chain[:,0],84))},
    'Om':    {'median': float(medians['Ω_m']),  '16th': float(np.percentile(chain[:,1],16)), '84th': float(np.percentile(chain[:,1],84))},
    'z_char':{'median': float(zc_med),           '16th': float(zc_lo), '84th': float(zc_hi)},
    'fraction_within_bounds': float(within),
    'chi2_BAO_fractal': float(cb_f), 'chi2_BAO_LCDM': float(cb_l), 'delta_chi2_BAO': float(cb_l-cb_f),
    'chi2_SN_fractal':  float(cs_f), 'chi2_SN_LCDM':  float(cs_l), 'delta_chi2_SN':  float(cs_l-cs_f),
    'chi2_joint_fractal': float(cb_f+cs_f), 'chi2_joint_LCDM': float(cb_l+cs_l),
    'delta_chi2_joint': float(cb_l+cs_l-cb_f-cs_f),
    'AIC_fractal': float(AIC_f), 'AIC_lcdm': float(AIC_l), 'delta_AIC': float(AIC_l-AIC_f),
    'BIC_fractal': float(BIC_f), 'BIC_lcdm': float(BIC_l), 'delta_BIC': float(BIC_l-BIC_f),
}
with open('results/step2_joint_mcmc.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Saved: results/step2_joint_mcmc.json")
