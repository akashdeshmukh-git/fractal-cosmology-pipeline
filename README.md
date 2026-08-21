FRACTAL HUBBLE MODEL — ANALYSIS PACKAGE
=========================================

WHAT THIS IS
------------
Complete analysis package for the paper:
"Redshift Evolution of the Fractal Cosmic Web Correction to the Hubble Constant"

This runs three analyses:
  1. LCDM benchmark — verifies the likelihood implementation against DESI's published results
  2. Joint MCMC     — fractal model vs DESI DR1 + Pantheon+ with full covariance matrices
  3. CMB check      — confirms the fractal correction is negligible at z=1100

REQUIREMENTS
------------
Python 3.13 

Install dependencies (one time only):
  pip install emcee corner numpy scipy pandas matplotlib

Or just run run_all.py — it installs missing packages automatically.

FILES
-----
run_all.py              — runs all three steps in sequence
step1_lcdm_benchmark.py — LCDM implementation verification
step2_joint_mcmc.py     — full 3-parameter joint MCMC (~20-40 min)
step3_cmb_check.py      — CMB consistency check (fast, ~1 min)
data/                   — all official data files included

HOW TO RUN
----------
Option A — run everything at once:
  python run_all.py

Option B — run steps individually (recommended):
  python step1_lcdm_benchmark.py    # ~5 min
  python step2_joint_mcmc.py        # ~20-40 min
  python step3_cmb_check.py         # ~1 min

OUTPUTS
-------
All results saved to results/ folder (created automatically):
  step1_lcdm_benchmark.json   — benchmark numbers
  step2_joint_mcmc.json       — full MCMC posteriors
  step3_cmb_check.json        — CMB check table
  corner_joint_full.png       — corner plot of all 3 parameters
  zchar_posterior_full.png    — z_char posterior with theoretical bounds

DATA SOURCES
------------
DESI DR1:   github.com/CobayaSampler/bao_data  (official DESI release)
Pantheon+:  github.com/PantheonPlusSH0ES/DataRelease  (official SH0ES release)

WHAT TO EXPECT
--------------
Step 1 should recover:
  H0 ≈ 68.5 ± 0.9  (DESI published: 68.52 ± 0.62)
  Om ≈ 0.294 ± 0.015  (DESI published: 0.2941 ± 0.0095)

Step 2 — three possible outcomes:
  If z_char lands in [0.3, 0.7]   → consistent with theoretical prediction ✓
  If z_char lands outside [0.3, 0.7] → discrepancy to investigate
  If z_char is unconstrained       → more data needed

Step 3 should confirm:
  f(z=1100) < 10^-1500 → correction completely negligible at CMB

NOTES
-----
- The Pantheon+ covariance matrix is 32MB. Loading it takes ~30 seconds.
- The MCMC progress bar shows estimated time remaining.
- Run on a laptop plugged in — CPU-intensive for ~30 minutes.
- If MCMC is too slow, reduce nsteps from 6000 to 3000 in step2_joint_mcmc.py
