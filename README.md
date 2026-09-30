# Fractal Hubble Model — Analysis Package

Complete analysis package for the paper:  
**"Redshift Evolution of the Fractal Cosmic Web Correction to the Hubble Constant"**

---

## 🔬 What This Is
This package runs three distinct analyses:
1. **LCDM Benchmark** – Verifies the likelihood implementation against DESI's published results.
2. **Joint MCMC** – Evaluates the fractal model vs DESI DR1 + Pantheon+ with full covariance matrices.
3. **CMB Check** – Confirms the fractal correction is negligible at $z=1100$.

---

## ✅ Independent checks (DESI DR2)
Two independent checks of the paper's DESI DR2 analysis live in `checks/`:
- [`checks/cobaya_reproduction/`](checks/cobaya_reproduction/) – the main DR2 fits re-run with [Cobaya](https://cobaya.readthedocs.io) and its official DESI DR2 BAO and Pantheon+ likelihoods. All 15 χ² minima of the paper agree to within 0.004. The check also found a half-bin offset in the grid quantiles (z_char = 0.52 → 0.54) and shows how the standard Pantheon+ cut (z > 0.01) changes each Δχ².
- [`checks/rd_robustness/`](checks/rd_robustness/) – tests whether the result depends on fixing the sound horizon r_d. The preference for the fractal model (Δχ² ≈ −6 to −7 with BAO+SN+CMB) survives with r_d fixed, free, or computed from ω_m.

Each folder has its own README with results, limitations and instructions to rerun.

---

## 🛠️ Requirements & Installation
- **Python 3.13** (or later)

### Install Dependencies (One-time setup):
```bash
pip install emcee corner numpy scipy pandas matplotlib
```
*Alternatively, running `python run_all.py` will automatically attempt to install missing packages.*

---

## 📁 Repository Files
- `run_all.py` – Runs all three analysis steps sequentially.
- `step1_lcdm_benchmark.py` – LCDM implementation verification pipeline.
- `step2_joint_mcmc.py` – Full 3-parameter joint MCMC script.
- `step3_cmb_check.py` – CMB consistency check pipeline.
- `/data/` – Local folder structure for official data inputs.
- `checks/` – Independent DESI DR2 checks (Cobaya reproduction, sound-horizon robustness).

---

## 🚀 How To Run

### Option A: Run everything at once
```bash
python run_all.py
```

### Option B: Run individual steps (Recommended)
```bash
python step1_lcdm_benchmark.py   # Runtime: ~5 min
python step2_joint_mcmc.py       # Runtime: ~20-40 min
python step3_cmb_check.py         # Runtime: ~1 min
```

---

## 📊 Outputs
All execution results save directly to an automatically generated `/results/` folder:
- `step1_lcdm_benchmark.json` – Benchmark numerical files.
- `step2_joint_mcmc.json` – Full MCMC posterior constraints.
- `step3_cmb_check.json` – CMB verification data table.
- `corner_joint_full.png` – Multi-parameter joint contour corner plot.
- `zchar_posterior_full.png` – Parameter posterior graph with theoretical bounds.

---

## 🛰️ Data Sources
- **DESI DR1/DR2 BAO:** [CobayaSampler/bao_data](https://github.com/CobayaSampler/bao_data) *(official DESI mean vectors and covariances)*
- **Pantheon+:** [PantheonPlusSH0ES/DataRelease](https://github.com/PantheonPlusSH0ES/DataRelease) *(official release)*

---

## 📈 What To Expect

### Step 1 Recovery:
- $H_0 = 68.5 \pm 0.9$ *(DESI published: $68.52 \pm 0.62$)*
- $\Omega_m = 0.294 \pm 0.015$ *(DESI published: $0.2941 \pm 0.0095$)*

### Step 2 Outcomes:
- If $z_{\text{char}}$ lands in $[0.3, 0.7]$ $\rightarrow$ Consistent with theoretical prediction ✓
- If $z_{\text{char}}$ lands outside $[0.3, 0.7]$ $\rightarrow$ Discrepancy to investigate
- If $z_{\text{char}}$ is unconstrained $\rightarrow$ More data points required

### Step 3 Confirmation:
- $f(z=1100) < 10^{-1500}$ $\rightarrow$ Correction completely negligible at CMB scales.

---

## 📝 Performance Notes
- **Memory Overhead:** The Pantheon+ covariance matrix takes roughly 32MB. Loading it into memory takes $\sim 30$ seconds.
- **Progress Tracking:** The active MCMC progress bar displays real-time estimated remaining time.
- **Hardware Warning:** Run this on a laptop plugged into wall power; joint MCMC execution is CPU-intensive for $\sim 30$ minutes.
- **Optimization:** If MCMC routines run too slowly, reduce `nsteps` from `6000` down to `3000` inside `step2_joint_mcmc.py`.

