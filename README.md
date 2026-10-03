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
- **DESI DR1:** [://github.com](https://://github.com) *(Official DESI release)*
- **Pantheon+:** [://github.com](https://://github.com) *(Official SH0ES release)*

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


---

## Independent Checks (DESI DR2)

The `checks/` directory contains three independent verifications of the DESI DR2 analysis reported in the paper. They are self-contained and do not import or modify the DR1 pipeline above.

- **[`checks/cobaya_reproduction/`](checks/cobaya_reproduction/)** — Reproduction of the main DR2 fits with [Cobaya](https://cobaya.readthedocs.io), using its official DESI DR2 BAO and Pantheon+ likelihoods. All fifteen χ² minima reported in the paper are recovered to within 0.004.
- **[`checks/rd_robustness/`](checks/rd_robustness/)** — Sensitivity of the results to the treatment of the BAO sound horizon r_d (fixed, free, or derived from ω_m).
- **[`checks/sn_cut_robustness/`](checks/sn_cut_robustness/)** — Sensitivity of the results to the minimum supernova redshift (z > 0.01, 0.02, 0.03), testing for contamination by peculiar velocities.

Each directory includes its own README describing the method, results, limitations, and instructions for reproduction.
