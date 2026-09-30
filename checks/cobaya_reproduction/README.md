# Cobaya reproduction of "Redshift Evolution of the Fractal Cosmic Web Correction to the Hubble Constant" (Annalen der Physik ms 9212113, revision R2)

This folder re-runs the DESI DR2 analysis of the paper with [Cobaya](https://cobaya.readthedocs.io) 3.6.2. It uses Cobaya's own DESI DR2 BAO and Pantheon+ likelihoods, a freshly written background theory for the fractal model and a freshly written compressed-CMB likelihood. None of the paper's code (`fraccosmo.py`) is imported, so every number below is an independent recomputation.

**Short version.**
- ΛCDM on DESI DR2 BAO reproduces the published DESI numbers.
- All 15 χ²_min values in Table 12 of the paper reproduce to within 0.004.
- Two things in the paper should be corrected before publication:
  1. The posterior medians and intervals quoted from the grid integration are biased low by half a grid step.
  2. The supernova sample is not the standard Pantheon+ cosmology cut, and the standard cut weakens every fractal-model preference by 0.4 to 0.9 in Δχ².

Both are described in detail below.

## 1. What was matched to the paper

| Choice | Paper | This run |
|---|---|---|
| BAO data | DESI DR2, 13 points, full covariance | `bao.desi_dr2.desi_bao_all` (Cobaya). Files are byte-identical to the paper's |
| SN data | Pantheon+, all 1624 non-calibrator SNe, STAT+SYS covariance, analytic M_B marginalisation | Cobaya `sn.pantheonplus` machinery with a subclass that applies the paper's cut (`PantheonPlusPaperCut`), plus Cobaya's default cut for comparison (Section 5) |
| CMB | compressed: D_M(z*) = r*/θ* = 13873 ± 25 Mpc and ω_m = 0.1430 ± 0.0011, model D_M(z*) calibrated so that the Planck ΛCDM best fit (67.36, 0.3153) reproduces r*/θ* | same numbers and same calibration (`CompressedCMB`) |
| r_d | fixed at 147.09 Mpc | fixed at 147.09 Mpc (`rdrag` provided as a constant) |
| Background | E² = Ωm(1+z)³ + Ωr(1+z)⁴ + ΩΛ + Ωm(1+z)³f/2, flatness imposed without the fractal term, Ωr h² = 2.469e-5 (1 + 0.2271·3.044), h = H_E/100 | same (`FractalBackground`) |
| Priors | H_E ∈ [50, 90], Ωm ∈ [0.1, 0.6], ln z_char ∈ [ln 0.02, ln 3], n ∈ [0.1, 10] | same |
| Best fits | Nelder–Mead from the best MCMC sample, checked against profiles | Cobaya `minimize` (BOBYQA), 16 random starts spread over the whole shape range for the free-shape models, 6 for the others |

The numerical integration differs on purpose: Simpson on a 16 001-point grid here, trapezoid with Δz = 10⁻³ in the paper. At four test points the two codes agree to better than 10⁻³ in every χ² term (`unit_check.py`).

## 2. ΛCDM validation on DESI DR2 BAO

| | this run (MCMC) | DESI DR2 published | paper (Table 3) |
|---|---|---|---|
| Ωm | 0.2973 ± 0.0084 | 0.2975 ± 0.0086 | 0.2976 ± 0.0086 |
| h r_d [Mpc] | 101.54 ± 0.72 | 101.54 ± 0.73 | 101.51 ± 0.73 |
| H_E (r_d = 147.09) | 69.03 ± 0.49 | – | 69.01 ± 0.50 |
| χ²_min (13 points) | 10.283 | – | 10.28 |

Best fit from the minimizer: Ωm = 0.2972, h r_d = 101.549. Plot: `plots/lcdm_validation.png` (the DESI marker sits under the Cobaya one because the two coincide).

## 3. Table 12: exponential versus power-law formation function

χ²_min, paper against Cobaya, with the paper's SN cut. Δχ² is relative to ΛCDM on the same data.

| Model | Data | Paper χ²_min | Cobaya χ²_min | Difference | Paper Δχ² | Cobaya Δχ² |
|---|---|---|---|---|---|---|
| ΛCDM | BAO | 10.28 | 10.283 | +0.003 | – | – |
| ΛCDM | BAO+SN | 1470.10 | 1470.098 | −0.002 | – | – |
| ΛCDM | BAO+SN+CMB | 1471.77 | 1471.774 | +0.004 | – | – |
| exp, z_char = 1/2.262 | BAO | 9.58 | 9.577 | −0.003 | −0.70 | −0.707 |
| exp, z_char = 1/2.262 | BAO+SN | 1466.19 | 1466.186 | −0.004 | −3.91 | −3.911 |
| exp, z_char = 1/2.262 | BAO+SN+CMB | 1467.88 | 1467.881 | +0.001 | −3.89 | −3.893 |
| power law, n = 2.262 | BAO | 9.29 | 9.289 | −0.001 | −0.99 | −0.994 |
| power law, n = 2.262 | BAO+SN | 1465.80 | 1465.797 | −0.003 | −4.30 | −4.301 |
| power law, n = 2.262 | BAO+SN+CMB | 1471.11 | 1471.114 | +0.004 | −0.66 | −0.661 |
| exp, z_char free | BAO | 8.30 | 8.300 | −0.000 | −1.98 | −1.984 |
| exp, z_char free | BAO+SN | 1463.96 | 1463.958 | −0.002 | −6.14 | −6.139 |
| exp, z_char free | BAO+SN+CMB | 1464.89 | 1464.889 | −0.001 | −6.88 | −6.886 |
| power law, n free | BAO | 9.23 | 9.234 | +0.004 | −1.05 | −1.049 |
| power law, n free | BAO+SN | 1465.12 | 1465.121 | +0.001 | −4.98 | −4.977 |
| power law, n free | BAO+SN+CMB | 1469.50 | 1469.504 | +0.004 | −2.27 | −2.270 |

Every row agrees within 0.004, which is below the paper's two-decimal rounding. The paper's main qualitative statements hold up in Cobaya:
- The two fixed forms fit BAO+SN equally well.
- The CMB separates them: exp −3.89, power law −0.66.
- The free exponential is the best model.

Best-fit shapes: z_char = 0.585 (BAO), 0.566 (BAO+SN), 0.542 (BAO+SN+CMB). n = 2.06, 1.81, 2.56. H_local = H_E√(1+X(0)) at the CMB-anchored exponential best fit is 67.72, as in the paper, so the Hubble-tension conclusion is unchanged.

Full table with the standard-cut column and best-fit parameters: `results/comparison_chi2.csv`. Plot: `plots/delta_chi2.png`.

### Bimodality

`plots/profiles.png` shows profile likelihoods from the Cobaya minimizer at fixed z_char (22 values) and fixed n (19 values), each profiled over H_E and Ωm, against the paper's grid profiles.
- **Exponential, BAO alone.** One dip below ΛCDM at z_char ≈ 0.57 (Δχ² = −2.0). The profile rises to +4 near z_char ≈ 1.9 and falls back towards ΛCDM at the z_char = 3 prior edge.
- **Exponential, BAO+SN.** The main dip deepens to −6.1, and a second, shallow dip opens at the prior edge (−1.3 at z_char = 3). This is the "secondary accumulation" of the paper.
- **Power law, BAO alone.** Nearly flat (−0.4 to −1.0) from n = 0.1 to 2.6, a barrier of +4.7 near n ≈ 8, then a slow return towards ΛCDM.
- **Power law, BAO+SN.** A broad trough (−3.7 to −4.9) for n between 0.3 and 2.3, heading back to ΛCDM as n → 0, where the model becomes ΛCDM with a rescaled matter density. The barrier is about +83 near n ≈ 20 to 25.
- **With the CMB,** both families have a single narrow minimum. The n → 0 end is excluded by Δχ² > 1000, and the n → ∞ side is again behind a barrier of about +85.

This confirms the explanation in Section 7.6 of the paper. The Cobaya profile sits slightly below the paper's grid profile in places because the grid minimises over a finite (H, Ωm) mesh.

## 4. Posteriors: an error in the paper's grid quantiles

Cobaya MCMC (single chain, R−1 < 0.003 on the means) against the numbers quoted in the paper:

| Posterior | Paper (grid) | Cobaya MCMC | Paper grid, recomputed with the edge fix |
|---|---|---|---|
| z_char, BAO+SN+CMB | 0.521 [0.463, 0.586] | 0.541 [0.481, 0.601] | 0.541 [0.480, 0.607] |
| z_char, BAO+SN | 0.566 [0.470, 0.730] | 0.590 [0.492, 0.764] | 0.587 [0.488, 0.757] |
| n, BAO+SN+CMB | 2.525 [2.278, 2.812] | 2.591 [2.347, 2.883] | 2.597 [2.350, 2.883] |
| fraction of z_char posterior in [0.3, 0.7], BAO+SN+CMB | 98% | 99% | – |
| same, BAO+SN | 74% | 76% | – |

The first pass disagreed by about 0.3σ, always in the same direction. The cause is in `grid_post.py` of the paper's code:

```python
cdf = np.cumsum(post_x) / post_x.sum()
q = np.interp(pc / 100, cdf, Xs)
```

The cumulative sum up to bin i includes the whole of bin i, so each CDF value belongs at the bin's upper edge X_i + ΔX/2, not at its centre X_i. Every quantile therefore comes out half a grid step too low: 0.036 in ln z_char and 0.072 in n. Shifting the CDF to the bin edges brings the paper's own grid into agreement with Cobaya to within 0.01 (last column). Your emcee chains, which don't use this code, already gave n = 2.60 and ln z_char = −0.626 (z_char = 0.535).

**What this changes in the paper.**
- z_char = 0.52 ± 0.06 becomes 0.54 ± 0.06. This value appears in:
  - the abstract
  - Section 7.3
  - Section 7.6
  - Section 9
  - the limitations section
  - the conclusions
  - the response letter
  - the graphical abstract label
- n = 2.53 (+0.29 −0.25) becomes 2.60 (+0.29 −0.25). Its distance from the Suhhonenko value 2.262 goes from 1.1σ to about 1.35σ.
- Other z_char or n values computed by the same function also shift, for example the DR2 prior-sensitivity table. The fractions inside [0.3, 0.7] and the χ² values are not affected.

## 5. Supernova sample: the paper's cut is not the standard one

The paper keeps every SN with `IS_CALIBRATOR == 0` (1624). The Pantheon+ cosmology analysis, and Cobaya's default likelihood, keep every SN with z_HD > 0.01 (1590). The two samples differ in two ways:
- The paper's sample includes 44 non-calibrator SNe at z ≤ 0.01. These are usually excluded because peculiar velocities dominate their redshifts.
- It leaves out 10 calibrator SNe at z > 0.01, which the standard sample uses as ordinary Hubble-flow SNe.

With the standard cut:

| Δχ² vs ΛCDM | BAO+SN, paper cut | BAO+SN, standard cut | BAO+SN+CMB, paper cut | BAO+SN+CMB, standard cut |
|---|---|---|---|---|
| exp, z_char = 1/2.262 | −3.91 | −3.42 | −3.89 | −3.46 |
| power law, n = 2.262 | −4.30 | −3.87 | −0.66 | +0.00 |
| exp, z_char free | −6.14 | −5.29 | −6.89 | −6.01 |
| power law, n free | −4.98 | −4.36 | −2.27 | −1.90 |

The ranking of the models does not change, and the best-fit z_char barely moves (0.542 to 0.535 with the CMB). Every improvement over ΛCDM becomes smaller, though, and the fixed power law with the CMB loses its preference entirely. A referee familiar with Pantheon+ could ask about the 1624 sample. Adopting the standard 1590 cut, or reporting both, would be the safer choice.

## 6. What this check does and does not cover

- **Independent here:** the BAO and SN likelihood code (Cobaya's own), the distance integration, the minimizer, and the sampler.
- **Not independent, and cannot be:** the model itself and the compressed-CMB numbers. These are written from the paper's equations, so a mistake in the physics would be copied, not caught.
- **Not reproduced:**
  - the growth-rate (fσ8) analysis
  - the DESI DR1 section
  - the δQ bounds
  - the next-release prediction
- **Sampling limits:**
  - The MCMC runs are one chain each. Cobaya's R−1 is computed from chain splits, and effective sample sizes are 2800 to 4500.
  - The BAO+SN exponential posterior has a long tail towards large z_char. MCMC captures it (76% inside the window against 74% for the grid) but converges slowly in that tail.
  - The ΛCDM chain is shorter (1736 rows after burn-in). That is enough for the quoted errors, not for the fourth digit.
- **The r_d assumption.** Fixing r_d = 147.09 Mpc follows the paper. Freeing r_d, as DESI does, would absorb part of the H_E difference between the models. That is a modelling choice this reproduction does not test.

## 7. Files

| File | Contents |
|---|---|
| `fractal_cobaya.py` | Cobaya theory (`FractalBackground`), compressed CMB likelihood, Pantheon+ subclass with the paper's cut |
| `common.py`, `runner.py` | parameter and likelihood set-up, minimizer and MCMC wrappers |
| `unit_check.py` | point-by-point comparison with the paper's code |
| `run_lcdm.py`, `run_fits.py` | the runs (`python run_fits.py A` and `B` run the two halves) |
| `analysis.py` | builds the tables and plots |
| `results/` | `min_*.json` for every minimization, `comparison_chi2.csv`, `summary.json` |
| `chains/` | Cobaya outputs (`*.minimum.txt`, MCMC chains) |
| `plots/` | `lcdm_validation.png`, `delta_chi2.png`, `profiles.png`, `posteriors.png` |

## 8. Reproducing

```bash
pip install cobaya getdist py-bobyqa
mkdir -p packages/data && cd packages/data
git clone --depth 1 --branch v2.6 https://github.com/CobayaSampler/bao_data.git && echo v2.6 > bao_data/version.dat
git clone --depth 1 --branch v1.8 https://github.com/CobayaSampler/sn_data.git  && echo v1.8 > sn_data/version.dat
cd ../..
python run_lcdm.py min mcmc
python run_fits.py A & python run_fits.py B
python analysis.py
```

`cobaya-install` would normally fetch the data. Here its tarball download was blocked, so the repositories were cloned at the exact tags Cobaya requires and the `version.dat` markers were written by hand. Cobaya packages are looked for in `packages/` next to the scripts (or set `COBAYA_PACKAGES`). `unit_check.py` and `analysis.py` also read the paper's own code and grid results; point `PAPER_CODE` at that folder. The full set of runs takes about 40 minutes on two cores.
