# Prediction for the next DESI BAO release

Recorded 9 November 2026, before any DESI release after DR2.

Source: A. Deshmukh, "Redshift Evolution of the Fractal Cosmic Web Correction to the Hubble Constant", *Annalen der Physik* **538**, e70321 (2026), [doi:10.1002/andp.70321](https://doi.org/10.1002/andp.70321). Received 13 August 2026, accepted 15 October 2026. Section 11 and Table 15.

The numbers below are copied from `prediction/table15.json`. That file and `prediction/posterior_samples.npz` were rebuilt with `make_prediction.py` and agree with the published Table 15 to every printed digit. Nothing in `prediction/` is changed after this date.

## What was predicted

The prediction is the posterior-predictive distribution of the fractal model (formation function exp(−z/z_char), z_char free) fitted to DESI DR2 BAO + Pantheon+ with a flat H0 prior, marginalised over H_E, Ωm and z_char. r_d is fixed at 147.09 Mpc. ΛCDM is fitted to the same data in the same way.

The test uses DESI DR2 only. DR2 already contains the DR1 galaxies, so DR1 adds nothing and the two cannot be combined without counting the same galaxies twice. The medians, 68% ranges and ΛCDM values below are exactly those of Table 15. The 95% column here is the DR2-based range. The paper's Table 15 printed a slightly wider 95% range that also covered a DR1-based fit, as a stability check. Both versions are kept in `prediction/table15.json`.

| z_eff | Observable | Fractal median | Fractal 68% | Fractal 95% | ΛCDM median | σ_next | Δ/σ_next |
|---|---|---|---|---|---|---|---|
| 0.295 | D_V/r_d | 7.938 | [7.895, 7.981] | [7.857, 8.024] | 7.931 | 0.059 | +0.1 |
| 0.510 | D_M/r_d | 13.295 | [13.230, 13.360] | [13.171, 13.426] | 13.297 | 0.128 | 0.0 |
| 0.510 | D_H/r_d | 22.323 | [22.222, 22.427] | [22.134, 22.527] | 22.480 | 0.339 | −0.5 |
| 0.706 | D_M/r_d | 17.425 | [17.349, 17.499] | [17.283, 17.574] | 17.454 | 0.140 | −0.2 |
| 0.706 | D_H/r_d | 19.868 | [19.761, 19.971] | [19.660, 20.068] | 19.978 | 0.264 | −0.4 |
| 0.934 | D_M/r_d | 21.666 | [21.578, 21.752] | [21.503, 21.835] | 21.712 | 0.126 | −0.4 |
| 0.934 | D_H/r_d | 17.399 | [17.299, 17.494] | [17.188, 17.581] | 17.434 | 0.154 | −0.2 |
| 1.321 | D_M/r_d | 27.717 | [27.609, 27.825] | [27.518, 27.921] | 27.759 | 0.252 | −0.2 |
| 1.321 | D_H/r_d | 14.039 | [13.959, 14.120] | [13.879, 14.197] | 13.985 | 0.172 | +0.3 |
| 1.484 | D_M/r_d | 29.910 | [29.796, 30.026] | [29.696, 30.129] | 29.941 | 0.580 | −0.1 |
| 1.484 | D_H/r_d | 12.891 | [12.812, 12.969] | [12.742, 13.044] | 12.814 | 0.404 | +0.2 |
| 2.330 | D_M/r_d | 38.860 | [38.712, 39.016] | [38.581, 39.153] | 38.810 | 0.410 | +0.1 |
| 2.330 | D_H/r_d | 8.686 | [8.603, 8.778] | [8.529, 8.867] | 8.587 | 0.079 | +1.3 |

σ_next is the DR2 error scaled by √(3/5), a rough forecast of the next-release precision. Summed over the 13 observables, the fractal and ΛCDM medians differ by Δχ² ≈ 3.2 with the forecast covariance.

The paper's three statements (Section 11):

1. D_H/r_d at z ≃ 0.5 to 0.7 comes out 0.5 to 0.7% below the ΛCDM fit.
2. At z = 2.33, D_H/r_d is 1.2% above ΛCDM. This is the largest single-bin signal, 1.3σ of the forecast error.
3. D_M/r_d stays within 0.3% of ΛCDM over 0.5 < z < 2.4.

The paper's falsification sentence: "Measurements falling outside the 95% ranges in several bins, or D_H/r_d at z = 2.33 below the ΛCDM value while D_H/r_d at z ≃ 0.5 to 0.7 lies above it, would disfavour the model in its present form."

## How the test will be scored

These choices are fixed now, before the new data exist. `check_next_release.py` implements them.

- **The test uses only the new data.** The next release will contain the DR2 galaxies, and the prediction was built from DR2, so part of the release is data the prediction has already seen. Comparing the whole release with the prediction would partly re-use DR2 and favour whichever model fitted DR2 better. The release is therefore split into DR2 plus independent new data. For Gaussian measurements of nested samples, the information (inverse covariance) adds, so the new data alone are x_new = F_new⁻¹ (F_release x_release − F_DR2 x_DR2) with covariance F_new⁻¹, where F_new = F_release − F_DR2. The test is applied to x_new, for every bin that matches a DR2 bin (same observable, redshift within 0.05). If a bin's effective redshift has moved, the DR2 value is carried to the new redshift with the ΛCDM shape. If F_new is not positive definite (the release is not more precise than DR2 in every combination of bins, which would point to an analysis change), the split is not possible and this is reported.
- **The whole release is reported too**, as a secondary result. Its simulations include the DR2 overlap, so its probabilities are honest, but its Δχ² is not an independent test. In mock releases where ΛCDM is true, the whole-release Δχ² has a median of −1.6, in favour of the fractal model, purely from the re-used DR2. The new-data-only Δχ² has a median of +0.5.
- **No refitting.** Both predictions are the frozen distributions above. If the release uses different effective redshifts, the same frozen posterior samples are evaluated at the new redshifts.
- **"Outside the 95% range"** means a measurement is inconsistent with the prediction once its own error is included: |pull| > 1.96, where the pull divides the difference from the predicted median by the measurement error and the prediction spread added in quadrature. Only the DR2-based prediction is used. Comparing the measured value with the prediction band alone, ignoring the measurement error, is not used: under that reading 5 of the 13 DR2 points, the data the prediction was built from, would already fall outside.
- **"Several bins"** means 3 or more, and only counts if that many is unlikely by chance (probability below 5% for 13 independent 95% tests).
- **The model fails** if the χ² of the new data against the fractal prediction has p < 0.01, or if several bins are inconsistent as defined above.
- **Fractal vs ΛCDM** is scored by Δχ² = χ²_fractal − χ²_ΛCDM, using the data covariance plus each prediction's own spread. Both are fixed predictions, so exp(−Δχ²/2) is their likelihood ratio. |Δχ²| < 2 counts as no preference, 2 to 6 as moderate and 6 or more as strong. Simulated releases under each model show how often a value as extreme as the observed one occurs.
- **The sign test** in the paper's sentence is reported, together with how often it fires in simulations. On new data at the forecast precision it fires about 9% of the time even when the fractal model is true, so it is not used on its own to reject the model.

## What one release can show

At the forecast precision (release covariance = 3/5 of DR2), about 40% of the information in each bin is new. On the new data alone, the median Δχ² is about −1.4 if the fractal model is true and +0.5 if ΛCDM is true, with a spread of several units. One release is therefore unlikely to be decisive on its own, and a clear verdict will probably need the final DESI data or other surveys.
