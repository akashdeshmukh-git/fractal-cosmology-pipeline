# Prediction for the next DESI BAO release

Recorded 9 November 2026, before any DESI release after DR2.

Source: A. Deshmukh, "Redshift Evolution of the Fractal Cosmic Web Correction to the Hubble Constant", *Annalen der Physik* **538**, e70321 (2026), [doi:10.1002/andp.70321](https://doi.org/10.1002/andp.70321). Received 13 August 2026, accepted 15 October 2026. Section 11 and Table 15.

The numbers below are copied from `prediction/table15.json`. That file and `prediction/posterior_samples.npz` were rebuilt with `make_prediction.py` and agree with the published Table 15 to every printed digit. Nothing in `prediction/` is changed after this date.

## What was predicted

The prediction is the posterior-predictive distribution of the fractal model (formation function exp(−z/z_char), z_char free) fitted to DESI DR2 BAO + Pantheon+ with a flat H0 prior, marginalised over H_E, Ωm and z_char. r_d is fixed at 147.09 Mpc. The 95% range is the union of the DR2-based and the DR1-based ranges. ΛCDM is fitted to the same data in the same way.

| z_eff | Observable | Fractal median | Fractal 68% | Fractal 95% (DR1 ∪ DR2) | ΛCDM median | σ_next | Δ/σ_next |
|---|---|---|---|---|---|---|---|
| 0.295 | D_V/r_d | 7.938 | [7.895, 7.981] | [7.807, 8.041] | 7.931 | 0.059 | +0.1 |
| 0.510 | D_M/r_d | 13.295 | [13.230, 13.360] | [13.093, 13.455] | 13.297 | 0.128 | 0.0 |
| 0.510 | D_H/r_d | 22.323 | [22.222, 22.427] | [22.027, 22.590] | 22.480 | 0.339 | −0.5 |
| 0.706 | D_M/r_d | 17.425 | [17.349, 17.499] | [17.180, 17.618] | 17.454 | 0.140 | −0.2 |
| 0.706 | D_H/r_d | 19.868 | [19.761, 19.971] | [19.569, 20.116] | 19.978 | 0.264 | −0.4 |
| 0.934 | D_M/r_d | 21.666 | [21.578, 21.752] | [21.381, 21.888] | 21.712 | 0.126 | −0.4 |
| 0.934 | D_H/r_d | 17.399 | [17.299, 17.494] | [17.109, 17.620] | 17.434 | 0.154 | −0.2 |
| 1.321 | D_M/r_d | 27.717 | [27.609, 27.825] | [27.379, 27.976] | 27.759 | 0.252 | −0.2 |
| 1.321 | D_H/r_d | 14.039 | [13.959, 14.120] | [13.797, 14.257] | 13.985 | 0.172 | +0.3 |
| 1.484 | D_M/r_d | 29.910 | [29.796, 30.026] | [29.552, 30.186] | 29.941 | 0.580 | −0.1 |
| 1.484 | D_H/r_d | 12.891 | [12.812, 12.969] | [12.653, 13.109] | 12.814 | 0.404 | +0.2 |
| 2.330 | D_M/r_d | 38.860 | [38.712, 39.016] | [38.396, 39.245] | 38.810 | 0.410 | +0.1 |
| 2.330 | D_H/r_d | 8.686 | [8.603, 8.778] | [8.444, 8.954] | 8.587 | 0.079 | +1.3 |

σ_next is the DR2 error scaled by √(3/5), a rough forecast of the next-release precision. Summed over the 13 observables, the fractal and ΛCDM medians differ by Δχ² ≈ 3.2 with the forecast covariance.

The paper's three statements (Section 11):

1. D_H/r_d at z ≃ 0.5 to 0.7 comes out 0.5 to 0.7% below the ΛCDM fit.
2. At z = 2.33, D_H/r_d is 1.2% above ΛCDM. This is the largest single-bin signal, 1.3σ of the forecast error.
3. D_M/r_d stays within 0.3% of ΛCDM over 0.5 < z < 2.4.

The paper's falsification sentence: "Measurements falling outside the 95% ranges in several bins, or D_H/r_d at z = 2.33 below the ΛCDM value while D_H/r_d at z ≃ 0.5 to 0.7 lies above it, would disfavour the model in its present form."

## How the test will be scored

These choices are fixed now, before the new data exist. `check_next_release.py` implements them.

- **No refitting.** The new release is compared with the frozen distributions above. If the release uses different effective redshifts, the same frozen posterior samples are evaluated at the new redshifts.
- **"Outside the 95% range"** means a measurement is inconsistent with the prediction once its own error is included: |pull| > 1.96, where the pull divides the difference from the predicted median by the measurement error and the prediction spread added in quadrature. The prediction must fail this against both the DR2-based and the DR1-based prediction. Comparing the measured value with the prediction band alone, ignoring the measurement error, is not used: under that reading 5 of the 13 DR2 points, the data the prediction was built from, would already fall outside.
- **"Several bins"** means 3 or more, and only counts if that many is unlikely by chance (probability below 5% for 13 independent 95% tests).
- **The model fails** if the χ² of the new data against the fractal prediction has p < 0.01, or if several bins are inconsistent as defined above.
- **Fractal vs ΛCDM** is scored by Δχ² = χ²_fractal − χ²_ΛCDM, using the new covariance plus each prediction's own spread. Both are fixed predictions, so exp(−Δχ²/2) is their likelihood ratio. |Δχ²| < 2 counts as no preference, 2 to 6 as moderate and 6 or more as strong.
- **The sign test** in the paper's sentence is reported, together with how often it fires by chance at the new precision. At DR2-like precision it fires in about 8% of simulated releases even when the fractal model is true, so it is not used on its own to reject the model.
