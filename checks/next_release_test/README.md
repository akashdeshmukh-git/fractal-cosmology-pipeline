# Next DESI release test

Tests the prediction published in Section 11 and Table 15 of the paper ([doi:10.1002/andp.70321](https://doi.org/10.1002/andp.70321)) against a new DESI BAO release. The prediction is frozen in `prediction/` and written out in [`PREDICTION.md`](PREDICTION.md), together with the scoring rules, which were fixed before the new data.

## Files

| File | Description |
|---|---|
| `PREDICTION.md` | The published prediction, the paper's three statements and the scoring rules |
| `prediction/table15.json` | Table 15 in machine-readable form |
| `prediction/desi_dr2_*.txt` | The DESI DR2 BAO data the prediction was built from (from [CobayaSampler/bao_data](https://github.com/CobayaSampler/bao_data)), used to separate the new data from DR2 |
| `prediction/posterior_samples.npz` | 6000 posterior samples of (H_E, Ωm, z_char) for the fractal model and ΛCDM fitted to DESI DR2 + Pantheon+, the frozen basis of the prediction (DR1-based samples are also stored, only to reproduce the paper's table) |
| `check_next_release.py` | Runs the test on a new release |
| `make_prediction.py` | Rebuilds `prediction/` and checks it against the paper's Table 15 (writes to `results/remake/`, never overwrites the frozen files) |
| `model.py` | Background model and likelihoods, the same numerics as the paper |
| `get_data.py` | Downloads DESI DR1/DR2 and Pantheon+ into `data/` |
| `make_mock_release.py` | Builds a mock next release (DR2 plus new data drawn from one model) to see how the test behaves |
| `mock/` | An example mock release with ΛCDM as the truth |
| `results/` | Dry run on DR2 and the example mock release |

## Requirements

Python 3 with NumPy, SciPy, pandas and Matplotlib.

## Usage

When DESI publishes a new release, get its BAO mean vector and covariance in the Cobaya `bao_data` format (a mean file with lines `z value quantity`, where the quantity is `DM_over_rs`, `DH_over_rs` or `DV_over_rs`, plus a covariance matrix file). Then:

```bash
pip install numpy scipy pandas matplotlib
python check_next_release.py NEW_mean.txt NEW_cov.txt --label DR3
```

This takes under a minute and needs no other data. It prints two sections, the new data only (the test) and the whole release (for reference). Each has a table of every measurement against the prediction, the χ² of both models, the paper's three statements and a plain-language verdict. It writes `results/DR3_report.json` and `results/DR3_figure.png`.

To also refit both models to the new BAO + Pantheon+ (same grid and priors as the paper):

```bash
python get_data.py
python check_next_release.py NEW_mean.txt NEW_cov.txt --label DR3 --refit
```

The refit uses the whole release, so it includes DR2. It updates the parameter constraints but is not an independent test. On the example ΛCDM-truth mock release, the refit still prefers the fractal model (Δχ²_min = −3.3), because the DR2 part and Pantheon+ are unchanged.

To check that the frozen prediction reproduces the paper:

```bash
python get_data.py
python make_prediction.py     # about 3 minutes on two cores
```

## Method

**Separating the new data.** The next release will contain the DR2 galaxies, and the prediction was built from DR2. Comparing the whole release with the prediction would therefore partly re-use DR2, and it would favour whichever model fitted DR2 better, which is the fractal model. The script splits the release into DR2 plus independent new data. For Gaussian measurements of nested samples the information (inverse covariance) adds, so

x_new = F_new⁻¹ (F_release x_release − F_DR2 x_DR2),  with covariance F_new⁻¹ and F_new = F_release − F_DR2.

x_new is independent of DR2, so the frozen prediction genuinely predicts it. This is done for every bin that matches a DR2 bin (same observable, redshift within 0.05). If a bin's effective redshift has moved, the DR2 value is carried to the new redshift with the ΛCDM shape. If F_new is not positive definite, the release cannot be treated as DR2 plus new data (probably an analysis change), and the script says so.

**The test.** Nothing is refitted. Both models' predictions are the frozen posterior-predictive distributions from DESI DR2 + Pantheon+, evaluated at the release's effective redshifts. DR1 is not used: DR2 already contains the DR1 galaxies. For each measurement the script computes the pull against the prediction, with the measurement error and the prediction spread added in quadrature. Over all measurements it computes χ² against each prediction, using the data covariance plus the prediction's own covariance, and Δχ² = χ²_fractal − χ²_ΛCDM. Because both predictions are fixed, exp(−Δχ²/2) is their likelihood ratio.

**Calibration.** The script simulates 4000 releases under each model to show how often a Δχ² as large as the observed one would occur under each, and how often the paper's sign test fires by chance. For the new data the simulations are a posterior-predictive draw plus noise. For the whole release they are the observed DR2 plus simulated new data, so they include the overlap.

**Verdict.** The model fails if the new data are inconsistent with it: χ² goodness of fit p < 0.01, or 3 or more measurements off at 95% when that many is unlikely by chance. Otherwise Δχ² decides which prediction the data prefer. The full rules are in [`PREDICTION.md`](PREDICTION.md).

## Validation

- `make_prediction.py` reproduces the published Table 15 to every printed digit, and its posterior samples are identical to those of the paper pipeline.
- Separating the new data recovers the injected new-data vector and covariance of a mock release exactly (to 10⁻¹⁴), including when bin redshifts are shifted.
- Mock releases at the forecast precision (release covariance = 3/5 of DR2, so about 40% of each bin's information is new) show why the separation matters:

| Truth | Median Δχ², whole release | Median Δχ², whole release if it were independent of DR2 | Median Δχ², new data only |
|---|---|---|---|
| ΛCDM | −1.6 | +0.8 | +0.5 |
| Fractal model | −3.3 | −2.8 | −1.4 |

When ΛCDM is true, the whole release still leans toward the fractal model because it re-uses DR2. The new data alone do not. The cost is that the new data carry less information, so one release on its own is unlikely to be decisive. The paper's sign test fires in about 9% of fractal-true and 24% of ΛCDM-true simulations of new data.

## Dry run and example

`results/DR2_dryrun_*` runs the script on DR2 itself. There is no new data, so the new-data test is reported as not available, and the whole-release comparison only checks the code (Δχ² = −2.27, 0 of 13 points off at 95%).

`results/mock_lcdm_*` runs it on the example mock release in `mock/`, made with `python make_mock_release.py lcdm`, where ΛCDM is the truth. The new data give Δχ² = +2.89 (moderate support for ΛCDM), while the whole release gives +0.90 (no clear preference), showing the DR2 overlap pulling the whole-release result toward the fractal model.

## Assumptions and caveats

- r_d is fixed at 147.09 Mpc, as in the paper. The prediction is for D/r_d with that r_d, so it is consistent, but a release-wide change in the r_d convention would need checking first.
- The integration grid ends at z = 2.5. A new bin above that stops the script with an error.
- The prediction uses BAO + Pantheon+ only, with no CMB, as in Section 11 of the paper.
- If the new release reports quantities other than D_M/r_d, D_H/r_d or D_V/r_d (for example D_M/D_H), the script stops and lists them. They would need to be added to `model.py`.
