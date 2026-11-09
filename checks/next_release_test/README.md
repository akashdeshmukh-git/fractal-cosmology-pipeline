# Next DESI release test

The paper ([doi:10.1002/andp.70321](https://doi.org/10.1002/andp.70321), Section 11 and Table 15) predicted what DESI's next BAO release will measure. This folder lets anyone check that prediction against the new data with one command, and compares it with ΛCDM. The prediction is frozen in [`prediction/`](prediction/) and written out, with the scoring rules, in [`PREDICTION.md`](PREDICTION.md). Both were fixed before the new data existed.

## Quick start

```bash
git clone https://github.com/akashdeshmukh-git/fractal-cosmology-pipeline.git
cd fractal-cosmology-pipeline/checks/next_release_test
pip install -r requirements.txt

python check_next_release.py --demo
```

`--demo` runs on a bundled mock release, so you can see what a result looks like right away. Its output is in [`results/demo_summary.md`](results/demo_summary.md).

When DESI publishes the next release, run it on the real data in one of two ways:

```bash
# if the release is in CobayaSampler/bao_data under the same file names as DR2:
python check_next_release.py --cobaya desi_bao_dr3

# or with the mean and covariance files from anywhere:
python check_next_release.py NEW_mean.txt NEW_cov.txt --label DR3
```

The files are in the Cobaya `bao_data` format: a mean file with lines `z value quantity` (quantity `DM_over_rs`, `DH_over_rs` or `DV_over_rs`) and a covariance matrix file. A run takes under a minute and needs nothing else.

## What you get

- **A summary at the end of the output**, in plain words. For example: "The data are consistent with the published prediction. They moderately prefer LCDM over the fractal model."
- **`results/<label>_summary.md`**: the result, a figure and a table of every measurement against both predictions. It is readable directly on GitHub.
- **`results/<label>_figure.png`**: the measurements against the predicted band, relative to ΛCDM.
- **`results/<label>_report.json`**: every number, for further analysis.

## How to read the result

The new release contains the DR2 galaxies the prediction was built from, so the script first removes the DR2 part and tests the prediction on the **new data only**. That is the real test. The **whole release** is shown too, for reference, but it partly re-uses DR2 and so leans toward the fractal model.

For each, the script answers two questions:

1. **Is the prediction still consistent with the data?** It fails if the data disagree with it clearly: the χ² goodness of fit has p < 0.01, or 3 or more measurements are off at 95% when that is unlikely by chance.
2. **Which prediction fits better, the fractal model or ΛCDM?** This is Δχ² = χ²_fractal − χ²_ΛCDM. Negative favours the fractal model. |Δχ²| < 2 means no clear preference, 2 to 6 moderate, 6 or more strong.

About 40% of each bin's information in a release like the forecast is new, so one release on its own is unlikely to be decisive. Expect "leaning one way" rather than a final answer.

## Other options

```bash
python get_data.py                                           # DESI DR1/DR2 and Pantheon+, needed only for the two below
python check_next_release.py NEW_mean.txt NEW_cov.txt --label DR3 --refit
python make_prediction.py                                    # rebuild the frozen prediction from scratch (about 3 min)
python make_mock_release.py frac                             # make your own mock release (truth: frac or lcdm)
```

`--refit` also refits both models to the new BAO + Pantheon+ with the paper's grid and priors. The refit uses the whole release, so it includes DR2. It updates the parameter constraints but is not an independent test: on the example ΛCDM-truth mock release it still prefers the fractal model (Δχ²_min = −3.3).

## Files

| File | Description |
|---|---|
| `check_next_release.py` | Runs the test |
| `PREDICTION.md` | The published prediction, the paper's three statements and the scoring rules |
| `prediction/` | The frozen prediction: Table 15 (`table15.json`), 6000 posterior samples of (H_E, Ωm, z_char) for both models fitted to DESI DR2 + Pantheon+ (`posterior_samples.npz`), and the DR2 BAO data they were built from (`desi_dr2_*.txt`, from [CobayaSampler/bao_data](https://github.com/CobayaSampler/bao_data)) |
| `model.py` | Background model and likelihoods, the same numerics as the paper |
| `make_prediction.py` | Rebuilds `prediction/` and checks it against the paper's Table 15 (writes to `results/remake/`, never overwrites the frozen files) |
| `get_data.py` | Downloads DESI DR1/DR2 and Pantheon+ into `data/` |
| `make_mock_release.py`, `mock/` | Mock releases (real DR2 plus new data drawn from one model), and the one used by `--demo` |
| `results/` | The demo and a dry run on DR2 |
| `requirements.txt` | Python packages (NumPy, SciPy, pandas, Matplotlib) |

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

## Dry run and demo

- `results/DR2_dryrun_*`: the script run on DR2 itself. There is no new data, so the new-data test is reported as not available, and the whole-release comparison only checks the code (Δχ² = −2.27, 0 of 13 points off at 95%).
- `results/demo_*`: the bundled mock release, where ΛCDM is the truth. The new data give Δχ² = +2.89 (moderate support for ΛCDM), while the whole release gives +0.90 (no clear preference). This shows the DR2 overlap pulling the whole-release result toward the fractal model.

## Assumptions and caveats

- r_d is fixed at 147.09 Mpc, as in the paper. The prediction is for D/r_d with that r_d, so it is consistent, but a release-wide change in the r_d convention would need checking first.
- The integration grid ends at z = 2.5. A new bin above that stops the script with an error.
- The prediction uses BAO + Pantheon+ only, with no CMB, as in Section 11 of the paper.
- If the new release reports quantities other than D_M/r_d, D_H/r_d or D_V/r_d (for example D_M/D_H), the script stops and lists them. They would need to be added to `model.py`.
