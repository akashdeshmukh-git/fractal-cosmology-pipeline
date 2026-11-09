# Next DESI release test

Tests the prediction published in Section 11 and Table 15 of the paper ([doi:10.1002/andp.70321](https://doi.org/10.1002/andp.70321)) against a new DESI BAO release. The prediction is frozen in `prediction/` and written out in [`PREDICTION.md`](PREDICTION.md), together with the scoring rules, which were fixed before the new data.

## Files

| File | Description |
|---|---|
| `PREDICTION.md` | The published prediction, the paper's three statements and the scoring rules |
| `prediction/table15.json` | Table 15 in machine-readable form |
| `prediction/posterior_samples.npz` | 6000 posterior samples of (H_E, Ωm, z_char) for the fractal model and ΛCDM fitted to DESI DR2 + Pantheon+, the frozen basis of the prediction (DR1-based samples are also stored, only to reproduce the paper's table) |
| `check_next_release.py` | Runs the test on a new release |
| `make_prediction.py` | Rebuilds `prediction/` and checks it against the paper's Table 15 (writes to `results/remake/`, never overwrites the frozen files) |
| `model.py` | Background model and likelihoods, the same numerics as the paper |
| `get_data.py` | Downloads DESI DR1/DR2 and Pantheon+ into `data/` |
| `results/` | Dry run on DR2 |

## Requirements

Python 3 with NumPy, SciPy, pandas and Matplotlib.

## Usage

When DESI publishes a new release, get its BAO mean vector and covariance in the Cobaya `bao_data` format (a mean file with lines `z value quantity`, where the quantity is `DM_over_rs`, `DH_over_rs` or `DV_over_rs`, plus a covariance matrix file). Then:

```bash
pip install numpy scipy pandas matplotlib
python check_next_release.py NEW_mean.txt NEW_cov.txt --label DR3
```

This takes about a minute and needs no other data. It prints a table of every measurement against the prediction, the χ² of both models, the paper's three statements and a plain-language verdict, and writes `results/DR3_report.json` and `results/DR3_figure.png`.

To also refit both models to the new BAO + Pantheon+ (same grid and priors as the paper):

```bash
python get_data.py
python check_next_release.py NEW_mean.txt NEW_cov.txt --label DR3 --refit
```

To check that the frozen prediction reproduces the paper:

```bash
python get_data.py
python make_prediction.py     # about 3 minutes on two cores
```

## Method

**The test.** Nothing is refitted. Both models' predictions are the frozen posterior-predictive distributions from DESI DR2 + Pantheon+, evaluated at the new release's effective redshifts. DR1 is not used: DR2 already contains the DR1 galaxies. For each measurement the script computes the pull against the prediction, with the measurement error and the prediction spread added in quadrature. Over all measurements it computes χ² against each prediction, using the new covariance plus the prediction's own covariance, and Δχ² = χ²_fractal − χ²_ΛCDM. Because both predictions are fixed, exp(−Δχ²/2) is their likelihood ratio.

**Calibration.** The script simulates 4000 releases under each model (a posterior-predictive draw plus noise with the new covariance) to show how often a Δχ² as large as the observed one would occur under each, and how often the paper's sign test fires by chance.

**Verdict.** The model fails if the new data are inconsistent with it: χ² goodness of fit p < 0.01, or 3 or more measurements off at 95% when that many is unlikely by chance. Otherwise Δχ² decides which prediction the data prefer. The full rules are in [`PREDICTION.md`](PREDICTION.md).

## Validation

`make_prediction.py` reproduces the published Table 15 to every printed digit, and its posterior samples are identical to those of the paper pipeline.

Mock releases drawn from each model at the forecast precision (DR2 covariance × 3/5) behave as expected. The median Δχ² is −2.8 when the fractal model is true and +0.7 when ΛCDM is true. The two are not symmetric because the fractal prediction has the wider spread, and χ² includes each prediction's spread. So one release on its own gives at most moderate evidence, as the paper says (about 1.8σ), and a ΛCDM universe would usually look inconclusive rather than clearly against the model. The paper's sign test fires in 6% of fractal-true mocks and 24% of ΛCDM-true mocks.

## Dry run

This uses DR2, the data the prediction was built from, so it checks the code, not the model.

| Data | χ² fractal | χ² ΛCDM | Δχ² | Points off at 95% |
|---|---|---|---|---|
| DESI DR2 (13 points) | 8.43 | 10.70 | −2.27 | 0 |

The script evaluates the frozen samples at whatever redshifts a release uses, so new or shifted bins are handled. This was checked with a release whose bins differ from Table 15.

## Assumptions and caveats

- r_d is fixed at 147.09 Mpc, as in the paper. The prediction is for D/r_d with that r_d, so it is consistent, but a release-wide change in the r_d convention would need checking first.
- The integration grid ends at z = 2.5. A new bin above that stops the script with an error.
- The prediction uses BAO + Pantheon+ only, with no CMB, as in Section 11 of the paper.
- If the new release reports quantities other than D_M/r_d, D_H/r_d or D_V/r_d (for example D_M/D_H), the script stops and lists them. They would need to be added to `model.py`.
