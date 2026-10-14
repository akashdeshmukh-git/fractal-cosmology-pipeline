# Low-redshift supernova cut

Refits ΛCDM and the fractal model with the Pantheon+ sample restricted to z > 0.01, 0.02 and 0.03, to test whether nearby supernovae, whose redshifts are most affected by peculiar velocities, drive the preference for the fractal correction.

## Files

| File | Description |
|---|---|
| `zcut.py` | Fits ΛCDM and the exponential model for a given minimum supernova redshift |
| `model.py`, `fit.py` | Likelihoods and minimiser, identical to those in `../rd_robustness/` |
| `results/` | Run logs for z > 0.01, 0.02 and 0.03 |

## Requirements

Python 3 with NumPy, SciPy and pandas. The BAO and supernova data come from Cobaya's data repositories (see Usage).

## Usage

```bash
pip install numpy scipy pandas
git clone --depth 1 --branch v2.6 https://github.com/CobayaSampler/bao_data.git
git clone --depth 1 --branch v1.8 https://github.com/CobayaSampler/sn_data.git
python zcut.py 0.02
```

The argument is the minimum supernova redshift. Each run takes about 10 minutes on one core.

## Method

- Data: DESI DR2 BAO (13 points, full covariance), Pantheon+ (STAT+SYS covariance, M_B marginalised analytically) and the paper's compressed CMB likelihood.
- r_d is fixed at 147.09 Mpc, as in the paper.
- Models: ΛCDM, the exponential formation function with z_char free, and the exponential with z_char = 1/2.262.
- Best fits from Nelder–Mead with multiple starting points.

The z > 0.01 case reproduces the standard-cut values of `../cobaya_reproduction/`.

## Results

Δχ² relative to ΛCDM.

| SN selection | N_SN | BAO+SN, z_char free | BAO+SN+CMB, z_char free | BAO+SN+CMB, z_char fixed | best-fit z_char (with CMB) | H_local (with CMB) |
|---|---|---|---|---|---|---|
| z > 0.01 | 1590 | −5.29 | −6.02 | −3.47 | 0.535 | 67.79 |
| z > 0.02 | 1436 | −4.44 | −5.15 | −2.95 | 0.530 | 67.83 |
| z > 0.03 | 1233 | −5.97 | −6.67 | −3.37 | 0.553 | 67.63 |

With the CMB included, the free-exponential Δχ² ranges from −5.15 to −6.67 across the three cuts and does not change monotonically with the cut. The best-fit z_char stays between 0.53 and 0.55 and H_local between 67.6 and 67.8 km s⁻¹ Mpc⁻¹.

## Assumptions and caveats

- Results are best fits only. No MCMC was run.
- r_d is fixed. Its treatment is tested separately in `../rd_robustness/`.

## Acknowledgement

This test was suggested by George Efstathiou (private communication).
