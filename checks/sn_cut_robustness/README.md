# Sensitivity to the low-redshift supernova cut

Peculiar velocities make the redshifts of nearby supernovae unreliable distance indicators, and z < 0.02 is often treated as problematic. This check refits the model with the Pantheon+ sample restricted to progressively higher minimum redshift.

## Setup

- Data: DESI DR2 BAO (13 points, full covariance), Pantheon+ (STAT+SYS covariance, analytic M_B marginalisation), and the paper's compressed CMB likelihood.
- r_d fixed at 147.09 Mpc, as in the paper.
- Models: ΛCDM; exponential formation function with z_char free; exponential with z_char = 1/2.262.
- The z > 0.01 row reproduces the standard-cut values of the Cobaya reproduction (`../cobaya_reproduction/`), which validates the set-up.

## Results (Δχ² relative to ΛCDM)

| SN selection | N_SN | BAO+SN, z_char free | BAO+SN+CMB, z_char free | BAO+SN+CMB, z_char fixed | best-fit z_char (with CMB) | H_local (with CMB) |
|---|---|---|---|---|---|---|
| z > 0.01 | 1590 | −5.29 | −6.02 | −3.47 | 0.535 | 67.79 |
| z > 0.02 | 1436 | −4.44 | −5.15 | −2.95 | 0.530 | 67.83 |
| z > 0.03 | 1233 | −5.97 | −6.67 | −3.37 | 0.553 | 67.63 |

## Conclusion

Removing the low-redshift supernovae changes the preference for the fractal correction by less than one unit of χ², with no monotonic trend with the cut. The best-fit z_char and H_local are essentially unchanged. The result is therefore not driven by peculiar-velocity-dominated supernovae.

## Running

```bash
git clone --depth 1 --branch v2.6 https://github.com/CobayaSampler/bao_data.git
git clone --depth 1 --branch v1.8 https://github.com/CobayaSampler/sn_data.git
python zcut.py 0.02    # minimum supernova redshift
```

`model.py` and `fit.py` are the same likelihood and minimiser as in `../rd_robustness/`. Each run takes about 10 minutes on one core. Results are best fits only; no MCMC was run.
