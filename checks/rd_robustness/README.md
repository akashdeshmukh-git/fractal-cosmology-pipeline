# Sensitivity to the treatment of the sound horizon r_d

A robustness test for *Redshift Evolution of the Fractal Cosmic Web Correction to the Hubble Constant* (Annalen der Physik ms 9212113). The paper fixes the BAO sound horizon at r_d = 147.09 Mpc, the Planck ΛCDM value. This folder refits every model with r_d treated four different ways and checks whether the preference for the fractal model survives.

Date of runs: 28 September 2026.

## Summary

- **BAO alone and BAO+SN: r_d cannot matter.** H_E is already a free parameter and BAO only measures distances in units of r_d, so the data constrain only the product H_E·r_d. Supernovae, with the absolute magnitude marginalised, do not depend on the absolute scale at all. Fixing or freeing r_d therefore gives the same χ²_min (checked numerically in `results/log_bao_only_fixed_vs_free.txt`). The BAO+SN preference (Δχ² = −6.1) is already independent of r_d.
- **BAO+SN+CMB: the main result is robust.** The free exponential beats ΛCDM in every treatment:
  - Δχ² from −5.9 to −7.0 with the 1624-SN sample
  - Δχ² from −5.2 to −6.2 with the standard SN cut
  - z_char stays at 0.52 to 0.57
  - H_local stays at 67.5 to 68.1 km/s/Mpc, so the Hubble-tension conclusion is unchanged
- **One claim depends on how r_d is treated: that the CMB favours the exponential over the power law.** If r_d is left completely free *and* the CMB's r* is allowed to move with it, the fixed power law (−4.5) beats the fixed exponential (−3.2). It does this by shrinking r_d to 145.6 Mpc. With pre-recombination physics unchanged (as in this model), that would need ω_m ≈ 0.149, about 5σ from the CMB value 0.1430 ± 0.0011. When r_d is computed from ω_m instead, the exponential–power-law gap gets wider (−4.6 against +3.4).

## Statement in the revised paper (R2)

Since the model leaves pre-recombination physics unchanged, r_d is set by ω_m. Computing r_d from ω_m instead of fixing it changes the Δχ² of the free exponential by about 0.2 (0.15 with the 1624-SN sample, 0.22 with the standard cut). Treating r_d as fully free leaves the preference for the fractal correction intact (Δχ² between −5.9 and −6.9), but removes the distinction between the exponential and power-law forms, and only for r_d values that the CMB excludes.

## Results (Δχ² relative to ΛCDM fitted the same way)

### BAO+SN+CMB, paper SN sample (1624 SNe)

| r_d treatment | exp fixed | power law fixed | exp free | power law free | z_char (exp free) | H_local (exp free) | r_d ΛCDM | r_d power law fixed |
|---|---|---|---|---|---|---|---|---|
| Fixed at 147.09 Mpc (paper) | −3.90 | −0.66 | −6.89 | −2.27 | 0.542 | 67.72 | 147.09 | 147.09 |
| Free, CMB r* fixed | −4.32 | +0.61 | −5.91 | −1.27 | 0.523 | 67.52 | 147.97 | 147.01 |
| Free, CMB r* scales with r_d | −3.17 | −4.53 | −6.92 | −4.55 | 0.565 | 68.06 | 146.43 | 145.58 |
| Computed from ω_m (physical) | −4.63 | +3.39 | −7.04 | −1.10 | 0.526 | 67.62 | 147.30 | 147.60 |

### BAO+SN+CMB, standard SN cut z > 0.01 (1590 SNe)

| r_d treatment | exp fixed | power law fixed | exp free | power law free | z_char (exp free) | H_local (exp free) | r_d ΛCDM | r_d power law fixed |
|---|---|---|---|---|---|---|---|---|
| Fixed at 147.09 Mpc (paper) | −3.47 | +0.01 | −6.02 | −1.89 | 0.535 | 67.79 | 147.09 | 147.09 |
| Free, CMB r* fixed | −3.88 | +1.15 | −5.16 | −1.01 | 0.515 | 67.58 | 147.93 | 146.97 |
| Free, CMB r* scales with r_d | −2.72 | −3.92 | −5.97 | −3.99 | 0.557 | 68.11 | 146.42 | 145.57 |
| Computed from ω_m (physical) | −4.22 | +4.11 | −6.24 | −0.84 | 0.519 | 67.69 | 147.31 | 147.61 |

Fixed shapes use z_char = 1/2.262 and n = 2.262. z_char values are best fits, not posterior medians.

## The four treatments

1. **Fixed** (`rdmode='fixed'`): r_d = 147.09 Mpc, as in the paper.
2. **Free, CMB r* fixed** (`'free'`): r_d is a free parameter for BAO, while the compressed CMB keeps D_M(z*) = 13873 ± 25 Mpc. This removes the absolute calibration from BAO, so BAO constrains only the shape of the expansion.
3. **Free, r* scales with r_d** (`'free_consistent'`): r_d is free and the CMB target becomes D_M(z*) = 13873 × (r_d/147.09), so θ* = r*/D_M(z*) is what is held fixed. This is the loosest test: it ignores that ω_m also sets r_d.
4. **Computed** (`fit_computed.py`): r_d is computed from ω_m = Ω_m h² with the Aubourg et al. (2015) fitting formula (ω_b = 0.02237, ω_ν = 0.0006), normalised so that ω_m = 0.1430 gives 147.09 Mpc. The CMB target scales as in (3). This is the physically consistent treatment for a model that leaves the early universe unchanged.

## Validation

`check_repro.py` evaluates this likelihood at the best-fit points of the Cobaya reproduction:

| Case | This code | Cobaya | Difference |
|---|---|---|---|
| ΛCDM, BAO | 10.284 | 10.283 | +0.001 |
| ΛCDM, BAO+SN | 1470.098 | 1470.098 | 0.000 |
| exp free, BAO+SN | 1463.960 | 1463.958 | +0.002 |
| ΛCDM, BAO+SN+CMB (minimised) | 1471.788 | 1471.774 | +0.014 |
| exp free, BAO+SN+CMB (minimised) | 1464.895 | 1464.889 | +0.006 |

The small CMB differences come from the integration grid and the choice z* = 1089.92 in the calibrated D_M(z*). They cancel in Δχ² to better than 0.01.

## Limitations

- Best fits only (Nelder–Mead, many starts over the shape range). No MCMC was run for the free-r_d cases, so there are no posterior intervals for them.
- The CMB is the paper's compressed two-number likelihood, not a full Boltzmann-code analysis.
- The r_d formula is a fitting function accurate to about 0.1%, which is enough for this test.

## Files

| File | Contents |
|---|---|
| `model.py` | Background, BAO, Pantheon+ (analytic M_B marginalisation) and compressed-CMB likelihoods |
| `fit.py` | Minimiser with the fixed, free and free-consistent r_d modes (running it directly reproduces the BAO and BAO+SN checks) |
| `fit_cmb.py` | BAO+SN+CMB fits in those three modes |
| `fit_computed.py` | BAO+SN+CMB fits with r_d computed from ω_m |
| `check_repro.py` | Comparison with the Cobaya numbers |
| `make_table.py` | Builds `results/summary_table.md` |
| `results/` | JSON results and run logs |

## Running

```bash
pip install numpy scipy pandas
git clone --depth 1 --branch v2.6 https://github.com/CobayaSampler/bao_data.git
git clone --depth 1 --branch v1.8 https://github.com/CobayaSampler/sn_data.git
python check_repro.py
python fit_cmb.py paper      # and: python fit_cmb.py standard
python fit_computed.py paper # and: python fit_computed.py standard
python make_table.py
```

Each `fit_cmb.py` run takes about 20 minutes on one core. `fit_computed.py` takes about 10 minutes.
