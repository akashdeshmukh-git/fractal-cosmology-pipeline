# Sound-horizon treatment

Refits ΛCDM and the fractal model to DESI DR2 BAO, Pantheon+ and the compressed CMB likelihood with the BAO sound horizon r_d handled in four ways, to test how much the results depend on fixing r_d = 147.09 Mpc as the paper does. Runs from 28 September 2026.

## Files

| File | Description |
|---|---|
| `model.py` | Background, BAO, Pantheon+ (analytic M_B marginalisation) and compressed-CMB likelihoods |
| `fit.py` | Minimiser with the fixed, free and free-consistent r_d modes. Run directly, it reproduces the BAO and BAO+SN comparison |
| `fit_cmb.py` | BAO+SN+CMB fits in the fixed, free and free-consistent modes |
| `fit_computed.py` | BAO+SN+CMB fits with r_d computed from ω_m |
| `check_repro.py` | Comparison with the Cobaya reproduction |
| `make_table.py` | Builds `results/summary_table.md` |
| `results/` | JSON results and run logs |

## Requirements

Python 3 with NumPy, SciPy and pandas. The BAO and supernova data come from Cobaya's data repositories (see Usage).

## Usage

```bash
pip install numpy scipy pandas
git clone --depth 1 --branch v2.6 https://github.com/CobayaSampler/bao_data.git
git clone --depth 1 --branch v1.8 https://github.com/CobayaSampler/sn_data.git
python check_repro.py
python fit_cmb.py paper
python fit_cmb.py standard
python fit_computed.py paper
python fit_computed.py standard
python make_table.py
```

The argument selects the supernova sample: `paper` uses all 1624 non-calibrator SNe, `standard` uses the 1590 SNe at z > 0.01. Each `fit_cmb.py` run takes about 20 minutes on one core and each `fit_computed.py` run about 10 minutes.

## Method

The four treatments of r_d are:

1. Fixed (`rdmode='fixed'`). r_d = 147.09 Mpc, as in the paper.
2. Free, with r* fixed (`'free'`). r_d is a free parameter for BAO, and the compressed CMB keeps D_M(z*) = 13873 ± 25 Mpc. BAO then constrains only the shape of the expansion history.
3. Free, with r* scaling (`'free_consistent'`). r_d is free and the CMB target becomes D_M(z*) = 13873 × (r_d / 147.09) Mpc, which holds θ* fixed. This ignores the dependence of r_d on ω_m and is the loosest of the four.
4. Computed (`fit_computed.py`). r_d is computed from ω_m = Ω_m h² with the fitting formula of Aubourg et al. (2015), using ω_b = 0.02237 and ω_ν = 0.0006, normalised so that ω_m = 0.1430 gives 147.09 Mpc. The CMB target scales as in treatment 3. This is the consistent choice for a model that leaves pre-recombination physics unchanged.

For BAO alone and BAO+SN the treatment makes no difference. H_E is already free, BAO constrains only the product H_E·r_d, and supernovae with M_B marginalised carry no absolute scale, so fixed and free r_d give the same χ²_min (`results/log_bao_only_fixed_vs_free.txt`). Only BAO+SN+CMB is shown below.

All fits are best fits from Nelder–Mead with multiple starting points spread over the shape parameter's range.

## Validation

`check_repro.py` evaluates the likelihood at the best-fit points of the Cobaya reproduction:

| Case | This code | Cobaya | Difference |
|---|---|---|---|
| ΛCDM, BAO | 10.284 | 10.283 | +0.001 |
| ΛCDM, BAO+SN | 1470.098 | 1470.098 | 0.000 |
| exp free, BAO+SN | 1463.960 | 1463.958 | +0.002 |
| ΛCDM, BAO+SN+CMB (minimised) | 1471.788 | 1471.774 | +0.014 |
| exp free, BAO+SN+CMB (minimised) | 1464.895 | 1464.889 | +0.006 |

The CMB differences come from the integration grid and the choice z* = 1089.92 in the calibrated D_M(z*). They cancel in Δχ² to better than 0.01.

## Results

Δχ² relative to ΛCDM fitted with the same r_d treatment, BAO+SN+CMB. Fixed shapes use z_char = 1/2.262 and n = 2.262. z_char values are best fits.

### 1624 SNe (all non-calibrators)

| r_d treatment | exp fixed | power law fixed | exp free | power law free | z_char (exp free) | H_local (exp free) | r_d ΛCDM | r_d power law fixed |
|---|---|---|---|---|---|---|---|---|
| Fixed at 147.09 Mpc | −3.90 | −0.66 | −6.89 | −2.27 | 0.542 | 67.72 | 147.09 | 147.09 |
| Free, r* fixed | −4.32 | +0.61 | −5.91 | −1.27 | 0.523 | 67.52 | 147.97 | 147.01 |
| Free, r* scaling | −3.17 | −4.53 | −6.92 | −4.55 | 0.565 | 68.06 | 146.43 | 145.58 |
| Computed from ω_m | −4.63 | +3.39 | −7.04 | −1.10 | 0.526 | 67.62 | 147.30 | 147.60 |

### 1590 SNe (z > 0.01)

| r_d treatment | exp fixed | power law fixed | exp free | power law free | z_char (exp free) | H_local (exp free) | r_d ΛCDM | r_d power law fixed |
|---|---|---|---|---|---|---|---|---|
| Fixed at 147.09 Mpc | −3.47 | +0.01 | −6.02 | −1.89 | 0.535 | 67.79 | 147.09 | 147.09 |
| Free, r* fixed | −3.88 | +1.15 | −5.16 | −1.01 | 0.515 | 67.58 | 147.93 | 146.97 |
| Free, r* scaling | −2.72 | −3.92 | −5.97 | −3.99 | 0.557 | 68.11 | 146.42 | 145.57 |
| Computed from ω_m | −4.22 | +4.11 | −6.24 | −0.84 | 0.519 | 67.69 | 147.31 | 147.61 |

The free exponential is preferred over ΛCDM in every treatment, with Δχ² between −5.9 and −7.0 for 1624 SNe and between −5.2 and −6.2 for 1590 SNe. Its best-fit z_char stays between 0.52 and 0.57 and H_local between 67.5 and 68.1 km s⁻¹ Mpc⁻¹. Computing r_d from ω_m rather than fixing it shifts the free-exponential Δχ² by 0.15 (1624 SNe) and 0.22 (1590 SNe).

The ordering of the two fixed shapes depends on the treatment. With r_d free and r* scaling, the fixed power law (−4.53) does better than the fixed exponential (−3.17) at r_d = 145.6 Mpc. With the early universe unchanged, that r_d requires ω_m ≈ 0.149, about 5σ from the CMB value of 0.1430 ± 0.0011. With r_d computed from ω_m, the exponential is preferred (−4.63 against +3.39).

## Assumptions and caveats

- Results are best fits only. No MCMC was run for the free-r_d treatments, so they have no posterior intervals.
- The CMB enters through the paper's compressed two-number likelihood, not a full Boltzmann-code analysis.
- The r_d fitting formula is accurate to about 0.1%.
