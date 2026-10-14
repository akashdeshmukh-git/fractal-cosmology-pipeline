# CMB secondary anisotropies

Computes the effect of the fractal late-time correction on the late integrated Sachs–Wolfe (ISW) signal and on CMB lensing, and estimates the bias this would introduce in early-universe parameters inferred under ΛCDM. The compressed CMB likelihood in the paper takes ω_m from a Planck ΛCDM fit, so the quantity of interest is the bias on ω_m.

## Files

| File | Description |
|---|---|
| `cmbsec.py` | Lensed and unlensed CMB spectra and lensing potential for ΛCDM and the fractal model, plus a check of the fractal H(z) against the model formula |
| `fisher.py` | Spectrum comparison and Fisher bias on ΛCDM parameters for a Planck-like experiment |
| `results/` | Output for z_char = 0.35, 0.542 and 0.75 |

## Requirements

Python 3 with CAMB (tested with 2.0.4), NumPy and SciPy.

## Usage

```bash
pip install camb numpy scipy
python cmbsec.py 0.542
python fisher.py 0.542
```

The argument is z_char. `cmbsec.py` writes `spec_<z_char>.npz` and `info_<z_char>.json`. `fisher.py` reads the spectra and writes `fisher_<z_char>.json`. All output goes to the working directory. Run times are about 10 s and 1 min.

## Method

The fractal correction enters CAMB as an effective dark-energy density

ρ_DE(z) ∝ Ω_Λ + ½ Ω_m (1+z)³ exp(−z/z_char),

with w(z) derived from it and passed as a table to the PPF dark-energy module, which allows w to cross −1. The resulting H(z) matches the model formula to within 5×10⁻⁵ for 0 ≤ z ≤ 2.

Both models use the Planck 2018 early-universe parameters: ω_b = 0.02237, ω_c = 0.1200, ln(10¹⁰A_s) = 3.044, n_s = 0.9649, τ = 0.0544 and Σm_ν = 0.06 eV. In the fractal model, H_E is set so that θ* equals its ΛCDM value.

For the parameter bias, the fractal-model TT, TE and EE spectra are treated as data and fitted with ΛCDM in (ω_b, ω_c, θ_MC, ln A_s, n_s, τ) using a Fisher approximation. The noise model is Planck-like: f_sky = 0.7, a 5′ beam, 33 and 70 μK·arcmin in temperature and polarization, ℓ = 2–2500, and a Gaussian prior σ(τ) = 0.007. With this set-up σ(ω_m) = 0.0011, the Planck 2018 value.

## Results

| z_char | ΔC_ℓ^φφ / C_ℓ^φφ, ℓ = 40–400 | max ΔC_ℓ^TT / C_ℓ^TT, unlensed, ℓ ≤ 30 | max ΔC_ℓ^TT / cosmic variance, ℓ ≤ 30 | bias on ω_m | Δχ² after ΛCDM refit |
|---|---|---|---|---|---|
| 0.35 | +0.53% | 0.6% | 0.017 | +0.08σ | 0.03 |
| 0.542 | +0.43% | 0.9% | 0.015 | +0.06σ | 0.01 |
| 0.75 | −0.02% | 1.6% | 0.036 | −0.01σ | 0.02 |

At z_char = 0.542, the paper's best fit, the bias on ω_m is 6.6×10⁻⁵ and no ΛCDM parameter shifts by more than 0.1σ. The change in lensing power is small next to the roughly 2.5% uncertainty on the Planck lensing amplitude. Within the assumptions below, the compressed CMB likelihood used in the paper is unaffected.

## Assumptions and caveats

- The correction is treated as smooth, with no clustering. A clustering correction would change the ISW and lensing results.
- The bias is a Fisher forecast with a simplified Planck-like noise model. It does not use the Planck likelihood or the lensing reconstruction.
- The early-universe parameters are fixed at the Planck 2018 values with θ* matched. They are not refitted jointly with BAO and supernovae.
