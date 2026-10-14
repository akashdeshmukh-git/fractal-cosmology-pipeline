# CMB secondary anisotropies (ISW and lensing)

A late-time change to H(z) leaves the sound horizon untouched but can still alter the late integrated Sachs–Wolfe (ISW) effect and CMB lensing. In principle that could bias the early-universe parameters, such as ω_m, that a ΛCDM analysis of the CMB returns, and that the compressed CMB likelihood in the paper relies on. This check measures the size of that effect.

## Method

- **Model in CAMB (v2.0.4).** The fractal correction is written as an effective dark-energy component, ρ_DE ∝ Ω_Λ + Ω_m(1+z)³ e^(−z/z_char)/2, with the equation of state w(z) computed from it. CAMB uses this w(z) table with its PPF dark-energy treatment, which allows w to cross −1. CAMB's H(z) agrees with the model formula to 5×10⁻⁵ at all redshifts tested.
- **Same early universe.** Both models use the Planck 2018 values ω_b = 0.02237, ω_c = 0.1200, ln(10¹⁰A_s) = 3.044, n_s = 0.9649, τ = 0.0544 and Σm_ν = 0.06 eV. The fractal model's H_E is chosen so that its θ* equals the ΛCDM value.
- **Spectra compared.** Lensed TT, TE and EE, the unlensed TT at ℓ ≤ 30 (where the late ISW appears), and the lensing potential C_ℓ^φφ.
- **Parameter bias.** A Fisher calculation for a Planck-like experiment (f_sky = 0.7, 5′ beam, 33 and 70 μK·arcmin noise in temperature and polarisation, ℓ = 2–2500, τ prior of 0.007) gives the shift in each ΛCDM parameter when the fractal-model spectra are fitted with ΛCDM. The forecast reproduces Planck's precision on ω_m (σ = 0.0011).

## Results

| z_char | change in C_ℓ^φφ (ℓ = 40–400) | largest low-ℓ TT change (ℓ ≤ 30) | largest TT change / cosmic variance | bias on ω_m | Δχ² after ΛCDM refit |
|---|---|---|---|---|---|
| 0.35 | +0.53% | 0.6% | 0.017 | +0.08σ | 0.03 |
| **0.542 (paper)** | **+0.43%** | **0.9%** | **0.015** | **+0.06σ** | **0.01** |
| 0.75 | −0.02% | 1.6% | 0.036 | −0.01σ | 0.02 |

At the paper's best fit, the shift in ω_m is 6.6×10⁻⁵, about 6% of its uncertainty. Every other ΛCDM parameter shifts by less than 0.1σ. The lensing power changes by under 0.6%, compared with the roughly 2.5% precision of the Planck lensing amplitude, and the ISW change stays below 4% of cosmic variance at every multipole.

## Conclusion

For a late-time correction that is consistent with the supernova data, the change in the CMB secondary anisotropies is far too small to bias the early-universe parameters. The compressed CMB likelihood used in the paper, D_M(z*) and ω_m from Planck ΛCDM, therefore remains valid for this model.

## Limitations

- The correction is treated as a smooth dark-energy component (CAMB's PPF prescription). A version in which the correction clusters would change the lensing and ISW results.
- This is a Fisher forecast for a Planck-like experiment, not a fit to the Planck likelihood. It includes TT, TE and EE but not the lensing reconstruction, which is compared through C_ℓ^φφ instead.
- The early-universe parameters are fixed at Planck values with θ* matched, not refitted jointly with BAO and supernovae.

## Running

```bash
pip install camb scipy
python cmbsec.py 0.542    # spectra for both models; argument is z_char
python fisher.py 0.542    # diagnostics and parameter bias
```

`cmbsec.py` takes about 10 seconds and `fisher.py` about a minute. Both write their output (`info_*.json`, `fisher_*.json`, `spec_*.npz`) to the current directory; the results for z_char = 0.35, 0.542 and 0.75 are in `results/`.
