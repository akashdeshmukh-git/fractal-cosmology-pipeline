"""
Test the published prediction against a new DESI BAO release.

    python check_next_release.py NEW_mean.txt NEW_cov.txt --label DR3
    python check_next_release.py NEW_mean.txt NEW_cov.txt --label DR3 --refit   # also refit (needs get_data.py)

Inputs are in the Cobaya bao_data format: a mean file with lines "z value quantity"
(quantity = DM_over_rs, DH_over_rs or DV_over_rs) and a covariance matrix file.

Part 1 (the test). Nothing is refitted. The fractal model and LCDM predictions are the frozen
posterior-predictive distributions published in Deshmukh, Ann. Phys. 538, e70321 (2026),
Section 11 and Table 15, built from DESI DR2 + Pantheon+ (DR1 + Pantheon+ for the 95% union).
They are evaluated at whatever effective redshifts the new release uses. The script reports
  - each measurement against the fractal 68%/95% ranges and the LCDM prediction,
  - chi^2 of the new data for both frozen predictions (data covariance + prediction spread),
  - the three statements made in the paper and the paper's own falsification criterion.

Part 2 (--refit, optional). Refits both models to the new BAO + Pantheon+ with the same grid and
priors as the paper, and reports z_char, H_local and the chi^2 difference.

Outputs: results/<label>_report.json and results/<label>_figure.png.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '2')
import sys
import json
import argparse
import numpy as np
import model as m

HERE = os.path.dirname(os.path.abspath(__file__))
FROZEN = os.path.join(HERE, 'prediction')
RESULTS = os.path.join(HERE, 'results')
DOI = '10.1002/andp.70321'

LOWZ_DH = (0.45, 0.75)      # "DH/r_d at z ~ 0.5-0.7 should come out 0.5-0.7% below the LCDM fit"
HIGHZ_DH = (2.2, 2.45)      # "At z = 2.33, DH/r_d is predicted 1.2% above LCDM"
DM_RANGE = (0.5, 2.4)       # "DM/r_d stays within 0.3% of LCDM over 0.5 < z < 2.4"
SEVERAL = 3                 # "outside the 95% ranges in several bins": taken as 3 or more
Z95 = 1.96                  # a measurement is inconsistent if |pull| > 1.96, pull including its own error


def load_samples():
    d = np.load(os.path.join(FROZEN, 'posterior_samples.npz'))
    return {k: (d[f'{k}_H'], d[f'{k}_Om'], d[f'{k}_zc']) for k in ('frac_DR2', 'frac_DR1', 'lcdm_DR2')}


def predictive(S, key, z, qty):
    H, Om, zc = S[key]
    mdl = 'frac' if key.startswith('frac') else 'lcdm'
    return m.bao_observables(z, qty, mdl, H, Om, zc if mdl == 'frac' else None)


def published_bin(table, z, q):
    for r in table:
        if r['quantity'] == q and abs(r['z'] - z) < 0.02:
            return r
    return None


def chi2(r, C):
    return float(r @ np.linalg.solve(C, r))


def weighted_offset(vals, errs):
    w = 1 / np.asarray(errs) ** 2
    return float(np.sum(w * vals) / np.sum(w)), float(1 / np.sqrt(np.sum(w)))


def sign_pattern(y, med_l, lowz_idx, highz_idx, err_pct):
    """paper's sign test on data vectors y (n_mock, n): True if low-z DH above LCDM and z=2.33 DH below"""
    off = 100 * (y / med_l - 1)
    w = 1 / err_pct[lowz_idx] ** 2
    low = (off[:, lowz_idx] * w).sum(axis=1) / w.sum()
    return (low > 0) & (off[:, highz_idx] < 0)


def calibrate(bao, Pf2, Pl2, med_f, med_l, Cf, Cl, lowz_idx, highz_idx, n=4000, seed=11):
    """Simulate the new release under each hypothesis (a posterior-predictive draw plus noise with the
    new covariance) and record Delta chi2 and the sign test, to show how decisive the observed values are."""
    rng = np.random.default_rng(seed)
    noise = rng.multivariate_normal(np.zeros(bao.n), bao.cov, size=(2, n))
    Wf, Wl = np.linalg.inv(Cf), np.linalg.inv(Cl)
    err_pct = 100 * bao.err / med_l
    out = {}
    for k, (name, P) in enumerate((('frac', Pf2), ('lcdm', Pl2))):
        y = P[rng.integers(0, len(P), n)] + noise[k]
        rf, rl = y - med_f, y - med_l
        d = np.einsum('ij,jk,ik->i', rf, Wf, rf) - np.einsum('ij,jk,ik->i', rl, Wl, rl)
        out[name] = dict(dchi2=d)
        if lowz_idx and highz_idx is not None:
            out[name]['sign'] = sign_pattern(y, med_l, lowz_idx, highz_idx, err_pct)
    return out


def test_prediction(bao, S, table):
    Pf2 = predictive(S, 'frac_DR2', bao.z, bao.qty)
    Pf1 = predictive(S, 'frac_DR1', bao.z, bao.qty)
    Pl2 = predictive(S, 'lcdm_DR2', bao.z, bao.qty)
    med_f, med_l = np.median(Pf2, axis=0), np.median(Pl2, axis=0)

    rows = []
    for i, (z, q, o, e) in enumerate(zip(bao.z, bao.qty, bao.obs, bao.err)):
        f2 = np.percentile(Pf2[:, i], [2.5, 16, 50, 84, 97.5])
        f1 = np.percentile(Pf1[:, i], [2.5, 97.5])
        lo, hi = min(f2[0], f1[0]), max(f2[4], f1[1])
        sf, sl, s1 = Pf2[:, i].std(), Pl2[:, i].std(), Pf1[:, i].std()
        pull_f2 = (o - f2[2]) / np.hypot(e, sf)
        pull_f1 = (o - np.median(Pf1[:, i])) / np.hypot(e, s1)
        pull_l = (o - med_l[i]) / np.hypot(e, sl)
        pub = published_bin(table, z, q)
        rows.append(dict(
            z=float(z), quantity=str(q), observed=float(o), error=float(e),
            frac_median=f2[2], frac_68=[f2[1], f2[3]], frac_95=[lo, hi],
            lcdm_median=float(med_l[i]),
            inside_model_band=bool(lo <= o <= hi),
            pull_frac=float(pull_f2), pull_frac_DR1=float(pull_f1), pull_lcdm=float(pull_l),
            consistent_frac=bool(min(abs(pull_f2), abs(pull_f1)) <= Z95),
            consistent_lcdm=bool(abs(pull_l) <= Z95),
            obs_vs_lcdm_pct=float(100 * (o / med_l[i] - 1)),
            obs_vs_lcdm_err_pct=float(100 * e / med_l[i]),
            frac_vs_lcdm_pct=float(100 * (f2[2] / med_l[i] - 1)),
            published_bin=None if pub is None else pub['z']))

    Cf = bao.cov + np.cov(Pf2, rowvar=False)
    Cl = bao.cov + np.cov(Pl2, rowvar=False)
    res = dict(
        n_points=bao.n,
        chi2_frac=chi2(bao.obs - med_f, Cf), chi2_lcdm=chi2(bao.obs - med_l, Cl),
        chi2_frac_datacov_only=chi2(bao.obs - med_f, bao.cov),
        chi2_lcdm_datacov_only=chi2(bao.obs - med_l, bao.cov))
    res['delta_chi2'] = res['chi2_frac'] - res['chi2_lcdm']
    res['likelihood_ratio_frac_over_lcdm'] = float(np.exp(-0.5 * res['delta_chi2']))

    # the three statements of Section 11
    def pick(q, zr):
        return [r for r in rows if r['quantity'] == q and zr[0] <= r['z'] <= zr[1]]

    st = {}
    low = pick('DH_over_rs', LOWZ_DH)
    if low:
        off, err = weighted_offset([r['obs_vs_lcdm_pct'] for r in low], [r['obs_vs_lcdm_err_pct'] for r in low])
        st['1_DH_lowz_below_LCDM'] = dict(
            claim='DH/r_d at z~0.5-0.7 is 0.5-0.7% below LCDM',
            predicted_pct=[round(r['frac_vs_lcdm_pct'], 2) for r in low],
            measured_pct=round(off, 2), measured_err_pct=round(err, 2),
            same_sign=bool(off < 0))
    high = pick('DH_over_rs', HIGHZ_DH)
    if high:
        r = high[0]
        st['2_DH_z233_above_LCDM'] = dict(
            claim='DH/r_d at z=2.33 is about 1.2% above LCDM',
            predicted_pct=round(r['frac_vs_lcdm_pct'], 2),
            measured_pct=round(r['obs_vs_lcdm_pct'], 2), measured_err_pct=round(r['obs_vs_lcdm_err_pct'], 2),
            same_sign=bool(r['obs_vs_lcdm_pct'] > 0))
    dm = pick('DM_over_rs', DM_RANGE)
    if dm:
        st['3_DM_close_to_LCDM'] = dict(
            claim='DM/r_d within 0.3% of LCDM over 0.5<z<2.4 (DM does not discriminate)',
            predicted_pct=[round(r['frac_vs_lcdm_pct'], 2) for r in dm],
            measured_pct=[round(r['obs_vs_lcdm_pct'], 2) for r in dm],
            pulls_frac=[round(r['pull_frac'], 2) for r in dm])
    res['statements'] = st

    from scipy.stats import chi2 as chi2dist
    res['goodness_of_fit_p_frac'] = float(chi2dist.sf(res['chi2_frac'], bao.n))
    res['goodness_of_fit_p_lcdm'] = float(chi2dist.sf(res['chi2_lcdm'], bao.n))
    lowz_idx = [i for i, r in enumerate(rows) if r['quantity'] == 'DH_over_rs' and LOWZ_DH[0] <= r['z'] <= LOWZ_DH[1]]
    hz = [i for i, r in enumerate(rows) if r['quantity'] == 'DH_over_rs' and HIGHZ_DH[0] <= r['z'] <= HIGHZ_DH[1]]
    cal = calibrate(bao, Pf2, Pl2, med_f, med_l, Cf, Cl, lowz_idx, hz[0] if hz else None)
    d = res['delta_chi2']
    res['calibration'] = dict(
        n_mocks=len(cal['frac']['dchi2']),
        median_dchi2_if_frac_true=float(np.median(cal['frac']['dchi2'])),
        median_dchi2_if_lcdm_true=float(np.median(cal['lcdm']['dchi2'])),
        prob_dchi2_this_high_if_frac_true=float(np.mean(cal['frac']['dchi2'] >= d)),
        prob_dchi2_this_low_if_lcdm_true=float(np.mean(cal['lcdm']['dchi2'] <= d)))
    if 'sign' in cal['frac']:
        res['calibration']['sign_test_fires_if_frac_true'] = float(np.mean(cal['frac']['sign']))
        res['calibration']['sign_test_fires_if_lcdm_true'] = float(np.mean(cal['lcdm']['sign']))

    # the paper's falsification criterion. A measurement counts as outside the predicted 95% range when
    # it is inconsistent with the prediction once its own error is included (|pull| > 1.96 against the
    # DR2- and the DR1-based prediction). Comparing the measured value with the prediction band alone,
    # ignoring the measurement error, would fail even the DR2 data the prediction was built from.
    from scipy.stats import binom
    n_out = sum(not r['consistent_frac'] for r in rows)
    n_out_lcdm = sum(not r['consistent_lcdm'] for r in rows)
    n_out_band = sum(not r['inside_model_band'] for r in rows)
    reversed_pattern = ('1_DH_lowz_below_LCDM' in st and '2_DH_z233_above_LCDM' in st
                        and not st['1_DH_lowz_below_LCDM']['same_sign']
                        and not st['2_DH_z233_above_LCDM']['same_sign'])
    res['falsification'] = dict(
        n_outside_95=n_out, n_outside_95_lcdm=n_out_lcdm,
        n_outside_model_band_ignoring_errors=n_out_band,
        expected_by_chance=round(0.05 * bao.n, 2),
        prob_this_many_by_chance=float(binom.sf(n_out - 1, bao.n, 0.05)) if n_out else 1.0,
        several_outside=bool(n_out >= SEVERAL),
        reversed_pattern=bool(reversed_pattern),
        disfavoured=bool(n_out >= SEVERAL or reversed_pattern))
    res['fails'] = bool(res['goodness_of_fit_p_frac'] < 0.01 or
                        (n_out >= SEVERAL and res['falsification']['prob_this_many_by_chance'] < 0.05))
    return rows, res


def verdict(res):
    """Plain-language verdict. The model fails outright if the new data are inconsistent with the
    published prediction (chi^2 goodness of fit p < 0.01, or 3+ points off at 95% when that is
    unlikely by chance). Otherwise the two frozen predictions are compared by Delta chi2, which is
    -2 ln of their likelihood ratio (both are fixed predictions, nothing is fitted)."""
    f, c, d = res['falsification'], res['calibration'], res['delta_chi2']
    lines = []
    if res['fails']:
        lines.append(f"FAILS: the new data are inconsistent with the published prediction "
                     f"(goodness of fit p = {res['goodness_of_fit_p_frac']:.3f}; "
                     f"{f['n_outside_95']} of {res['n_points']} points off at 95%).")
    else:
        lines.append(f"SURVIVES: the new data are consistent with the published prediction "
                     f"(goodness of fit p = {res['goodness_of_fit_p_frac']:.2f}; "
                     f"{f['n_outside_95']} of {res['n_points']} points off at 95%).")
    if d <= -6:
        tail = 'strong support for the fractal prediction over LCDM'
    elif d <= -2:
        tail = 'moderate support for the fractal prediction over LCDM'
    elif d < 2:
        tail = 'no clear preference; the release cannot yet separate the two'
    elif d < 6:
        tail = 'moderate support for LCDM over the fractal prediction'
    else:
        tail = 'strong support for LCDM over the fractal prediction'
    lines.append(f"Fractal vs LCDM: Delta chi2 = {d:+.2f} (likelihood ratio {res['likelihood_ratio_frac_over_lcdm']:.2g}), "
                 f"{tail}. If the fractal model were true a value this high or higher would occur "
                 f"{100 * c['prob_dchi2_this_high_if_frac_true']:.0f}% of the time; if LCDM were true, this low or lower "
                 f"{100 * c['prob_dchi2_this_low_if_lcdm_true']:.0f}% of the time.")
    if 'sign_test_fires_if_frac_true' in c:
        lines.append(f"Paper's sign test {'FIRED' if f['reversed_pattern'] else 'did not fire'} "
                     f"(at this precision it fires {100 * c['sign_test_fires_if_frac_true']:.0f}% of the time even if the "
                     f"fractal model is right, and {100 * c['sign_test_fires_if_lcdm_true']:.0f}% if LCDM is right, "
                     f"so on its own it is weak evidence).")
    return '\n'.join(lines)


def refit(bao, data_dir):
    sn = m.SN(os.path.join(data_dir, 'Pantheon+SH0ES.dat'), os.path.join(data_dir, 'Pantheon+SH0ES_STAT+SYS.cov'))
    rng = np.random.default_rng(7)
    out = {}
    for mdl in ('frac', 'lcdm'):
        print(f'refit {mdl} on new BAO + Pantheon+ ...', flush=True)
        H, Om, zc, cmin = m.sample_grid(*m.grid_chi2(bao, sn, mdl), rng)
        out[f'{mdl}_chi2_min'] = cmin
        if mdl == 'frac':
            out['zchar_16_50_84'] = np.percentile(zc, [16, 50, 84]).tolist()
            out['zchar_fraction_in_0.3_0.7'] = float(np.mean((zc > 0.3) & (zc < 0.7)))
            out['Hlocal_16_50_84'] = np.percentile(m.H_local(H, Om), [16, 50, 84]).tolist()
            out['HE_16_50_84'] = np.percentile(H, [16, 50, 84]).tolist()
    out['delta_chi2_min'] = out['frac_chi2_min'] - out['lcdm_chi2_min']
    return out


def figure(bao, rows, S, path, label):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    zg = np.linspace(0.25, 2.45, 120)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, q in zip(axes, ('DH_over_rs', 'DM_over_rs')):
        qq = np.full(zg.size, q)
        l2 = np.median(predictive(S, 'lcdm_DR2', zg, qq), axis=0)
        f2 = 100 * (predictive(S, 'frac_DR2', zg, qq) / l2 - 1)
        f1 = 100 * (predictive(S, 'frac_DR1', zg, qq) / l2 - 1)
        p2 = np.percentile(f2, [2.5, 16, 50, 84, 97.5], axis=0)
        p1 = np.percentile(f1, [2.5, 97.5], axis=0)
        ax.fill_between(zg, p2[0], p2[4], color='C0', alpha=0.18, lw=0, label='Fractal, published 95%')
        ax.fill_between(zg, p2[1], p2[3], color='C0', alpha=0.35, lw=0, label='Fractal, published 68%')
        ax.plot(zg, p1[0], 'C1--', lw=1, label='Fractal, DR1-based 95%')
        ax.plot(zg, p1[1], 'C1--', lw=1)
        ax.axhline(0, color='k', lw=0.8, label='LCDM (DR2 + Pantheon+)')
        sel = [r for r in rows if r['quantity'] == q]
        ax.errorbar([r['z'] for r in sel], [r['obs_vs_lcdm_pct'] for r in sel],
                    yerr=[r['obs_vs_lcdm_err_pct'] for r in sel], fmt='o', color='C3', ms=5,
                    capsize=2, label=f'DESI {label}')
        ax.set_xlabel('z')
        ax.set_ylabel(f"{q.split('_')[0]}/r_d relative to LCDM [%]")
        ax.set_ylim(-6, 6)
    axes[0].legend(fontsize=8, loc='upper left')
    fig.suptitle(f'Published prediction (doi:{DOI}) vs DESI {label}', fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mean')
    ap.add_argument('cov')
    ap.add_argument('--label', default='new')
    ap.add_argument('--refit', action='store_true')
    ap.add_argument('--data-dir', default=os.path.join(HERE, 'data'))
    a = ap.parse_args()

    bao = m.BAO(a.mean, a.cov)
    S = load_samples()
    table = json.load(open(os.path.join(FROZEN, 'table15.json')))['rows']
    rows, res = test_prediction(bao, S, table)

    print(f'\nPublished prediction (doi:{DOI}) vs DESI {a.label}: {bao.n} measurements\n')
    print(f"{'z':>6} {'obs':>4} {'measured':>17} {'fractal [95% band]':>27} {'LCDM':>8} {'pull':>6} {'pullL':>6}  ok   bin")
    for r in rows:
        print(f"{r['z']:6.3f} {r['quantity'][:2]:>4} {r['observed']:9.3f}±{r['error']:<7.3f}"
              f"{r['frac_median']:9.3f} [{r['frac_95'][0]:.3f},{r['frac_95'][1]:.3f}]"
              f"{r['lcdm_median']:9.3f} {r['pull_frac']:+6.2f} {r['pull_lcdm']:+6.2f}  {'yes' if r['consistent_frac'] else 'NO ':4} "
              f"{'published' if r['published_bin'] is not None else 'new z'}")
    print(f"\nchi2 fractal = {res['chi2_frac']:.2f}   chi2 LCDM = {res['chi2_lcdm']:.2f}   "
          f"Delta chi2 = {res['delta_chi2']:+.2f}   (data-cov only: "
          f"{res['chi2_frac_datacov_only'] - res['chi2_lcdm_datacov_only']:+.2f})")
    for k, s in res['statements'].items():
        if 'measured_err_pct' in s:
            print(f"statement {k[0]}: {s['claim']}: predicted {s['predicted_pct']}%, measured "
                  f"{s['measured_pct']:+.2f} ± {s['measured_err_pct']:.2f}% -> {'same sign' if s['same_sign'] else 'OPPOSITE sign'}")
        else:
            print(f"statement {k[0]}: {s['claim']}: pulls {s['pulls_frac']}")
    f = res['falsification']
    print(f"inconsistent at 95% (|pull| > 1.96): fractal {f['n_outside_95']}, LCDM {f['n_outside_95_lcdm']} "
          f"(expected by chance {f['expected_by_chance']}; probability of {f['n_outside_95']} or more "
          f"{f['prob_this_many_by_chance']:.2f})")
    print(f"paper's sign test (low-z DH above LCDM and z=2.33 DH below): "
          f"{'TRIGGERED' if f['reversed_pattern'] else 'not triggered'}")
    res['verdict'] = verdict(res)
    print('\nVERDICT: ' + res['verdict'])

    if a.refit:
        res['refit'] = refit(bao, a.data_dir)
        rf = res['refit']
        print(f"\nRefit to {a.label} + Pantheon+: z_char = {rf['zchar_16_50_84'][1]:.2f} "
              f"[{rf['zchar_16_50_84'][0]:.2f}, {rf['zchar_16_50_84'][2]:.2f}], "
              f"{100 * rf['zchar_fraction_in_0.3_0.7']:.0f}% in [0.3, 0.7]; H_local = {rf['Hlocal_16_50_84'][1]:.1f}; "
              f"Delta chi2_min (fractal - LCDM) = {rf['delta_chi2_min']:+.2f}")

    os.makedirs(RESULTS, exist_ok=True)
    out = dict(label=a.label, mean_file=os.path.basename(a.mean), cov_file=os.path.basename(a.cov),
               prediction_doi=DOI, rows=rows, **res)
    json.dump(out, open(os.path.join(RESULTS, f'{a.label}_report.json'), 'w'), indent=1, default=float)
    figure(bao, rows, S, os.path.join(RESULTS, f'{a.label}_figure.png'), a.label)
    print(f"\nwrote results/{a.label}_report.json and results/{a.label}_figure.png")


if __name__ == '__main__':
    main()
