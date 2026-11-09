"""
Test the published prediction against a new DESI BAO release.

    python check_next_release.py NEW_mean.txt NEW_cov.txt --label DR3
    python check_next_release.py NEW_mean.txt NEW_cov.txt --label DR3 --refit   # also refit (needs get_data.py)

Inputs are in the Cobaya bao_data format: a mean file with lines "z value quantity"
(quantity = DM_over_rs, DH_over_rs or DV_over_rs) and a covariance matrix file.

The prediction. The fractal model and LCDM predictions are the frozen posterior-predictive
distributions published in Deshmukh, Ann. Phys. 538, e70321 (2026), Section 11 and Table 15,
built from DESI DR2 + Pantheon+. Nothing is refitted.

The overlap problem. A new DESI release contains the DR2 galaxies. Part of it is therefore data the
prediction was built from, and comparing the whole release with the prediction would partly re-use
DR2. The primary test removes this overlap. Treating the new release as DR2 plus independent new
data, the information (inverse covariance) adds: F_new = F_release - F_DR2, and the new-data-only
measurement is x_new = F_new^-1 (F_release x_release - F_DR2 x_DR2), with covariance F_new^-1. This
x_new is independent of DR2, so the frozen prediction genuinely predicts it. It is computed for the
bins that match a DR2 bin (same observable, |dz| <= 0.05).

The script reports, for the new-data-only measurement (primary) and for the whole release:
  - each measurement against the fractal prediction and LCDM,
  - chi^2 for both frozen predictions (data covariance + prediction spread) and Delta chi2,
  - simulated releases under each model, including the DR2 overlap, to calibrate Delta chi2,
  - the paper's three statements and its falsification criterion.

--refit additionally refits both models to the new BAO + Pantheon+ with the paper's grid and priors.

Outputs: results/<label>_report.json and results/<label>_figure.png.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '2')
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
MATCH_DZ = 0.05             # a new bin matches a DR2 bin if same observable and |dz| <= 0.05


class Data:
    def __init__(self, z, qty, obs, cov):
        self.z = np.asarray(z, float)
        self.qty = np.asarray(qty)
        self.obs = np.asarray(obs, float)
        self.cov = np.atleast_2d(np.asarray(cov, float))
        self.err = np.sqrt(np.diag(self.cov))
        self.n = len(self.z)

    @classmethod
    def from_files(cls, mean_file, cov_file):
        b = m.BAO(mean_file, cov_file)
        return cls(b.z, b.qty, b.obs, b.cov)


def load_samples():
    d = np.load(os.path.join(FROZEN, 'posterior_samples.npz'))
    return {k: (d[f'{k}_H'], d[f'{k}_Om'], d[f'{k}_zc']) for k in ('frac_DR2', 'lcdm_DR2')}


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


# ---------------------------------------------------------------- overlap with DR2
def match_bins(new, dr2):
    """one-to-one matches (i_new, j_dr2): same observable, nearest redshift within MATCH_DZ"""
    pairs, used = [], set()
    for i, (z, q) in enumerate(zip(new.z, new.qty)):
        cand = [(abs(z - z2), j) for j, (z2, q2) in enumerate(zip(dr2.z, dr2.qty))
                if q2 == q and abs(z - z2) <= MATCH_DZ and j not in used]
        if cand:
            j = min(cand)[1]
            pairs.append((i, j)); used.add(j)
    return pairs


def separate_new_data(new, dr2, S):
    """Split the new release into DR2 + independent new data (Gaussian, nested samples).
    Returns a dict with the new-data-only Data object and the pieces needed to simulate the release,
    or a dict with 'reason' if the split is not possible."""
    pairs = match_bins(new, dr2)
    if not pairs:
        return dict(reason='no bin of the new release matches a DR2 bin')
    i_n = np.array([p[0] for p in pairs]); j_2 = np.array([p[1] for p in pairs])
    x3 = new.obs[i_n]; C3 = new.cov[np.ix_(i_n, i_n)]
    x2 = dr2.obs[j_2].copy(); C2 = dr2.cov[np.ix_(j_2, j_2)].copy()
    # if a bin's effective redshift moved slightly, carry the DR2 value to the new redshift with the
    # LCDM shape (a small correction; the same factor scales its errors)
    l_new = np.median(predictive(S, 'lcdm_DR2', new.z[i_n], new.qty[i_n]), axis=0)
    l_old = np.median(predictive(S, 'lcdm_DR2', dr2.z[j_2], dr2.qty[j_2]), axis=0)
    r = l_new / l_old
    x2 = x2 * r; C2 = C2 * np.outer(r, r)
    F3, F2 = np.linalg.inv(C3), np.linalg.inv(C2)
    Fn = F3 - F2
    Fn = 0.5 * (Fn + Fn.T)
    ev = np.linalg.eigvalsh(Fn)
    frac_new = 1 - np.diag(C3) / np.diag(C2)
    if np.abs(ev).max() <= 1e-6 * np.abs(F3).max():
        return dict(reason='the release carries no information beyond DR2 (is it DR2 itself?)',
                    frac_new_per_bin=frac_new.tolist())
    if ev.min() <= 1e-8 * ev.max():
        return dict(reason='the new release is not more precise than DR2 in every combination of bins, '
                           'so it cannot be treated as DR2 plus independent new data '
                           '(probably an analysis change)', frac_new_per_bin=frac_new.tolist())
    Cn = np.linalg.inv(Fn); Cn = 0.5 * (Cn + Cn.T)
    xn = Cn @ (F3 @ x3 - F2 @ x2)
    return dict(data=Data(new.z[i_n], new.qty[i_n], xn, Cn), i_new=i_n, j_dr2=j_2,
                A=C3 @ F2, B=C3 @ Fn, x2=x2, Cn=Cn, frac_new_per_bin=frac_new.tolist(),
                dz_max=float(np.max(np.abs(new.z[i_n] - dr2.z[j_2]))),
                all_matched=bool(len(pairs) == new.n))


# ---------------------------------------------------------------- the test
def sign_pattern(y, med_l, lowz_idx, highz_idx, err_pct):
    """paper's sign test on data vectors y (n_mock, n): True if low-z DH above LCDM and z=2.33 DH below"""
    off = 100 * (y / med_l - 1)
    w = 1 / err_pct[lowz_idx] ** 2
    low = (off[:, lowz_idx] * w).sum(axis=1) / w.sum()
    return (low > 0) & (off[:, highz_idx] < 0)


def independent_sim(data):
    def sim(P, rng, n):
        return P[rng.integers(0, len(P), n)] + rng.multivariate_normal(np.zeros(data.n), data.cov, size=n)
    return sim


def overlap_sim(sep):
    """simulate the whole release as DR2 (as observed) plus independent new data drawn from the model"""
    i_n, Cn, A, B, x2 = sep['i_new'], sep['Cn'], sep['A'], sep['B'], sep['x2']

    def sim(P, rng, n):
        xn = P[rng.integers(0, len(P), n)][:, i_n] + rng.multivariate_normal(np.zeros(len(i_n)), Cn, size=n)
        y = np.empty((n, P.shape[1]))
        y[:, i_n] = (A @ x2)[None, :] + xn @ B.T
        return y
    return sim


def test_prediction(data, S, table, simulate, n_mocks=4000, seed=11):
    Pf = predictive(S, 'frac_DR2', data.z, data.qty)
    Pl = predictive(S, 'lcdm_DR2', data.z, data.qty)
    med_f, med_l = np.median(Pf, axis=0), np.median(Pl, axis=0)

    rows = []
    for i, (z, q, o, e) in enumerate(zip(data.z, data.qty, data.obs, data.err)):
        f = np.percentile(Pf[:, i], [2.5, 16, 50, 84, 97.5])
        sf, sl = Pf[:, i].std(), Pl[:, i].std()
        pf = (o - f[2]) / np.hypot(e, sf)
        pl = (o - med_l[i]) / np.hypot(e, sl)
        pub = published_bin(table, z, q)
        rows.append(dict(
            z=float(z), quantity=str(q), observed=float(o), error=float(e),
            frac_median=f[2], frac_68=[f[1], f[3]], frac_95=[f[0], f[4]],
            lcdm_median=float(med_l[i]),
            inside_model_band=bool(f[0] <= o <= f[4]),
            pull_frac=float(pf), pull_lcdm=float(pl),
            consistent_frac=bool(abs(pf) <= Z95), consistent_lcdm=bool(abs(pl) <= Z95),
            obs_vs_lcdm_pct=float(100 * (o / med_l[i] - 1)),
            obs_vs_lcdm_err_pct=float(100 * e / med_l[i]),
            frac_vs_lcdm_pct=float(100 * (f[2] / med_l[i] - 1)),
            published_bin=None if pub is None else pub['z']))

    Cf = data.cov + np.cov(Pf, rowvar=False)
    Cl = data.cov + np.cov(Pl, rowvar=False)
    res = dict(n_points=data.n,
               chi2_frac=chi2(data.obs - med_f, Cf), chi2_lcdm=chi2(data.obs - med_l, Cl),
               chi2_frac_datacov_only=chi2(data.obs - med_f, data.cov),
               chi2_lcdm_datacov_only=chi2(data.obs - med_l, data.cov))
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
            measured_pct=round(off, 2), measured_err_pct=round(err, 2), same_sign=bool(off < 0))
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

    from scipy.stats import chi2 as chi2dist, binom
    res['goodness_of_fit_p_frac'] = float(chi2dist.sf(res['chi2_frac'], data.n))
    res['goodness_of_fit_p_lcdm'] = float(chi2dist.sf(res['chi2_lcdm'], data.n))

    # calibration: simulated releases under each model
    rng = np.random.default_rng(seed)
    Wf, Wl = np.linalg.inv(Cf), np.linalg.inv(Cl)
    lowz_idx = [i for i, r in enumerate(rows) if r['quantity'] == 'DH_over_rs' and LOWZ_DH[0] <= r['z'] <= LOWZ_DH[1]]
    hz = [i for i, r in enumerate(rows) if r['quantity'] == 'DH_over_rs' and HIGHZ_DH[0] <= r['z'] <= HIGHZ_DH[1]]
    err_pct = 100 * data.err / med_l
    cal = {}
    for name, P in (('frac', Pf), ('lcdm', Pl)):
        y = simulate(P, rng, n_mocks)
        rf, rl = y - med_f, y - med_l
        cal[name] = dict(dchi2=np.einsum('ij,jk,ik->i', rf, Wf, rf) - np.einsum('ij,jk,ik->i', rl, Wl, rl))
        if lowz_idx and hz:
            cal[name]['sign'] = sign_pattern(y, med_l, lowz_idx, hz[0], err_pct)
    d = res['delta_chi2']
    res['calibration'] = dict(
        n_mocks=n_mocks,
        median_dchi2_if_frac_true=float(np.median(cal['frac']['dchi2'])),
        median_dchi2_if_lcdm_true=float(np.median(cal['lcdm']['dchi2'])),
        prob_dchi2_this_high_if_frac_true=float(np.mean(cal['frac']['dchi2'] >= d)),
        prob_dchi2_this_low_if_lcdm_true=float(np.mean(cal['lcdm']['dchi2'] <= d)))
    if 'sign' in cal['frac']:
        res['calibration']['sign_test_fires_if_frac_true'] = float(np.mean(cal['frac']['sign']))
        res['calibration']['sign_test_fires_if_lcdm_true'] = float(np.mean(cal['lcdm']['sign']))

    # the paper's falsification criterion. "Outside the 95% range" counts a measurement that is
    # inconsistent with the prediction once its own error is included (|pull| > 1.96). Comparing the
    # value with the prediction band alone, ignoring the error, would fail even the DR2 data.
    n_out = sum(not r['consistent_frac'] for r in rows)
    reversed_pattern = ('1_DH_lowz_below_LCDM' in st and '2_DH_z233_above_LCDM' in st
                        and not st['1_DH_lowz_below_LCDM']['same_sign']
                        and not st['2_DH_z233_above_LCDM']['same_sign'])
    res['falsification'] = dict(
        n_outside_95=n_out, n_outside_95_lcdm=sum(not r['consistent_lcdm'] for r in rows),
        n_outside_model_band_ignoring_errors=sum(not r['inside_model_band'] for r in rows),
        expected_by_chance=round(0.05 * data.n, 2),
        prob_this_many_by_chance=float(binom.sf(n_out - 1, data.n, 0.05)) if n_out else 1.0,
        several_outside=bool(n_out >= SEVERAL), reversed_pattern=bool(reversed_pattern))
    res['fails'] = bool(res['goodness_of_fit_p_frac'] < 0.01 or
                        (n_out >= SEVERAL and res['falsification']['prob_this_many_by_chance'] < 0.05))
    return rows, res


def pct(p):
    return '<1%' if p < 0.01 else ('>99%' if p > 0.99 else f'{100 * p:.0f}%')


def verdict(res):
    """Plain-language verdict. The model fails if the data are inconsistent with the published
    prediction (chi^2 goodness of fit p < 0.01, or 3+ points off at 95% when that is unlikely by
    chance). Otherwise the two frozen predictions are compared by Delta chi2 = -2 ln(likelihood ratio)."""
    f, c, d = res['falsification'], res['calibration'], res['delta_chi2']
    word = 'FAILS: inconsistent with' if res['fails'] else 'SURVIVES: consistent with'
    lines = [f"{word} the published prediction (goodness of fit p = {res['goodness_of_fit_p_frac']:.3f}; "
             f"{f['n_outside_95']} of {res['n_points']} points off at 95%)."]
    if d <= -6:
        tail = 'strong support for the fractal prediction over LCDM'
    elif d <= -2:
        tail = 'moderate support for the fractal prediction over LCDM'
    elif d < 2:
        tail = 'no clear preference; these data cannot yet separate the two'
    elif d < 6:
        tail = 'moderate support for LCDM over the fractal prediction'
    else:
        tail = 'strong support for LCDM over the fractal prediction'
    lines.append(f"Fractal vs LCDM: Delta chi2 = {d:+.2f} (likelihood ratio {res['likelihood_ratio_frac_over_lcdm']:.2g}), "
                 f"{tail}. In simulations, a value this high or higher occurs "
                 f"{pct(c['prob_dchi2_this_high_if_frac_true'])} of the time if the fractal model is true; "
                 f"this low or lower {pct(c['prob_dchi2_this_low_if_lcdm_true'])} of the time if LCDM is true.")
    if 'sign_test_fires_if_frac_true' in c:
        weak = ', so on its own it is weak evidence' if c['sign_test_fires_if_frac_true'] >= 0.05 else ''
        lines.append(f"Paper's sign test {'FIRED' if f['reversed_pattern'] else 'did not fire'} "
                     f"(in simulations it fires {pct(c['sign_test_fires_if_frac_true'])} of the time if the fractal "
                     f"model is true and {pct(c['sign_test_fires_if_lcdm_true'])} if LCDM is true{weak}).")
    return '\n'.join(lines)


def refit(bao_data, data_dir):
    sn = m.SN(os.path.join(data_dir, 'Pantheon+SH0ES.dat'), os.path.join(data_dir, 'Pantheon+SH0ES_STAT+SYS.cov'))
    bao = m.BAO.__new__(m.BAO)
    bao.z, bao.qty, bao.obs, bao.cov = bao_data.z, bao_data.qty, bao_data.obs, bao_data.cov
    from scipy.linalg import cholesky, solve_triangular
    bao.W = solve_triangular(cholesky(bao.cov, lower=True), np.eye(bao_data.n), lower=True)
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


# ---------------------------------------------------------------- output
def print_section(title, rows, res):
    print(f'\n=== {title}: {res["n_points"]} measurements ===\n')
    print(f"{'z':>6} {'obs':>4} {'measured':>17} {'fractal [95%]':>27} {'LCDM':>8} {'pull':>6} {'pullL':>6}  ok")
    for r in rows:
        print(f"{r['z']:6.3f} {r['quantity'][:2]:>4} {r['observed']:9.3f}±{r['error']:<7.3f}"
              f"{r['frac_median']:9.3f} [{r['frac_95'][0]:.3f},{r['frac_95'][1]:.3f}]"
              f"{r['lcdm_median']:9.3f} {r['pull_frac']:+6.2f} {r['pull_lcdm']:+6.2f}  "
              f"{'yes' if r['consistent_frac'] else 'NO'}")
    print(f"\nchi2 fractal = {res['chi2_frac']:.2f}   chi2 LCDM = {res['chi2_lcdm']:.2f}   "
          f"Delta chi2 = {res['delta_chi2']:+.2f}")
    for k, s in res['statements'].items():
        if 'measured_err_pct' in s:
            print(f"statement {k[0]}: {s['claim']}: predicted {s['predicted_pct']}%, measured "
                  f"{s['measured_pct']:+.2f} ± {s['measured_err_pct']:.2f}% -> "
                  f"{'same sign' if s['same_sign'] else 'OPPOSITE sign'}")
        else:
            print(f"statement {k[0]}: {s['claim']}: pulls {s['pulls_frac']}")
    print('\n' + verdict(res))


def figure(full_rows, new_rows, S, path, label):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    zg = np.linspace(0.25, 2.45, 120)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    ext = [abs(r['obs_vs_lcdm_pct']) + r['obs_vs_lcdm_err_pct'] for r in (full_rows or []) + (new_rows or [])]
    lim = min(15.0, max(6.0, max(ext, default=0) + 0.5))
    for ax, q in zip(axes, ('DH_over_rs', 'DM_over_rs')):
        qq = np.full(zg.size, q)
        l2 = np.median(predictive(S, 'lcdm_DR2', zg, qq), axis=0)
        p2 = np.percentile(100 * (predictive(S, 'frac_DR2', zg, qq) / l2 - 1), [2.5, 16, 84, 97.5], axis=0)
        ax.fill_between(zg, p2[0], p2[3], color='C0', alpha=0.18, lw=0, label='Fractal prediction (DR2), 95%')
        ax.fill_between(zg, p2[1], p2[2], color='C0', alpha=0.35, lw=0, label='Fractal prediction (DR2), 68%')
        ax.axhline(0, color='k', lw=0.8, label='LCDM (DR2 + Pantheon+)')
        for rows, kw in ((full_rows, dict(fmt='o', color='C3', label=f'DESI {label}, whole release')),
                         (new_rows, dict(fmt='s', color='C2', mfc='none', label=f'DESI {label}, new data only'))):
            sel = [r for r in (rows or []) if r['quantity'] == q]
            if sel:
                dz = 0.012 if kw['fmt'] == 's' else 0.0
                ax.errorbar([r['z'] + dz for r in sel], [r['obs_vs_lcdm_pct'] for r in sel],
                            yerr=[r['obs_vs_lcdm_err_pct'] for r in sel], ms=5, capsize=2, **kw)
        ax.set_xlabel('z')
        ax.set_ylabel(f"{q.split('_')[0]}/r_d relative to LCDM [%]")
        ax.set_ylim(-lim, lim)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=3, fontsize=8, frameon=False)
    fig.suptitle(f'Published prediction (doi:{DOI}) vs DESI {label}', fontsize=10)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
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

    new = Data.from_files(a.mean, a.cov)
    dr2 = Data.from_files(os.path.join(FROZEN, 'desi_dr2_mean.txt'), os.path.join(FROZEN, 'desi_dr2_cov.txt'))
    S = load_samples()
    table = json.load(open(os.path.join(FROZEN, 'table15.json')))['rows']
    print(f'Published prediction (doi:{DOI}) vs DESI {a.label}')

    sep = separate_new_data(new, dr2, S)
    out = dict(label=a.label, mean_file=os.path.basename(a.mean), cov_file=os.path.basename(a.cov),
               prediction_doi=DOI)
    new_rows = None
    if 'data' in sep:
        new_rows, new_res = test_prediction(sep['data'], S, table, independent_sim(sep['data']))
        new_res['verdict'] = verdict(new_res)
        fn = np.array(sep['frac_new_per_bin'])
        rng_txt = (f"{100 * fn.min():.0f}%" if round(100 * fn.min()) == round(100 * fn.max())
                   else f"{100 * fn.min():.0f}-{100 * fn.max():.0f}%")
        print(f"\nNew data only: DR2 removed from {sep['data'].n} matched bins "
              f"(new information per bin {rng_txt} of the release; "
              f"largest redshift shift {sep['dz_max']:.3f}).")
        print_section('NEW DATA ONLY (primary test)', new_rows, new_res)
        out['new_data_only'] = dict(rows=new_rows, frac_new_per_bin=sep['frac_new_per_bin'], **new_res)
        full_sim = overlap_sim(sep) if sep['all_matched'] else independent_sim(new)
        overlap_note = ('simulations include the DR2 overlap' if sep['all_matched'] else
                        'some bins have no DR2 counterpart, so the simulations ignore the overlap')
    else:
        print(f"\nNew data only: NOT AVAILABLE: {sep['reason']}.")
        out['new_data_only'] = dict(available=False, reason=sep['reason'],
                                    frac_new_per_bin=sep.get('frac_new_per_bin'))
        full_sim = independent_sim(new)
        overlap_note = ('the DR2 overlap could not be modelled, so this comparison partly re-uses DR2 '
                        'and favours whichever model fitted DR2 better')

    full_rows, full_res = test_prediction(new, S, table, full_sim)
    full_res['verdict'] = verdict(full_res)
    full_res['overlap_note'] = overlap_note
    print_section(f'WHOLE RELEASE (contains DR2; {overlap_note})', full_rows, full_res)
    out['whole_release'] = dict(rows=full_rows, **full_res)

    if a.refit:
        rf = refit(new, a.data_dir)
        out['refit'] = rf
        print(f"\nRefit to {a.label} + Pantheon+: z_char = {rf['zchar_16_50_84'][1]:.2f} "
              f"[{rf['zchar_16_50_84'][0]:.2f}, {rf['zchar_16_50_84'][2]:.2f}], "
              f"{100 * rf['zchar_fraction_in_0.3_0.7']:.0f}% in [0.3, 0.7]; H_local = {rf['Hlocal_16_50_84'][1]:.1f}; "
              f"Delta chi2_min (fractal - LCDM) = {rf['delta_chi2_min']:+.2f}")

    os.makedirs(RESULTS, exist_ok=True)
    json.dump(out, open(os.path.join(RESULTS, f'{a.label}_report.json'), 'w'), indent=1, default=float)
    figure(full_rows, new_rows, S, os.path.join(RESULTS, f'{a.label}_figure.png'), a.label)
    print(f"\nwrote results/{a.label}_report.json and results/{a.label}_figure.png")


if __name__ == '__main__':
    main()
