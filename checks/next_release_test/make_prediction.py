"""
Rebuild the published next-release prediction (paper Table 15) from scratch.

Fits the fractal model (z_char free) and LCDM to DESI DR2 + Pantheon+ and DESI DR1 + Pantheon+
on a 3-D grid with flat priors, draws 6000 posterior samples from each, and computes the
posterior-predictive BAO observables at the DR2 effective redshifts.

The frozen files in prediction/ were made with this script. Re-running it writes to results/remake/
and compares the result with the frozen files; it never overwrites them.

    python get_data.py
    python make_prediction.py          # about 3 minutes on two cores
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '2')
import json
import numpy as np
import model as m

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
FROZEN = os.path.join(HERE, 'prediction')
OUT = os.path.join(HERE, 'results', 'remake')


def predictive(samples, z, qty, mdl):
    H, Om, zc = samples
    return m.bao_observables(z, qty, mdl, H, Om, zc if mdl == 'frac' else None)


def build_table(S, bao2):
    """Table 15: fractal DR2 median/68%, DR1-or-DR2 95% union, LCDM DR2 median/68%, sigma_next."""
    P = {k: predictive(v, bao2.z, bao2.qty, 'frac' if k.startswith('frac') else 'lcdm') for k, v in S.items()}
    rows = []
    for i, (z, q) in enumerate(zip(bao2.z, bao2.qty)):
        f2 = np.percentile(P['frac_DR2'][:, i], [2.5, 16, 50, 84, 97.5])
        f1 = np.percentile(P['frac_DR1'][:, i], [2.5, 97.5])
        l2 = np.percentile(P['lcdm_DR2'][:, i], [16, 50, 84])
        sig = bao2.err[i] / bao2.obs[i] * np.sqrt(3 / 5) * f2[2]
        rows.append(dict(z=round(float(z), 3), quantity=str(q),
                         frac_median=f2[2], frac_68=[f2[1], f2[3]],
                         frac_95_union=[min(f2[0], f1[0]), max(f2[4], f1[1])],
                         frac_95_DR2=[f2[0], f2[4]], frac_95_DR1=list(f1),
                         lcdm_median=l2[1], lcdm_68=[l2[0], l2[2]],
                         sigma_next=sig, delta_over_sigma_next=(f2[2] - l2[1]) / sig))
    med_f = np.median(P['frac_DR2'], axis=0)
    med_l = np.median(P['lcdm_DR2'], axis=0)
    d = med_f - med_l
    dchi2 = float(d @ np.linalg.inv(bao2.cov * 0.6) @ d)
    return rows, dchi2


def main():
    os.makedirs(OUT, exist_ok=True)
    bao2 = m.BAO(os.path.join(DATA, 'desi_dr2_mean.txt'), os.path.join(DATA, 'desi_dr2_cov.txt'))
    bao1 = m.BAO(os.path.join(DATA, 'desi_dr1_mean.txt'), os.path.join(DATA, 'desi_dr1_cov.txt'))
    sn = m.SN(os.path.join(DATA, 'Pantheon+SH0ES.dat'), os.path.join(DATA, 'Pantheon+SH0ES_STAT+SYS.cov'))

    # same seed and draw order as the paper pipeline
    rng = np.random.default_rng(7)
    S, chimin = {}, {}
    for key, bao, mdl in [('frac_DR2', bao2, 'frac'), ('frac_DR1', bao1, 'frac'),
                          ('lcdm_DR2', bao2, 'lcdm'), ('lcdm_DR1', bao1, 'lcdm')]:
        print(f'grid {key}', flush=True)
        H, Om, zc, cmin = m.sample_grid(*m.grid_chi2(bao, sn, mdl), rng)
        S[key] = (H, Om, zc); chimin[key] = cmin

    rows, dchi2 = build_table(S, bao2)
    summary = dict(
        rows=rows, dchi2_next_forecast=dchi2, chi2_min=chimin,
        zchar_DR2_16_50_84=np.percentile(S['frac_DR2'][2], [16, 50, 84]).tolist(),
        zchar_DR1_16_50_84=np.percentile(S['frac_DR1'][2], [16, 50, 84]).tolist(),
        Hlocal_DR2_16_50_84=np.percentile(m.H_local(S['frac_DR2'][0], S['frac_DR2'][1]), [16, 50, 84]).tolist(),
        Hlocal_DR1_16_50_84=np.percentile(m.H_local(S['frac_DR1'][0], S['frac_DR1'][1]), [16, 50, 84]).tolist())
    json.dump(summary, open(os.path.join(OUT, 'table15.json'), 'w'), indent=1, default=float)
    np.savez_compressed(os.path.join(OUT, 'posterior_samples.npz'),
                        **{f'{k}_{n}': v for k, t in S.items() for n, v in zip(('H', 'Om', 'zc'), t)})

    for r in rows:
        print(f"{r['z']:.3f} {r['quantity'][:2]}  frac {r['frac_median']:.3f} "
              f"95% [{r['frac_95_union'][0]:.3f}, {r['frac_95_union'][1]:.3f}]  LCDM {r['lcdm_median']:.3f}  "
              f"sig_next {r['sigma_next']:.3f}  {r['delta_over_sigma_next']:+.1f}")
    print(f'Delta chi2 (forecast) = {dchi2:.2f}')

    frozen = os.path.join(FROZEN, 'table15.json')
    if os.path.exists(frozen):
        F = json.load(open(frozen))['rows']
        diff = max(abs(a['frac_median'] - b['frac_median']) for a, b in zip(F, rows))
        print(f'max |frac median - frozen| = {diff:.2e}')


if __name__ == '__main__':
    main()
