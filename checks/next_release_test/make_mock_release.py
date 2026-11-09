"""
Make a mock next release, to see how the test behaves before real data exist.

The mock is built the way a real release is: the real DR2 measurement plus independent new data drawn
from one model, combined with inverse-covariance weights, with the release covariance equal to
FACTOR x the DR2 covariance (default 3/5, the precision forecast in the paper).

    python make_mock_release.py lcdm           # writes mock/mock_lcdm_mean.txt and mock/mock_lcdm_cov.txt
    python make_mock_release.py frac --seed 3
    python check_next_release.py mock/mock_lcdm_mean.txt mock/mock_lcdm_cov.txt --label mock_lcdm
"""
import os
import argparse
import numpy as np
import check_next_release as c

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('truth', choices=('lcdm', 'frac'))
    ap.add_argument('--seed', type=int, default=21)
    ap.add_argument('--factor', type=float, default=0.6)
    a = ap.parse_args()
    S = c.load_samples()
    dr2 = c.Data.from_files(os.path.join(c.FROZEN, 'desi_dr2_mean.txt'), os.path.join(c.FROZEN, 'desi_dr2_cov.txt'))
    F2 = np.linalg.inv(dr2.cov)
    C3 = a.factor * dr2.cov
    Fn = np.linalg.inv(C3) - F2
    Cn = np.linalg.inv(Fn)
    rng = np.random.default_rng(a.seed)
    P = c.predictive(S, f'{a.truth}_DR2', dr2.z, dr2.qty)
    xn = P[rng.integers(0, len(P))] + rng.multivariate_normal(np.zeros(dr2.n), Cn)
    x3 = C3 @ (F2 @ dr2.obs + Fn @ xn)
    out = os.path.join(HERE, 'mock')
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, f'mock_{a.truth}_mean.txt'), 'w') as f:
        f.write('# z value quantity  (MOCK release: DR2 + new data drawn from the %s prediction)\n' % a.truth)
        for z, v, q in zip(dr2.z, x3, dr2.qty):
            f.write(f'{z:.3f} {v:.8f} {q}\n')
    np.savetxt(os.path.join(out, f'mock_{a.truth}_cov.txt'), C3)
    print(f'wrote mock/mock_{a.truth}_mean.txt and mock/mock_{a.truth}_cov.txt')


if __name__ == '__main__':
    main()
