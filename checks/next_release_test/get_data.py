"""Download the data used to build the published prediction into ./data/.

DESI DR1 and DR2 BAO (CobayaSampler/bao_data) and Pantheon+ (PantheonPlusSH0ES/DataRelease).
Only needed for make_prediction.py and for check_next_release.py --refit.
"""
import os
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')

BAO = 'https://raw.githubusercontent.com/CobayaSampler/bao_data/master/'
PP = ('https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/'
      'Pantheon%2B_Data/4_DISTANCES_AND_COVAR/')

FILES = {
    'desi_dr2_mean.txt': BAO + 'desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt',
    'desi_dr2_cov.txt': BAO + 'desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt',
    'desi_dr1_mean.txt': BAO + 'desi_2024_gaussian_bao_ALL_GCcomb_mean.txt',
    'desi_dr1_cov.txt': BAO + 'desi_2024_gaussian_bao_ALL_GCcomb_cov.txt',
    'Pantheon+SH0ES.dat': PP + 'Pantheon%2BSH0ES.dat',
    'Pantheon+SH0ES_STAT+SYS.cov': PP + 'Pantheon%2BSH0ES_STAT%2BSYS.cov',
}


def main(force=False):
    os.makedirs(DATA, exist_ok=True)
    for name, url in FILES.items():
        path = os.path.join(DATA, name)
        if os.path.exists(path) and not force:
            print(f'have  {name}')
            continue
        print(f'get   {name}', flush=True)
        urllib.request.urlretrieve(url, path)


if __name__ == '__main__':
    main(force='--force' in sys.argv)
