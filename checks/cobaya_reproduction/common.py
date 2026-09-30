import os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.environ.get('COBAYA_PACKAGES', os.path.join(HERE, 'packages'))
LNZC = [np.log(0.02), np.log(3.0)]
def info(family='lcdm', kind='exp', shape_fixed=None, data=('bao',), sn_cut='paper'):
    th = {'fractal_cobaya.FractalBackground': {'python_path': HERE,
          'family': family, 'kind': kind, 'shape_fixed': shape_fixed}}
    params = {'H0': {'prior': {'min': 50, 'max': 90}, 'ref': {'dist': 'norm', 'loc': 66, 'scale': 2}, 'proposal': 0.5, 'latex': 'H_E'},
              'omegam': {'prior': {'min': 0.1, 'max': 0.6}, 'ref': {'dist': 'norm', 'loc': 0.31, 'scale': 0.02}, 'proposal': 0.01, 'latex': r'\Omega_m'},
              'rdrag': {'latex': 'r_d'}, 'omegamh2': {'latex': r'\omega_m'}, 'H0_local': {'latex': r'H_{\rm local}'},
              'hrdrag': {'derived': 'lambda H0: H0/100*147.09', 'latex': r'h r_d'}}
    if family != 'lcdm' and shape_fixed is None:
        if kind == 'exp':
            params['lnzc'] = {'prior': {'min': LNZC[0], 'max': LNZC[1]}, 'ref': {'dist': 'uniform', 'min': -1.5, 'max': 0.0}, 'proposal': 0.2, 'latex': r'\ln z_{\rm char}'}
            params['zchar'] = {'derived': 'lambda lnzc: np.exp(lnzc)', 'latex': r'z_{\rm char}'}
        else:
            params['n'] = {'prior': {'min': 0.1, 'max': 10.0}, 'ref': {'dist': 'uniform', 'min': 1.0, 'max': 4.0}, 'proposal': 0.3, 'latex': 'n'}
    lik = {}
    if 'bao' in data:
        lik['bao.desi_dr2.desi_bao_all'] = None
    if 'sn' in data:
        if sn_cut == 'paper':
            lik['fractal_cobaya.PantheonPlusPaperCut'] = {'python_path': HERE,
                'path': PKG + '/data/sn_data', 'dataset_file': 'PantheonPlus/config.dataset', 'use_abs_mag': False}
        else:
            lik['sn.pantheonplus'] = None
    if 'cmb' in data:
        lik['fractal_cobaya.CompressedCMB'] = {'python_path': HERE}
    return {'theory': th, 'likelihood': lik, 'params': params, 'packages_path': PKG}
