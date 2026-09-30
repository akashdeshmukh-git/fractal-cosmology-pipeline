import sys, os, json, copy, numpy as np
from cobaya.run import run
from common import info, LNZC

def minimize(name, seed_ranges=None, best_of=10, **kw):
    """multi-start minimizer; seed_ranges overrides 'ref' of the shape param to cover all modes"""
    inf = info(**kw)
    if seed_ranges:
        for p, (a, b) in seed_ranges.items():
            inf['params'][p]['ref'] = {'dist': 'uniform', 'min': a, 'max': b}
    inf['params']['H0']['ref'] = {'dist': 'uniform', 'min': 55, 'max': 75}
    inf['params']['omegam']['ref'] = {'dist': 'uniform', 'min': 0.25, 'max': 0.40}
    inf['sampler'] = {'minimize': {'method': 'bobyqa', 'best_of': best_of, 'seed': 1,
                                   'override_bobyqa': {'rhoend': 1e-9}}}
    inf['output'] = f'chains/{name}'; inf['force'] = True
    upd, s = run(inf)
    hdr = open(f'chains/{name}.minimum.txt').readline().lstrip('#').split()
    val = np.loadtxt(f'chains/{name}.minimum.txt', ndmin=2)[0]
    m = dict(zip(hdr, val))
    res = {'name': name, 'minus_logpost': m['minuslogpost'], 'chi2': m['chi2']}
    res['chi2_parts'] = {k: m[k] for k in m if k.startswith('chi2__') and '.' in k}
    res['params'] = {k: m[k] for k in ['H0', 'omegam', 'lnzc', 'n', 'zchar', 'H0_local', 'hrdrag', 'omegamh2'] if k in m}
    json.dump(res, open(f'results/min_{name}.json', 'w'), indent=1)
    return res

def mcmc(name, Rstop=0.003, **kw):
    inf = info(**kw)
    inf['sampler'] = {'mcmc': {'Rminus1_stop': Rstop, 'Rminus1_cl_stop': 0.2, 'max_tries': 100000, 'learn_proposal': True, 'seed': 3}}
    inf['output'] = f'chains/{name}'; inf['force'] = True
    upd, s = run(inf)
    return s
