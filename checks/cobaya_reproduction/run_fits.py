import sys, json, numpy as np
from runner import minimize, mcmc, LNZC
DS = {'BAO': ('bao',), 'BAO+SN': ('bao', 'sn'), 'BAO+SN+CMB': ('bao', 'sn', 'cmb')}
MODELS = {'LCDM': dict(family='lcdm'), 'exp_fix': dict(family='A', kind='exp', shape_fixed=1/2.262),
          'pl_fix': dict(family='A', kind='pl', shape_fixed=2.262),
          'exp_free': dict(family='A', kind='exp'), 'pl_free': dict(family='A', kind='pl')}
SEED = {'exp_free': {'lnzc': LNZC}, 'pl_free': {'n': (0.1, 10.0)}}
def table(cut):
    for mname, mk in MODELS.items():
        for dname, d in DS.items():
            if cut == 'cobaya' and 'sn' not in d: continue
            tag = f'{mname}_{dname.replace("+","")}_{cut}'
            r = minimize(tag, seed_ranges=SEED.get(mname), best_of=16 if mname in SEED else 6, data=d, sn_cut=cut, **mk)
            print('RES', tag, round(r['chi2'], 3), r['params'], flush=True)
def profile(kind):
    vals = np.log(np.geomspace(0.02, 3.0, 22)) if kind == 'exp' else np.array([0.1,0.3,0.6,1,1.5,2,2.3,2.6,3,3.5,4,5,6,8,10,15,25,50,100])
    for dname, d in DS.items():
        for v in vals:
            shape = float(np.exp(v)) if kind == 'exp' else float(v)
            tag = f'prof_{kind}_{dname.replace("+","")}_{v:.4f}'
            r = minimize(tag, best_of=3, data=d, sn_cut='paper', family='A', kind=kind, shape_fixed=shape)
            print('PROF', kind, dname, shape, round(r['chi2'], 3), flush=True)
if __name__ == '__main__':
    job = sys.argv[1]
    if job == 'A':
        table('paper'); profile('exp')
        mcmc('mc_exp_free_BAOSNCMB', data=DS['BAO+SN+CMB'], family='A', kind='exp')
        mcmc('mc_exp_free_BAOSN', data=DS['BAO+SN'], family='A', kind='exp')
    else:
        table('cobaya'); profile('pl')
        mcmc('mc_pl_free_BAOSNCMB', data=DS['BAO+SN+CMB'], family='A', kind='pl')
        mcmc('lcdm_bao_mcmc', data=('bao',))
