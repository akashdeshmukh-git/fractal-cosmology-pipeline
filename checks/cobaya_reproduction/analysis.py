import os
"""Collect Cobaya results, compare with the paper, write tables and plots."""
import json, glob, os, re, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from getdist import loadMCSamples

R = 'results'; os.makedirs('plots', exist_ok=True)
PAPER_T12 = {  # Table 12 of the revised paper: chi2_min (BAO, BAO+SN, BAO+SN+CMB)
    'LCDM': (10.28, 1470.10, 1471.77), 'exp_fix': (9.58, 1466.19, 1467.88), 'pl_fix': (9.29, 1465.80, 1471.11),
    'exp_free': (8.30, 1463.96, 1464.89), 'pl_free': (9.23, 1465.12, 1469.50)}
LABEL = {'LCDM': 'ΛCDM', 'exp_fix': 'exp, z_char = 1/2.262 (fixed)', 'pl_fix': 'power law, n = 2.262 (fixed)',
         'exp_free': 'exp, z_char free', 'pl_free': 'power law, n free'}
DS = ['BAO', 'BAOSN', 'BAOSNCMB']; DSL = {'BAO': 'BAO', 'BAOSN': 'BAO+SN', 'BAOSNCMB': 'BAO+SN+CMB'}

def mn(tag):
    f = f'{R}/min_{tag}.json'
    return json.load(open(f)) if os.path.exists(f) else None

# ---------------- chi2 table
rows = []
for m in PAPER_T12:
    for i, d in enumerate(DS):
        p = PAPER_T12[m][i]; cp = mn(f'{m}_{d}_paper'); cc = mn(f'{m}_{d}_cobaya') if d != 'BAO' else cp
        pl = PAPER_T12['LCDM'][i]; lp = mn(f'LCDM_{d}_paper'); lc = mn(f'LCDM_{d}_cobaya') if d != 'BAO' else lp
        rows.append(dict(model=m, data=DSL[d], paper=p, cobaya=cp['chi2'], diff=cp['chi2'] - p,
                         dpaper=p - pl, dcob=cp['chi2'] - lp['chi2'], dcob_std=cc['chi2'] - lc['chi2'],
                         shape=cp['params'].get('zchar', cp['params'].get('n')), shape_std=cc['params'].get('zchar', cc['params'].get('n')),
                         H0=cp['params']['H0'], Om=cp['params']['omegam'], Hloc=cp['params']['H0_local']))
with open(f'{R}/comparison_chi2.csv', 'w') as f:
    f.write('model,data,chi2_paper,chi2_cobaya_papercut,difference,dchi2_paper,dchi2_cobaya_papercut,dchi2_cobaya_standardcut,bestfit_shape_papercut,bestfit_shape_standardcut,H_E,Omega_m,H_local\n')
    for r in rows:
        f.write(f"{r['model']},{r['data']},{r['paper']:.2f},{r['cobaya']:.3f},{r['diff']:+.3f},{r['dpaper']:+.2f},{r['dcob']:+.3f},{r['dcob_std']:+.3f},"
                f"{'' if r['shape'] is None else round(r['shape'],3)},{'' if r['shape_std'] is None else round(r['shape_std'],3)},{r['H0']:.2f},{r['Om']:.4f},{r['Hloc']:.2f}\n")

# ---------------- profiles
def prof(kind):
    out = {}
    for f in glob.glob(f'{R}/min_prof_{kind}_*.json'):
        _, _, _, d, v = os.path.basename(f)[:-5].split('_')
        r = json.load(open(f)); v = float(v)
        out.setdefault(d, []).append((np.exp(v) if kind == 'exp' else v, r['chi2']))
    return {d: np.array(sorted(v)) for d, v in out.items()}
PE, PP = prof('exp'), prof('pl')
G = json.load(open(os.path.join(os.environ.get('PAPER_CODE', '../../'), 'results/grid_post.json'))); GC = json.load(open(os.path.join(os.environ.get('PAPER_CODE', '../../'), 'results/grid_post_cmb.json')))
GK = {'BAO': ('A-{}_bao_flat', G), 'BAOSN': ('A-{}_bao+sn_flat', G), 'BAOSNCMB': ('A-{}_bao+sn+cmb_flat', GC)}

plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
COL = {'BAO': '#8a8f98', 'BAOSN': '#2a6fdb', 'BAOSNCMB': '#d9731a'}
fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.6))
for ax, kind, P in [(axs[0], 'exp', PE), (axs[1], 'pl', PP)]:
    for d in DS:
        a = P[d]; lc = mn(f'LCDM_{d}_paper')['chi2']
        ax.plot(a[:, 0], a[:, 1] - lc, 'o-', ms=3, lw=1.3, color=COL[d], label=f'{DSL[d]} (Cobaya)')
        key, src = GK[d]; g = src[key.format(kind)]
        X = np.array(g['X']); x = np.exp(X) if kind == 'exp' else X
        ax.plot(x, np.array(g['prof_X']) - PAPER_T12['LCDM'][DS.index(d)], '--', lw=0.9, color=COL[d], alpha=0.8)
    ax.axhline(0, color='k', lw=0.5)
    ax.set_xscale('log'); ax.set_ylim(-8, 6)
    ax.set_xlabel(r'$z_{\rm char}$' if kind == 'exp' else r'$n$')
    ax.set_ylabel(r'$\chi^2_{\rm min}(\mathrm{shape})-\chi^2_{\Lambda{\rm CDM}}$')
    if kind == 'exp': ax.axvspan(0.3, 0.7, color='#2fa36b', alpha=0.08, lw=0)
    else: ax.axvline(2.262, color='#2fa36b', lw=0.8, ls=':')
    ax.set_title('exponential $f=e^{-z/z_{\\rm char}}$' if kind == 'exp' else 'power law $f=(1+z)^{-n}$', fontsize=9)
axs[0].legend(frameon=False, fontsize=7.5, loc='lower left')
axs[1].text(0.98, 0.04, 'solid: Cobaya profile (minimizer)\ndashed: paper grid profile', transform=axs[1].transAxes, ha='right', fontsize=7, color='#555')
fig.tight_layout(); fig.savefig('plots/profiles.png', dpi=200); plt.close(fig)

# ---------------- posteriors
def load(name):
    s = loadMCSamples(f'chains/{name}', settings={'ignore_rows': 0.3})
    return s
post = {}
fig, axs = plt.subplots(1, 3, figsize=(11, 3.3))
for ax, name, par, key, src, xl in [
        (axs[0], 'mc_exp_free_BAOSNCMB', 'zchar', 'A-exp_bao+sn+cmb_flat', GC, r'$z_{\rm char}$ (BAO+SN+CMB)'),
        (axs[1], 'mc_exp_free_BAOSN', 'zchar', 'A-exp_bao+sn_flat', G, r'$z_{\rm char}$ (BAO+SN)'),
        (axs[2], 'mc_pl_free_BAOSNCMB', 'n', 'A-pl_bao+sn+cmb_flat', GC, r'$n$ (BAO+SN+CMB)')]:
    if not glob.glob(f'chains/{name}.1.txt'):
        continue
    s = load(name); x = s[par]; w = s.weights
    lo, md, hi = [float(np.interp(q, np.cumsum(w[np.argsort(x)]) / w.sum(), np.sort(x))) for q in (0.16, 0.5, 0.84)]
    R1 = None
    post[name] = dict(par=par, median=md, lo=lo, hi=hi, N=int(s.numrows), neff=float(w.sum() ** 2 / (w ** 2).sum()))
    if par == 'zchar':
        post[name]['frac_0.3_0.7'] = float(w[(x > 0.3) & (x < 0.7)].sum() / w.sum())
    g = src[key]; X = np.array(g['X']); gx = np.exp(X) if par == 'zchar' else X
    gp = np.array(g['post_X'])
    if par == 'zchar':
        bins = np.geomspace(0.02, 3, 60); dens_g = gp / gx  # posterior per unit z from per unit ln z
    else:
        bins = np.linspace(0.1, 10, 100); dens_g = gp
    ax.hist(x, bins=bins, weights=w, density=True, histtype='stepfilled', color='#2a6fdb', alpha=0.35, label='Cobaya MCMC')
    ax.plot(gx, dens_g / np.trapezoid(dens_g, gx), color='#d9731a', lw=1.3, label='paper (grid)')
    ax.set_xscale('log'); ax.set_xlabel(xl); ax.set_yticks([])
    if par == 'zchar': ax.axvspan(0.3, 0.7, color='#2fa36b', alpha=0.08, lw=0); ax.set_xlim(0.02, 3)
    else: ax.axvline(2.262, color='#2fa36b', lw=0.8, ls=':'); ax.set_xlim(0.5, 10)
axs[0].legend(frameon=False, fontsize=7.5)
fig.tight_layout(); fig.savefig('plots/posteriors.png', dpi=200); plt.close(fig)

# LCDM validation
s = load('lcdm_bao_mcmc')
om, hr = s['omegam'], s['hrdrag']; w = s.weights
stat = lambda a: (float(np.average(a, weights=w)), float(np.sqrt(np.cov(a, aweights=w))))
post['lcdm'] = {'omegam': stat(om), 'hrdrag': stat(hr), 'H0': stat(s['H0']), 'N': int(s.numrows)}
fig, ax = plt.subplots(figsize=(3.8, 3.4))
ax.scatter(om, hr, s=2, c='#2a6fdb', alpha=0.25, lw=0)
ax.errorbar([0.2975], [101.54], xerr=[0.0086], yerr=[0.73], fmt='o', color='#d9731a', ms=4, capsize=2, label='DESI DR2 published')
ax.errorbar([post['lcdm']['omegam'][0]], [post['lcdm']['hrdrag'][0]], xerr=[post['lcdm']['omegam'][1]], yerr=[post['lcdm']['hrdrag'][1]],
            fmt='s', color='k', ms=3.5, capsize=2, label='this Cobaya run')
ax.set_xlabel(r'$\Omega_m$'); ax.set_ylabel(r'$h\,r_d$ [Mpc]'); ax.legend(frameon=False, fontsize=7.5)
ax.set_title('ΛCDM, DESI DR2 BAO only', fontsize=9)
fig.tight_layout(); fig.savefig('plots/lcdm_validation.png', dpi=200); plt.close(fig)

# delta chi2 comparison plot
fig, ax = plt.subplots(figsize=(7.5, 3.4))
ms = ['exp_fix', 'pl_fix', 'exp_free', 'pl_free']
for j, d in enumerate(DSL.values()):
    for i, m in enumerate(ms):
        r = [r for r in rows if r['model'] == m and r['data'] == d][0]
        x0 = j * 5 + i
        ax.plot([x0 - 0.2], [r['dpaper']], 'o', color='#d9731a', ms=4.5, label='paper' if (i == 0 and j == 0) else None)
        ax.plot([x0], [r['dcob']], 's', color='#2a6fdb', ms=4, label='Cobaya, paper SN cut (1624)' if (i == 0 and j == 0) else None)
        ax.plot([x0 + 0.2], [r['dcob_std']], '^', color='#555', ms=4, label='Cobaya, standard SN cut (1590)' if (i == 0 and j == 0) else None)
ax.axhline(0, color='k', lw=0.5)
ax.set_xticks([j * 5 + i for j in range(3) for i in range(4)])
ax.set_xticklabels(['exp fix', 'pl fix', 'exp free', 'pl free'] * 3, rotation=60, fontsize=7.5)
for j, d in enumerate(DSL.values()): ax.text(j * 5 + 1.5, 1.6, d, ha='center', fontsize=8.5)
ax.set_ylabel(r'$\Delta\chi^2$ vs ΛCDM'); ax.set_ylim(-8, 2.4); ax.legend(frameon=False, fontsize=7.5, loc='upper center', bbox_to_anchor=(0.5, -0.32), ncol=3)
fig.tight_layout(); fig.savefig('plots/delta_chi2.png', dpi=200); plt.close(fig)

json.dump({'rows': rows, 'posteriors': post}, open(f'{R}/summary.json', 'w'), indent=1, default=float)
for r in rows: print(r['model'], r['data'], r['paper'], round(r['cobaya'], 3), f"{r['diff']:+.3f}", f"{r['dpaper']:+.2f}", f"{r['dcob']:+.3f}", f"{r['dcob_std']:+.3f}")
print(json.dumps(post, indent=1, default=float))
