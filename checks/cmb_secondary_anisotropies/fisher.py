import numpy as np, json, sys
from cmbsec import base, spectra, P0, LMAX
import camb

zc = sys.argv[1]
S = np.load(f'spec_{float(zc):.3f}.npz')
ell = np.arange(LMAX+1)
FSKY, FWHM, DT, DP = 0.7, 5.0, 33.0, 70.0   # Planck-like: arcmin, uK-arcmin
th = np.radians(FWHM/60); b2 = np.exp(ell*(ell+1)*th**2/(8*np.log(2)))
NT = (np.radians(DT/60))**2*b2; NP = (np.radians(DP/60))**2*b2

def cl(sp):  # TT, EE, TE (raw Cl)
    return sp[:, 0], sp[:, 1], sp[:, 3]

def inv_cov(TT, EE, TE):
    a, d, b = TT+NT, EE+NP, TE
    det = a*d-b*b
    return a, d, b, det

def trprod(C, X, Y):
    # sum_l fsky(2l+1)/2 Tr[C^-1 X C^-1 Y], X,Y = (TT,EE,TE)
    a, d, b, det = C
    out = 0.0
    for l in range(2, LMAX+1):
        Ci = np.array([[d[l], -b[l]], [-b[l], a[l]]])/det[l]
        Xm = np.array([[X[0][l], X[2][l]], [X[2][l], X[1][l]]])
        Ym = np.array([[Y[0][l], Y[2][l]], [Y[2][l], Y[1][l]]])
        out += FSKY*(2*l+1)/2*np.trace(Ci@Xm@Ci@Ym)
    return out

L = cl(S['L']); F = cl(S['F'])
C = inv_cov(*L)
names = ['ombh2', 'omch2', 'theta', 'lnAs', 'ns', 'tau']
step = dict(ombh2=1e-4, omch2=6e-4, theta=1e-6, lnAs=0.01, ns=0.003, tau=0.004)
D = []
for n in names:
    s = []
    for sg in (+1, -1):
        p = dict(P0); p[n] += sg*step[n]
        s.append(cl(spectra(base(p))[1]))
    D.append(tuple((s[0][k]-s[1][k])/(2*step[n]) for k in range(3)))
dC = tuple(F[k]-L[k] for k in range(3))
Fm = np.array([[trprod(C, D[i], D[j]) for j in range(6)] for i in range(6)])
Fm[5, 5] += 1/0.007**2          # low-l EE tau prior
B = np.array([trprod(C, D[i], dC) for i in range(6)])
cov = np.linalg.inv(Fm); bias = cov@B
sig = np.sqrt(np.diag(cov))
J = np.array([1, 1, 0, 0, 0, 0])  # omega_m = ombh2 + omch2 (+ fixed omnu)
dwm = J@bias; swm = np.sqrt(J@cov@J)
chi2_diff = trprod(C, dC, dC); chi2_refit = chi2_diff-B@cov@B
# lensing and low-l ISW diagnostics
rphi = S['phiF']/S['phiL']
lw = slice(40, 401)
cv = np.sqrt(2/(2*ell+1))
isw = (S['Fu'][:, 0]-S['Lu'][:, 0])/S['Lu'][:, 0]
res = dict(zc=float(zc),
    lensing_ratio_mean_l40_400=float(rphi[lw].mean()-1),
    lensing_ratio_range_l10_1000=[float(rphi[10:1001].min()-1), float(rphi[10:1001].max()-1)],
    lowl_TT_unlensed_frac_change_l2_30=[float(isw[2:31].min()), float(isw[2:31].max())],
    lowl_TT_change_over_cosmic_variance_max=float(np.max(np.abs(isw[2:31])/cv[2:31])),
    bias_in_sigma=dict(zip(names, (bias/sig).round(3).tolist())),
    sigma=dict(zip(names, sig.tolist())),
    omega_m=dict(bias=float(dwm), sigma=float(swm), bias_over_sigma=float(dwm/swm)),
    chi2_fractal_vs_LCDM_same_params=float(chi2_diff), chi2_after_LCDM_refit=float(chi2_refit))
json.dump(res, open(f'fisher_{float(zc):.3f}.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
