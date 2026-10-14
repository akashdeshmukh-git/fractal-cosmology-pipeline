"""Effect of the fractal late-time correction on CMB secondary anisotropies (ISW, lensing),
and the bias it would induce on early-universe parameters inferred with LCDM."""
import numpy as np, camb, json, sys
from scipy.optimize import brentq

LMAX = 2500
P0 = dict(ombh2=0.02237, omch2=0.1200, theta=0.0104092, lnAs=3.044, ns=0.9649, tau=0.0544)
MNU = 0.06

def base(p, H0=None):
    pars = camb.CAMBparams()
    kw = dict(ombh2=p['ombh2'], omch2=p['omch2'], mnu=MNU, tau=p['tau'])
    if H0 is None: kw['cosmomc_theta'] = p['theta']
    else: kw['H0'] = H0
    pars.set_cosmology(**kw)
    pars.InitPower.set_params(As=np.exp(p['lnAs'])*1e-10, ns=p['ns'])
    pars.set_for_lmax(LMAX+500, lens_potential_accuracy=1)
    pars.set_matter_power(redshifts=[0.], kmax=2.0)
    return pars

def omnuh2(): 
    return MNU/93.14

def fractal_pars(p, HE, zc, A=1.0):
    wm = p['ombh2']+p['omch2']+omnuh2()
    OmE = wm/(HE/100)**2
    Hloc = HE*np.sqrt(1+A*OmE/2)
    pars = base(p, H0=Hloc)
    a = np.logspace(-5, 0, 3000); x = 1/a
    OrE = 2.469e-5*(1+0.2271*3.044)/(HE/100)**2
    OL = 1-OmE-OrE
    corr = A*OmE/2*x**3*np.exp(-(x-1)/zc)
    rho = OL+corr
    dlnrho = A*OmE/2*x**3*np.exp(-(x-1)/zc)*(3-x/zc)/rho
    w = -1+dlnrho/3
    pars.set_dark_energy_w_a(a, w, dark_energy_model='ppf')
    return pars, dict(HE=HE, OmE=OmE, Hloc=Hloc, OL=OL, wm=wm)

def thetastar(pars):
    return camb.get_background(pars).get_derived_params()['thetastar']

def spectra(pars):
    r = camb.get_results(pars)
    tot = r.get_cmb_power_spectra(pars, CMB_unit='muK', raw_cl=True)
    return r, tot['total'][:LMAX+1], tot['unlensed_scalar'][:LMAX+1], tot['lens_potential'][:LMAX+1, 0]

if __name__ == '__main__':
    zc = float(sys.argv[1]) if len(sys.argv) > 1 else 0.542
    lp = base(P0); ts0 = thetastar(lp)
    HE = brentq(lambda h: thetastar(fractal_pars(P0, h, zc)[0])-ts0, 55, 72, xtol=1e-6)
    fp, info = fractal_pars(P0, HE, zc)
    rL, L, Lu, phiL = spectra(lp); rF, F, Fu, phiF = spectra(fp)
    # background check against the paper's formula
    z = np.array([0, .3, .5, 1, 2])
    OrE = 2.469e-5*(1+0.2271*3.044)/(HE/100)**2
    E2 = info['OmE']*(1+z)**3*(1+np.exp(-z/zc)/2)+OrE*(1+z)**4+info['OL']
    Hform = HE*np.sqrt(E2); Hcamb = rF.hubble_parameter(z)
    out = dict(zc=zc, info=info, H0_LCDM=lp.H0, thetastar=ts0,
               H_check=dict(z=z.tolist(), rel_diff=((Hcamb-Hform)/Hform).tolist()),
               sigma8=dict(LCDM=float(rL.get_sigma8_0()), fractal=float(rF.get_sigma8_0())))
    np.savez(f'spec_{zc:.3f}.npz', L=L, F=F, Lu=Lu, Fu=Fu, phiL=phiL, phiF=phiF)
    json.dump(out, open(f'info_{zc:.3f}.json', 'w'), indent=1)
    print(json.dumps(out, indent=1))
