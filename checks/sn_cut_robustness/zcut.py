import sys, numpy as np, warnings; warnings.filterwarnings('ignore')
import model
from model import _sn, _cov_full
zmin=float(sys.argv[1])
m=(_sn.zHD>zmin).values
d=_sn[m]; ic=np.linalg.inv(_cov_full[np.ix_(m,m)])
model.SN['z']=dict(zhd=d.zHD.values,zhel=d.zHEL.values,mb=d.m_b_corr.values,icov=ic,s=ic.sum(),N=m.sum())
from fit import fit
for data in ['BAO+SN','BAO+SN+CMB']:
    L=fit('lcdm',data,'fixed','z'); E=fit('exp',data,'fixed','z'); F=fit('exp',data,'fixed','z',fixed_shape=1/2.262)
    print(f"zmin={zmin} N={m.sum()} {data}: dchi2 exp_free={E['chi2']-L['chi2']:+.2f} zc={E['shape']:.3f} Hloc={E['H_local']:.2f} | exp_fix={F['chi2']-L['chi2']:+.2f}",flush=True)
