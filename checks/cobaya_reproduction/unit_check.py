import sys, numpy as np
import os; sys.path.insert(0, os.environ.get('PAPER_CODE', '../../')); import fraccosmo as fc
from cobaya.model import get_model
from common import info
bao, sn, cmb = fc.BAO(), fc.SN(), fc.CMB()
pts=[('lcdm',None,None,68.5,0.30),('A','exp',None,62.3,0.366,np.log(0.52)),('A','pl',None,62.7,0.36,2.53),('A','exp',1/2.262,64.0,0.34)]
for p in pts:
    fam,kind,fix=p[:3]; H,Om=p[3],p[4]
    m=get_model(info(fam,kind or 'exp',fix,data=('bao','sn','cmb')))
    pv={'H0':H,'omegam':Om}
    shape=None
    if fam!='lcdm':
        if fix is None:
            if kind=='exp': pv['lnzc']=p[5]; shape=np.exp(p[5])
            else: pv['n']=p[5]; shape=p[5]
        else: shape=fix
    ll=m.loglikes(pv,as_dict=True)[0]
    c={k:-2*v for k,v in ll.items()}
    ref={'bao':bao.chi2(fam,kind,np.array([H]),np.array([Om]),None if shape is None else np.array([shape]))[0],
         'sn':sn.chi2(fam,kind,np.array([H]),np.array([Om]),None if shape is None else np.array([shape]))[0],
         'cmb':cmb.chi2(fam,kind,np.array([H]),np.array([Om]),None if shape is None else np.array([shape]))[0]}
    print(fam,kind,fix, {k:round(v,4) for k,v in c.items()}); print('   paper code:', {k:round(v,4) for k,v in ref.items()})
