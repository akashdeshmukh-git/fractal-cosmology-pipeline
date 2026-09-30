from runner import *
import sys
if 'min' in sys.argv: print('MIN', minimize('lcdm_bao', data=('bao',)))
if 'mcmc' in sys.argv: mcmc('lcdm_bao_mcmc', data=('bao',))
