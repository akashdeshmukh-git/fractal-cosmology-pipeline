"""
FRACTAL HUBBLE MODEL — COMPLETE ANALYSIS PACKAGE
=================================================
Runs all three remaining analyses in sequence:
  1. Full 3-parameter joint MCMC (DESI + Pantheon+)
  2. LCDM benchmark verification
  3. CMB consistency check (analytic)

Requirements: pip install emcee corner numpy scipy pandas matplotlib

Estimated runtime: 20-40 minutes on a modern laptop.
Results saved to: results/ folder
"""

import subprocess, sys

def install(pkg):
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', pkg, '-q'])

print("Checking dependencies...")
for pkg in ['emcee', 'corner', 'numpy', 'scipy', 'pandas', 'matplotlib']:
    try:
        __import__(pkg)
    except ImportError:
        print(f"Installing {pkg}...")
        install(pkg)

import os
os.makedirs('results', exist_ok=True)

print("\n" + "="*60)
print("Step 1: LCDM Benchmark")
print("="*60)
exec(open('step1_lcdm_benchmark.py').read())

print("\n" + "="*60)
print("Step 2: Full Joint MCMC (fractal model)")
print("="*60)
exec(open('step2_joint_mcmc.py').read())

print("\n" + "="*60)
print("Step 3: CMB Consistency Check")
print("="*60)
exec(open('step3_cmb_check.py').read())

print("\n" + "="*60)
print("ALL DONE. Check the results/ folder.")
print("="*60)
