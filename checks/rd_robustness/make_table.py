"""Build results/summary_table.md from the JSON results."""
import json

MODES = [('fixed', 'r_d fixed at 147.09 Mpc (paper)'),
         ('free', 'r_d free, CMB r* fixed'),
         ('free_consistent', 'r_d free, CMB r* scales with r_d'),
         ('computed', 'r_d computed from omega_m (physical)')]
MODELS = ['exp_fix', 'pl_fix', 'exp_free', 'pl_free']
lines = []
for cut, label in [('paper', 'paper SN cut (1624 SNe)'), ('standard', 'standard SN cut, z > 0.01 (1590 SNe)')]:
    a = json.load(open(f'results/cmb_{cut}.json'))
    c = json.load(open(f'results/computed_{cut}.json'))
    lines += [f'### BAO+SN+CMB, {label}', '',
              '| r_d treatment | exp fixed | power law fixed | exp free | power law free | z_char (exp free) | H_local (exp free) | r_d LCDM | r_d power law fixed |',
              '|---|---|---|---|---|---|---|---|---|']
    for key, name in MODES:
        d = c if key == 'computed' else {k.split('|')[0]: v for k, v in a.items() if k.endswith('|' + key)}
        L = d['LCDM']['chi2']
        cells = [f"{d[m]['chi2'] - L:+.2f}" for m in MODELS]
        e = d['exp_free']
        lines.append(f"| {name} | " + ' | '.join(cells) +
                     f" | {e['shape']:.3f} | {e['H_local']:.2f} | {d['LCDM']['rd']:.2f} | {d['pl_fix']['rd']:.2f} |")
    lines.append('')
open('results/summary_table.md', 'w').write('\n'.join(lines))
print('\n'.join(lines))
