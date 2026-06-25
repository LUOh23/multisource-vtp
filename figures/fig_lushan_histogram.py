"""
图9 — 庐山12.5m DEM 测地/欧氏距离比值分布直方图
加载已计算的lushan_results.npz
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

print('Loading Lushan results...')
data = np.load('e:/vtp_geodesic/lushan_results.npz')
ratio_pct = data['ratio_pct']
valid_count = len(ratio_pct)

print(f'Valid: {valid_count:,}, mean={np.mean(ratio_pct):.2f}%, max={np.max(ratio_pct):.2f}%')

# ── Plot ──
fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))

# Bins: 0-20% range covers most data, long tail beyond
xmax = min(np.percentile(ratio_pct, 99.5), 25)
bins = np.linspace(0, xmax, 100)
counts, edges = np.histogram(ratio_pct, bins=bins)
centers = (edges[:-1] + edges[1:]) / 2
width = (edges[1] - edges[0]) * 0.9

for i in range(len(counts)):
    c = centers[i]
    if c < 2.0:
        ax.bar(c, counts[i], width=width, color='#2E8B57', alpha=0.82,
               edgecolor='white', lw=0.3, zorder=3)
    elif c < 5.0:
        ax.bar(c, counts[i], width=width, color='#DAA520', alpha=0.82,
               edgecolor='white', lw=0.3, zorder=3)
    else:
        ax.bar(c, counts[i], width=width, color='#B22222', alpha=0.82,
               edgecolor='white', lw=0.3, zorder=3)

# Threshold lines
ax.axvline(x=2.0, color='#DAA520', ls='--', lw=1.8, alpha=0.7)
ax.axvline(x=5.0, color='#B22222', ls=':', lw=1.5, alpha=0.5)

# Zone labels at top
ymax = ax.get_ylim()[1]
ax.text(1.0, ymax * 0.94, '< 2%\n(Minor)', ha='center', fontsize=9,
        color='#2E8B57', fontweight='bold')
ax.text(3.5, ymax * 0.94, '2%–5%\n(Moderate)', ha='center', fontsize=9,
        color='#B8860B', fontweight='bold')
ax.text(8.0, ymax * 0.94, '> 5%\n(Significant)', ha='center', fontsize=9,
        color='#B22222', fontweight='bold')

# Stats box
f_lt2 = 100 * (ratio_pct < 2).sum() / len(ratio_pct)
f_2to5 = 100 * ((ratio_pct >= 2) & (ratio_pct < 5)).sum() / len(ratio_pct)
f_gt5 = 100 * (ratio_pct > 5).sum() / len(ratio_pct)
f_gt10 = 100 * (ratio_pct > 10).sum() / len(ratio_pct)
stats = (f'Mean: {np.mean(ratio_pct):.2f}%\n'
         f'Max:  {np.max(ratio_pct):.2f}%\n'
         f'Std:  {np.std(ratio_pct):.2f}%\n'
         f'<2%:  {f_lt2:.1f}%\n'
         f'2–5%: {f_2to5:.1f}%\n'
         f'>5%:  {f_gt5:.1f}%\n'
         f'>10%: {f_gt10:.1f}%')
ax.text(0.97, 0.95, stats, transform=ax.transAxes, fontsize=8.5, va='top', ha='right',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#FFF8DC', edgecolor='#CCAA66', alpha=0.92),
        fontfamily='monospace')

# Legend
legend_el = [
    Patch(facecolor='#2E8B57', alpha=0.82, label='< 2%   Minor deviation'),
    Patch(facecolor='#DAA520', alpha=0.82, label='2%–5%  Moderate'),
    Patch(facecolor='#B22222', alpha=0.82, label='> 5%   Significant'),
]
ax.legend(handles=legend_el, loc='upper center', fontsize=8.5, ncol=3,
          framealpha=0.9, edgecolor='gray')

ax.set_xlabel(r'$(d_g / d_e - 1)$ / %', fontsize=12, fontweight='bold')
ax.set_ylabel('Vertex Count', fontsize=12, fontweight='bold')
ax.set_xlim(0, xmax)
ax.grid(True, alpha=0.2, color='gray', lw=0.3, axis='y')

plt.tight_layout()
plt.savefig('e:/vtp_geodesic/论文截图/fig9_histogram.png', dpi=400,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('Fig 9 (histogram) saved.')
