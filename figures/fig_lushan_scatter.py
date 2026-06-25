"""
图8 — 庐山12.5m DEM 测地vs欧氏距离散点
加载已计算的lushan_results.npz，不重复跑VTP
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

print('Loading Lushan results...')
data = np.load('e:/vtp_geodesic/lushan_results.npz')
d_geo = data['d_geo']
d_euc = data['d_euc']
valid = data['valid']
ratio_pct = data['ratio_pct']
n_verts = len(d_geo)

d_v = d_geo[valid]
e_v = d_euc[valid]
rp = ratio_pct

print(f'Vertices: {n_verts:,}, Valid: {valid.sum():,}')
print(f'Mean r: {np.mean(d_v/(e_v+1e-10)):.4f}, Max r: {np.max(d_v/(e_v+1e-10)):.4f}')

# ── Plot ──
fig, ax = plt.subplots(1, 1, figsize=(7.5, 6.5))

# Random sample for scatter (3000 points, enough to show density)
np.random.seed(42)
samp = np.random.choice(len(d_v), min(3000, len(d_v)), replace=False)

# Color by deviation
colors = rp[samp]
sc = ax.scatter(e_v[samp], d_v[samp], c=colors, cmap='YlOrRd',
                s=8, alpha=0.55, edgecolors='none', zorder=3,
                vmin=0, vmax=15)

lim = max(e_v.max(), d_v.max())
ax.plot([0, lim], [0, lim], '--', color='#333333', lw=2.0,
        label=r'$d_g = d_e$', zorder=4)

ax.set_xlabel('Euclidean Distance / m', fontsize=12, fontweight='bold')
ax.set_ylabel('Geodesic Distance / m', fontsize=12, fontweight='bold')
ax.legend(loc='upper left', fontsize=10, framealpha=0.9, edgecolor='gray')
ax.grid(True, alpha=0.2, color='gray', lw=0.4)
ax.set_xlim(0, lim * 1.02)
ax.set_ylim(0, lim * 1.02)

# Colorbar
cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
cbar.set_label(r'$(d_g/d_e - 1)$ / %', fontsize=11, fontweight='bold')

# Stats box
mean_r = np.mean(d_v / (e_v + 1e-10))
max_r = np.max(d_v / (e_v + 1e-10))
gt2 = 100 * (rp > 2).sum() / len(rp)
gt5 = 100 * (rp > 5).sum() / len(rp)
stats = (f'Vertices: {n_verts:,}\n'
         f'Mean r: {mean_r:.3f}\n'
         f'Max r: {max_r:.3f}\n'
         f'>2%: {gt2:.1f}%\n'
         f'>5%: {gt5:.1f}%')
ax.text(0.97, 0.06, stats, transform=ax.transAxes, fontsize=9, va='bottom', ha='right',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#AAAAAA', alpha=0.9),
        fontfamily='monospace')

plt.tight_layout()
plt.savefig('e:/vtp_geodesic/论文截图/fig8_scatter.png', dpi=400,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('Fig 8 (scatter) saved.')
