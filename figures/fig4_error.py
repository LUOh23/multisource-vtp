"""
图4 — Heat Method误差分析 (彩色恢复版)
橙虚线Max Error + 红实线Mean Error + 淡黄高亮收敛区
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

vertices   = np.array([50, 226, 530, 962, 2210, 3970, 6242, 9802, 19322, 32042])
mean_error = np.array([3.35, 0.86, 0.39, 0.22, 0.10, 0.07, 0.05, 0.05, 0.08, 0.12])
max_error  = np.array([7.22, 1.90, 0.85, 0.48, 0.22, 0.13, 0.10, 0.13, 0.17, 0.21])

fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))

# Convergence zone highlight
ax.axvspan(1000, 40000, alpha=0.12, color='#E8C820', zorder=0)
ax.axvline(x=1000, color='#B8960F', linestyle='--', lw=0.8, alpha=0.5, zorder=1)
ax.text(7000, 6.8, 'Convergence Zone', fontsize=9, color='#B8960F', fontstyle='italic', ha='center')

# Max Error
ax.semilogx(vertices, max_error, 's--', color='#D98736', lw=2.0, ms=8,
            markerfacecolor='#F0A050', markeredgecolor='#B8651A', markeredgewidth=1.0,
            label=r'$\varepsilon_{\max}$', zorder=5)
# Mean Error
ax.semilogx(vertices, mean_error, 'o-', color='#C00000', lw=2.5, ms=8,
            markerfacecolor='#E04040', markeredgecolor='#8B0000', markeredgewidth=1.0,
            label=r'$\bar{\varepsilon}$', zorder=5)

# Annotations
ax.annotate('Min: 0.05% @ 6,242 vertices', xy=(6242, 0.05), xytext=(3500, 2.5),
            arrowprops=dict(arrowstyle='->', color='#C00000', lw=1.2),
            fontsize=9, color='#8B0000', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#C00000', alpha=0.9), zorder=10)
ax.annotate('Coarse mesh: 3.35%', xy=(50, 3.35), xytext=(250, 5.8),
            arrowprops=dict(arrowstyle='->', color='#D98736', lw=1.0),
            fontsize=8.5, color='#B8651A', zorder=10)
ax.annotate('Finest: 0.12%', xy=(32042, 0.12), xytext=(18000, 2.5),
            arrowprops=dict(arrowstyle='->', color='#C00000', lw=1.0),
            fontsize=8.5, color='#8B0000', zorder=10)

ax.set_xlabel('Number of Vertices', fontsize=12, fontweight='bold')
ax.set_ylabel('Relative Error / %', fontsize=12, fontweight='bold')
ax.legend(loc='upper right', fontsize=10, framealpha=0.9, edgecolor='gray')
ax.grid(True, which='major', alpha=0.25, color='gray', lw=0.4)
ax.set_xlim(40, 40000); ax.set_ylim(-0.3, 7.5)
ax.set_xticks([50, 100, 500, 1000, 5000, 10000, 32000])
ax.set_xticklabels(['50', '100', '500', '1k', '5k', '10k', '32k'])

plt.tight_layout()
plt.savefig('e:/vtp_geodesic/论文截图/fig3_heat_error.png', dpi=300,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('Fig.4 (error) saved')
