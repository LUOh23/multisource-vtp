"""
图4 重绘: 多源加速比曲线
《测绘科学》标准: 红实线+数据点, 标签对齐, y=x灰色参考线
Data from Table 3
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 10
plt.rcParams['mathtext.fontset'] = 'stix'

# ── Data from Table 3 ──
k_values = np.array([1, 2, 5, 10, 20, 50])
eta_values = np.array([1.1, 1.9, 4.6, 9.6, 23.8, 50.7])
labels = ['1.1×', '1.9×', '4.6×', '9.6×', '23.8×', '50.7×']

fig, ax = plt.subplots(1, 1, figsize=(7.5, 5.5))

# ── y=x ideal reference line ──
k_ideal = np.linspace(0, 55, 100)
ax.plot(k_ideal, k_ideal, '--', color='#AAAAAA', lw=1.5, alpha=0.7,
        label=r'$\eta = k$ (ideal linear)', zorder=2)
ax.fill_between(k_ideal, k_ideal*0.9, k_ideal*1.1, alpha=0.05, color='gray', zorder=0)
ax.text(42, 39, 'Ideal $\pm$10%', fontsize=8, color='#AAAAAA', style='italic')

# ── Speedup curve: red solid with filled markers ──
ax.plot(k_values, eta_values, 'o-', color='#8B0000', lw=2.5, ms=10,
        markerfacecolor='#C00000', markeredgecolor='#600000', markeredgewidth=1.2,
        label=r'$\eta(k)$ measured', zorder=5)

# ── Labels above each data point ──
offsets = [(-0.5, 0.8), (0.5, 2.0), (0.5, 1.2), (0.5, 1.5), (1.0, 1.5), (1.5, 0.5)]
for i, (k, eta, lbl) in enumerate(zip(k_values, eta_values, labels)):
    ox, oy = offsets[i]
    ax.annotate(lbl, (k, eta), xytext=(k + ox, eta + oy),
                fontsize=10, fontweight='bold', color='#8B0000',
                ha='center', va='bottom',
                arrowprops=dict(arrowstyle='->', color='#CC6666', lw=0.8),
                zorder=10)

# ── Highlight k=50 result ──
ax.annotate(r'$\eta(50) = 50.7$' + '\n' + r'$E = 101\%$' + '\n' + r'$T_{\rm multi}=0.032$ s',
            xy=(50, 50.7), xytext=(35, 52),
            arrowprops=dict(arrowstyle='->', color='#8B0000', lw=1.5),
            fontsize=9, color='#8B0000', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#8B0000', alpha=0.9),
            zorder=10)

# ── Labels ──
ax.set_xlabel('Number of Sources', fontsize=12, fontweight='bold')
ax.set_ylabel('Speedup Ratio ' + r'$\eta$', fontsize=12, fontweight='bold')
ax.legend(loc='upper left', fontsize=9, framealpha=0.9, edgecolor='gray')

# Grid
ax.grid(True, alpha=0.3, color='gray', lw=0.5)
ax.set_xlim(-1, 58)
ax.set_ylim(-2, 60)
ax.set_xticks([0, 10, 20, 30, 40, 50])
ax.set_yticks([0, 10, 20, 30, 40, 50])

# ── Inset: Sh10 mesh info ──
ax.text(0.98, 0.08, 'Sh10 mesh\n8,800 vertices | 17,222 faces\n$T_{\\rm clear}$ + $T_{\\rm prop}$ = 0.032 s',
        transform=ax.transAxes, fontsize=8.5, color='#555555',
        verticalalignment='bottom', horizontalalignment='right',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#F5F5F5', edgecolor='#CCCCCC', alpha=0.9))

plt.tight_layout()
plt.savefig('e:/vtp_geodesic/论文截图/fig4_speedup.png', dpi=400,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('Fig 4 saved')
