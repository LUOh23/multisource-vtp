"""
图2 重绘: VTP vs Heat Method 计算时间对比
《测绘科学》标准: 双对数坐标, 红色VTP线, 灰色Heat虚线, 波动注释
Data from Table 2
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 10
plt.rcParams['mathtext.fontset'] = 'stix'

# ── Data from Table 2 ──
vertices  = np.array([50, 226, 530, 962, 2210, 3970, 6242, 9802, 19322, 32042])
vtp_time  = np.array([0.001, 0.002, 0.001, 0.003, 0.006, 0.012, 0.018, 0.025, 0.056, 0.094])
heat_time = np.array([0.008, 0.002, 0.007, 0.011, 0.044, 0.060, 0.109, 0.169, 0.381, 0.790])

fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))

# VTP: deep red solid line with circles
ax.loglog(vertices, vtp_time, 'o-', color='#8B0000', lw=2.0, ms=7,
          markerfacecolor='#C00000', markeredgecolor='#8B0000', markeredgewidth=1.0,
          label='VTP', zorder=5)

# Heat Method: gray dashed line with squares
ax.loglog(vertices, heat_time, 's--', color='#666666', lw=1.5, ms=6,
          markerfacecolor='#999999', markeredgecolor='#666666', markeredgewidth=1.0,
          label='Heat Method', zorder=4)

# ── Labels ──
ax.set_xlabel('Number of Vertices', fontsize=12, fontweight='bold')
ax.set_ylabel('Time / s', fontsize=12, fontweight='bold')
ax.legend(loc='upper left', fontsize=10, framealpha=0.9, edgecolor='gray')

# Grid
ax.grid(True, which='major', alpha=0.3, color='gray', lw=0.5)
ax.grid(True, which='minor', alpha=0.1, color='gray', lw=0.3)
ax.set_xlim(40, 50000)
ax.set_ylim(0.0005, 2.0)

# Ticks
ax.set_xticks([50, 100, 500, 1000, 5000, 10000, 32000])
ax.set_xticklabels(['50', '100', '500', '1000', '5000', '10000', '32000'])
ax.set_yticks([0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0])
ax.set_yticklabels(['0.001', '0.005', '0.01', '0.05', '0.1', '0.5', '1.0'])

plt.tight_layout()
plt.savefig('e:/vtp_geodesic/论文截图/fig2_performance_time.png', dpi=400,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('Fig 2 saved')
