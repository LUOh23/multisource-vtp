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
          label='VTP (exact)', zorder=5)

# Heat Method: gray dashed line with squares
ax.loglog(vertices, heat_time, 's--', color='#666666', lw=1.5, ms=6,
          markerfacecolor='#999999', markeredgecolor='#666666', markeredgewidth=1.0,
          label='Heat Method (approx.)', zorder=4)

# ── Annotate the spike at 530 vertices ──
# The VTP time drops from 0.002→0.001→0.003 (anomalous dip at 530)
ax.annotate('Measurement artifact:\n timer resolution limit\n at small mesh size',
            xy=(530, 0.001), xytext=(150, 0.00025),
            arrowprops=dict(arrowstyle='->', color='#555555', lw=1.2, connectionstyle='arc3,rad=0.3'),
            fontsize=8.5, color='#555555', style='italic',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF8DC', edgecolor='#CCCCCC', alpha=0.9),
            zorder=10)

# ── Annotate the key result at 32k ──
ax.annotate('VTP 8.4× faster\nat 32k vertices',
            xy=(32042, 0.094), xytext=(8000, 0.25),
            arrowprops=dict(arrowstyle='->', color='#8B0000', lw=1.5),
            fontsize=9, color='#8B0000', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#8B0000', alpha=0.9),
            zorder=10)

# ── Labels ──
ax.set_xlabel('Number of Vertices |V|', fontsize=12, fontweight='bold')
ax.set_ylabel('Computation Time / s', fontsize=12, fontweight='bold')
ax.legend(loc='upper left', fontsize=10, framealpha=0.9, edgecolor='gray')

# Grid
ax.grid(True, which='major', alpha=0.3, color='gray', lw=0.5)
ax.grid(True, which='minor', alpha=0.1, color='gray', lw=0.3)
ax.set_xlim(40, 50000)
ax.set_ylim(0.0005, 2.0)

# Ticks
ax.set_xticks([50, 100, 500, 1000, 5000, 10000, 32000])
ax.set_xticklabels(['50', '100', '500', '1k', '5k', '10k', '32k'])
ax.set_yticks([0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0])
ax.set_yticklabels(['0.001', '0.005', '0.01', '0.05', '0.1', '0.5', '1.0'])

plt.tight_layout()
plt.savefig('output/fig2_performance_time.png', dpi=400,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('Fig 2 saved')
