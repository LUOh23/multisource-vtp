"""
图: 多源波前互抑效应 — k vs 窗口创建数/互抑率
Times New Roman, clean academic style
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 10
plt.rcParams['mathtext.fontset'] = 'stix'

# Data from experiment
k_vals = np.array([2, 5, 10, 20, 50])
multi_created = np.array([3356, 3263, 3442, 3749, 4242])
sum_single = np.array([6857, 17388, 34912, 70024, 174797])
suppression = np.array([51.1, 81.2, 90.1, 94.6, 97.6])
speedup = np.array([2.4, 5.7, 14.8, 30.0, 69.4])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

# === Left panel ===
ax1.plot(k_vals, multi_created, 'o-', color='#003366', lw=2.5, ms=10,
         label='Multi-source (single call)')
ax1.plot(k_vals, sum_single, 's--', color='#C00000', lw=2.0, ms=10,
         label='Sum of k single-source calls')

ax1.annotate(f'{multi_created[-1]:,}',
             (k_vals[-1], multi_created[-1]),
             textcoords="offset points", xytext=(-10, 10),
             fontsize=10, fontweight='bold', color='#003366')
ax1.annotate(f'{sum_single[-1]:,}',
             (k_vals[-1], sum_single[-1]),
             textcoords="offset points", xytext=(-65, -10),
             fontsize=10, fontweight='bold', color='#C00000')

ax1.set_xlabel('Number of sources', fontsize=12, fontweight='bold')
ax1.set_ylabel('Total windows', fontsize=12, fontweight='bold')
ax1.set_title('(a) Window Creation Count', fontsize=13, fontweight='bold')
ax1.legend(fontsize=10, framealpha=0.9, loc='upper left')
ax1.grid(True, alpha=0.3)
ax1.set_xscale('log')
ax1.set_xticks(k_vals)
ax1.set_xticklabels([str(k) for k in k_vals])

# === Right panel ===
color_supp = '#003366'
color_spd = '#C00000'

ax2b = ax2.twinx()

line1, = ax2.plot(k_vals, suppression, 'o-', color=color_supp, lw=2.5, ms=10,
                  label='Suppression rate S (%)')
line2, = ax2b.plot(k_vals, speedup, 's--', color=color_spd, lw=2.0, ms=10,
                   label=r'Speedup ratio $\eta$')

ax2b.plot(k_vals, k_vals, ':', color='gray', lw=1.0, alpha=0.6)

# Annotations — each point has its own offset to avoid overlap
supp_offsets = [(20, 36), (0, 14), (0, 14), (0, 14), (0, 14)]
spd_offsets  = [(10, 10), (0, -14), (10, -12), (15, -10), (-8, -30)]
for i in range(len(k_vals)):
    ax2.annotate(f'{suppression[i]:.1f}%',
                 (k_vals[i], suppression[i]),
                 textcoords="offset points", xytext=supp_offsets[i],
                 fontsize=9, fontweight='bold', color=color_supp, ha='center')
    ax2b.annotate(f'{speedup[i]:.1f}x',
                  (k_vals[i], speedup[i]),
                  textcoords="offset points", xytext=spd_offsets[i],
                  fontsize=9, fontweight='bold', color=color_spd)

ax2.set_xlabel('Number of sources', fontsize=12, fontweight='bold')
ax2.set_ylabel('Suppression rate S (%)', fontsize=12, fontweight='bold', color=color_supp)
ax2b.set_ylabel(r'Speedup ratio $\eta$', fontsize=12, fontweight='bold', color=color_spd)
ax2.set_title('(b) Suppression Rate & Speedup', fontsize=13, fontweight='bold')
ax2.set_ylim(0, 108)
ax2b.set_ylim(0, 82)
ax2.grid(True, alpha=0.3)
ax2.set_xscale('log')
ax2.set_xticks(k_vals)
ax2.set_xticklabels([str(k) for k in k_vals])

# Legend moved to upper left
lines = [line1, line2]
labels = [r'Suppression rate $S$ (%)', r'Speedup ratio $\eta$']
ax2.legend(lines, labels, fontsize=10, framealpha=0.9, loc='upper left')

ax2.tick_params(axis='y', labelcolor=color_supp)
ax2b.tick_params(axis='y', labelcolor=color_spd)

plt.tight_layout()
plt.savefig('e:/vtp_geodesic/论文截图/fig_suppression_effect.png', dpi=300,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('Figure saved')
