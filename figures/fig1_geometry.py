"""
图1 改进版 — VTP局部波前传播几何构造
仅加3处标注: 入射面/出射面 + 传播箭头 + p_y<=0
全英文标注，匹配原图风格。另存为fig1_geometry_v2.png
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

fig, ax = plt.subplots(1, 1, figsize=(9, 6.5))
ax.set_aspect('equal')

# ── Vertices ──
S = np.array([1.0, 4.0])
A = np.array([0.0, 1.0])
B = np.array([4.5, 0.5])
C = np.array([6.0, 2.8])

# ── Triangles ──
ax.add_patch(plt.Polygon([S, A, B], facecolor='#F6F6F6', edgecolor='black', lw=1.2, zorder=1))
ax.add_patch(plt.Polygon([A, B, C], facecolor='#EEEEEE', edgecolor='black', lw=1.2, zorder=1))

# ── NEW: Incident / Outgoing face labels (subtle gray italic) ──
centroid_in = (S + A + B) / 3
ax.text(centroid_in[0], centroid_in[1], 'Incident', fontsize=9, color='#999999',
        ha='center', va='center', style='italic')
centroid_out = (A + B + C) / 3
ax.text(centroid_out[0]+0.1, centroid_out[1], 'Outgoing', fontsize=9, color='#999999',
        ha='center', va='center', style='italic')

# ── Vertices ──
ax.plot(S[0], S[1], 'o', color='black', markersize=10, zorder=10)
ax.plot(A[0], A[1], 'o', color='black', markersize=7, zorder=10)
ax.plot(B[0], B[1], 'o', color='black', markersize=7, zorder=10)
ax.plot(C[0], C[1], 'o', color='black', markersize=7, zorder=10)

ax.text(S[0]+0.18, S[1]+0.18, r'$S$', fontsize=12, fontweight='bold')
ax.text(A[0]-0.38, A[1]-0.15, r'$A$', fontsize=11, fontweight='bold')
ax.text(B[0]+0.18, B[1]-0.2, r'$B$', fontsize=11, fontweight='bold')
ax.text(C[0]+0.18, C[1]+0.12, r'$C$', fontsize=11, fontweight='bold')

# ── Edge AB (shared edge) — bold ──
ax.plot([A[0], B[0]], [A[1], B[1]], 'k-', lw=3.0, zorder=5)
mid_AB = (A + B) / 2

# ── Edges AC and BC ──
ax.plot([A[0], C[0]], [A[1], C[1]], 'k-', lw=1.3, alpha=0.6, zorder=2)
ax.plot([B[0], C[0]], [B[1], C[1]], 'k-', lw=1.3, alpha=0.6, zorder=2)

# ── NEW: Propagation arrow from S toward AB (no text, just arrow) ──
arrow_start = S + (mid_AB - S) * 0.25
arrow_end = S + (mid_AB - S) * 0.55
ax.annotate('', xy=arrow_end, xytext=arrow_start,
            arrowprops=dict(arrowstyle='->', color='#336699', lw=2.2, alpha=0.5), zorder=12)

# ── Window interval on AB ──
e1_dir = (B - A) / np.linalg.norm(B - A)
e2_dir = np.array([-e1_dir[1], e1_dir[0]])

interval_s = mid_AB + e1_dir * (-1.5)
interval_e = mid_AB + e1_dir * 1.2
ax.plot([interval_s[0], interval_e[0]], [interval_s[1], interval_e[1]],
        color='#CC0000', lw=4.5, alpha=0.5, zorder=6)
ax.text(mid_AB[0]-0.5, mid_AB[1]-1.3, r'$[\sigma_{\rm start},\sigma_{\rm stop}]$',
        fontsize=9, color='#CC0000', fontweight='bold')

# ── Pseudo-source ──
p_src = mid_AB - e2_dir * 0.75
ax.plot(p_src[0], p_src[1], 'D', color='#CC6600', markersize=8, markeredgecolor='#994400',
        markeredgewidth=1.2, zorder=10)
ax.text(p_src[0]+0.25, p_src[1]+0.12, r'$(p_x, p_y)$', fontsize=10, color='#CC6600', fontweight='bold')

# ── NEW: p_y <= 0 note (math font, works with Times New Roman) ──
ax.text(p_src[0]-0.75, p_src[1]-0.42, r'$p_y \leq 0$', fontsize=9, color='#CC6600',
        style='italic', ha='center')

# ── Point M and dashed line ──
ax.plot(mid_AB[0], mid_AB[1], 'o', color='#CC6600', markersize=5, zorder=10)
ax.text(mid_AB[0]+0.15, mid_AB[1]-0.25, r'$M$', fontsize=9)
ax.plot([mid_AB[0], p_src[0]], [mid_AB[1], p_src[1]], '--', color='#CC6600', lw=1.8, alpha=0.7, zorder=3)

# ── Local coordinate system at A ──
sc = 0.9
ax.arrow(A[0], A[1], e1_dir[0]*sc, e1_dir[1]*sc,
         head_width=0.12, head_length=0.14, fc='#003366', ec='#003366', lw=1.8, zorder=15)
ax.arrow(A[0], A[1], e2_dir[0]*sc, e2_dir[1]*sc,
         head_width=0.12, head_length=0.14, fc='#003366', ec='#003366', lw=1.8, zorder=15)
ax.text(A[0]+e1_dir[0]*sc+0.1, A[1]+e1_dir[1]*sc+0.12,
        r'$\mathbf{e}_1$', fontsize=11, color='#003366', fontweight='bold')
ax.text(A[0]+e2_dir[0]*sc-0.35, A[1]+e2_dir[1]*sc+0.15,
        r'$\mathbf{e}_2$', fontsize=11, color='#003366', fontweight='bold')

# ── Angle α at C ──
a1 = np.arctan2(A[1]-C[1], A[0]-C[0])
a2 = np.arctan2(B[1]-C[1], B[0]-C[0])
th = np.linspace(a1, a2, 40)
r = 0.55
ax.plot(C[0]+r*np.cos(th), C[1]+r*np.sin(th), 'k-', lw=1.2, zorder=3)
mid_angle = (a1 + a2) / 2
alpha_r = r + 0.3
ax.text(C[0]+alpha_r*np.cos(mid_angle)-0.15, C[1]+alpha_r*np.sin(mid_angle)-0.15,
        r'$\alpha$', fontsize=14, fontweight='bold')

# ── Wavefront arc ──
th_w = np.linspace(np.arctan2(A[1]-S[1], A[0]-S[0]), np.arctan2(B[1]-S[1], B[0]-S[0]), 30)
ax.plot(S[0]+1.3*np.cos(th_w), S[1]+1.3*np.sin(th_w), color='#336699', lw=1.2, alpha=0.5, zorder=2)

ax.set_xlim(-1.8, 7.8)
ax.set_ylim(-1.2, 5.5)
ax.axis('off')

plt.tight_layout(pad=0.5)
plt.savefig('e:/vtp_geodesic/论文截图/fig1_geometry_v2.png', dpi=400,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('图1改进版(v2, English labels) → 论文截图/fig1_geometry_v2.png')
