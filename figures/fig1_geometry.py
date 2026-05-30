"""
新Fig.1 — VTP局部波前传播几何构造 (极简版)
去掉左上角文字框，仅保留图面核心几何标注
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.sans-serif'] = ['SimSun', 'SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False

fig, ax = plt.subplots(1, 1, figsize=(9, 6.5))
ax.set_aspect('equal')

# ── Vertices ──
S = np.array([1.0, 4.0])  # source
A = np.array([0.0, 1.0])
B = np.array([4.5, 0.5])
C = np.array([6.0, 2.8])

# ── Triangles ──
ax.add_patch(plt.Polygon([S, A, B], facecolor='#F6F6F6', edgecolor='black', lw=1.2, zorder=1))
ax.add_patch(plt.Polygon([A, B, C], facecolor='#EEEEEE', edgecolor='black', lw=1.2, zorder=1))

ax.text((S[0]+A[0]+B[0])/3, (S[1]+A[1]+B[1])/3, r'$\triangle_1$', fontsize=9, ha='center', color='#888888')
ax.text((A[0]+B[0]+C[0])/3, (A[1]+B[1]+C[1])/3, r'$\triangle_2$', fontsize=9, ha='center', color='#888888')

# ── Vertices ──
ax.plot(S[0], S[1], 'o', color='black', markersize=10, zorder=10)
ax.plot(A[0], A[1], 'o', color='black', markersize=7, zorder=10)
ax.plot(B[0], B[1], 'o', color='black', markersize=7, zorder=10)
ax.plot(C[0], C[1], 'o', color='black', markersize=7, zorder=10)

ax.text(S[0]+0.18, S[1]+0.18, r'$S$', fontsize=12, fontweight='bold')
ax.text(A[0]-0.38, A[1]-0.15, r'$A$', fontsize=11, fontweight='bold')
ax.text(B[0]+0.18, B[1]-0.2, r'$B$', fontsize=11, fontweight='bold')
ax.text(C[0]+0.18, C[1]+0.12, r'$C$', fontsize=11, fontweight='bold')

# ── Edge AB (bottom edge) — bold ──
ax.plot([A[0], B[0]], [A[1], B[1]], 'k-', lw=3.0, zorder=5)
mid_AB = (A + B) / 2
ax.text(mid_AB[0]-0.6, mid_AB[1]+0.3, r'$e_{\rm bottom}$', fontsize=9, fontweight='bold')

# ── Edges AC and BC ──
ax.plot([A[0], C[0]], [A[1], C[1]], 'k-', lw=1.3, alpha=0.6, zorder=2)
ax.plot([B[0], C[0]], [B[1], C[1]], 'k-', lw=1.3, alpha=0.6, zorder=2)

# ── Window interval on AB ──
e1_dir = (B - A) / np.linalg.norm(B - A)
e2_dir = np.array([-e1_dir[1], e1_dir[0]])

interval_s = mid_AB + e1_dir * (-1.5)
interval_e = mid_AB + e1_dir * 1.2
ax.plot([interval_s[0], interval_e[0]], [interval_s[1], interval_e[1]],
        color='#CC0000', lw=4.5, alpha=0.5, zorder=6)
ax.text(mid_AB[0]-0.5, mid_AB[1]-0.7, r'$[\sigma_{\rm start},\sigma_{\rm stop}]$',
        fontsize=9, color='#CC0000', fontweight='bold')

# ── Pseudo-source (p_x, p_y) ──
p_src = mid_AB - e2_dir * 1.0
ax.plot(p_src[0], p_src[1], 'D', color='#CC6600', markersize=9, zorder=10)
ax.text(p_src[0]+0.2, p_src[1]-0.3, r'$(p_x, p_y)$', fontsize=10, color='#CC6600', fontweight='bold')
ax.text(p_src[0]+0.2, p_src[1]-0.65, r'$p_y \leq 0$', fontsize=8, color='#CC6600', style='italic')

# ── Perpendicular bisector ──
perp_pt = mid_AB + e2_dir * 2.0
ax.plot([mid_AB[0], perp_pt[0]], [mid_AB[1], perp_pt[1]], 'k--', lw=0.8, alpha=0.5, zorder=3)
ax.plot(mid_AB[0], mid_AB[1], 'o', color='black', markersize=4, zorder=10)
ax.text(mid_AB[0]+0.15, mid_AB[1]-0.2, r'$M$', fontsize=9)

# ── Local coordinate system at A ──
sc = 0.9
ax.arrow(A[0], A[1], e1_dir[0]*sc, e1_dir[1]*sc,
         head_width=0.12, head_length=0.14, fc='#003366', ec='#003366', lw=1.8, zorder=15)
ax.arrow(A[0], A[1], e2_dir[0]*sc, e2_dir[1]*sc,
         head_width=0.12, head_length=0.14, fc='#003366', ec='#003366', lw=1.8, zorder=15)
ax.text(A[0]+e1_dir[0]*sc+0.1, A[1]+e1_dir[1]*sc+0.12,
        r'$\mathbf{e}_1$', fontsize=11, color='#003366', fontweight='bold')
ax.text(A[0]+e2_dir[0]*sc+0.1, A[1]+e2_dir[1]*sc-0.15,
        r'$\mathbf{e}_2$', fontsize=11, color='#003366', fontweight='bold')

# ── Angle α at C ──
a1 = np.arctan2(A[1]-C[1], A[0]-C[0])
a2 = np.arctan2(B[1]-C[1], B[0]-C[0])
th = np.linspace(a1, a2, 40)
r = 0.55
ax.plot(C[0]+r*np.cos(th), C[1]+r*np.sin(th), 'k-', lw=1.2, zorder=3)
ax.text(C[0]-0.35, C[1]-0.5, r'$\alpha$', fontsize=12, fontweight='bold')

# ── Wavefront arc ──
th_w = np.linspace(np.arctan2(A[1]-S[1], A[0]-S[0]), np.arctan2(B[1]-S[1], B[0]-S[0]), 30)
ax.plot(S[0]+1.3*np.cos(th_w), S[1]+1.3*np.sin(th_w), color='#336699', lw=1.2, alpha=0.5, zorder=2)
ax.text(S[0]+1.0, S[1]-0.4, 'wavefront', fontsize=8, color='#336699', style='italic')

# ── Clean legend at bottom ──
legend_items = [
    (r'$S$: source vertex  |  $\triangle_1$: incident face  |  $\triangle_2$: target face', 'black'),
    (r'$e_{\rm bottom}=AB$: shared edge  |  $(p_x,p_y)$: pseudo-source ($p_y\leq 0$, Property 1)', '#CC6600'),
    (r'$[\sigma_{\rm start},\sigma_{\rm stop}]$: window interval on $AB$', '#CC0000'),
    (r'$\mathbf{e}_1,\mathbf{e}_2$: local frame at $A$  |  $\alpha$: interior angle at $C$ (Eq. 5 rotation $R(\alpha)$)', '#003366'),
]
for i, (text, color) in enumerate(legend_items):
    ax.text(0.5, -0.12 - i*0.08, text, transform=ax.transAxes, fontsize=8.5,
            color=color, ha='center', va='top')

ax.set_xlim(-1.8, 7.8)
ax.set_ylim(-1.2, 5.5)
ax.axis('off')

plt.tight_layout(pad=0.5)
plt.savefig('output/fig1_geometry.png', dpi=600,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('New Fig.1 (geometry) saved')
