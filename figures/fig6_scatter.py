"""
图6 — 测地vs欧氏距离散点 (彩色恢复版)
半透明红点 + 蓝三角离群点
"""
import numpy as np
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import vtp_geodesic
from pyproj import Transformer

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False

print('Loading DEM for scatter...')
dem = rasterio.open(r'E:\vtp_geodesic\高程原始数据\ASTGTM2_N29E116_dem.tif')
elev_full = dem.read(1).astype(float)
elev_full[elev_full < -1000] = np.nan
step = 10
elev = elev_full[::step, ::step]
rows, cols = elev.shape
lons = np.linspace(dem.bounds.left, dem.bounds.right, dem.width)[::step][:cols]
lats = np.linspace(dem.bounds.top, dem.bounds.bottom, dem.height)[::step][:rows]
lon_grid, lat_grid = np.meshgrid(lons, lats)
transformer = Transformer.from_crs("EPSG:4326", "EPSG:32650", always_xy=True)
x_m, y_m = transformer.transform(lon_grid.ravel(), lat_grid.ravel())
x_m = np.array(x_m).reshape(rows, cols); y_m = np.array(y_m).reshape(rows, cols)
x_m -= x_m.mean(); y_m -= y_m.mean()
points = np.column_stack([x_m.ravel(), y_m.ravel(), elev.ravel()]).astype(np.float64)
valid_mask = ~np.isnan(points[:, 2])
faces_list = []
for i in range(rows-1):
    for j in range(cols-1):
        a=i*cols+j; b=a+1; c=a+cols; d=c+1
        if all(valid_mask[[a,b,c]]): faces_list.append([a,b,c])
        if all(valid_mask[[b,d,c]]): faces_list.append([b,d,c])
faces_vtp = np.array(faces_list, dtype=np.int32)
peak_idx = np.where(valid_mask)[0][np.nanargmax(points[np.where(valid_mask)[0], 2])]
vtp = vtp_geodesic.PyVTPAlgorithmExact(points, faces_vtp)
d_peak = vtp.compute_distances([peak_idx])
euc = np.sqrt(np.sum((points - points[peak_idx])**2, axis=1))
valid = np.isfinite(d_peak) & (d_peak > 1.0)
ratio_pct = (d_peak[valid]/(euc[valid]+1e-10) - 1)*100
d_v = d_peak[valid]; e_v = euc[valid]
print(f'Done. ratio mean={np.mean(d_v/(e_v+1e-10)):.3f}')

# ── Plot ──
fig, ax = plt.subplots(1, 1, figsize=(7.5, 6.5))
samp = np.random.choice(len(d_v), min(3000, len(d_v)), replace=False)
out = ratio_pct[samp] > 5
norm = ~out

ax.scatter(e_v[samp][norm], d_v[samp][norm], c='#C03030', s=6, alpha=0.45,
           edgecolors='none', label='Normal (< 5%)', zorder=3)
if out.sum() > 0:
    ax.scatter(e_v[samp][out], d_v[samp][out], c='#0055CC', s=55, alpha=0.9,
               marker='^', edgecolors='#003388', linewidths=0.8,
               label=f'Outlier > 5% (n={out.sum()})', zorder=5)

lim = max(e_v.max(), d_v.max())
ax.plot([0, lim], [0, lim], '--', color='#333333', lw=2.0, label=r'$d_g = d_e$', zorder=4)
ax.set_xlabel('Euclidean Distance / m', fontsize=12, fontweight='bold')
ax.set_ylabel('Geodesic Distance / m', fontsize=12, fontweight='bold')
ax.legend(loc='upper left', fontsize=10, framealpha=0.9, edgecolor='gray')
ax.grid(True, alpha=0.2, color='gray', lw=0.4)
ax.set_xlim(0, lim*1.02); ax.set_ylim(0, lim*1.02)

stats = 'Vertices: 130,321\nMean r: 1.008\nMax r: 1.080'
ax.text(0.97, 0.06, stats, transform=ax.transAxes, fontsize=9, va='bottom', ha='right',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#AAAAAA', alpha=0.9))

plt.tight_layout()
plt.savefig('output/fig5_scatter.png', dpi=600,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('Fig.6 (scatter) saved')
