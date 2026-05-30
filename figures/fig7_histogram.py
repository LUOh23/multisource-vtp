"""
图7 — 距离比值分布直方图 (彩色恢复版)
绿(<1%) / 黄(1%-2%) / 红(>2%)
"""
import numpy as np
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import vtp_geodesic
from pyproj import Transformer

plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False

print('Loading DEM for histogram...')
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
print(f'Done. ratio mean={np.mean(ratio_pct/100+1):.3f} max={np.max(ratio_pct/100+1):.3f}')

# ── Plot ──
fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))
bins = np.concatenate([np.arange(0, 5.05, 0.15), [max(ratio_pct.max(), 5.1)+0.01]])
counts, edges = np.histogram(ratio_pct, bins=bins)
centers = (edges[:-1]+edges[1:])/2

for i in range(len(counts)):
    c = centers[i]; w = (edges[i+1]-edges[i])*0.9
    if c < 1.0:
        ax.bar(c, counts[i], width=w, color='#2E8B57', alpha=0.80, edgecolor='white', lw=0.3, zorder=3)
    elif c < 2.0:
        ax.bar(c, counts[i], width=w, color='#DAA520', alpha=0.80, edgecolor='white', lw=0.3, zorder=3)
    else:
        ax.bar(c, counts[i], width=w, color='#B22222', alpha=0.80, edgecolor='white', lw=0.3, zorder=3)

ax.axvline(x=1.0, color='#2E8B57', ls='--', lw=1.5, alpha=0.7)
ax.axvline(x=2.0, color='#DAA520', ls='--', lw=1.5, alpha=0.7)
ax.axvline(x=5.0, color='#B22222', ls=':', lw=1.2, alpha=0.5)

# Zone labels
ym = ax.get_ylim()[1]*1.1
ax.text(0.5, ym*0.92, 'Excellent\n(< 1%)', ha='center', fontsize=9, color='#2E8B57', fontweight='bold')
ax.text(1.5, ym*0.92, 'Good\n(1%–2%)', ha='center', fontsize=9, color='#B8860B', fontweight='bold')
ax.text(3.5, ym*0.92, 'Poor\n(> 2%)', ha='center', fontsize=9, color='#B22222', fontweight='bold')

# Stats box
f1 = (ratio_pct<1).sum()/len(ratio_pct)*100
f2 = (ratio_pct<2).sum()/len(ratio_pct)*100
f5 = (ratio_pct>5).sum()/len(ratio_pct)*100
stats = (f'Mean: {np.mean(ratio_pct):.2f}%\nMax:  {np.max(ratio_pct):.2f}%\n'
         f'Std:  {np.std(ratio_pct):.2f}%\n<1%:  {f1:.1f}%\n<2%:  {f2:.1f}%\n>5%:  {f5:.2f}%')
ax.text(0.97, 0.95, stats, transform=ax.transAxes, fontsize=9, va='top', ha='right',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#FFF8DC', edgecolor='#CCAA66', alpha=0.92),
        fontfamily='monospace')

# Legend
legend_el = [
    Patch(facecolor='#2E8B57', alpha=0.8, label='[0, 1%)  Excellent'),
    Patch(facecolor='#DAA520', alpha=0.8, label='[1%, 2%)  Good'),
    Patch(facecolor='#B22222', alpha=0.8, label='[2%, ...)  Poor'),
]
ax.legend(handles=legend_el, loc='upper center', fontsize=8.5, ncol=3, framealpha=0.9, edgecolor='gray')

ax.set_xlabel(r'$(d_g / d_e - 1)$ / %', fontsize=12, fontweight='bold')
ax.set_ylabel('Vertex Count', fontsize=12, fontweight='bold')
ax.set_xlim(-0.1, max(ratio_pct.max(), 5.5)+0.3)
ax.grid(True, alpha=0.2, color='gray', lw=0.3, axis='y')

plt.tight_layout()
plt.savefig('output/fig6_histogram.png', dpi=600,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('Fig.7 (histogram) saved')
