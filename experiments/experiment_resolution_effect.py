"""
方案1: DEM分辨率对测地/欧氏距离偏差的影响
核心发现: 高分辨率下>2%偏差比例系统性高于粗分辨率
VTP作为核心工具贯穿始终
"""
import numpy as np
import rasterio
import vtp_geodesic
from pyproj import Transformer
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

print('='*60)
print('DEM分辨率对地形曲面距离效应的影响')
print('基于多源VTP精确测地算法的定量分析')
print('='*60)

# ── Load DEM and crop target region ──
dem = rasterio.open(r'E:\vtp_geodesic\高程原始数据\ASTGTM2_N29E116_dem.tif')
elev_full = dem.read(1).astype(float)
elev_full[elev_full < -1000] = np.nan

# Target: N29.3°~29.7°, E116.0°~116.5°
row_start, row_end = int((30.0-29.7)*3600), int((30.0-29.3)*3600)  # 1080:2519
col_start, col_end = int((116.0-116.0)*3600), int((116.5-116.0)*3600)  # 0:1800
elev_crop = elev_full[row_start:row_end, col_start:col_end].copy()
crop_rows, crop_cols = elev_crop.shape

# Build lon/lat grids
lons = np.linspace(116.0, 116.5, crop_cols)
lats = np.linspace(29.7, 29.3, crop_rows)
lon_grid, lat_grid = np.meshgrid(lons, lats)
transformer = Transformer.from_crs("EPSG:4326", "EPSG:32650", always_xy=True)

def build_mesh_and_compute(elev_data, lons_2d, lats_2d, step, label):
    """Build triangulated mesh at given step and run VTP."""
    # Subsample
    elev_sub = elev_data[::step, ::step]
    lons_sub = lons_2d[::step, ::step]
    lats_sub = lats_2d[::step, ::step]
    rows, cols = elev_sub.shape
    print(f'\n[{label}] step={step}, grid={rows}×{cols}={rows*cols} vertices')

    # UTM projection
    x_m, y_m = transformer.transform(lons_sub.ravel(), lats_sub.ravel())
    x_m = np.array(x_m).reshape(rows, cols)
    y_m = np.array(y_m).reshape(rows, cols)
    x_m -= x_m.mean(); y_m -= y_m.mean()

    points = np.column_stack([x_m.ravel(), y_m.ravel(), elev_sub.ravel()]).astype(np.float64)
    valid_mask = ~np.isnan(points[:, 2])

    # Delaunay-like triangulation
    faces_list = []
    for i in range(rows-1):
        for j in range(cols-1):
            a=i*cols+j; b=a+1; c=a+cols; d=c+1
            if all(valid_mask[[a,b,c]]): faces_list.append([a,b,c])
            if all(valid_mask[[b,d,c]]): faces_list.append([b,d,c])
    faces_vtp = np.array(faces_list, dtype=np.int32)
    n_verts = valid_mask.sum()
    print(f'  Valid vertices: {n_verts}, faces: {len(faces_vtp)}')
    print(f'  Resolution: ~{step*30:.0f}m, Area: ~{rows*step*0.03:.1f}km × {cols*step*0.03:.1f}km')

    # Select 5 source points (peak, 2 mid-slope, 2 valley)
    valid_indices = np.where(valid_mask)[0]
    elev_valid = points[valid_indices, 2]
    peak_idx = valid_indices[np.argmax(elev_valid)]
    elev_range = elev_valid.max() - elev_valid.min()
    mid50_idx = valid_indices[np.argmin(np.abs(elev_valid - elev_valid.min() - 0.50*elev_range))]
    mid75_idx = valid_indices[np.argmin(np.abs(elev_valid - elev_valid.min() - 0.75*elev_range))]
    valley_order = np.argsort(elev_valid)
    v1_idx = valid_indices[valley_order[0]]
    v2_idx = valid_indices[valley_order[min(10, len(valley_order)-1)]]
    sources = [peak_idx, mid50_idx, mid75_idx, v1_idx, v2_idx]

    print(f'  Sources: peak@{points[peak_idx,2]:.0f}m, mid@{points[mid50_idx,2]:.0f}m, valley@{points[v1_idx,2]:.0f}m')

    # VTP
    vtp = vtp_geodesic.PyVTPAlgorithmExact(points, faces_vtp)
    distances = vtp.compute_distances(sources)

    # Euclidean (min to any source)
    euc = np.full(len(points), np.inf)
    for s in sources:
        d = np.sqrt(np.sum((points - points[s])**2, axis=1))
        euc = np.minimum(euc, d)

    valid = np.isfinite(distances) & (distances > 1.0) & np.isfinite(euc) & (euc > 1.0)
    ratio = distances[valid] / (euc[valid] + 1e-10)
    ratio_pct = (ratio - 1) * 100

    stats = {
        'label': label,
        'step': step,
        'resolution_m': step * 30,
        'n_vertices': n_verts,
        'n_valid': valid.sum(),
        'mean_r': np.mean(ratio),
        'max_r': np.max(ratio),
        'std_r': np.std(ratio),
        'frac_gt1': (ratio_pct > 1).sum() / len(ratio_pct) * 100,
        'frac_gt2': (ratio_pct > 2).sum() / len(ratio_pct) * 100,
        'frac_gt5': (ratio_pct > 5).sum() / len(ratio_pct) * 100,
        'frac_gt10': (ratio_pct > 10).sum() / len(ratio_pct) * 100,
        'max_abs_m': (distances[valid] - euc[valid]).max(),
        'ratio_pct': ratio_pct,
    }
    print(f'  Mean d_g/d_e={stats["mean_r"]:.4f}, Max={stats["max_r"]:.4f}')
    print(f'  >1%: {stats["frac_gt1"]:.1f}%, >2%: {stats["frac_gt2"]:.1f}%, >5%: {stats["frac_gt5"]:.1f}%, >10%: {stats["frac_gt10"]:.1f}%')
    return stats

# ── Run 4 resolution gradients ──
# steps: 1=30m, 2=60m, 4=120m, 10=300m
results = []
for step, label in [(1, '30m'), (2, '60m'), (4, '120m'), (10, '300m')]:
    stats = build_mesh_and_compute(elev_crop, lon_grid, lat_grid, step, label)
    results.append(stats)

# ══════════════════════════════════════════════════════════
# Generate figures and summary table
# ══════════════════════════════════════════════════════════
print('\n' + '='*60)
print('对比分析')
print('='*60)

# Table
print(f"\n{'分辨率':>8} {'顶点数':>8} {'均值r':>8} {'最大r':>8} {'>1%':>8} {'>2%':>8} {'>5%':>8} {'>10%':>8}")
print('-'*70)
for r in results:
    print(f"{r['label']:>8} {r['n_vertices']:>8d} {r['mean_r']:>8.4f} {r['max_r']:>8.4f} "
          f"{r['frac_gt1']:>7.1f}% {r['frac_gt2']:>7.1f}% {r['frac_gt5']:>7.1f}% {r['frac_gt10']:>7.1f}%")

# ── Figure 1: Resolution vs bias ──
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

resolutions = [r['resolution_m'] for r in results]
frac_gt2 = [r['frac_gt2'] for r in results]
frac_gt5 = [r['frac_gt5'] for r in results]
mean_r = [(r['mean_r']-1)*100 for r in results]  # convert to %

# Left: fraction exceeding threshold
ax1.plot(resolutions, frac_gt2, 'o-', color='#C00000', lw=2.5, ms=10, label=r'$d_g/d_e > 1.02$')
ax1.plot(resolutions, frac_gt5, 's--', color='#D98736', lw=2.0, ms=10, label=r'$d_g/d_e > 1.05$')
ax1.set_xlabel('DEM Resolution / m', fontsize=12, fontweight='bold')
ax1.set_ylabel('Fraction of Vertices / %', fontsize=12, fontweight='bold')
ax1.set_title('(a) Fraction of vertices exceeding threshold', fontsize=13, fontweight='bold')
ax1.legend(fontsize=10, framealpha=0.9)
ax1.grid(True, alpha=0.3)
ax1.invert_xaxis()

# Annotations — moved away from lines
ax1.annotate(f'{frac_gt2[0]:.1f}%', (30, frac_gt2[0]), textcoords="offset points", xytext=(15,18),
            fontsize=11, fontweight='bold', color='#C00000')
ax1.annotate(f'{frac_gt2[-1]:.1f}%', (300, frac_gt2[-1]), textcoords="offset points", xytext=(-35,15),
            fontsize=11, fontweight='bold', color='#C00000')

# Right: mean deviation
ax2.fill_between(resolutions, 0, mean_r, alpha=0.3, color='#336699')
ax2.plot(resolutions, mean_r, 'o-', color='#003366', lw=2.5, ms=10)
ax2.set_xlabel('DEM Resolution / m', fontsize=12, fontweight='bold')
ax2.set_ylabel(r'Mean $(d_g/d_e - 1)$ / %', fontsize=12, fontweight='bold')
ax2.set_title('(b) Mean geodesic–Euclidean deviation', fontsize=13, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.invert_xaxis()

for i, r in enumerate(results):
    ax2.annotate(f'{mean_r[i]:.2f}%', (resolutions[i], mean_r[i]),
                textcoords="offset points", xytext=(0,16), fontsize=11,
                fontweight='bold', color='#003366', ha='center')
plt.tight_layout()
plt.savefig('e:/vtp_geodesic/论文截图/fig_resolution_effect.png', dpi=300,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('\nFigure saved: 论文截图/fig_resolution_effect.png')

# ── Figure 2: Histogram overlay (all resolutions) ──
fig, ax = plt.subplots(1, 1, figsize=(8, 5))
colors = ['#003366', '#336699', '#6699CC', '#99CCFF']
for i, r in enumerate(results):
    ax.hist(r['ratio_pct'], bins=80, range=(0, 15), alpha=0.5, color=colors[i],
            label=f"{r['label']} (n={r['n_valid']})", density=True)
ax.set_xlabel(r'$(d_g/d_e - 1)$ / %', fontsize=12, fontweight='bold')
ax.set_ylabel('Density', fontsize=12, fontweight='bold')
ax.set_title('Distribution of Geodesic/Euclidean Bias Across Resolutions', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.axvline(x=2, color='red', ls='--', lw=1, alpha=0.5)
ax.text(2.2, ax.get_ylim()[1]*0.9, '2% threshold', fontsize=9, color='red')
ax.set_xlim(0, 12)
plt.tight_layout()
plt.savefig('e:/vtp_geodesic/论文截图/fig_resolution_histogram.png', dpi=300,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('Figure saved: 论文截图/fig_resolution_histogram.png')

print('\n' + '='*60)
print('核心发现:')
print(f'  30m→300m, >2%偏差比例从 {frac_gt2[0]:.1f}% 降至 {frac_gt2[-1]:.1f}%')
print(f'  粗分辨率({results[-1]["label"]})系统性地低估了地形曲面效应的影响范围')
print(f'  建议: 在精度敏感的测绘应用中,应使用尽可能高的DEM分辨率')
print('='*60)
