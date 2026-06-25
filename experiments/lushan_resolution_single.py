"""
庐山DEM — 单源分辨率效应 (仅山顶源点)
4级分辨率: 12.5m / 25m / 50m / 100m
"""
import numpy as np
import rasterio
import vtp_geodesic
import time

print('=' * 60)
print('庐山 DEM — 单源分辨率效应')
print('=' * 60)

dem = rasterio.open(r'E:\vtp_geodesic\data\dem_clip.tif')
elev_full = dem.read(1).astype(float)
rows_full, cols_full = elev_full.shape
res_x, res_y = dem.res
x_coords_full = np.linspace(dem.bounds.left + res_x/2, dem.bounds.right - res_x/2, cols_full)
y_coords_full = np.linspace(dem.bounds.top - res_y/2, dem.bounds.bottom + res_y/2, rows_full)
dem.close()

def run_single_source(elev, step, label):
    elev_sub = elev[::step, ::step]
    x_sub = x_coords_full[::step]
    y_sub = y_coords_full[::step]
    xg, yg = np.meshgrid(x_sub, y_sub)
    xg -= xg.mean(); yg -= yg.mean()
    rows, cols = elev_sub.shape
    n_all = rows * cols
    points = np.column_stack([xg.ravel(), yg.ravel(), elev_sub.ravel()]).astype(np.float64)

    faces_list = []
    for i in range(rows-1):
        for j in range(cols-1):
            a=i*cols+j; b=a+1; c=a+cols; d=c+1
            faces_list.append([a,b,c]); faces_list.append([b,d,c])
    faces_vtp = np.array(faces_list, dtype=np.int32)

    peak_idx = int(np.argmax(elev_sub.ravel()))
    print(f'\n[{label}] {n_all:,} verts, {len(faces_vtp):,} faces, peak={elev_sub.ravel()[peak_idx]:.0f}m')

    t0 = time.perf_counter()
    vtp = vtp_geodesic.PyVTPAlgorithmExact(points, faces_vtp)
    d_geo = vtp.compute_distances([peak_idx])
    t_vtp = time.perf_counter() - t0

    d_euc = np.sqrt(np.sum((points - points[peak_idx])**2, axis=1))
    valid = np.isfinite(d_geo) & (d_geo > 1.0) & (d_euc > 1.0)
    ratio = d_geo[valid] / (d_euc[valid] + 1e-10)
    ratio_pct = (ratio - 1) * 100

    stats = {
        'label': label, 'step': step,
        'resolution_m': step * res_x,
        'n_vertices': n_all, 'n_faces': len(faces_vtp),
        'n_valid': int(valid.sum()),
        'mean_r': float(np.mean(ratio)),
        'max_r': float(np.max(ratio)),
        'std_r': float(np.std(ratio)),
        'mean_pct': float(np.mean(ratio_pct)),
        'frac_gt1': float(100*(ratio_pct>1).sum()/len(ratio_pct)),
        'frac_gt2': float(100*(ratio_pct>2).sum()/len(ratio_pct)),
        'frac_gt5': float(100*(ratio_pct>5).sum()/len(ratio_pct)),
        'frac_gt10': float(100*(ratio_pct>10).sum()/len(ratio_pct)),
        't_vtp': t_vtp,
    }
    print(f'  mean={stats["mean_r"]:.4f}, max={stats["max_r"]:.4f}, >2%={stats["frac_gt2"]:.1f}%, '
          f'>5%={stats["frac_gt5"]:.1f}%, >10%={stats["frac_gt10"]:.1f}%, VTP={t_vtp:.1f}s')
    return stats

results = []
for step, label in [(1, '12.5m'), (2, '25m'), (4, '50m'), (8, '100m')]:
    results.append(run_single_source(elev_full, step, label))

print(f'\n{"="*80}')
print(f'Table 2 单源 — DEM分辨率对测地/欧氏距离偏差的影响（庐山，山顶单源）')
print(f'{"="*80}')
print(f'{"分辨率":>8} {"顶点数":>10} {"均值dg/de":>10} {"最大dg/de":>10} {">1%":>8} {">2%":>8} {">5%":>8} {">10%":>8} {"VTP/s":>8}')
print('-'*82)
for r in results:
    print(f'{r["label"]:>8} {r["n_vertices"]:>10,} {r["mean_r"]:>10.4f} {r["max_r"]:>10.4f} '
          f'{r["frac_gt1"]:>7.1f}% {r["frac_gt2"]:>7.1f}% {r["frac_gt5"]:>7.1f}% {r["frac_gt10"]:>7.1f}% {r["t_vtp"]:>7.1f}')

print(f'\n关键对比:')
for r in results:
    print(f'  {r["label"]}: >2%={r["frac_gt2"]:.1f}%, >5%={r["frac_gt5"]:.1f}%, >10%={r["frac_gt10"]:.1f}%')
