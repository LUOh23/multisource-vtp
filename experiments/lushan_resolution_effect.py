"""
庐山12.5m DEM — DEM分辨率对测地/欧氏距离偏差的影响
替代原Table 2 (Poyang 30m/60m/120m/300m)
4级分辨率: 12.5m / 25m / 50m / 100m
5源点多源VTP
"""
import numpy as np
import rasterio
import vtp_geodesic
import time

print('=' * 60)
print('庐山 DEM — 分辨率效应实验 (Table 2)')
print('=' * 60)

# ── 1. 加载DEM ──
dem = rasterio.open(r'E:\vtp_geodesic\data\dem_clip.tif')
elev_full = dem.read(1).astype(float)
rows_full, cols_full = elev_full.shape
res_x, res_y = dem.res
x_coords_full = np.linspace(dem.bounds.left + res_x/2, dem.bounds.right - res_x/2, cols_full)
y_coords_full = np.linspace(dem.bounds.top - res_y/2, dem.bounds.bottom + res_y/2, rows_full)
dem.close()

def build_and_compute(elev, step, label):
    """在给定步长下建网并运行5源VTP"""
    elev_sub = elev[::step, ::step]
    x_sub = x_coords_full[::step]
    y_sub = y_coords_full[::step]
    xg, yg = np.meshgrid(x_sub, y_sub)
    xg -= xg.mean(); yg -= yg.mean()
    rows, cols = elev_sub.shape
    n_all = rows * cols

    points = np.column_stack([xg.ravel(), yg.ravel(), elev_sub.ravel()]).astype(np.float64)

    # 三角剖分
    faces_list = []
    for i in range(rows - 1):
        for j in range(cols - 1):
            a = i * cols + j
            b = a + 1
            c = a + cols
            d = c + 1
            faces_list.append([a, b, c])
            faces_list.append([b, d, c])
    faces_vtp = np.array(faces_list, dtype=np.int32)

    print(f'\n[{label}] step={step}, grid={rows}x{cols}={n_all:,} verts, {len(faces_vtp):,} faces')

    # 选5个源点: peak + 2 mid-slope + 2 valley
    elev_flat = elev_sub.ravel()
    valid_all = np.ones(n_all, dtype=bool)

    # Peak
    peak_idx = int(np.argmax(elev_flat))

    # Valley (lowest 2 distinct locations)
    # Pick the lowest point, and another low point far from it
    sorted_idx = np.argsort(elev_flat)
    v1_idx = int(sorted_idx[0])
    # Find another low point at least 2km away from v1
    v1_pt = points[v1_idx]
    for idx in sorted_idx[1:]:
        dist = np.sqrt(np.sum((points[idx] - v1_pt)**2))
        if dist > 2000:
            v2_idx = int(idx)
            break
    else:
        v2_idx = int(sorted_idx[min(10, len(sorted_idx)-1)])

    # Mid-slope: P50 and P75 elevation
    e_range = elev_flat.max() - elev_flat.min()
    mid50_idx = int(np.argmin(np.abs(elev_flat - elev_flat.min() - 0.50*e_range)))
    mid75_idx = int(np.argmin(np.abs(elev_flat - elev_flat.min() - 0.75*e_range)))

    sources = [peak_idx, mid50_idx, mid75_idx, v1_idx, v2_idx]
    src_elevs = [elev_flat[s] for s in sources]
    print(f'  源点高程: peak={src_elevs[0]:.0f}, mid50={src_elevs[1]:.0f}, mid75={src_elevs[2]:.0f}, v1={src_elevs[3]:.0f}, v2={src_elevs[4]:.0f}')

    # VTP multi-source
    t0 = time.perf_counter()
    vtp = vtp_geodesic.PyVTPAlgorithmExact(points, faces_vtp)
    d_geo = vtp.compute_distances(sources)
    t_vtp = time.perf_counter() - t0
    print(f'  VTP: {t_vtp:.1f}s')

    # Euclidean to nearest source
    d_euc = np.full(len(points), np.inf)
    for s in sources:
        d = np.sqrt(np.sum((points - points[s])**2, axis=1))
        d_euc = np.minimum(d_euc, d)

    valid = np.isfinite(d_geo) & (d_geo > 1.0) & (d_euc > 1.0)
    ratio = d_geo[valid] / (d_euc[valid] + 1e-10)
    ratio_pct = (ratio - 1) * 100

    stats = {
        'label': label,
        'step': step,
        'resolution_m': step * res_x,
        'n_vertices': n_all,
        'n_faces': len(faces_vtp),
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

    print(f'  mean dg/de={stats["mean_r"]:.4f}, max={stats["max_r"]:.4f}')
    print(f'  >1%: {stats["frac_gt1"]:.1f}%, >2%: {stats["frac_gt2"]:.1f}%, >5%: {stats["frac_gt5"]:.1f}%, >10%: {stats["frac_gt10"]:.1f}%')

    return stats, points, faces_vtp, d_geo, d_euc, valid, ratio_pct

# ── 2. 跑四级分辨率 ──
results = []
all_data = []
for step, label in [(1, '12.5m'), (2, '25m'), (4, '50m'), (8, '100m')]:
    stats, pts, faces, dg, de, v, rp = build_and_compute(elev_full, step, label)
    results.append(stats)
    all_data.append((pts, faces, dg, de, v, rp))

# ── 3. 输出表格 ──
print(f'\n{"="*80}')
print('Table 2 — DEM分辨率对测地/欧氏距离偏差的影响（庐山，5源点）')
print(f'{"="*80}')
print(f'{"分辨率":>8} {"顶点数":>10} {"均值dg/de":>10} {"最大dg/de":>10} {">1%":>8} {">2%":>8} {">5%":>8} {">10%":>8} {"VTP/s":>8}')
print('-'*82)
for r in results:
    print(f'{r["label"]:>8} {r["n_vertices"]:>10,} {r["mean_r"]:>10.4f} {r["max_r"]:>10.4f} '
          f'{r["frac_gt1"]:>7.1f}% {r["frac_gt2"]:>7.1f}% {r["frac_gt5"]:>7.1f}% {r["frac_gt10"]:>7.1f}% {r["t_vtp"]:>7.1f}')

# ── 4. 核心发现 ──
print(f'\n核心发现:')
print(f'  12.5m→100m, >2%偏差比例从 {results[0]["frac_gt2"]:.1f}% 降至 {results[-1]["frac_gt2"]:.1f}%')
print(f'  12.5m→100m, 均值偏差从 {results[0]["mean_pct"]:.2f}% 降至 {results[-1]["mean_pct"]:.2f}%')
print(f'  >5%偏差比例从 {results[0]["frac_gt5"]:.1f}% 降至 {results[-1]["frac_gt5"]:.1f}%')

# Save
np.savez_compressed('e:/vtp_geodesic/lushan_resolution_results.npz', results=results)
print('\nSaved: lushan_resolution_results.npz')
print('Done.')
