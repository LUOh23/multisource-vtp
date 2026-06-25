"""
庐山12.5m DEM — 测地/欧氏距离偏差分析
替代原Table 6 (Poyang 30m)
"""
import numpy as np
import rasterio
import vtp_geodesic
import time

print('=' * 60)
print('庐山 12.5m DEM — 地形曲面偏差分析 (Table 6)')
print('=' * 60)

# ── 1. 加载DEM ──
dem = rasterio.open(r'E:\vtp_geodesic\data\dem_clip.tif')
elev = dem.read(1).astype(float)
rows, cols = elev.shape
res_x, res_y = dem.res

x_coords = np.linspace(dem.bounds.left + res_x/2, dem.bounds.right - res_x/2, cols)
y_coords = np.linspace(dem.bounds.top - res_y/2, dem.bounds.bottom + res_y/2, rows)
x_grid, y_grid = np.meshgrid(x_coords, y_coords)
x_offset = x_grid.mean()
y_offset = y_grid.mean()
x_grid -= x_offset
y_grid -= y_offset
dem.close()

# ── 2. 构建三角网格 ──
print('\n构建三角网格...')
points = np.column_stack([x_grid.ravel(), y_grid.ravel(), elev.ravel()]).astype(np.float64)

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
n_verts = len(points)
n_faces = len(faces_vtp)
print(f'顶点: {n_verts:,}  三角面: {n_faces:,}')
print(f'范围: {rows*res_x/1000:.1f}×{cols*res_y/1000:.1f} km')
print(f'高差: {elev.min():.0f}~{elev.max():.0f} m')

# ── 3. 单源VTP ──
peak_idx = int(np.argmax(elev.ravel()))
print(f'\n源点: 山顶 #{peak_idx}, 高程 {elev.ravel()[peak_idx]:.0f}m')

print('运行 VTP...')
t0 = time.perf_counter()
vtp = vtp_geodesic.PyVTPAlgorithmExact(points, faces_vtp)
d_geo = vtp.compute_distances([peak_idx])
t_vtp = time.perf_counter() - t0
print(f'VTP 耗时: {t_vtp:.1f}s')

# ── 4. 欧氏距离 ──
src_pt = points[peak_idx]
d_euc = np.sqrt(np.sum((points - src_pt)**2, axis=1))

# ── 5. 统计 ──
valid = np.isfinite(d_geo) & (d_geo > 1.0) & (d_euc > 1.0)
ratio = d_geo[valid] / (d_euc[valid] + 1e-10)
ratio_pct = (ratio - 1) * 100

print(f'\n{"="*60}')
print(f'RESULTS — 庐山12.5m DEM 地形偏差')
print(f'{"="*60}')
print(f'  网格: {n_verts:,} 顶点, {n_faces:,} 面')
print(f'  高差: {int(elev.min())} ~ {int(elev.max())} m')
print(f'  源点高程: {elev.ravel()[peak_idx]:.0f} m')
print(f'  VTP 时间: {t_vtp:.1f} s')
print(f'  有效点数: {valid.sum():,}')
print(f'  r_mean = d_g/d_e 均值: {np.mean(ratio):.4f}')
print(f'  r_max = d_g/d_e 最大值: {np.max(ratio):.4f}')
print(f'  r_std = d_g/d_e 标准差: {np.std(ratio):.4f}')
print(f'  (d_g/d_e - 1) 均值: {np.mean(ratio_pct):.2f}%')
print(f'  r > 1.01 比例: {100*(ratio_pct>1).sum()/len(ratio_pct):.1f}%')
print(f'  r > 1.02 比例: {100*(ratio_pct>2).sum()/len(ratio_pct):.1f}%')
print(f'  r > 1.05 比例: {100*(ratio_pct>5).sum()/len(ratio_pct):.1f}%')
print(f'  r > 1.10 比例: {100*(ratio_pct>10).sum()/len(ratio_pct):.1f}%')
print(f'  最大绝对偏差: {(d_geo[valid] - d_euc[valid]).max():.0f} m')

# ── 6. 保存数据供后续绘图使用 ──
np.savez_compressed('e:/vtp_geodesic/lushan_results.npz',
    points=points, faces=faces_vtp,
    d_geo=d_geo, d_euc=d_euc, elev=elev.ravel(),
    peak_idx=peak_idx, ratio=ratio, ratio_pct=ratio_pct,
    valid=valid, x_offset=x_offset, y_offset=y_offset,
    rows=rows, cols=cols, res_x=res_x, res_y=res_y,
    x_grid=x_grid.ravel(), y_grid=y_grid.ravel())
print('\nSaved: lushan_results.npz')
print('Done.')
