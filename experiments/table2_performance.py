"""
实验2: 大规模网格性能测试
对比 VTP Exact vs Heat Method 在不同网格规模下的精度和速度
"""
import numpy as np
import pyvista as pv
import vtp_geodesic
import potpourri3d as p3d
import matplotlib.pyplot as plt
import time

# ============================================================
# 1. 生成多尺度测试网格 (球面 — 标准测试体)
# ============================================================
resolutions = [8, 12, 16, 24, 32, 48, 64, 80, 100, 140, 180]
results = []

for res in resolutions:
    print(f"\n{'='*50}")
    print(f"分辨率={res} → ", end='', flush=True)

    # 生成球面网格
    sphere = pv.Sphere(radius=100.0, theta_resolution=res, phi_resolution=res)
    points = sphere.points.astype(np.float64)
    faces_raw = sphere.faces.reshape(-1, 4)[:, 1:].astype(np.int32)

    n_verts = points.shape[0]
    n_faces = faces_raw.shape[0]

    # --- VTP Exact ---
    t0 = time.perf_counter()
    vtp = vtp_geodesic.PyVTPAlgorithmExact(points, faces_raw)
    dt_init = time.perf_counter() - t0

    # 用顶点0作为源点
    t0 = time.perf_counter()
    d_vtp = vtp.compute_distances([0])
    dt_vtp = time.perf_counter() - t0

    # --- Heat Method ---
    t0 = time.perf_counter()
    sol = p3d.MeshHeatMethodDistanceSolver(points, faces_raw)
    d_heat = sol.compute_distance(0)
    dt_heat = time.perf_counter() - t0

    # --- 误差统计 ---
    valid = np.isfinite(d_vtp) & (d_vtp > 0) & np.isfinite(d_heat)
    rel_error = np.abs(d_heat[valid] - d_vtp[valid]) / (d_vtp[valid] + 1e-10)
    mean_err = np.mean(rel_error) * 100
    max_err = np.max(rel_error) * 100

    results.append({
        'n': n_verts,
        'n_faces': n_faces,
        't_init': dt_init,
        't_vtp': dt_vtp,
        't_heat': dt_heat,
        'err_mean': mean_err,
        'err_max': max_err,
    })

    print(f"{n_verts}顶点 {n_faces}面 | "
          f"VTP={dt_vtp:.3f}s Heat={dt_heat:.4f}s "
          f"误差均值={mean_err:.2f}% 最大={max_err:.2f}%")

# ============================================================
# 2. 画图
# ============================================================
ns = [r['n'] for r in results]
ns_faces = [r['n_faces'] for r in results]
tv = [r['t_vtp'] for r in results]
th = [r['t_heat'] for r in results]
er = [r['err_mean'] for r in results]
emax = [r['err_max'] for r in results]

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# 左: 时间 vs 网格规模
ax = axes[0]
ax.plot(ns, tv, 'o-', color='#d62728', linewidth=2, markersize=8, label='VTP Exact')
ax.plot(ns, th, 's-', color='#1f77b4', linewidth=2, markersize=8, label='Heat Method')
ax.set_xlabel('Vertices', fontsize=12)
ax.set_ylabel('Time (s)', fontsize=12)
ax.set_title('Computation Time vs Mesh Size', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_xscale('log')
ax.set_yscale('log')

# 中: 误差 vs 网格规模
ax = axes[1]
ax.fill_between(ns, 0, emax, alpha=0.2, color='orange', label='Max Error')
ax.plot(ns, er, 'o-', color='#d62728', linewidth=2, markersize=8, label='Mean Error')
ax.set_xlabel('Vertices', fontsize=12)
ax.set_ylabel('Relative Error (%)', fontsize=12)
ax.set_title('Heat Method vs VTP (ground truth)', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_xscale('log')

# 右: 精度-速度 Trade-off (每个点是一个网格规模)
ax = axes[2]
sc = ax.scatter(tv, er, c=ns, cmap='viridis', s=120, edgecolors='k', linewidth=0.5)
for i, n in enumerate(ns):
    ax.annotate(f'{n//1000}k' if n>=1000 else str(n),
                (tv[i], er[i]), textcoords="offset points", xytext=(8, 5), fontsize=8)
ax.set_xlabel('VTP Time (s)', fontsize=12)
ax.set_ylabel('Heat Method Mean Error (%)', fontsize=12)
ax.set_title('Precision–Speed Trade-off', fontsize=14)
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('Vertices', fontsize=11)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('experiment2_performance.png', dpi=300, bbox_inches='tight')
print(f"\n图片已保存: experiment2_performance.png")
plt.show()

# ============================================================
# 3. 关键数据
# ============================================================
print("\n" + "="*60)
print("实验2 关键数据汇总")
print("="*60)
print(f"{'顶点':<10} {'VTP(s)':<10} {'Heat(s)':<10} {'误差均值%':<12} {'误差最大%':<12}")
print("-"*54)
for r in results:
    print(f"{r['n']:<10} {r['t_vtp']:<10.4f} {r['t_heat']:<10.4f} "
          f"{r['err_mean']:<12.3f} {r['err_max']:<12.3f}")

# 论文核心句
r200k = results[-1]
print(f"\n论文可用结论:")
print(f"- 在{r200k['n']}顶点网格上，VTP 耗时{r200k['t_vtp']:.2f}s，Heat Method 耗时{r200k['t_heat']:.4f}s")
print(f"- Heat Method 平均误差 {np.mean(er):.2f}%，最大误差 {np.max(emax):.2f}%")
print(f"- VTP 提供精确结果，Heat Method 在速度上有 {np.mean(tv)/np.mean(th):.0f}× 优势但以精度为代价")
