"""
实验4: 多源加速比测试
对比: k次单源 propagate vs 1次多源 propagate
"""
import numpy as np
import pyvista as pv
import vtp_geodesic
import matplotlib.pyplot as plt
import time

# 1. 加载网格
mesh = pv.read('sh10.vtk')
points = mesh.points.astype(np.float64)
faces = mesh.faces.reshape(-1, 4)[:, 1:].astype(np.int32)
n_verts = points.shape[0]
print(f"网格: {n_verts} 顶点, {faces.shape[0]} 面")

# 2. 测试不同源点数量
source_counts = [1, 2, 5, 10, 20, 50]
results = []

for k in source_counts:
    sources = list(range(0, k))  # 前k个顶点

    # --- 方式A: k次单源（模拟原论文）---
    # 每次创建新 solver，模拟独立的单源计算
    t_single_total = 0
    for s in sources:
        solver = vtp_geodesic.PyVTPAlgorithmExact(points, faces)
        t0 = time.perf_counter()
        solver.compute_distances([s])
        t_single_total += time.perf_counter() - t0

    # --- 方式B: 1次多源（你的改进）---
    solver_multi = vtp_geodesic.PyVTPAlgorithmExact(points, faces)
    t0 = time.perf_counter()
    d_multi = solver_multi.compute_distances(sources)
    t_multi = time.perf_counter() - t0

    # 验证正确性：多源结果 = 每个顶点到最近源点的距离
    # 用单源结果逐个比较，取 min
    d_min_single = np.full(n_verts, np.inf)
    for s in sources:
        solver = vtp_geodesic.PyVTPAlgorithmExact(points, faces)
        d_single = solver.compute_distances([s])
        d_min_single = np.minimum(d_min_single, d_single)

    valid = np.isfinite(d_multi) & np.isfinite(d_min_single)
    consistency = np.allclose(d_multi[valid], d_min_single[valid])

    speedup = t_single_total / t_multi

    results.append({
        'k': k,
        't_single': t_single_total,
        't_multi': t_multi,
        'speedup': speedup,
        'correct': consistency
    })
    print(f"k={k:3d} | 单源总时间={t_single_total:.3f}s | 多源时间={t_multi:.3f}s | "
          f"加速比={speedup:.1f}× | 结果一致={consistency}")

# 3. 画图
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

ks = [r['k'] for r in results]
ts = [r['t_single'] for r in results]
tm = [r['t_multi'] for r in results]
sp = [r['speedup'] for r in results]

# 左: 绝对时间对比
axes[0].bar(np.arange(len(ks)) - 0.15, ts, 0.3, label='Single-source (k calls)')
axes[0].bar(np.arange(len(ks)) + 0.15, tm, 0.3, label=f'Multi-source (1 call)')
axes[0].set_xticks(range(len(ks)))
axes[0].set_xticklabels([f'k={k}' for k in ks])
axes[0].set_ylabel('Time (s)')
axes[0].set_title('Computation Time')
axes[0].legend()

# 中: 加速比
axes[1].plot(ks, sp, 'o-', linewidth=2, markersize=8)
axes[1].plot([1, max(ks)], [1, max(ks)], 'k--', alpha=0.3, label='Linear (ideal)')
axes[1].set_xlabel('Number of Sources (k)')
axes[1].set_ylabel('Speedup (×)')
axes[1].set_title('Multi-Source Speedup')
axes[1].legend()

# 右: 效率
eff = [s / k for s, k in zip(sp, ks)]
axes[2].plot(ks, eff, 's-', linewidth=2, color='green')
axes[2].axhline(y=1.0, color='k', linestyle='--', alpha=0.3)
axes[2].set_xlabel('Number of Sources (k)')
axes[2].set_ylabel('Efficiency (speedup/k)')
axes[2].set_title('Parallel Efficiency')

plt.tight_layout()
plt.savefig('multi_source_speedup.png', dpi=300)
print("\n图片已保存: multi_source_speedup.png")

# 4. 输出论文可用的关键数字
print("\n=== 论文关键数据 ===")
print(f"k=1 基准时间: {ts[0]:.3f}s")
print(f"k=10 加速比: {sp[ks.index(10) if 10 in ks else 2]:.1f}×")
print(f"k=50 加速比: {sp[ks.index(50) if 50 in ks else -1]:.1f}×")
print(f"全部结果与单源取min一致: {all(r['correct'] for r in results)}")
