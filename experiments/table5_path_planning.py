"""
交互式测地路径规划: 点击选起终点, 自动避开水体, VTP算法计算最短路径
"""
import numpy as np
import pyvista as pv
import vtp_geodesic

# ============================================================
# 1. 加载数据 & 过滤水体
# ============================================================
print("加载数据...")
mesh_raw = pv.read(r'E:\vtp_geodesic\高程原始数据\final_3d_terrain_converted.vtk')
all_points = mesh_raw.points.astype(np.float64)
all_faces_raw = mesh_raw.faces.reshape(-1, 4)[:, 1:].astype(np.int32)  # Mx3
in_water = mesh_raw.cell_data['in_water']

# 分离陆地/水体
land_mask = (in_water == 0)
water_mask = (in_water == 1)
land_faces = all_faces_raw[land_mask]   # 陆地三角面
water_faces = all_faces_raw[water_mask] # 水体三角面

# 找陆地面连通分量 (避免断开的水体面导致半边结构错误)
print("查找陆地最大连通分量...")
n_faces = len(land_faces)
face_adj = [[] for _ in range(n_faces)]
edge_to_faces = {}
for fi, f in enumerate(land_faces):
    for i in range(3):
        a, b = sorted([f[i], f[(i+1)%3]])
        key = (a, b)
        if key in edge_to_faces:
            edge_to_faces[key].append(fi)
        else:
            edge_to_faces[key] = [fi]

for edge_faces in edge_to_faces.values():
    if len(edge_faces) == 2:
        fi, fj = edge_faces
        face_adj[fi].append(fj)
        face_adj[fj].append(fi)

# BFS 找最大分量
visited = np.zeros(n_faces, dtype=bool)
largest_comp = []
for fi in range(n_faces):
    if not visited[fi]:
        comp = []
        stack = [fi]
        visited[fi] = True
        while stack:
            f = stack.pop()
            comp.append(f)
            for nf in face_adj[f]:
                if not visited[nf]:
                    visited[nf] = True
                    stack.append(nf)
        if len(comp) > len(largest_comp):
            largest_comp = comp

print(f"陆地: {n_faces}面, {len(edge_to_faces)}边, 最大连通分量={len(largest_comp)}面")

# 只用最大分量
land_faces_main = land_faces[largest_comp]

# 重新映射顶点
used_verts = np.unique(land_faces_main.ravel())
old_to_new = np.full(all_points.shape[0], -1, dtype=np.int32)
old_to_new[used_verts] = np.arange(len(used_verts))

points = all_points[used_verts].copy()
faces_vtp = old_to_new[land_faces_main]

# 构建PyVista格式
faces_pv = np.column_stack([np.full((len(faces_vtp), 1), 3, dtype=np.int32), faces_vtp]).ravel()

print(f"原始: {all_points.shape[0]}顶点, {len(all_faces_raw)}面")
print(f"陆地: {len(points)}顶点, {len(faces_vtp)}面 (去除{water_mask.sum()}个水面, {n_faces - len(largest_comp)}个断片)")
print(f"XY: {points[:,0].min():.0f}~{points[:,0].max():.0f}, {points[:,1].min():.0f}~{points[:,1].max():.0f}")
print(f"Z: {points[:,2].min():.1f}~{points[:,2].max():.1f}")

# ============================================================
# 1.5 三角网细分 (细化网格，路径更平滑)
# ============================================================
def subdivide_mesh(pts, faces, iterations=2):
    """中点细分: 每个三角形裂成4个，迭代iterations次"""
    cur_pts = pts.copy()
    cur_faces = faces.copy()
    for it in range(iterations):
        edge_mid = {}  # (a,b)排序key → 中点索引
        new_faces = []
        for f in cur_faces:
            mids = []
            for k in range(3):
                e = tuple(sorted([int(f[k]), int(f[(k+1)%3])]))
                if e not in edge_mid:
                    mp = (cur_pts[e[0]] + cur_pts[e[1]]) * 0.5
                    edge_mid[e] = len(cur_pts)
                    cur_pts = np.vstack([cur_pts, mp])
                mids.append(edge_mid[e])
            # 4个子三角形
            a, b, c = f[0], f[1], f[2]
            m0, m1, m2 = mids  # m0=mid(a,b), m1=mid(b,c), m2=mid(c,a)
            new_faces.extend([
                [a, m0, m2], [m0, b, m1],
                [m0, m1, m2], [m2, m1, c]
            ])
        cur_faces = np.array(new_faces, dtype=np.int32)
        print(f"  细分{it+1}: {len(cur_pts)}顶点, {len(cur_faces)}面")
    return cur_pts, cur_faces

print("三角网细分...")
points, faces_vtp = subdivide_mesh(points, faces_vtp, iterations=2)
faces_pv = np.column_stack([np.full((len(faces_vtp), 1), 3, dtype=np.int32), faces_vtp]).ravel()

# ============================================================
# 2. 预构建邻接表 (用于路径回溯)
# ============================================================
print("构建邻接表 & 面邻接...")
# 顶点邻接
adj = [[] for _ in range(len(points))]
for f in faces_vtp:
    for i in range(3):
        a, b = f[i], f[(i+1)%3]
        if b not in adj[a]: adj[a].append(b)
        if a not in adj[b]: adj[b].append(a)

# 面邻接 + 边→面 映射 (用于连续梯度追踪)
edge2faces = {}
for fi, f in enumerate(faces_vtp):
    for i in range(3):
        e = tuple(sorted([int(f[i]), int(f[(i+1)%3])]))
        edge2faces.setdefault(e, []).append(fi)

face_neighbors = [[] for _ in range(len(faces_vtp))]
for e, flist in edge2faces.items():
    if len(flist) == 2:
        face_neighbors[flist[0]].append(flist[1])
        face_neighbors[flist[1]].append(flist[0])

# 找顶点所在的面列表
vert2faces = [[] for _ in range(len(points))]
for fi, f in enumerate(faces_vtp):
    for vi in f:
        vert2faces[vi].append(fi)

def trace_continuous_path(distances, start_v, end_v, step=1.5):
    """连续梯度追踪: 在三角形内沿距离梯度下降, 生成平滑路径"""
    pts = [points[end_v].copy()]
    cur_v = end_v

    for _ in range(10000):
        if cur_v == start_v:
            break

        # 在当前顶点的邻域找最近的下一个顶点 (用于初始化面)
        best_nbr = min(adj[cur_v], key=lambda n: distances[n])
        if distances[best_nbr] >= distances[cur_v]:
            break

        # 用当前顶点所在的面进行一步梯度追踪
        cand_faces = vert2faces[cur_v]

        # 在共享 best_nbr 的面中追踪
        shared_faces = [fi for fi in cand_faces
                       if best_nbr in faces_vtp[fi]]

        if not shared_faces:
            cur_v = best_nbr
            pts.append(points[cur_v].copy())
            continue

        # 取第一个共享面, 进行亚网格梯度追踪
        fi = shared_faces[0]
        f = faces_vtp[fi]
        # 三角形三个顶点的距离
        d0, d1, d2 = distances[f[0]], distances[f[1]], distances[f[2]]
        p0, p1, p2 = points[f[0]], points[f[1]], points[f[2]]

        # 当前点到对边的梯度追踪
        # 从当前顶点沿三角形向内追踪直到到达best_nbr
        local_cur = p0 if f[0] == cur_v else (p1 if f[1] == cur_v else p2)
        local_target = p0 if f[0] == best_nbr else (p1 if f[1] == best_nbr else p2)
        local_other = [p for p in [p0, p1, p2]
                       if not np.array_equal(p, local_cur) and not np.array_equal(p, local_target)][0]

        # 在三角形内插值追踪
        dist_to_target = np.linalg.norm(local_target - local_cur)
        n_steps = max(2, int(dist_to_target / step))
        for s in range(1, n_steps + 1):
            t = s / n_steps
            pt = local_cur + t * (local_target - local_cur)
            pts.append(pt.copy())

        cur_v = best_nbr

    # 反转得到起点→终点
    pts = pts[::-1]
    # 去掉重复的起点/终点附近的点
    if len(pts) > 2:
        pts = [p for i, p in enumerate(pts)
               if i == 0 or i == len(pts)-1 or
               np.linalg.norm(p - pts[i-1]) > step * 0.3]

    # 三次样条平滑
    pts_arr = np.array(pts)
    if len(pts_arr) >= 4:
        # 对每段做简单移动平均平滑
        smooth = pts_arr.copy()
        for _ in range(3):
            for i in range(1, len(smooth)-1):
                smooth[i] = (pts_arr[i-1] + pts_arr[i] + pts_arr[i+1]) / 3
            pts_arr = smooth.copy()

    return pts_arr

# ============================================================
# 3. VTP 初始化
# ============================================================
print("初始化VTP求解器...")
solver = vtp_geodesic.PyVTPAlgorithmExact(points, faces_vtp)
print("完成!\n")

# ============================================================
# 4. 可视化 & 交互
# ============================================================
# 陆地网格
land_mesh = pv.PolyData(points, faces_pv)
land_mesh.point_data['elev'] = points[:, 2]

# 水体网格 — 用原始顶点, 不重新映射
water_faces_pv_raw = np.column_stack([
    np.full((len(water_faces), 1), 3, dtype=np.int32), water_faces
]).ravel()
water_mesh = pv.PolyData(all_points, water_faces_pv_raw)

pl = pv.Plotter(window_size=(1500, 900))

# 水体 (蓝色半透明)
pl.add_mesh(water_mesh, color='#4a90d9', opacity=0.55,
            show_edges=False, smooth_shading=True, label='Water (non-walkable)')

# 陆地 (高程着色 + 浅色三角边线)
pl.add_mesh(land_mesh, scalars='elev', cmap='terrain', opacity=1.0,
            show_edges=True, edge_color='#666633', line_width=0.3,
            smooth_shading=True, label='Land',
            scalar_bar_args={'title': 'Elev (m)', 'vertical': True,
                            'position_x': 0.02, 'position_y': 0.3,
                            'height': 0.4})

pl.add_legend()

# 全局状态
state = {
    'start_idx': None,
    'end_idx': None,
    'start_pt_actor': None,
    'end_pt_actor': None,
    'path_actor': None,
    'profile_actor': None,
}

def clear_overlay():
    for k in ('start_pt_actor', 'end_pt_actor', 'path_actor', 'profile_actor'):
        if state[k] is not None:
            pl.remove_actor(state[k])
            state[k] = None

def add_start_marker(idx):
    state['start_pt_actor'] = pl.add_points(
        points[idx:idx+1], color='#00ff00', point_size=40,
        render_points_as_spheres=True)

def add_end_marker(idx):
    state['end_pt_actor'] = pl.add_points(
        points[idx:idx+1], color='#ff00ff', point_size=40,
        render_points_as_spheres=True)

def compute_and_show(sidx, eidx):
    clear_overlay()

    # ---- VTP ----
    distances = solver.compute_distances([sidx])
    geo_dist = distances[eidx]
    print(f"路径: {sidx}→{eidx} | 起点({points[sidx,0]:.0f},{points[sidx,1]:.0f},{points[sidx,2]:.1f}) "
          f"终点({points[eidx,0]:.0f},{points[eidx,1]:.0f},{points[eidx,2]:.1f}) | 距离={geo_dist:.1f}m")

    if not np.isfinite(geo_dist):
        pl.add_text("Unreachable! Start/end separated by water.",
                    position='upper_center', font_size=16, color='red')
        return

    # ---- 连续梯度追踪平滑路径 ----
    path_pts = trace_continuous_path(distances, sidx, eidx, step=1.5)
    n_pts = len(path_pts)

    # 稍微抬高让管不陷入面
    offset = np.array([0, 0, 0.5])
    path_pts_raised = path_pts + offset

    path_poly = pv.PolyData(path_pts_raised)
    path_poly.lines = np.column_stack([
        np.full(n_pts-1, 2),
        np.arange(n_pts-1),
        np.arange(1, n_pts)
    ]).ravel()
    tube = path_poly.tube(radius=2.5, n_sides=8)
    state['path_actor'] = pl.add_mesh(tube, color='#ff2222', smooth_shading=True)

    add_start_marker(sidx)
    add_end_marker(eidx)

    pl.add_text(f'Dist: {geo_dist:.1f}m  |  Points: {n_pts}',
                position='upper_left', font_size=13, color='white')

def pick_callback(point):
    # 最近顶点
    vidx = np.argmin(np.sum((points - point)**2, axis=1))

    # 状态机: None → start → end → reset
    if state['start_idx'] is None:
        # 第一次: 设起点
        state['start_idx'] = vidx
        clear_overlay()
        add_start_marker(vidx)
        pl.add_text('Start selected — click to place END',
                    position='upper_right', font_size=13, color='lime')
        print(f"S: {vidx} ({points[vidx,0]:.0f},{points[vidx,1]:.0f},{points[vidx,2]:.1f})")

    elif state['end_idx'] is None:
        # 第二次: 设终点并计算
        state['end_idx'] = vidx
        compute_and_show(state['start_idx'], state['end_idx'])
        pl.add_text('Click to start new path | Close window to exit',
                    position='lower_right', font_size=11, color=(0.7,0.7,0.7))

    else:
        # 已有路径: 重置, 新起点
        clear_overlay()
        state['start_idx'] = vidx
        state['end_idx'] = None
        add_start_marker(vidx)
        pl.add_text('Start selected — click to place END',
                    position='upper_right', font_size=13, color='lime')
        print(f"\nS: {vidx} ({points[vidx,0]:.0f},{points[vidx,1]:.0f},{points[vidx,2]:.1f})")

pl.enable_surface_point_picking(callback=pick_callback, show_point=False,
                                left_clicking=True, pickable_window=False)

pl.add_text('CLICK terrain to pick START point',
            position='upper_left', font_size=15, color='yellow')
pl.add_text('Land (terrain) | Water (blue) | Green=Start | Magenta=End | Red=Path',
            position='lower_left', font_size=10, color=(0.8,0.8,0.8))

# 更好的观察角度
pl.camera_position = [(points[:,0].mean(), points[:,1].mean() - 2000, points[:,2].max() * 2),
                      (points[:,0].mean(), points[:,1].mean(), points[:,2].mean()),
                      (0, 0, 1)]
pl.show()
