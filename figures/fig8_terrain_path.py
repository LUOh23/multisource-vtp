"""
Fig.7 — 三维地形与测地路径 (Sepia地形 + Teal水体 + Gold路径 + 黑边)
"""
import numpy as np
import pyvista as pv
import vtp_geodesic

print('Loading water mesh...')
mesh_raw = pv.read(r'E:\vtp_geodesic\高程原始数据\final_3d_terrain_converted.vtk')
all_points = mesh_raw.points.astype(np.float64)
all_faces_raw = mesh_raw.faces.reshape(-1, 4)[:, 1:].astype(np.int32)
in_water = mesh_raw.cell_data['in_water']

land_faces = all_faces_raw[in_water == 0]
water_faces_raw = all_faces_raw[in_water == 1]
print(f'Original: {len(all_points)} verts, {len(all_faces_raw)} faces ({len(water_faces_raw)} water)')

# Largest land component
print('Finding largest land component...')
edge_to_faces = {}
for fi, f in enumerate(land_faces):
    for i in range(3):
        a, b = sorted([int(f[i]), int(f[(i+1)%3])])
        edge_to_faces.setdefault((a, b), []).append(fi)
face_adj = [[] for _ in range(len(land_faces))]
for flist in edge_to_faces.values():
    if len(flist) == 2:
        face_adj[flist[0]].append(flist[1])
        face_adj[flist[1]].append(flist[0])
visited = np.zeros(len(land_faces), dtype=bool)
largest_comp = []
for fi in range(len(land_faces)):
    if not visited[fi]:
        comp = []; stack = [fi]; visited[fi] = True
        while stack:
            f = stack.pop(); comp.append(f)
            for nf in face_adj[f]:
                if not visited[nf]: visited[nf] = True; stack.append(nf)
        if len(comp) > len(largest_comp): largest_comp = comp
land_faces_main = land_faces[largest_comp]
used_verts = np.unique(land_faces_main.ravel())
old_to_new = np.full(all_points.shape[0], -1, dtype=np.int32)
old_to_new[used_verts] = np.arange(len(used_verts))
points = all_points[used_verts].copy()
faces_vtp = old_to_new[land_faces_main]
print(f'Land: {len(points)} verts, {len(faces_vtp)} faces')

# Subdivide
def subdivide(pts, faces, iterations=2):
    cur_pts, cur_faces = pts.copy(), faces.copy()
    for it in range(iterations):
        edge_mid = {}; new_faces = []
        for f in cur_faces:
            mids = []
            for k in range(3):
                e = tuple(sorted([int(f[k]), int(f[(k+1)%3])]))
                if e not in edge_mid:
                    edge_mid[e] = len(cur_pts)
                    cur_pts = np.vstack([cur_pts, (cur_pts[e[0]]+cur_pts[e[1]])*0.5])
                mids.append(edge_mid[e])
            a,b,c = f[0],f[1],f[2]; m0,m1,m2 = mids
            new_faces.extend([[a,m0,m2],[m0,b,m1],[m0,m1,m2],[m2,m1,c]])
        cur_faces = np.array(new_faces, dtype=np.int32)
        print(f'  Subdiv {it+1}: {len(cur_pts)} verts, {len(cur_faces)} faces')
    return cur_pts, cur_faces

print('Subdividing...')
points, faces_vtp = subdivide(points, faces_vtp, iterations=2)
faces_pv = np.column_stack([np.full((len(faces_vtp),1),3,dtype=np.int32), faces_vtp]).ravel()

# Build adjacency
adj = [[] for _ in range(len(points))]
for f in faces_vtp:
    for i in range(3):
        a,b = f[i], f[(i+1)%3]
        if b not in adj[a]: adj[a].append(b)
        if a not in adj[b]: adj[b].append(a)

# Pick start/end
target_start = np.array([380560.0, 2828296.0, 228.0])
target_end   = np.array([381882.0, 2828266.0, 214.0])
start_idx = int(np.argmin(np.sum((points-target_start)**2, axis=1)))
end_idx   = int(np.argmin(np.sum((points-target_end)**2, axis=1)))
print(f'Start: {start_idx}, End: {end_idx}')

# VTP
print('Running VTP...')
solver = vtp_geodesic.PyVTPAlgorithmExact(points, faces_vtp)
distances = solver.compute_distances([start_idx])
geo_dist = distances[end_idx]
print(f'Distance: {geo_dist:.1f} m')

# Backtrack path
print('Backtracking...')
cur_v = end_idx; path_raw = [points[cur_v].copy()]
for _ in range(20000):
    if cur_v == start_idx: break
    best = min(adj[cur_v], key=lambda n: distances[n])
    if distances[best] >= distances[cur_v]: break
    ns = max(2, int(np.linalg.norm(points[best]-points[cur_v])/0.8))
    for s in range(1, ns+1):
        path_raw.append(points[cur_v]+(s/ns)*(points[best]-points[cur_v]))
    cur_v = best
path_raw.append(points[start_idx].copy())
path_pts = np.array(path_raw[::-1])
for _ in range(2):
    for i in range(1, len(path_pts)-1):
        path_pts[i] = (path_pts[i-1]+path_pts[i]+path_pts[i+1])/3
print(f'Path: {len(path_pts)} points')

# ══════════════════════════════════════════════════════════
# RENDER — Sepia + Teal + Gold scheme
# ══════════════════════════════════════════════════════════
print('Rendering...')
pl = pv.Plotter(window_size=(2600, 2000), off_screen=True)

# Water: teal (#006666)
water_faces_pv = np.column_stack([np.full((len(water_faces_raw),1),3,dtype=np.int32), water_faces_raw]).ravel()
water_mesh = pv.PolyData(all_points, water_faces_pv)
pl.add_mesh(water_mesh, color='#006666', opacity=0.70, show_edges=False,
            smooth_shading=True, label='Water (non-walkable)')

# Terrain: Sepia tones via gist_earth colormap (brownish)
land_mesh = pv.PolyData(points, faces_pv)
land_mesh.point_data['elev'] = points[:, 2]
pl.add_mesh(land_mesh, scalars='elev', cmap='gist_earth', opacity=1.0,
            show_edges=False, smooth_shading=True, lighting=True,
            specular=0.12, specular_power=8,
            scalar_bar_args={'title': 'Elev / m', 'vertical': True,
                            'position_x': 0.02, 'position_y': 0.25,
                            'height': 0.35, 'font_family': 'times',
                            'n_labels': 6, 'color': 'black'})

# Path: gold (#FFD700) thick tube + black edge tube (slightly wider)
path_raised = path_pts + np.array([0, 0, 2.0])
n_pts = len(path_pts)
path_poly = pv.PolyData(path_raised)
path_poly.lines = np.column_stack([np.full(n_pts-1,2), np.arange(n_pts-1), np.arange(1,n_pts)]).ravel()

# Black "edge" tube behind gold for definition
black_tube = path_poly.tube(radius=6.5, n_sides=14)
pl.add_mesh(black_tube, color='black', opacity=1.0, smooth_shading=True)

# Gold main tube
gold_tube = path_poly.tube(radius=4.5, n_sides=14)
pl.add_mesh(gold_tube, color='#FFD700', opacity=1.0, smooth_shading=True,
            specular=0.6, specular_power=30)

# Start: green sphere
pl.add_points(points[start_idx:start_idx+1]+np.array([0,0,3.5]),
              color='#00AA00', point_size=60, render_points_as_spheres=True)
# End: red sphere
pl.add_points(points[end_idx:end_idx+1]+np.array([0,0,3.5]),
              color='#DD2222', point_size=60, render_points_as_spheres=True)

# Camera: 45-degree overhead oblique view from southeast
center = np.array([points[:,0].mean(), points[:,1].mean(), points[:,2].mean()])
bounds = points[:, 0].max()-points[:, 0].min(), points[:, 1].max()-points[:, 1].min()
pl.camera_position = [(center[0] + max(bounds)*0.7,
                        center[1] - max(bounds)*0.7,
                        center[2] + (points[:,2].max()-points[:,2].min())*3.0),
                       center,
                       (0, 0, 1)]

# Labels
pl.add_text(f'Geodesic Path with Water Obstacle Avoidance\n'
            f'Distance: {geo_dist:.0f} m  |  Path points: {n_pts}',
            position='upper_edge', font_size=15, color='black', font='times')
pl.add_text('o Start   o End   Gold line = Geodesic path   Teal = Water',
            position='lower_left', font_size=12, color='#333333', font='times')
pl.set_background('white')
pl.enable_anti_aliasing()

pl.screenshot('output/fig7_path_planning.png', return_img=False)
pl.close()
print('Fig.7 saved')
