"""
图 — 庐山研究区地形概况（双面板）
(a) 区域位置: ASTGTM2 30m DEM 渲染，框出研究区
(b) 研究区地形: 12.5m DEM hillshade + 高程渲染，标源点
"""
import numpy as np
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'

# ═══════════════════════════════════════════════
# Load data
# ═══════════════════════════════════════════════
print('Loading DEMs...')

# Panel (a): Regional context from ASTGTM2
dem_30m = rasterio.open(r'E:\vtp_geodesic\高程原始数据\ASTGTM2_N29E116_dem.tif')
elev_30m = dem_30m.read(1).astype(float)
elev_30m[elev_30m < -1000] = np.nan
rows30, cols30 = elev_30m.shape  # 3600 x 3600

# Crop to show Lushan + Poyang Lake context (N29.0°~N29.9°, E115.7°~E116.7°)
# ASTGTM2 tile: N29°-N30°, E116°-E117°. Row 0 = N30°, Col 0 = E116°
# N29.9° = row at (30.0-29.9)*3600 = 360
# N29.0° = row at (30.0-29.0)*3600 = 3600
# E115.7°: outside tile. Let me crop within tile: E116.0°-E116.7°
# Col 0 = E116°, Col 0.7*3600 = 2520
r1, r2 = int((30.0-29.9)*3600), int((30.0-29.0)*3600)  # 360:3600
c1, c2 = 0, int((116.7-116.0)*3600)  # 0:2520
elev_regional = elev_30m[r1:r2, c1:c2]
reg_rows, reg_cols = elev_regional.shape
print(f'Regional: {reg_rows}x{reg_cols}')

# Lon/lat for regional
lons_reg = np.linspace(116.0, 116.7, reg_cols)
lats_reg = np.linspace(29.9, 29.0, reg_rows)
lon_grid, lat_grid = np.meshgrid(lons_reg, lats_reg)

# Panel (b): Study area detail from 12.5m DEM
dem_12m = rasterio.open(r'E:\vtp_geodesic\data\dem_clip.tif')
elev_12m = dem_12m.read(1).astype(float)
rows12, cols12 = elev_12m.shape
res_x = dem_12m.res[0]  # ~12.5m

# UTM coordinates for study area
x_coords = np.linspace(dem_12m.bounds.left, dem_12m.bounds.right, cols12)
y_coords = np.linspace(dem_12m.bounds.bottom, dem_12m.bounds.top, rows12)
# Use km for display
x_km = (x_coords - x_coords.mean()) / 1000
y_km = (y_coords - y_coords.mean()) / 1000
extent_km = [x_km.min(), x_km.max(), y_km.min(), y_km.max()]

peak_idx = np.argmax(elev_12m.ravel())
peak_row, peak_col = peak_idx // cols12, peak_idx % cols12
peak_x_km = x_km[peak_col]
peak_y_km = y_km[peak_row]
peak_elev = elev_12m[peak_row, peak_col]
print(f'Peak: ({peak_x_km:.2f}, {peak_y_km:.2f}) km, {peak_elev:.0f}m')

# ═══════════════════════════════════════════════
# Plot
# ═══════════════════════════════════════════════
fig = plt.figure(figsize=(14, 6.5))

# ── Panel (a): Regional context ──
ax1 = fig.add_axes([0.04, 0.10, 0.38, 0.82])
ls1 = LightSource(azdeg=315, altdeg=45)
elev_filled = elev_regional.copy()
elev_filled[np.isnan(elev_filled)] = np.nanmin(elev_regional)
hs1 = ls1.hillshade(elev_filled, vert_exag=3)

ax1.imshow(hs1, cmap='gray', extent=[lons_reg.min(), lons_reg.max(), lats_reg.min(), lats_reg.max()],
           origin='upper', aspect='auto')
# Elevation overlay with transparency
elev_masked = np.ma.masked_where(np.isnan(elev_regional), elev_regional)
im1 = ax1.imshow(elev_masked, cmap='terrain', alpha=0.35,
                 extent=[lons_reg.min(), lons_reg.max(), lats_reg.min(), lats_reg.max()],
                 origin='upper', aspect='auto', vmin=-50, vmax=1500)

# Mark study area (approximate — the Lushan DEM bounds in lat/lon)
# UTM bounds: X[392143, 404447], Y[3255631, 3266687] ≈ lon[115.97, 116.09], lat[29.42, 29.52]
study_lon = [115.97, 116.09, 116.09, 115.97, 115.97]
study_lat = [29.42, 29.42, 29.52, 29.52, 29.42]
ax1.plot(study_lon, study_lat, 'r-', lw=2.0, zorder=10)
# Place text clearly inside the rectangle (center of box)
center_lon = (115.97 + 116.09) / 2
center_lat = (29.42 + 29.52) / 2
ax1.text(center_lon, center_lat, 'Study\narea', fontsize=9, color='red', fontweight='bold',
         ha='center', va='center',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='red', alpha=0.85))

# Label key features
ax1.text(116.45, 29.62, 'Poyang\nLake', fontsize=9, color='#336699', fontstyle='italic', ha='center', alpha=0.8)
ax1.text(116.05, 29.55, 'Lushan', fontsize=10, color='#8B0000', fontweight='bold', ha='center')

ax1.set_xlabel('Longitude / °E', fontsize=11, fontweight='bold')
ax1.set_ylabel('Latitude / °N', fontsize=11, fontweight='bold')
ax1.set_title('(a) Regional location', fontsize=13, fontweight='bold', loc='left')
ax1.tick_params(labelsize=9)

# ── Panel (b): Study area detail ──
ax2 = fig.add_axes([0.48, 0.10, 0.50, 0.82])
ls2 = LightSource(azdeg=315, altdeg=45)
elev_filled2 = elev_12m.copy()
hs2 = ls2.hillshade(elev_filled2, vert_exag=2.5)

# Hillshade base
ax2.imshow(hs2, cmap='gray', extent=extent_km, origin='upper', aspect='equal', alpha=0.55)

# Elevation color overlay
im2 = ax2.imshow(elev_12m, cmap='terrain', alpha=0.55,
                 extent=extent_km, origin='upper', aspect='equal',
                 vmin=0, vmax=1500)

# Colorbar
cbar_ax = fig.add_axes([0.99, 0.10, 0.012, 0.82])
cbar = fig.colorbar(im2, cax=cbar_ax)
cbar.set_label('Elevation / m', fontsize=11, fontweight='bold')
cbar.ax.tick_params(labelsize=9)

# Mark source point (peak)
ax2.plot(peak_x_km, peak_y_km, '^', color='#CC0000', markersize=14,
         markeredgecolor='white', markeredgewidth=1.5, zorder=15)
ax2.text(peak_x_km-0.8, peak_y_km-1.2, 'Source\n(1460 m)', fontsize=10,
         fontweight='bold', color='#CC0000', ha='center',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#CC0000', alpha=0.85))

# Scale bar (2 km)
bar_x = extent_km[0] + 0.5
bar_y = extent_km[2] + 0.6
ax2.plot([bar_x, bar_x + 2.0], [bar_y, bar_y], 'k-', lw=3.0)
ax2.plot([bar_x, bar_x], [bar_y-0.15, bar_y+0.15], 'k-', lw=2.0)
ax2.plot([bar_x+2.0, bar_x+2.0], [bar_y-0.15, bar_y+0.15], 'k-', lw=2.0)
ax2.text(bar_x+1.0, bar_y+0.35, '2 km', fontsize=10, fontweight='bold', ha='center')

# North arrow
north_x = extent_km[1] - 0.5
north_y = extent_km[2] + 0.8
ax2.annotate('N', xy=(north_x, north_y+0.6), xytext=(north_x, north_y),
             arrowprops=dict(arrowstyle='->', color='black', lw=2.0),
             fontsize=12, fontweight='bold', ha='center')

ax2.set_xlabel('UTM Easting (relative) / km', fontsize=11, fontweight='bold')
ax2.set_ylabel('UTM Northing (relative) / km', fontsize=11, fontweight='bold')
ax2.set_title('(b) Lushan study area (12.5 m DEM)', fontsize=13, fontweight='bold', loc='left')
ax2.tick_params(labelsize=9)

plt.savefig('e:/vtp_geodesic/论文截图/fig_study_area.png', dpi=200,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print('Study area figure saved → 论文截图/fig_study_area.png')
