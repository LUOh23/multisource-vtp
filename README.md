# Multi-Source VTP Exact Geodesic Algorithm

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Code and experiments accompanying the paper:

> **罗昊, 段新桥.** 多源VTP波前互抑效应的发现与地形影响因素分析[J]. 测绘科学, 2026.

A Cython Python wrapper for the C++ VTP (Virtual Times Parallel) exact geodesic algorithm [Qin et al., ACM TOG 2016], extended with multi-source synchronous wavefront propagation and window-counting instrumentation.

---

## Quick Start

```bash
git clone https://github.com/LUOh23/multisource-vtp.git
cd multisource-vtp
pip install -r requirements.txt
python setup.py build_ext --inplace
```

---

## Reproducing Paper Results

### Tables

| Table | Script | Description |
|-------|--------|-------------|
| Table 2 | `experiments/experiment_resolution_effect.py` | DEM resolution effect (Poyang Lake, 30-300 m) |
| Table 3 | Window-counting instrumentation (C++ source) | Multi-source vs. single-source window creation |
| Table 4 | `experiments/table2_performance.py` | VTP vs. Heat Method across 11 sphere meshes |
| Table 5 | `experiments/table3_speedup.py` | Multi-source speedup on Sh10 mesh (k=1-50) |
| Table 6 | `experiments/lushan_terrain_deviation.py` | Geodesic/Euclidean deviation on Lushan 12.5 m DEM |
| Table 7 | Synthesis table in paper | Multi-source VTP vs. baseline methods |

### Figures

| Figure | Script | Description |
|--------|--------|-------------|
| Fig. 1 | `figures/fig1_geometry.py` | VTP wavefront propagation geometry |
| Fig. 2 | `figures/fig2_flowchart.py` | Algorithm comparison (single-source vs. multi-source) |
| Fig. 3 | `figures/fig_study_area.py` | Lushan study area overview (location + hillshade) |
| Fig. 4 | `experiments/experiment_resolution_effect.py` | DEM resolution vs. geodesic/Euclidean deviation |
| Fig. 5 | `figures/fig_suppression_effect.py` | Multi-source wavefront mutual suppression effect |
| Fig. 6 | `figures/fig3_performance.py` | VTP vs. Heat Method computation time |
| Fig. 7 | `figures/fig4_error.py` | Heat Method approximation error vs. mesh size |
| Fig. 8 | `figures/fig5_speedup.py` | Multi-source speedup curve (Sh10 mesh) |
| Fig. 9 | `figures/fig_lushan_scatter.py` | Geodesic vs. Euclidean distance scatter (Lushan) |
| Fig. 10 | `figures/fig_lushan_histogram.py` | Distance ratio distribution histogram (Lushan) |

---

## Data

### Lushan 12.5 m DEM (Table 6, Figs. 3, 9-10)

Requires ALOS World 3D DSM data (`dem_clip.tif`) placed in `data/`. The exact tile used in the paper covers UTM 50N coordinates X: 392.1-404.4 km, Y: 3255.6-3266.7 km (~136 km², 885x985 pixels at 12.5 m resolution).

### Poyang Lake 30 m DEM (Tables 2, 4; Fig. 4)

Requires ASTER GDEM v2 tile `ASTGTM2_N29E116_dem.tif`:

1. Register at [NASA Earthdata](https://urs.earthdata.nasa.gov/) (free)
2. Search for `ASTGTM2_N29E116` at [Earthdata Search](https://search.earthdata.nasa.gov/)
3. Download and place in `data/ASTER_GDEM/`

### Sphere & Sh10 meshes (Tables 3-5, Figs. 5-8)

No external data required. Sphere meshes are generated programmatically via PyVista.

---

## Citation

```bibtex
@article{luo2026multisource,
  author  = {罗昊 and 段新桥},
  title   = {多源VTP波前互抑效应的发现与地形影响因素分析},
  journal = {测绘科学},
  year    = {2026},
}
```

## License

MIT
