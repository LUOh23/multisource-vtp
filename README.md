# Multi-Source VTP Exact Geodesic Algorithm

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Code and experiments accompanying the paper:

> **罗昊, 段新桥.** 多源VTP精确测地算法的精度与效率分析[J]. 测绘科学, 2026.

A Cython Python wrapper for the C++ VTP (Virtual Times Parallel) exact geodesic algorithm [Qin et al., ACM TOG 2016], extended with multi-source synchronous wavefront propagation.

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
| Table 2 | `experiments/table2_performance.py` | VTP vs. Heat Method across 11 sphere meshes |
| Table 3 | `experiments/table3_speedup.py` | Multi-source speedup on Sh10 mesh ($k=1$–$50$) |
| Table 4 | `experiments/table4_terrain_ratio.py` | Geodesic/Euclidean ratio on Poyang Lake DEM |
| Table 5 | `experiments/table5_path_planning.py` | Water-obstacle path planning |

Tables 2–3 run immediately. Tables 4–5 require ASTER GDEM data (see [Data](#data)).

### Figures

| Figure | Script | Description |
|--------|--------|-------------|
| Fig. 1 | `figures/fig1_geometry.py` | VTP wavefront propagation geometry |
| Fig. 2 | `figures/fig2_flowchart.py` | Algorithm comparison (traditional vs. multi-source) |
| Fig. 3 | `figures/fig3_performance.py` | VTP vs. Heat Method computation time |
| Fig. 4 | `figures/fig4_error.py` | Heat Method approximation error |
| Fig. 5 | `figures/fig5_speedup.py` | Multi-source speedup curve |
| Fig. 6 | `figures/fig6_scatter.py` | Geodesic vs. Euclidean distance scatter |
| Fig. 7 | `figures/fig7_histogram.py` | Distance ratio distribution |
| Fig. 8 | `figures/fig8_terrain_path.py` | 3D terrain path with water avoidance |

Generated figures are saved to `output/`.

---

## Data

The real terrain experiments (Tables 4–5, Figs. 6–8) use **ASTER GDEM v2** tile `ASTGTM2_N29E116_dem.tif`:

1. Register at [NASA Earthdata](https://urs.earthdata.nasa.gov/) (free)
2. Search for `ASTGTM2_N29E116` at [Earthdata Search](https://search.earthdata.nasa.gov/)
3. Download the `.tif` file
4. Place it in `data/ASTER_GDEM/` (rename `ASTGTM2_N29E116_dem.tif`)

---

## Directory Structure

```
multisource-vtp/
├── setup.py                    # Build script
├── vtp_geodesic.pyx            # Cython wrapper
├── vtp_src/                    # C++ VTP source
│   ├── geodesic_algorithm_exact.h
│   ├── geodesic_mesh.h
│   └── ...
├── experiments/                # Experiment scripts
│   ├── table2_performance.py
│   ├── table3_speedup.py
│   ├── table4_terrain_ratio.py
│   └── table5_path_planning.py
├── figures/                    # Figure generation scripts
│   ├── fig1_geometry.py
│   ├── fig2_flowchart.py
│   ├── fig3_performance.py
│   ├── fig4_error.py
│   ├── fig5_speedup.py
│   ├── fig6_scatter.py
│   ├── fig7_histogram.py
│   └── fig8_terrain_path.py
├── data/                       # Place downloaded DEM here
│   └── ASTER_GDEM/
├── output/                     # Generated figures
├── requirements.txt
└── README.md
```

---

## Citation

```bibtex
@article{luo2026multisource,
  author  = {罗昊 and 段新桥},
  title   = {多源VTP精确测地算法的精度与效率分析},
  journal = {测绘科学},
  year    = {2026},
}
```

## License

MIT
