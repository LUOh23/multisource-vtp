import sys
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

ext_modules = [
    Extension(
        name="vtp_geodesic",
        sources=[
            "vtp_geodesic.pyx", 
        ],
        include_dirs=[
            np.get_include(),
            "vtp_src" 
        ],
        language="c++",
        extra_compile_args=["/std:c++17"] if sys.platform == 'win32' else ["-std=c++17"],
    )
]

setup(
    name="vtp_geodesic",
    version="0.1.0",
    description="Python wrapper for VTP Geodesic library",
    ext_modules=cythonize(
        ext_modules, 
        compiler_directives={'language_level': "3"}
    ),
    install_requires=['numpy'],
    zip_safe=False,
)