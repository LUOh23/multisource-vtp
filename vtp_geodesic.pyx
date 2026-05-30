# -*- coding: utf-8 -*-
#!python
#cython: language_level=3

import numpy
cimport numpy as np
from libcpp.vector cimport vector

# 初始化NumPy的C-API
np.import_array()
# -----------------------------------------------------------------------------
# C++ Class Declarations
# -----------------------------------------------------------------------------

cdef extern from "geodesic_mesh_elements.h" namespace "geodesic":
    cdef cppclass Vertex:
        Vertex()
        double& geodesic_distance()

cdef extern from "geodesic_mesh_elements.h" namespace "geodesic":
    cdef cppclass Face:
        Face()

cdef extern from "geodesic_mesh.h" namespace "geodesic":
    cdef cppclass Mesh:
        Mesh()
        void initialize_mesh_data(vector[double]&, vector[unsigned]&) nogil
        #
        vector[Vertex]& vertices() nogil
        vector[Face]& faces() nogil
        size_t vertices_size() nogil

cdef extern from "geodesic_algorithm_exact.h" namespace "geodesic":
    cdef cppclass GeodesicAlgorithmExact:
        GeodesicAlgorithmExact(Mesh*)
        void propagate(const vector[unsigned]&) nogil

cdef extern from "geodesic_constants_and_simple_functions.h" namespace "geodesic":
    double GEODESIC_INF

# -----------------------------------------------------------------------------
# Python Wrapper Class
# -----------------------------------------------------------------------------

cdef class PyVTPAlgorithmExact:
    cdef GeodesicAlgorithmExact *algorithm
    cdef Mesh mesh
    cdef vector[double] points
    cdef vector[unsigned] faces
    cdef public int num_vertices

    def __cinit__(self, _points, _faces):
        try:
            _points = numpy.asarray(_points, dtype=numpy.float64)
            _faces = numpy.asarray(_faces, dtype=numpy.int32)
            
            if len(_points.shape) != 2 or _points.shape[-1] != 3:
                raise ValueError("'points' must be an Nx3 array")
            if len(_faces.shape) != 2 or _faces.shape[-1] != 3:
                raise ValueError("'faces' must be an Mx3 array")
            if _faces.min() < 0 or _faces.max() >= _points.shape[0]:
                raise ValueError("Face indices out of bounds")
                
            self.num_vertices = _points.shape[0]
        except Exception as e:
            print(f'Error in initialization: {e}')
            self.algorithm = NULL
            self.num_vertices = 0
            return

        cdef double coord
        for coord in _points.flatten():
            self.points.push_back(coord)

        cdef unsigned indx
        for indx in _faces.flatten():
            self.faces.push_back(indx)

        self.mesh.initialize_mesh_data(self.points, self.faces)
        self.algorithm = new GeodesicAlgorithmExact(&self.mesh)

    def compute_distances(self, source_indices=None):
        """
        Compute geodesic distances from sources to all vertices.
        Returns a numpy array of distances.
        """
        cdef vector[unsigned] src_vec
        cdef np.ndarray[np.float64_t, ndim=1] distances
        cdef Py_ssize_t i
        cdef double dist_val
        
        if self.algorithm == NULL:
            raise RuntimeError("Algorithm not initialized.")

        # 1. Prepare Sources
        if source_indices is None:
            source_indices = [0]
        
        try:
            src_array = numpy.asarray(source_indices, dtype=numpy.int32)
            if len(src_array.shape) != 1:
                raise ValueError("source_indices must be 1D")
            if src_array.min() < 0 or src_array.max() >= self.num_vertices:
                raise ValueError("Source indices out of bounds")
            
            for i in range(len(src_array)):
                src_vec.push_back(<unsigned>src_array[i])
        except Exception as e:
            raise ValueError(f"Invalid source_indices: {e}")

        if src_vec.empty():
            raise ValueError("Source list cannot be empty")

        # 2. Run Propagation
        with nogil:
            self.algorithm.propagate(src_vec)

        # 3. Extract Results
        distances = numpy.empty(self.num_vertices, dtype=numpy.float64)
        cdef double[:] dist_view = distances
        
        #cdef vector[Vertex]& verts = self.mesh.vertices()
        cdef vector[Vertex]* verts_ptr = &self.mesh.vertices() 
        
        for i in range(self.num_vertices):
            #dist_val = verts[i].geodesic_distance()
            dist_val = self.mesh.vertices()[i].geodesic_distance()
            
            if dist_val >= GEODESIC_INF:
                dist_view[i] = numpy.inf
            else:
                dist_view[i] = dist_val

        return distances

    def __dealloc__(self):
        if self.algorithm != NULL:
            del self.algorithm