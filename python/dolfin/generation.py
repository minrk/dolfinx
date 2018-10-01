# -*- coding: utf-8 -*-
# Copyright (C) 2017 Chris N. Richardson
#
# This file is part of DOLFIN (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

"""Simple mesh generation module"""

from dolfin import fem
from dolfin import cpp
import numpy as np

__all__ = [
    "UnitIntervalMesh",
    "UnitSquareMesh",
    "UnitCubeMesh"
]

# FIXME: Remove, and use 'create' method?


def UnitIntervalMesh(comm, nx,
                     ghost_mode=cpp.mesh.GhostMode.none):
    """Create a mesh on the unit interval"""
    mesh = cpp.generation.IntervalMesh.create(comm, nx, [0.0, 1.0], ghost_mode)
    mesh.geometry.coord_mapping = fem.create_coordinate_map(mesh)
    return mesh


def UnitSquareMesh(comm, nx, ny, cell_type=cpp.mesh.CellType.Type.triangle,
                   ghost_mode=cpp.mesh.GhostMode.none, diagonal="right"):
    """Create a mesh of a unit square"""
    mesh = cpp.generation.RectangleMesh.create(
        comm, [np.array((0.0, 0.0, 0.0), dtype=np.float64), np.array((1.0, 1.0, 0.0), dtype=np.float64)],
        [nx, ny], cell_type, ghost_mode, diagonal)
    mesh.geometry.coord_mapping = fem.create_coordinate_map(mesh)
    return mesh


def UnitCubeMesh(comm, nx, ny, nz, cell_type=cpp.mesh.CellType.Type.tetrahedron,
                 ghost_mode=cpp.mesh.GhostMode.none):
    """Create a mesh of a unit cube"""
    from dolfin.geometry import Point
    mesh = cpp.generation.BoxMesh.create(
        comm, [np.array((0.0, 0.0, 0.0), dtype=np.float64), np.array((1.0, 1.0, 1.0), dtype=np.float64)],
        [nx, ny, nz], cell_type, ghost_mode)
    mesh.geometry.coord_mapping = fem.create_coordinate_map(mesh)
    return mesh
