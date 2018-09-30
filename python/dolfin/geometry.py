# -*- coding: utf-8 -*-
# Copyright (C) 2018 Michal Habera
#
# This file is part of DOLFIN (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

from dolfin import cpp


class BoundingBoxTree:
    def __init__(self, gdim=None):
        """Initialise from geometric dimension"""
        if gdim is not None:
            self._cpp_object = cpp.geometry.BoundingBoxTree(gdim)

    def build_points(self, points: list):
        """Build from cloud of points"""
        # Unpack to cpp points
        points_cpp = (point._cpp_object for point in points)
        self._cpp_object.build(points_cpp)

    def build_mesh(self, mesh, tdim: int):
        """Build from mesh entities of given topological dimension"""
        self._cpp_object.build(mesh, tdim)

    def compute_collisions_point(self, point: "Point"):
        """Compute collisions with the point"""
        return self._cpp_object.compute_collisions(point._cpp_object)

    def compute_collisions_bb(self, bb: "BoundingBoxTree"):
        """Compute collisions with the bounding box"""
        return self._cpp_object.compute_collisions(bb._cpp_object)

    def compute_entity_collisions_mesh(self, point: "Point", mesh):
        """Compute collisions between the point and entities of the mesh"""
        return self._cpp_object.compute_entity_collisions(point._cpp_object,
                                                          mesh)

    def compute_entity_collisions_bb_mesh(self, bb: "BoundingBoxTree",
                                          mesh1, mesh2):
        """Compute collisions between the bounding box and entities of meshes"""
        return self._cpp_object.compute_entity_collisions(bb._cpp_object,
                                                          mesh1, mesh2)

    def compute_first_collision(self, point: "Point"):
        """Compute first collision with the point"""
        return self._cpp_object.compute_first_collision(point._cpp_object)

    def compute_first_entity_collision(self, point: "Point", mesh):
        """Compute fist collision between entities of mesh and the point"""
        return self._cpp_object.compute_first_entity_collision(point._cpp_object, mesh)

    def compute_closest_entity(self, point: "Point", mesh):
        """Compute closest entity of the mesh to the point"""
        return self._cpp_object.compute_closest_entity(point._cpp_object, mesh)

    def str(self):
        """Print for debbuging"""
        return self._cpp_object.str()
