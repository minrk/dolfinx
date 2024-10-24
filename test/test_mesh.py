# Create a mesh consisting of points only
# Author: Jørgen S. Dokken
# SPDX-License-Identifier: MIT

from mpi4py import MPI
import basix.ufl

import dolfinx
import numpy as np
import numpy.typing as npt
import ufl
from dolfinx import log

log.set_log_level(log.LogLevel.DEBUG)


def create_point_mesh(comm: MPI.Intracomm, points: npt.NDArray[np.float32] | npt.NDArray[np.float64]) -> dolfinx.mesh.Mesh:
    """
    Create a mesh consisting of points only.

    Note:
        No nodes are shared between processes.

    Args:
        comm: MPI communicator to create the mesh on.
        points: Points local to the process in the mesh.
    """
    print("start")
    # Create mesh topology
    cells = np.arange(points.shape[0], dtype=np.int32).reshape(-1, 1)
    print("topology")
    topology = dolfinx.cpp.mesh.Topology(
        comm, dolfinx.mesh.CellType.point)
    num_nodes_local = cells.shape[0]
    print("imap")
    imap = dolfinx.common.IndexMap(comm, num_nodes_local)
    local_range = imap.local_range[0]
    igi = np.arange(num_nodes_local, dtype=np.int64)+local_range
    print("set_index")
    topology.set_index_map(0, imap)
    print("set_conn")
    topology.set_connectivity(dolfinx.graph.adjacencylist(cells), 0, 0)

    # Create mesh geometry
    e = basix.ufl.element("Lagrange", "point", 0, shape=(points.shape[1],))
    c_el = dolfinx.fem.coordinate_element(e.basix_element)
    print("create_geometry")
    geometry = dolfinx.mesh.create_geometry(
        imap, cells, c_el._cpp_object, points,  igi)

    # Create DOLFINx mesh
    print("mesh", points.dtype)
    if points.dtype == np.float64:
        cpp_mesh = dolfinx.cpp.mesh.Mesh_float64(
            comm, topology, geometry._cpp_object)
    elif points.dtype == np.float32:
        cpp_mesh = dolfinx.cpp.mesh.Mesh_float32(
            comm, topology, geometry._cpp_object)
    else:
        raise RuntimeError(f"Unsupported dtype for mesh {points.dtype}")
    # Wrap as Python object
    print("mesh.Mesh")
    return dolfinx.mesh.Mesh(cpp_mesh, domain=ufl.Mesh(e))


if __name__ == "__main__":
    p_m = create_point_mesh(MPI.COMM_WORLD, np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]))
    print(p_m.geometry.x)
    print(p_m.topology.index_map(0).size_local)
