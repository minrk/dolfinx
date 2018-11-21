// Copyright (C) 2018 Chris N. Richardson
//
// This file is part of DOLFIN (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include <dolfin/mesh/CellType.h>
#include <vector>

#pragma once

namespace dolfin
{

namespace fem
{

/// Data needed to lay out dofs on a cell

class CellDofLayout
{
public:
  CellDofLayout(mesh::CellType::Type cell_type,
                const std::vector<std::vector<std::vector<int>>>& entity_dofs)
  {

    if (cell_type == mesh::CellType::Type::tetrahedron)
    {
      // nfacet_dofs should be a triangular number
      int nfacet_dofs = entity_dofs[2][0].size();
      int n = (std::sqrt(1 + 8 * nfacet_dofs) - 1) / 2;
      if (n * (n + 1) != 2 * nfacet_dofs)
        throw std::runtime_error("Tetrahedron facet dofs not triangular");

      unsigned int c = 0;
      _facet_dof_coords.resize(n);
      for (unsigned int j = 0; j < n; ++j)
      {
        for (unsigned int i = 0; i < n - j; ++i)
        {
          _facet_dof_coords[j].push_back(c);
          ++c;
        }
      }
      assert(c == nfacet_dofs);
    }
    else if (cell_type == mesh::CellType::Type::hexahedron)
    {
      // nfacet_dofs should be a square number
      int nfacet_dofs = entity_dofs[2][0].size();
      int n = std::sqrt(nfacet_dofs);
      if (n * n != nfacet_dofs)
        throw std::runtime_error("Hexahedron facet dofs not square");

      unsigned int c = 0;
      _facet_dof_coords.resize(n);
      for (unsigned int j = 0; j < n; ++j)
      {
        for (unsigned int i = 0; i < n; ++i)
        {
          _facet_dof_coords[j].push_back(c);
          ++c;
        }
      }
      assert(c == nfacet_dofs);
    }
  }

  int facet_dof_coords(int i, int j) { return _facet_dof_coords[i][j]; }

  /// Move constructor
  CellDofLayout(CellDofLayout&& layout) = default;

  /// Destructor
  ~CellDofLayout() = default;

private:
  std::vector<std::vector<unsigned int>> _facet_dof_coords;
};
} // namespace fem
} // namespace dolfin
