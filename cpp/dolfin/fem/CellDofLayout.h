// Copyright (C) 2018 Chris N. Richardson
//
// This file is part of DOLFIN (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include <dolfin/mesh/CellType.h>
#include <string>
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
  /// Constructor
  CellDofLayout(mesh::CellType::Type cell_type,
                const std::vector<std::vector<std::vector<int>>>& entity_dofs);

  /// Move constructor
  CellDofLayout(CellDofLayout&& layout) = default;

  /// Destructor
  ~CellDofLayout() = default;

  /// Get the permutation of the dofs on this cell, in order to match with a
  /// global reference orientation for a given set of input global vertex
  /// indices. Only remaps the dofs which are on edges and facets.
  /// @param perm
  ///      Permutation vector required
  /// @param vertex_indices
  ///      Global vertex indices of the cell
  void permutation(std::vector<int>& perm, const int64_t* vertex_indices);

private:
  // Cell type
  mesh::CellType::Type _cell_type;

  // Layout of dofs by entity, ordered by [dim][entity][dof], e.g. for a
  // tetrahedron _entity_dofs[2][1][9] is the tenth dof of facet 1.
  std::vector<std::vector<std::vector<int>>> _entity_dofs;
};
} // namespace fem
} // namespace dolfin
