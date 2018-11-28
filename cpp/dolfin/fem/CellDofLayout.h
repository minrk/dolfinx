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

  void permutation(std::vector<int>& perm, const int64_t* vertex_indices);

private:
  mesh::CellType::Type _cell_type;
  std::vector<std::vector<std::vector<int>>> _entity_dofs;
  std::vector<std::vector<unsigned int>> _facet_dof_coords;
};
} // namespace fem
} // namespace dolfin
