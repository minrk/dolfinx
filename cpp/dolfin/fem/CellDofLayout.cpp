// Copyright (C) 2018 Chris N. Richardson
//
// This file is part of DOLFIN (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include <iostream>

#include "CellDofLayout.h"
#include <dolfin/log/log.h>

using namespace dolfin;
using namespace dolfin::fem;

//-----------------------------------------------------------------------------
CellDofLayout::CellDofLayout(
    mesh::CellType::Type cell_type,
    const std::vector<std::vector<std::vector<int>>>& entity_dofs)
    : _cell_type(cell_type), _entity_dofs(entity_dofs)
{
}
//-----------------------------------------------------------------------------
// Helper function to get layout of dofs on a triangular facet
int tri_layout_ij(int i, int j, int n) { return i + ((2 * n + 1 - j) * j) / 2; }
//-----------------------------------------------------------------------------
void CellDofLayout::permutation(std::vector<int>& perm,
                                const int64_t* vertex_indices)
{
  // Reset to identity
  for (unsigned int i = 0; i < perm.size(); ++i)
    perm[i] = i;

  if (_cell_type == mesh::CellType::Type::tetrahedron)
  {
    // nfacet_dofs should be a triangular number
    unsigned int nfacet_dofs = _entity_dofs[2][0].size();
    unsigned int n = (std::sqrt(1 + 8 * nfacet_dofs) - 1) / 2;
    if (n * (n + 1) != 2 * nfacet_dofs)
    {
      log::warning("Tetrahedron facet dofs not triangular: "
                   + std::to_string(nfacet_dofs))
          + " - not permuting...";
      return;
    }

    // Get ordering of edges
    bool edge_ordering[6];
    edge_ordering[0] = vertex_indices[2] > vertex_indices[3];
    edge_ordering[1] = vertex_indices[1] > vertex_indices[3];
    edge_ordering[2] = vertex_indices[1] > vertex_indices[2];
    edge_ordering[3] = vertex_indices[0] > vertex_indices[3];
    edge_ordering[4] = vertex_indices[0] > vertex_indices[2];
    edge_ordering[5] = vertex_indices[0] > vertex_indices[1];

    for (unsigned int j = 0; j < 6; ++j)
    {
      if (edge_ordering[j])
      {
        // Reverse dofs along this edge
        const std::vector<int>& edge_dofs = _entity_dofs[1][j];
        unsigned int ne = edge_dofs.size();
        for (unsigned int i = 0; i < ne; ++i)
          perm[edge_dofs[i]] = edge_dofs[ne - i - 1];
      }
    }

    // Edges on each facet
    static unsigned int facet_edges[4][3]
        = {{0, 1, 2}, {0, 3, 4}, {1, 3, 5}, {2, 4, 5}};

    // Facet ordering
    for (unsigned int m = 0; m < 4; ++m)
    {
      const std::vector<int>& facet_dofs = _entity_dofs[2][m];
      unsigned int c = 0;

      if (facet_dofs.size() > 1)
      {
        const unsigned int* fe = facet_edges[m];
        int facet_ordering
            = edge_ordering[fe[0]]
              + 2 * (edge_ordering[fe[1]] + edge_ordering[fe[2]]);

        std::cout << "Ordering = " << facet_ordering << "\n";

        // Do stuff based on value of facet_ordering (0-5)
        switch (facet_ordering)
        {
        case 0:
          break;
        case 1:
          for (unsigned int j = 0; j < n; ++j)
          {
            for (unsigned int i = 0; i < n - j; ++i)
            {
              perm[facet_dofs[c]] = facet_dofs[tri_layout_ij(j, i, n)];
              ++c;
            }
          }
          break;
        case 2:
          for (unsigned int j = 0; j < n; ++j)
          {
            for (unsigned int i = 0; i < n - j; ++i)
            {
              unsigned int k = n - i - j - 1;
              perm[facet_dofs[c]] = facet_dofs[tri_layout_ij(k, j, n)];
              ++c;
            }
          }
          break;
        case 3:
          for (unsigned int j = 0; j < n; ++j)
          {
            for (unsigned int i = 0; i < n - j; ++i)
            {
              unsigned int k = n - i - j - 1;
              perm[facet_dofs[c]] = facet_dofs[tri_layout_ij(k, i, n)];
              ++c;
            }
          }
          break;
        case 4:
          for (unsigned int j = 0; j < n; ++j)
          {
            for (unsigned int i = 0; i < n - j; ++i)
            {
              unsigned int k = n - i - j - 1;
              perm[facet_dofs[c]] = facet_dofs[tri_layout_ij(j, k, n)];
              ++c;
            }
          }
          break;
        case 5:
          for (unsigned int j = 0; j < n; ++j)
          {
            for (unsigned int i = 0; i < n - j; ++i)
            {
              unsigned int k = n - i - j - 1;
              perm[facet_dofs[c]] = facet_dofs[tri_layout_ij(i, k, n)];
              ++c;
            }
          }
          break;
        }
      }
    }
  }
  else if (_cell_type == mesh::CellType::Type::triangle)
  {
    // Get ordering of edges
    bool edge_ordering[3];
    edge_ordering[0] = vertex_indices[1] > vertex_indices[2];
    edge_ordering[1] = vertex_indices[0] > vertex_indices[2];
    edge_ordering[2] = vertex_indices[0] > vertex_indices[1];

    for (unsigned int j = 0; j < 3; ++j)
    {
      if (edge_ordering[j])
      {
        // Reverse dofs along this edge
        const std::vector<int>& edge_dofs = _entity_dofs[1][j];
        const unsigned int ne = edge_dofs.size();
        for (unsigned int i = 0; i < ne; ++i)
          perm[edge_dofs[i]] = edge_dofs[ne - i - 1];
      }
    }
  }

  // Debug printout
  for (auto& q : perm)
    std::cout << q << " ";
  std::cout << "\n";
}
//-----------------------------------------------------------------------------
