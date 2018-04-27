// Copyright (C) 2018 Garth N. Wells
//
// This file is part of DOLFIN (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#include "Assembler.h"
#include "DirichletBC.h"
#include "Form.h"
#include "GenericDofMap.h"
#include "SparsityPatternBuilder.h"
#include "UFC.h"
#include "utils.h"
#include <dolfin/common/types.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/Scalar.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshIterator.h>
#include <string>

using namespace dolfin;
using namespace dolfin::fem;

//-----------------------------------------------------------------------------
Assembler::Assembler(std::vector<std::vector<std::shared_ptr<const Form>>> a,
                     std::vector<std::shared_ptr<const Form>> L,
                     std::vector<std::shared_ptr<const DirichletBC>> bcs)
    : _a(a), _l(L), _bcs(bcs)
{
  assert(!a.empty());
  assert(!a[0].empty());

  // FIXME: check that a is square

  // Check shape of a and L

  // Check rank of forms
  /*
  if (a and a->rank() != 2)
  {
    throw std::runtime_error(
        "Expecting bilinear form (rank 2), but form has rank \'"
        + std::to_string(a->rank()) + "\'");
  }
  if (L and L->rank() != 1)
  {
    throw std::runtime_error(
        "Expecting linear form (rank 1), but form has rank \'"
        + std::to_string(L->rank()) + "\'");
  }
  */
}
//-----------------------------------------------------------------------------
void Assembler::assemble(la::PETScMatrix& A, BlockType block_type)
{
  // Check if matrix should be nested
  assert(!_a.empty());
  const bool block_matrix = _a.size() > 1 or _a[0].size() > 1;

  if (A.empty())
  {
    // Initialise matrix if empty

    // Build array of pointers to forms
    std::vector<std::vector<const Form*>> forms(
        _a.size(), std::vector<const Form*>(_a[0].size()));
    for (std::size_t i = 0; i < _a.size(); ++i)
      for (std::size_t j = 0; j < _a[i].size(); ++j)
        forms[i][j] = _a[i][j].get();

    // Initialise matrix
    if (block_type == BlockType::nested)
    {
      fem::init_nest(A, forms);
    }
    else if (block_matrix and block_type == BlockType::monolithic)
    {
      // Initialise matrix
      fem::init_monolithic(A, forms);

      // Create local-to-global maps and attach to matrix
      std::vector<PetscInt> _map0, _map1;

      // Build list of index maps
      std::vector<const common::IndexMap*> maps0, maps1;
      for (std::size_t i = 0; i < _a.size(); ++i)
      {
        auto map = _a[i][0]->function_space(0)->dofmap()->index_map();
        maps0.push_back(map.get());
      }
      for (std::size_t i = 0; i < _a[0].size(); ++i)
      {
        auto map = _a[0][i]->function_space(1)->dofmap()->index_map();
        maps1.push_back(map.get());
      }

      // Row local-to-global map
      for (std::size_t i = 0; i < _a.size(); ++i)
      {
        auto map = _a[i][0]->function_space(0)->dofmap()->index_map();
        for (std::size_t k = 0; k < map->size(common::IndexMap::MapSize::ALL);
             ++k)
        {
          auto index_k = map->local_to_global(k);
          std::size_t index = get_global_index(maps0, i, index_k);
          _map0.push_back(index);
        }
      }

      // Column local-to-global map
      for (std::size_t i = 0; i < _a[0].size(); ++i)
      {
        auto map = _a[0][i]->function_space(1)->dofmap()->index_map();
        for (std::size_t k = 0; k < map->size(common::IndexMap::MapSize::ALL);
             ++k)
        {
          auto index_k = map->local_to_global(k);
          std::size_t index = get_global_index(maps1, i, index_k);
          _map1.push_back(index);
        }
      }

      // Create PETSc local-to-global map/index sets and attach to matrix
      ISLocalToGlobalMapping petsc_local_to_global0, petsc_local_to_global1;
      ISLocalToGlobalMappingCreate(MPI_COMM_SELF, 1, _map0.size(), _map0.data(),
                                   PETSC_COPY_VALUES, &petsc_local_to_global0);
      ISLocalToGlobalMappingCreate(MPI_COMM_SELF, 1, _map1.size(), _map1.data(),
                                   PETSC_COPY_VALUES, &petsc_local_to_global1);
      MatSetLocalToGlobalMapping(A.mat(), petsc_local_to_global0,
                                 petsc_local_to_global1);

      // Clean up local-to-global maps
      ISLocalToGlobalMappingDestroy(&petsc_local_to_global0);
      ISLocalToGlobalMappingDestroy(&petsc_local_to_global1);
    }
    else
    {
      init(A, *_a[0][0]);
    }
  }

  // Get matrix type
  MatType mat_type;
  MatGetType(A.mat(), &mat_type);
  const bool is_matnest = strcmp(mat_type, MATNEST) == 0 ? true : false;

  // Assemble matrix

  if (is_matnest)
  {
    for (std::size_t i = 0; i < _a.size(); ++i)
    {
      for (std::size_t j = 0; j < _a[i].size(); ++j)
      {
        // Get submatrix
        Mat subA;
        MatNestGetSubMat(A.mat(), i, j, &subA);
        if (_a[i][j])
        {
          la::PETScMatrix mat(subA);
          this->assemble(mat, *_a[i][j], _bcs);
        }
        else
        {
          // FIXME: Figure out how to check that matrix block is null
          // Null block, do nothing
        }
      }
    }
  }
  else if (block_matrix)
  {
    std::int64_t offset_row = 0;
    for (std::size_t i = 0; i < _a.size(); ++i)
    {
      // Loop over columns
      std::int64_t offset_col = 0;
      for (std::size_t j = 0; j < _a[i].size(); ++j)
      {
        if (_a[i][j])
        {
          // Build index set for block
          auto map0 = _a[i][j]->function_space(0)->dofmap()->index_map();
          auto map1 = _a[i][j]->function_space(1)->dofmap()->index_map();
          auto map0_size = map0->size(common::IndexMap::MapSize::ALL);
          auto map1_size = map1->size(common::IndexMap::MapSize::ALL);

          std::vector<PetscInt> index0(map0_size);
          std::vector<PetscInt> index1(map1_size);
          std::iota(index0.begin(), index0.end(), offset_row);
          std::iota(index1.begin(), index1.end(), offset_col);

          IS is0, is1;
          ISCreateBlock(MPI_COMM_SELF, map0->block_size(), index0.size(),
                        index0.data(), PETSC_COPY_VALUES, &is0);
          ISCreateBlock(MPI_COMM_SELF, map1->block_size(), index1.size(),
                        index1.data(), PETSC_COPY_VALUES, &is1);

          // Get sub-matrix
          Mat subA;
          MatGetLocalSubMatrix(A.mat(), is0, is1, &subA);

          // Assemble block
          la::PETScMatrix mat(subA);
          this->assemble(mat, *_a[i][j], _bcs);

          // Restore sub-matrix and destroy index sets
          MatRestoreLocalSubMatrix(A.mat(), is0, is1, &subA);
          ISDestroy(&is0);
          ISDestroy(&is1);

          offset_col += map1_size;
        }
        else
        {
          // FIXME: Figure out how to check that matrix block is null
          // Null block, do nothing
        }
      }
      auto map0 = _a[i][0]->function_space(0)->dofmap()->index_map();
      auto map0_size = map0->size(common::IndexMap::MapSize::ALL);
      offset_row += map0_size;
    }
  }
  else
  {
    this->assemble(A, *_a[0][0], _bcs);
  }

  A.apply(la::PETScMatrix::AssemblyType::FINAL);
}
//-----------------------------------------------------------------------------
void Assembler::assemble(la::PETScVector& b, BlockType block_type)
{
  // Check if matrix should be nested
  assert(!_l.empty());
  const bool block_vector = _l.size() > 1;

  if (b.empty())
  {
    // Initialise matrix if empty

    // Build array of pointers to forms
    std::vector<const Form*> forms(_a.size());
    for (std::size_t i = 0; i < _l.size(); ++i)
      forms[i] = _l[i].get();

    // Initialise vector
    if (block_type == BlockType::nested)
      fem::init_nest(b, forms);
    else if (block_vector and block_type == BlockType::monolithic)
      fem::init_monolithic(b, forms);
    else
      init(b, *_l[0]);
  }

  // Get vector type
  VecType vec_type;
  VecGetType(b.vec(), &vec_type);
  const bool is_vecnest = strcmp(vec_type, VECNEST) == 0;

  if (is_vecnest)
  {
    for (std::size_t i = 0; i < _l.size(); ++i)
    {
      // Get subvector
      Vec sub_b;
      VecNestGetSubVec(b.vec(), i, &sub_b);
      if (_l[i])
      {
        la::PETScVector vec(sub_b);
        this->assemble(vec, *_l[i]);
      }
      else
      {
        // FIXME: Figure out how to check that vector block is null
        // Null block, do nothing
      }
    }
  }
  else if (block_vector)
  {
    std::int64_t offset = 0;
    for (std::size_t i = 0; i < _l.size(); ++i)
    {
      if (_l[i])
      {
        auto map = _l[i]->function_space(0)->dofmap()->index_map();
        auto map_size = map->size(common::IndexMap::MapSize::ALL);

        std::vector<PetscInt> index(map_size);
        std::iota(index.begin(), index.end(), offset);

        IS is;
        ISCreateBlock(b.mpi_comm(), map->block_size(), index.size(),
                      index.data(), PETSC_COPY_VALUES, &is);

        Vec sub_b;
        VecGetSubVector(b.vec(), is, &sub_b);
        la::PETScVector vec(sub_b);

        // FIXME: Does it pick up the block size?

        // FIXME: Update for parallel
        // Attach local-to-global map

        // Fill vector with [i0 + 0, i0 + 1, i0 +2, . . .]
        std::vector<PetscInt> local_to_global_map(vec.size());
        std::iota(local_to_global_map.begin(), local_to_global_map.end(), 1);

        // Create PETSc local-to-global map
        ISLocalToGlobalMapping petsc_local_to_global;
        // PetscErrorCode ierr = 0;
        ISLocalToGlobalMappingCreate(PETSC_COMM_SELF, 1,
                                     local_to_global_map.size(),
                                     local_to_global_map.data(),
                                     PETSC_COPY_VALUES, &petsc_local_to_global);
        // CHECK_ERROR("ISLocalToGlobalMappingCreate");

        // Apply local-to-global map to vector
        VecSetLocalToGlobalMapping(sub_b, petsc_local_to_global);
        // CHECK_ERROR("VecSetLocalToGlobalMapping");

        // Clean-up PETSc local-to-global map
        ISLocalToGlobalMappingDestroy(&petsc_local_to_global);
        // CHECK_ERROR("ISLocalToGlobalMappingDestroy");

        this->assemble(vec, *_l[i]);

        VecRestoreSubVector(b.vec(), is, &sub_b);
        ISDestroy(&is);

        offset += map_size;
      }
    }
  }
  else
  {
    this->assemble(b, *_l[0]);
  }

  /*
  // Assemble vector
  this->assemble(b, *_l[0]);

  // Apply bcs to RHS of vector
  for (std::size_t i = 0; i < _l.size(); ++i)
    for (std::size_t j = 0; j < _a[i].size(); ++j)
      apply_bc(b, *_a[i][j], _bcs);

  // Set bc values
  set_bc(b, *_l[0], _bcs);
  `*/

  // // Assemble blocks (b)
  // for (auto row : _l)
  // {
  //   this->assemble(b, *row);
  // }
}
//-----------------------------------------------------------------------------
void Assembler::assemble(la::PETScMatrix& A, la::PETScVector& b)
{
  // TODO: pre common boundary condition data

  // Assemble matrix
  assemble(A);

  // Assemble vector
  assemble(b);
}
//-----------------------------------------------------------------------------
void Assembler::assemble(la::PETScMatrix& A, const Form& a,
                         std::vector<std::shared_ptr<const DirichletBC>> bcs)
{
  assert(!A.empty());

  // Get mesh from form
  assert(a.mesh());
  const mesh::Mesh& mesh = *a.mesh();

  // FIXME: Remove UFC
  // Create data structures for local assembly data
  UFC ufc(a);

  const std::size_t gdim = mesh.geometry().dim();
  const std::size_t tdim = mesh.topology().dim();
  mesh.init(tdim);

  // Function spaces for each axis
  std::array<const function::FunctionSpace*, 2> spaces
      = {{a.function_space(0).get(), a.function_space(1).get()}};

  // Collect pointers to dof maps
  std::array<const GenericDofMap*, 2> dofmaps
      = {{spaces[0]->dofmap().get(), spaces[1]->dofmap().get()}};

  // FIXME: Move out of this function
  // FIXME: For the matrix, we only need to know if there is a boundary
  // condition on the entry. The value is not required.
  // FIXME: Avoid duplication when spaces[0] == spaces[1]
  // Collect boundary conditions by matrix axis
  std::array<DirichletBC::Map, 2> boundary_values;
  for (std::size_t i = 0; i < bcs.size(); ++i)
  {
    assert(bcs[i]);
    assert(bcs[i]->function_space());
    for (std::size_t axis = 0; axis < 2; ++axis)
    {
      if (spaces[axis]->contains(*bcs[i]->function_space()))
      {
        // FIXME: find way to avoid gather, or perform with a single
        // gather
        bcs[i]->get_boundary_values(boundary_values[axis]);
        if (MPI::size(mesh.mpi_comm()) > 1
            and bcs[i]->method() != DirichletBC::Method::pointwise)
        {
          bcs[i]->gather(boundary_values[axis]);
        }
      }
    }
  }

  // FIXME: Better way of applying bcs to local matrix?
  const auto apply_bc_to_local_matrix =
          [](dolfin::EigenRowMatrixXd& Ae,
             const std::array<DirichletBC::Map, 2>& boundary_values,
             const Eigen::Map<const Eigen::Array<dolfin::la_index_t, Eigen::Dynamic, 1>>& dmap0,
             const Eigen::Map<const Eigen::Array<dolfin::la_index_t, Eigen::Dynamic, 1>>& dmap1) -> void {
    // FIXME: Pass in list  of cells, and list of local dofs, with
    // Dirichlet conditions
    // Note: could use zero dof indices to have PETSc do this
    // Zero rows/columns for Dirichlet bcs
    for (int i = 0; i < Ae.rows(); ++i)
    {
      const std::size_t ii = dmap0[i];
      DirichletBC::Map::const_iterator bc_value = boundary_values[0].find(ii);
      if (bc_value != boundary_values[0].end())
        Ae.row(i).setZero();
    }
    // Loop over columns
    for (int j = 0; j < Ae.cols(); ++j)
    {
      const std::size_t jj = dmap1[j];
      DirichletBC::Map::const_iterator bc_value = boundary_values[1].find(jj);
      if (bc_value != boundary_values[1].end())
        Ae.col(j).setZero();
    }
  };

  // Data structures used in assembly
  EigenRowArrayXXd coordinate_dofs;
  EigenRowMatrixXd Ae;

  // Check whether integral is domain-dependent
  auto cell_domains = a.cell_domains();
  bool use_cell_domains = cell_domains && cell_domains->size() > 0;

  if (a.integrals().num_cell_integrals() > 0)
  {
    // FIXME: This will not be valid if num_element_dofs can vary by cell in the future
    const auto init_element_vector = [&dofmaps](dolfin::EigenRowMatrixXd& Ae)
            -> void { Ae.resize(dofmaps[0]->num_element_dofs(0),
                                dofmaps[1]->num_element_dofs(0)); };

    const auto add_local_matrix_to_global = [&A, &dofmaps, &boundary_values, &apply_bc_to_local_matrix]
            (dolfin::EigenRowMatrixXd& Ae, dolfin::mesh::MeshEntity& cell)
            -> void {
      const auto dmap0 = dofmaps[0]->cell_dofs(cell.index());
      const auto dmap1 = dofmaps[1]->cell_dofs(cell.index());

      apply_bc_to_local_matrix(Ae, boundary_values, dmap0, dmap1);

      A.add_local(Ae.data(), dmap0.size(), dmap0.data(), dmap1.size(),
                  dmap1.data());
    };

    assemble_over_cells(a, init_element_vector, add_local_matrix_to_global);
  }


  if (a.integrals().num_exterior_facet_integrals() > 0)
  {
    // Exterior facet assembly
    mesh.init(tdim - 1);
    mesh.init(tdim - 1, tdim);

    // Get exterior facet integral
    auto exterior_facet_integral = a.integrals().exterior_facet_integral();

    // Check whether facet integral is domain dependent
    auto exterior_facet_domains = a.exterior_facet_domains();
    auto use_exterior_facet_domains = exterior_facet_domains && exterior_facet_domains->size() > 0;

    // Iterate over exterior facets
    for (const auto &facet : mesh::MeshRange<mesh::Facet>(mesh)) {
      if (!facet.exterior())
        continue;

      if (use_exterior_facet_domains)
        exterior_facet_integral = a.integrals().exterior_facet_integral((*exterior_facet_domains)[facet]);

      if (!exterior_facet_integral)
        continue;

      assert(facet.num_entities(tdim) == 1);
      const mesh::Cell mesh_cell(mesh, facet.entities(tdim)[0]);

      assert(!mesh_cell.is_ghost());
      const std::size_t local_facet = mesh_cell.index(facet);

      coordinate_dofs.resize(mesh_cell.num_vertices(), gdim);
      mesh_cell.get_coordinate_dofs(coordinate_dofs);

      ufc.update(mesh_cell, coordinate_dofs, exterior_facet_integral->enabled_coefficients);

      // Get dof maps for cell
      auto dmap0 = dofmaps[0]->cell_dofs(mesh_cell.index());
      auto dmap1 = dofmaps[1]->cell_dofs(mesh_cell.index());

      // Size data structure for assembly
      Ae.resize(dmap0.size(), dmap1.size());
      Ae.setZero();

      // TODO: the cell orientation
      exterior_facet_integral->tabulate_tensor(Ae.data(), ufc.w(), coordinate_dofs.data(), local_facet, 1);

      apply_bc_to_local_matrix(Ae, boundary_values, dmap0, dmap1);

      // Add to matrix
      A.add_local(Ae.data(), dmap0.size(), dmap0.data(), dmap1.size(),
                  dmap1.data());
    }
  }


  // Interior facet assembly
  if (a.integrals().num_interior_facet_integrals() > 0)
  {
    // Sanity check of ghost mode (proper check in AssemblerBase::check)
    assert(mesh.get_ghost_mode() == "shared_vertex"
           or mesh.get_ghost_mode() == "shared_facet"
           or MPI::size(mesh.mpi_comm()) == 1);

    mesh.init(tdim - 1);
    mesh.init(tdim - 1, tdim);

    auto interior_facet_integral = a.integrals().interior_facet_integral();

    // Check whether facet integral is domain dependent
    auto interior_facet_domains = a.interior_facet_domains();
    auto use_interior_facet_domains = interior_facet_domains && interior_facet_domains->size() > 0;

    std::vector<EigenRowArrayXXd, Eigen::aligned_allocator<EigenRowArrayXXd>> neighbour_coordinate_dofs(2);

    const std::size_t mpi_rank = MPI::rank(mesh.mpi_comm());

    for (const auto &facet : mesh::MeshRange<mesh::Facet>(mesh)) {
      if (facet.exterior())
        continue;

      assert(!facet.is_ghost());

      if (use_interior_facet_domains)
        interior_facet_integral = a.integrals().interior_facet_integral((*interior_facet_domains)[facet]);

      if (!interior_facet_integral)
        continue;

      assert(facet.num_entities(tdim) == 2);
      auto cell_index_plus = facet.entities(tdim)[0];
      auto cell_index_minus = facet.entities(tdim)[1];

      if (use_cell_domains and (*cell_domains)[cell_index_plus] < (*cell_domains)[cell_index_minus])
        std::swap(cell_index_plus, cell_index_minus);

      // The convention '+' = 0, '-' = 1 is from ffc
      const mesh::Cell cell_plus(mesh, cell_index_plus);
      const mesh::Cell cell_minus(mesh, cell_index_minus);

      // Get local index of facet with respect to each cell
      std::size_t local_facet0 = cell_plus.index(facet);
      std::size_t local_facet1 = cell_minus.index(facet);

      // Update to current pair of cells
      neighbour_coordinate_dofs[0].resize(cell_plus.num_vertices(), gdim);
      cell_plus.get_coordinate_dofs(neighbour_coordinate_dofs[0]);

      neighbour_coordinate_dofs[1].resize(cell_minus.num_vertices(), gdim);
      cell_minus.get_coordinate_dofs(neighbour_coordinate_dofs[1]);

      ufc.update(cell_plus, neighbour_coordinate_dofs[0],
                 cell_minus, neighbour_coordinate_dofs[1],
                 interior_facet_integral->enabled_coefficients);

      // Get dof maps for cells
      auto dmap0_plus = dofmaps[0]->cell_dofs(cell_plus.index());
      auto dmap0_minus = dofmaps[0]->cell_dofs(cell_minus.index());
      auto dmap1_plus = dofmaps[1]->cell_dofs(cell_plus.index());
      auto dmap1_minus = dofmaps[1]->cell_dofs(cell_minus.index());

      const std::size_t macro_dmap0_size = dmap0_plus.size() + dmap0_minus.size();
      const std::size_t macro_dmap1_size = dmap1_plus.size() + dmap1_minus.size();

      // Prepare for assembly -- Ae is now a macro local matrix
      Ae.resize(macro_dmap0_size, macro_dmap1_size);
      Ae.setZero();

      // TODO: sort out the orientation
      interior_facet_integral->tabulate_tensor(Ae.data(), ufc.macro_w(),
                                               neighbour_coordinate_dofs[0].data(),
                                               neighbour_coordinate_dofs[1].data(),
                                               local_facet0, local_facet1,
                                               1, 1);

      if (cell_plus.is_ghost() != cell_minus.is_ghost()) {
        const std::size_t ghost_rank = cell_plus.is_ghost() ? cell_plus.owner() : cell_minus.owner();
        assert(mpi_rank != ghost_rank);
        if (ghost_rank < mpi_rank)
          continue;
      }

      // FIXME: Can this be replaced with a Ref or Map somehow? This code is horrendous
      Eigen::Array<dolfin::la_index_t, Eigen::Dynamic, 1> macro_dofs0(macro_dmap0_size);
      macro_dofs0 << dmap0_plus, dmap0_minus;
      const Eigen::Map<const Eigen::Array<dolfin::la_index_t, Eigen::Dynamic, 1>>
              macro_dofs_map0(macro_dofs0.data(), macro_dmap0_size);
      Eigen::Array<dolfin::la_index_t, Eigen::Dynamic, 1> macro_dofs1(macro_dmap1_size);
      macro_dofs1 << dmap1_plus, dmap1_minus;
      const Eigen::Map<const Eigen::Array<dolfin::la_index_t, Eigen::Dynamic, 1>>
              macro_dofs_map1(macro_dofs1.data(), macro_dmap1_size);

      // FIXME: is this BC application correct?
      apply_bc_to_local_matrix(Ae, boundary_values, macro_dofs_map0, macro_dofs_map1);

      // Add to matrix
      A.add_local(Ae.data(),
                  macro_dmap0_size, macro_dofs0.data(),
                  macro_dmap1_size, macro_dofs1.data());
    }
  }

  // Flush matrix
  A.apply(la::PETScMatrix::AssemblyType::FLUSH);

  // FIXME: Move this outside of function?
  // Place '1' on diagonal for bc entries
  if (spaces[0] == spaces[1])
  {
    // Note: set diagonal using PETScMatrix::set_local since other functions,
    // e.g. PETScMatrix::set_local, do not work for all PETSc Mat types
    for (auto bc : boundary_values[0])
    {
      la_index_t row = bc.first;
      double one = 1.0;
      A.set_local(&one, 1, &row, 1, &row);
    }
  }

  // Finalise matrix
  A.apply(la::PETScMatrix::AssemblyType::FINAL);
}
//-----------------------------------------------------------------------------
void Assembler::assemble(la::PETScVector& b, const Form& L)
{
  // if (b.empty())
  //  init(b, L);
  // FIXME: Apply BCs

  // Get mesh from form
  assert(L.mesh());
  const mesh::Mesh& mesh = *L.mesh();

  // FIXME: Remove UFC
  // Create data structures for local assembly data
  UFC ufc(L);

  const std::size_t gdim = mesh.geometry().dim();
  const std::size_t tdim = mesh.topology().dim();
  mesh.init(tdim);

  // Collect pointers to dof maps
  auto dofmap = L.function_space(0)->dofmap();

  // Data structures used in assembly
  EigenRowArrayXXd coordinate_dofs;
  EigenVectorXd be;

  // Check whether integral is domain-dependent
  auto cell_domains = L.cell_domains();
  bool use_cell_domains = cell_domains && cell_domains->size() > 0;

  if (L.integrals().num_cell_integrals() > 0)
  {
    // FIXME: This will not be valid if num_element_dofs can vary by cell in the future
    const auto init_element_vector = [&dofmap](dolfin::EigenRowMatrixXd& be)
            -> void { be.resize(dofmap->num_element_dofs(0), 1); };

    const auto add_local_vector_to_global = [&b, &dofmap]
            (dolfin::EigenRowMatrixXd& be, dolfin::mesh::MeshEntity& cell)
            -> void {
      const auto dmap = dofmap->cell_dofs(cell.index());
      b.add_local(be.data(), dmap.size(), dmap.data());
    };

    assemble_over_cells(L, init_element_vector, add_local_vector_to_global);
  }

  // Exterior facet assembly
  mesh.init(tdim - 1);
  mesh.init(tdim - 1, tdim);

  // Get exterior facet integral
  auto exterior_facet_integral = L.integrals().exterior_facet_integral();

  // Check whether facet integral is domain dependent
  auto exterior_facet_domains = L.exterior_facet_domains();
  auto use_exterior_facet_domains = exterior_facet_domains && exterior_facet_domains->size() > 0;

  // Iterate over exterior facets
  for (const auto& facet : mesh::MeshRange<mesh::Facet>(mesh))
  {
    if (!facet.exterior())
      continue;

    if (use_exterior_facet_domains)
      exterior_facet_integral = L.integrals().exterior_facet_integral((*exterior_facet_domains)[facet]);

    if (!exterior_facet_integral)
      continue;

    assert(facet.num_entities(tdim) == 1);
    const mesh::Cell mesh_cell(mesh, facet.entities(tdim)[0]);

    // Get dof maps for cell
    auto dmap = dofmap->cell_dofs(mesh_cell.index());

    assert(!mesh_cell.is_ghost());
    const std::size_t local_facet = mesh_cell.index(facet);

    coordinate_dofs.resize(mesh_cell.num_vertices(), gdim);
    mesh_cell.get_coordinate_dofs(coordinate_dofs);

    ufc.update(mesh_cell, coordinate_dofs, exterior_facet_integral->enabled_coefficients);

    // Size data structure for assembly
    be.resize(dmap.size());
    be.setZero();

    // TODO: the cell orientation
    exterior_facet_integral->tabulate_tensor(be.data(), ufc.w(), coordinate_dofs.data(), local_facet, 1);

    // Add to vector
    b.add_local(be.data(), dmap.size(), dmap.data());
  }


  // Interior facet assembly
  // Sanity check of ghost mode (proper check in AssemblerBase::check)
  assert(mesh.get_ghost_mode() == "shared_vertex"
         or mesh.get_ghost_mode() == "shared_facet"
         or MPI::size(mesh.mpi_comm()) == 1);

  mesh.init(tdim - 1);
  mesh.init(tdim - 1, tdim);

  auto interior_facet_integral = L.integrals().interior_facet_integral();

  // Check whether facet integral is domain dependent
  auto interior_facet_domains = L.interior_facet_domains();
  auto use_interior_facet_domains = interior_facet_domains && interior_facet_domains->size() > 0;

  std::vector<EigenRowArrayXXd, Eigen::aligned_allocator<EigenRowArrayXXd>> neighbour_coordinate_dofs(2);

  const std::size_t mpi_rank = MPI::rank(mesh.mpi_comm());

  for (const auto& facet : mesh::MeshRange<mesh::Facet>(mesh))
  {
    if (facet.exterior())
      continue;

    assert(!facet.is_ghost());

    if (use_interior_facet_domains)
      interior_facet_integral = L.integrals().interior_facet_integral((*interior_facet_domains)[facet]);

    if (!interior_facet_integral)
      continue;

    assert(facet.num_entities(tdim) == 2);
    auto cell_index_plus = facet.entities(tdim)[0];
    auto cell_index_minus = facet.entities(tdim)[1];

    if (use_cell_domains and (*cell_domains)[cell_index_plus] < (*cell_domains)[cell_index_minus])
      std::swap(cell_index_plus, cell_index_minus);

    // The convention '+' = 0, '-' = 1 is from ffc
    const mesh::Cell cell0(mesh, cell_index_plus);
    const mesh::Cell cell1(mesh, cell_index_minus);

    // Get local index of facet with respect to each cell
    std::size_t local_facet0 = cell0.index(facet);
    std::size_t local_facet1 = cell1.index(facet);

    // Update to current pair of cells
    neighbour_coordinate_dofs[0].resize(cell0.num_vertices(), gdim);
    cell0.get_coordinate_dofs(neighbour_coordinate_dofs[0]);

    neighbour_coordinate_dofs[1].resize(cell1.num_vertices(), gdim);
    cell1.get_coordinate_dofs(neighbour_coordinate_dofs[1]);

    ufc.update(cell0, neighbour_coordinate_dofs[0],
               cell1, neighbour_coordinate_dofs[1],
               interior_facet_integral->enabled_coefficients);

    // Get dof maps for cells
    auto dmap0 = dofmap->cell_dofs(cell0.index());
    auto dmap1 = dofmap->cell_dofs(cell1.index());
    const std::size_t macro_dmap_size = dmap0.size() + dmap1.size();

    // Prepare for assembly
    be.resize(macro_dmap_size);
    be.setZero();

    // TODO: sort out the orientation
    interior_facet_integral->tabulate_tensor(be.data(), ufc.macro_w(),
                                             neighbour_coordinate_dofs[0].data(),
                                             neighbour_coordinate_dofs[1].data(),
                                             local_facet0, local_facet1,
                                             1, 1);

    if (cell0.is_ghost() != cell1.is_ghost())
    {
      const std::size_t ghost_rank = cell0.is_ghost() ? cell0.owner() : cell1.owner();
      assert(mpi_rank != ghost_rank);
      if (ghost_rank < mpi_rank)
        continue;
    }

    // FIXME: Can this be replaced with a Ref or Map somehow?
    Eigen::Array<dolfin::la_index_t, Eigen::Dynamic, 1> macro_dofs(macro_dmap_size);
    macro_dofs << dmap0, dmap1;
    b.add_local(be.data(), macro_dmap_size, macro_dofs.data());
  }

  // FIXME: Put this elsewhere?
  // Finalise matrix
  b.apply();
}
//-----------------------------------------------------------------------------
void Assembler::assemble(la::Scalar& m, const Form& M)
{
  // Get mesh from form
  assert(M.mesh());
  const mesh::Mesh& mesh = *M.mesh();

  if (M.rank() != 0)
  {
    throw std::runtime_error("Expected form of rank 0 for scalar assembly, "
                             "form has rank " + std::to_string(M.rank()));
  }

  // FIXME: Remove UFC
  // Create data structures for local assembly data
  UFC ufc(M);

  const std::size_t gdim = mesh.geometry().dim();
  const std::size_t tdim = mesh.topology().dim();
  mesh.init(tdim);

//  // Data structures used in assembly
  EigenRowArrayXXd coordinate_dofs;
  EigenVectorXd me;
  me.resize(1);

  // Check whether integral is domain-dependent
  auto cell_domains = M.cell_domains();
  bool use_cell_domains = cell_domains && cell_domains->size() > 0;

  // Iterate over all cells
  if (M.integrals().num_cell_integrals() > 0)
  {
    const auto init_scalar = [](dolfin::EigenRowMatrixXd& Ae) -> void { Ae.resize(1, 1); };
    const auto add_scalar_to_global_scalar = [&m]
            (dolfin::EigenRowMatrixXd& Ae, dolfin::mesh::MeshEntity& cell)
            -> void { m.add(Ae.data()[0]); };

    assemble_over_cells(M, init_scalar, add_scalar_to_global_scalar);
  }

  // Exterior facet assembly
  mesh.init(tdim - 1);
  mesh.init(tdim - 1, tdim);

  // Get exterior facet integral
  auto exterior_facet_integral = M.integrals().exterior_facet_integral();

  // Check whether facet integral is domain dependent
  auto exterior_facet_domains = M.exterior_facet_domains();
  auto use_exterior_facet_domains = exterior_facet_domains && exterior_facet_domains->size() > 0;

  // Iterate over exterior facets
  for (const auto& facet : mesh::MeshRange<mesh::Facet>(mesh))
  {
    if (!facet.exterior())
      continue;

    if (use_exterior_facet_domains)
      exterior_facet_integral = M.integrals().exterior_facet_integral((*exterior_facet_domains)[facet]);

    if (!exterior_facet_integral)
      continue;

    assert(facet.num_entities(tdim) == 1);
    const mesh::Cell mesh_cell(mesh, facet.entities(tdim)[0]);

    assert(!mesh_cell.is_ghost());
    const std::size_t local_facet = mesh_cell.index(facet);

    coordinate_dofs.resize(mesh_cell.num_vertices(), gdim);
    mesh_cell.get_coordinate_dofs(coordinate_dofs);

    ufc.update(mesh_cell, coordinate_dofs, exterior_facet_integral->enabled_coefficients);

    me.setZero();

    // TODO: the cell orientation
    exterior_facet_integral->tabulate_tensor(me.data(), ufc.w(), coordinate_dofs.data(), local_facet, 1);

    m.add(me.data()[0]);
  }


  // Interior facet assembly
  // Sanity check of ghost mode (proper check in AssemblerBase::check)
  assert(mesh.get_ghost_mode() == "shared_vertex"
         or mesh.get_ghost_mode() == "shared_facet"
         or MPI::size(mesh.mpi_comm()) == 1);

  mesh.init(tdim - 1);
  mesh.init(tdim - 1, tdim);

  auto interior_facet_integral = M.integrals().interior_facet_integral();

  // Check whether facet integral is domain dependent
  auto interior_facet_domains = M.interior_facet_domains();
  auto use_interior_facet_domains = interior_facet_domains && interior_facet_domains->size() > 0;

  std::vector<EigenRowArrayXXd, Eigen::aligned_allocator<EigenRowArrayXXd>> neighbour_coordinate_dofs(2);

  const std::size_t mpi_rank = MPI::rank(mesh.mpi_comm());

  for (const auto& facet : mesh::MeshRange<mesh::Facet>(mesh))
  {
    if (facet.exterior())
      continue;

    assert(!facet.is_ghost());

    if (use_interior_facet_domains)
      interior_facet_integral = M.integrals().interior_facet_integral((*interior_facet_domains)[facet]);

    if (!interior_facet_integral)
      continue;

    assert(facet.num_entities(tdim) == 2);
    auto cell_index_plus = facet.entities(tdim)[0];
    auto cell_index_minus = facet.entities(tdim)[1];

    if (use_cell_domains and (*cell_domains)[cell_index_plus] < (*cell_domains)[cell_index_minus])
      std::swap(cell_index_plus, cell_index_minus);

    // The convention '+' = 0, '-' = 1 is from ffc
    const mesh::Cell cell0(mesh, cell_index_plus);
    const mesh::Cell cell1(mesh, cell_index_minus);

    // Get local index of facet with respect to each cell
    std::size_t local_facet0 = cell0.index(facet);
    std::size_t local_facet1 = cell1.index(facet);

    // Update to current pair of cells
    neighbour_coordinate_dofs[0].resize(cell0.num_vertices(), gdim);
    cell0.get_coordinate_dofs(neighbour_coordinate_dofs[0]);

    neighbour_coordinate_dofs[1].resize(cell1.num_vertices(), gdim);
    cell1.get_coordinate_dofs(neighbour_coordinate_dofs[1]);

    ufc.update(cell0, neighbour_coordinate_dofs[0],
               cell1, neighbour_coordinate_dofs[1],
               interior_facet_integral->enabled_coefficients);

    // Prepare for assembly
    me.setZero();

    // TODO: sort out the orientation
    interior_facet_integral->tabulate_tensor(me.data(), ufc.macro_w(),
                                             neighbour_coordinate_dofs[0].data(),
                                             neighbour_coordinate_dofs[1].data(),
                                             local_facet0, local_facet1,
                                             1, 1);

    if (cell0.is_ghost() != cell1.is_ghost())
    {
      const std::size_t ghost_rank = cell0.is_ghost() ? cell0.owner() : cell1.owner();
      assert(mpi_rank != ghost_rank);
      if (ghost_rank < mpi_rank)
        continue;
    }

    m.add(me.data()[0]);
  }

  // FIXME: Put this elsewhere?
  // Finalise matrix
  m.apply();
}
//-----------------------------------------------------------------------------
void Assembler::apply_bc(la::PETScVector& b, const Form& a,
                         std::vector<std::shared_ptr<const DirichletBC>> bcs)
{
  // Get mesh from form
  assert(a.mesh());
  const mesh::Mesh& mesh = *a.mesh();

  const std::size_t gdim = mesh.geometry().dim();

  // Get bcs
  DirichletBC::Map boundary_values;
  for (std::size_t i = 0; i < bcs.size(); ++i)
  {
    assert(bcs[i]);
    assert(bcs[i]->function_space());
    if (a.function_space(1)->contains(*bcs[i]->function_space()))
    {
      bcs[i]->get_boundary_values(boundary_values);
      if (MPI::size(mesh.mpi_comm()) > 1
          and bcs[i]->method() != DirichletBC::Method::pointwise)
      {
        bcs[i]->gather(boundary_values);
      }
    }
  }

  // std::array<const function::FunctionSpace*, 2> spaces
  //    = {{a.function_space(0).get(), a.function_space(1).get()}};

  // Get dofmap for columns a a[i]
  auto dofmap0 = a.function_space(0)->dofmap();
  auto dofmap1 = a.function_space(1)->dofmap();

  EigenRowMatrixXd Ae;
  EigenVectorXd be;
  EigenRowArrayXXd coordinate_dofs;

  // Create data structures for local assembly data
  UFC ufc(a);

  // Get cell integral
  auto cell_integral = a.integrals().cell_integral();

  // Iterate over all cells
  for (auto& cell : mesh::MeshRange<mesh::Cell>(mesh))
  {
    // Check that cell is not a ghost
    assert(!cell.is_ghost());

    // Get dof maps for cell
    auto dmap1 = dofmap1->cell_dofs(cell.index());

    // Check if bc is applied to cell
    bool has_bc = false;
    for (int i = 0; i < dmap1.size(); ++i)
    {
      const std::size_t ii = dmap1[i];
      if (boundary_values.find(ii) != boundary_values.end())
      {
        has_bc = true;
        break;
      }
    }

    // std::cout << "Applying bcs" << std::endl;
    if (!has_bc)
      continue;
    // std::cout << "  has bc" << std::endl;

    // Get cell vertex coordinates
    coordinate_dofs.resize(cell.num_vertices(), gdim);
    cell.get_coordinate_dofs(coordinate_dofs);

    // Update UFC data to current cell
    ufc.update(cell, coordinate_dofs, cell_integral->enabled_coefficients);

    // Size data structure for assembly
    auto dmap0 = dofmap1->cell_dofs(cell.index());
    Ae.resize(dmap0.size(), dmap1.size());
    Ae.setZero();
    cell_integral->tabulate_tensor(Ae.data(), ufc.w(), coordinate_dofs.data(),
                                   1);

    // FIXME: Is this required?
    // Zero Dirichlet rows in Ae
    /*
    if (spaces[0] == spaces[1])
    {
      for (int i = 0; i < dmap0.size(); ++i)
      {
        const std::size_t ii = dmap0[i];
        auto bc = boundary_values.find(ii);
        if (bc != boundary_values.end())
          Ae.row(i).setZero();
      }
    }
    */

    // Size data structure for assembly
    be.resize(dmap0.size());
    be.setZero();

    for (int j = 0; j < dmap1.size(); ++j)
    {
      const std::size_t jj = dmap1[j];
      auto bc = boundary_values.find(jj);
      if (bc != boundary_values.end())
      {
        be -= Ae.col(j) * bc->second;
      }
    }

    // Add to vector
    b.add_local(be.data(), dmap0.size(), dmap0.data());
  }

  // FIXME: Put this elsewhere?
  // Finalise matrix
  b.apply();
}
//-----------------------------------------------------------------------------
void Assembler::set_bc(la::PETScVector& b, const Form& L,
                       std::vector<std::shared_ptr<const DirichletBC>> bcs)
{
  // Get mesh from form
  assert(L.mesh());
  const mesh::Mesh& mesh = *L.mesh();

  auto V = L.function_space(0);

  // Get bcs
  DirichletBC::Map boundary_values;
  for (std::size_t i = 0; i < bcs.size(); ++i)
  {
    assert(bcs[i]);
    assert(bcs[i]->function_space());
    if (V->contains(*bcs[i]->function_space()))
    {
      bcs[i]->get_boundary_values(boundary_values);
      if (MPI::size(mesh.mpi_comm()) > 1
          and bcs[i]->method() != DirichletBC::Method::pointwise)
      {
        bcs[i]->gather(boundary_values);
      }
    }
  }

  std::vector<double> values;
  values.reserve(boundary_values.size());
  std::vector<la_index_t> rows;
  rows.reserve(boundary_values.size());
  for (auto bc : boundary_values)
  {
    rows.push_back(bc.first);
    values.push_back(bc.second);
  }

  b.set_local(values.data(), values.size(), rows.data());
  b.apply();
}
//-----------------------------------------------------------------------------
void Assembler::assemble_over_cells(const Form &form,
                                    const std::function<void(EigenRowMatrixXd& Ae)>& initialise_element_tensor,
                                    const std::function<void(EigenRowMatrixXd& Ae, mesh::MeshEntity& cell)>& insert_element_to_global_tensor)
{
  assert(form.mesh());
  const mesh::Mesh& mesh = *form.mesh();

  // Create data structures for local assembly data
  UFC ufc(form);

  const std::size_t gdim = mesh.geometry().dim();
  const std::size_t tdim = mesh.topology().dim();
  mesh.init(tdim);

  // Get cell integral
  auto cell_integral = form.integrals().cell_integral();

  // Check whether integral is domain-dependent
  auto cell_domains = form.cell_domains();
  bool use_cell_domains = cell_domains && cell_domains->size() > 0;

  // DoF coords
  EigenRowArrayXXd coordinate_dofs;

  // Element tensor
  EigenRowMatrixXd Ae;
  initialise_element_tensor(Ae);

  // Iterate over all cells
  for (auto &cell : mesh::MeshRange<mesh::Cell>(mesh))
  {
    // Get integral for sub domain (if any)
    if (use_cell_domains)
      cell_integral = form.integrals().cell_integral((*cell_domains)[cell]);

    // Skip if no integral on current domain
    if (!cell_integral)
      continue;

    // Check that cell is not a ghost
    assert(!cell.is_ghost());

    // Get cell vertex coordinates
    coordinate_dofs.resize(cell.num_vertices(), gdim);
    cell.get_coordinate_dofs(coordinate_dofs);

    // Update UFC data to current cell
    ufc.update(cell, coordinate_dofs, cell_integral->enabled_coefficients);

    // Zero Ae
    Ae.setZero();

    // Compute cell matrix
    cell_integral->tabulate_tensor(Ae.data(), ufc.w(), coordinate_dofs.data(),
                                   1);
    insert_element_to_global_tensor(Ae, cell);
  }
}
//-----------------------------------------------------------------------------
