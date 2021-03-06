# Copyright (C) 2018 Garth N. Wells
#
# This file is part of DOLFIN (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Unit tests for Newton solver assembly"""

import numpy as np
from petsc4py import PETSc

import dolfin
import dolfin.fem as fem
import dolfin.function as function
import ufl
from ufl import derivative, dx, grad, inner


class NonlinearPDEProblem(dolfin.cpp.nls.NonlinearProblem):
    """Nonlinear problem class for a PDE problem."""

    def __init__(self, F, u, bc):
        super().__init__()
        V = u.function_space()
        du = function.TrialFunction(V)
        self.L = F
        self.a = derivative(F, u, du)
        self.bc = bc
        self._F, self._J = None, None

    def form(self, x):
        x.update_ghosts()

    def F(self, x):
        """Assemble residual vector."""
        if self._F is None:
            self._F = fem.assemble_vector([self.L], [[self.a]], [self.bc],
                                          dolfin.cpp.fem.BlockType.monolithic,
                                          x)
        else:
            self._F = fem.assemble(self._F, self.L, [self.a], [self.bc], x)
        return self._F

    def J(self, x):
        """Assemble Jacobian matrix."""
        if self._J is None:
            self._J = fem.assemble_matrix([[self.a]], [self.bc],
                                          dolfin.cpp.fem.BlockType.monolithic)
        else:
            self._J = fem.assemble(self._J, self.a, [self.bc])
        return self._J


class NonlinearPDE_SNESProblem():
    def __init__(self, F, u, bc):
        super().__init__()
        V = u.function_space()
        du = function.TrialFunction(V)
        self.L = F
        self.a = derivative(F, u, du)
        self.a_comp = dolfin.fem.Form(self.a)
        self.bc = bc
        self._F, self._J = None, None
        self.u = u

    def F(self, snes, x, F):
        """Assemble residual vector."""
        _F = dolfin.cpp.la.PETScVector(F)
        _x = dolfin.cpp.la.PETScVector(x)
        _x.update_ghosts()
        x.copy(self.u.vector().vec())
        self.u.vector().update_ghosts()
        fem.assemble(_F, self.L, [self.a], [self.bc], _x)

    def J(self, snes, x, J, P):
        """Assemble Jacobian matrix."""
        _J = dolfin.cpp.la.PETScMatrix(J)
        # _x = dolfin.cpp.la.PETScVector(x)
        fem.assemble(_J, self.a, [self.bc])


def test_linear_pde():
    """Test Newton solver for a linear PDE"""
    # Create mesh and function space
    mesh = dolfin.generation.UnitSquareMesh(dolfin.MPI.comm_world, 12, 12)
    V = dolfin.function.FunctionSpace(mesh, ("Lagrange", 1))
    u = dolfin.function.Function(V)
    v = function.TestFunction(V)
    F = inner(10.0, v) * dx - inner(grad(u), grad(v)) * dx

    def boundary(x):
        """Define Dirichlet boundary (x = 0 or x = 1)."""
        return np.logical_or(x[:, 0] < 1.0e-8, x[:, 0] > 1.0 - 1.0e-8)

    u_bc = function.Function(V)
    u_bc.vector().set(1.0)
    u_bc.vector().update_ghosts()
    bc = fem.DirichletBC(V, u_bc, boundary)

    # Create nonlinear problem
    problem = NonlinearPDEProblem(F, u, bc)

    # Create Newton solver and solve
    solver = dolfin.cpp.nls.NewtonSolver(dolfin.MPI.comm_world)
    n, converged = solver.solve(problem, u.vector())
    assert converged
    assert n == 1


def test_nonlinear_pde():
    """Test Newton solver for a simple nonlinear PDE"""
    # Create mesh and function space
    mesh = dolfin.generation.UnitSquareMesh(dolfin.MPI.comm_world, 12, 15)
    V = dolfin.function.FunctionSpace(mesh, ("Lagrange", 1))
    u = dolfin.function.Function(V)
    v = function.TestFunction(V)
    F = inner(2.0, v) * dx - ufl.sqrt(u * u) * inner(
        grad(u), grad(v)) * dx - inner(u, v) * dx

    def boundary(x):
        """Define Dirichlet boundary (x = 0 or x = 1)."""
        return np.logical_or(x[:, 0] < 1.0e-8, x[:, 0] > 1.0 - 1.0e-8)

    u_bc = function.Function(V)
    u_bc.vector().set(1.0)
    u_bc.vector().update_ghosts()
    bc = fem.DirichletBC(V, u_bc, boundary)

    # Create nonlinear problem
    problem = NonlinearPDEProblem(F, u, bc)

    # Create Newton solver and solve
    u.vector().set(0.9)
    u.vector().update_ghosts()
    solver = dolfin.cpp.nls.NewtonSolver(dolfin.MPI.comm_world)
    n, converged = solver.solve(problem, u.vector())
    assert converged
    assert n < 6


def test_nonlinear_pde_snes():
    """Test Newton solver for a simple nonlinear PDE"""
    # Create mesh and function space
    mesh = dolfin.generation.UnitSquareMesh(dolfin.MPI.comm_world, 12, 15)
    V = dolfin.function.FunctionSpace(mesh, ("Lagrange", 1))
    u = dolfin.function.Function(V)
    v = function.TestFunction(V)
    F = inner(2.0, v) * dx - ufl.sqrt(u * u) * inner(
        grad(u), grad(v)) * dx - inner(u, v) * dx

    def boundary(x):
        """Define Dirichlet boundary (x = 0 or x = 1)."""
        return np.logical_or(x[:, 0] < 1.0e-8, x[:, 0] > 1.0 - 1.0e-8)

    u_bc = function.Function(V)
    u_bc.vector().set(1.0)
    u_bc.vector().update_ghosts()
    bc = fem.DirichletBC(V, u_bc, boundary)

    # Create nonlinear problem
    problem = NonlinearPDE_SNESProblem(F, u, bc)

    u.vector().set(0.9)
    u.vector().update_ghosts()

    b = dolfin.cpp.la.PETScVector(V.dofmap().index_map())
    J = dolfin.cpp.fem.init_matrix(problem.a_comp._cpp_object)

    # Create Newton solver and solve
    snes = PETSc.SNES().create()
    snes.setFunction(problem.F, b.vec())
    snes.setJacobian(problem.J, J.mat())

    snes.setTolerances(rtol=1.0e-9, max_it=10)
    snes.setFromOptions()

    snes.getKSP().setTolerances(rtol=1.0e-9)
    snes.solve(None, u.vector().vec())

    assert snes.getConvergedReason() > 0
    assert snes.getIterationNumber() < 6
    # print(snes.getIterationNumber())
    # print(snes.getFunctionNorm())
