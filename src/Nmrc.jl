module Nmrc

using Gmsh
using Ferrite, FerriteGmsh
using LinearAlgebra
using SparseArrays
using KrylovKit: svdsolve

# PML parameters
include("pml.jl")
export PML, get_width, coord_transform

# FEM routines
include("fem/setups.jl")
export setup_grid, setup_fevs, setup_dofs
export setup_bdcs, allocate_matrices
export assemble_A2, assemble_A1, assemble_A0, apply_all_bds!

# Quadratic eigenvalue problem solver
include("quadeig.jl")
export get_scaling_factors, compute_scaling

end
