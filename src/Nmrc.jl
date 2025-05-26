module Nmrc

using Gmsh
using Ferrite
using FerriteGmsh: togrid
using LinearAlgebra
using SparseArrays
using KrylovKit: svdsolve

# PML parameters
include("pml.jl")
export PML, get_width, coord_transform

# Incident plane wave
include("incident.jl")
export Incident
export get_wavenumber, get_alpha, get_beta
export beta_n

# FEM routines
include("fem/fem_pml.jl")
export setup_grid, setup_fevs, setup_dofs
export setup_bdcs, allocate_matrices
export assemble_pml_A2, assemble_pml_A1, assemble_pml_A0, apply_all_bds!

include("fem/fem_dtn.jl")
export periodic_cell
export setup_vals, setup_dh
export dofs_on_dtn, setup_bcs, allocate_stiff_matrix
export assemble_A, assemble_load, assemble_tbc
export sub_preserve_structure

# Quadratic eigenvalue problem solver
include("quadeig.jl")
export get_scaling_factors, compute_scaling

end