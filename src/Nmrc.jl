module Nmrc

using Gmsh
using Ferrite, FerriteGmsh
using SparseArrays
using KrylovKit: svdsolve

# PML parameters
include("pml.jl")
export PML, get_width, coord_transform

# FEM routines
include("fem/setups.jl")
export setup_grid

# Quadratic eigenvalue problem solver
include("quadeig.jl")
export get_scaling_factors, compute_scaling

end
