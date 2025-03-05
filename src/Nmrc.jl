module Nmrc

using Gmsh
using Ferrite, FerriteGmsh
using SparseArrays
using KrylovKit: svdsolve

export get_scaling_factors, compute_scaling

# PML parameters
include("pml.jl")

# FEM routines
include("fem/setups.jl")

# Quadratic eigenvalue problem solver
include("quadeig.jl")

end
