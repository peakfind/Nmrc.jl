module Nmrc

using Gmsh
using Ferrite, FerriteGmsh
using SparseArrays

# PML parameters
include("pml.jl")
# FEM routines
include("fem/setups.jl")

end
