using Nmrc
using Test

@testset "Square root" begin
    include("test_squareroot.jl")
end

@testset "Computing integrals" begin
    include("test_computing_integrals.jl")
end
