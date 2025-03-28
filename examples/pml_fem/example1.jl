using Nmrc
using Ferrite, FerriteGmsh
using MAT

# Periodic layer with height 1.0
function medium(x)
    n = 1.0

    if x[2] < 1.0
        n = 3.0 + sin(2.0 * x[1])
    end

    return n
end

# Wavenumber
k = 8.1

# PML parameters
p = PML(10, 1.0 + 2.0im, 2, 1.1, 2.0)

# Set up the grid

# Generate mesh by using `Gmsh.jl`
# grid = setup_grid(;d = 2π, ĥ=1.1, δ=2.0, lc=0.02, lp=0.02)

# Read `.msh` file directly
# grid = togrid("mesh1.msh")
grid = togrid("mesh2.msh")

# Define the interpolation: linear lagrange
ip = Lagrange{RefTriangle, 1}()

# Set up FE values, dofs, and boundary conditions
cv = setup_fevs(ip)
dh = setup_dofs(grid, ip)
cst = setup_bdcs(dh, 2π)

# Allocate A₀, A₁, A₂
A₀ = allocate_matrices(dh, cst)
A₁ = allocate_matrices(dh, cst)
A₂ = allocate_matrices(dh, cst)

# Assemble A₀, A₁, A₂
A₀ = assemble_A0(cv, dh, A₀, medium, p, k)
A₁ = assemble_A1(cv, dh, A₁, p)
A₂ = assemble_A2(cv, dh, A₂, p)

# Impose the boundary conditions
apply_all_bds!(A₀, cst)
apply_all_bds!(A₁, cst)
apply_all_bds!(A₂, cst)

# Write three matrices into `.mat` file and 
# compute the quadratic eigenvalue problem by MATLAB
matwrite("A0.mat", Dict("A0" => A₀))
matwrite("A1.mat", Dict("A1" => A₁))
matwrite("A2.mat", Dict("A2" => A₂))