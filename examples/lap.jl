using Ferrite
using Nmrc

# functions from Nmrc we need:
#     periodic_cell, setup_vals, setup_dh, setup_bcs
#     allocate_stiff_matrix, assemble_A0, assemble_A1, assemble_A2
#     compute_coef!, assemble_tbc TODO: We have modified these two functions
#     assemble_load
#     dofs_on_dtn, Incident

using SparseArrays
using Random, LinearAlgebra, Cim


# ---------------------------------------------------------
# Functions 
# ---------------------------------------------------------

function medium(x)
    if x[2] < 1.0
        return 3.9
    else
        return 1.0
    end
end

# Assemble TBC part for nonlinear eigenvalue problems
#
# TODO: Since we use the new compute_coef!, we also need to 
#.      modify this function
# Difference between this function and the assmeble_tbc() in Nmrc.jl:
#     in this function, we do not have inc::Incident as our argument
function assemble_tbc(fv::FacetValues, dh::DofHandler, F::SparseMatrixCSC, facetset, dofsDtN, N, k, α; period = 2π)
    # Allocate the vector Θ
    Θ = sparsevec(dofsDtN, zeros(ComplexF64, length(dofsDtN)), ndofs(dh))
    
    # Loop over truncated terms
    for n in -N:N 
        # Reset the vector Θ to zero
        fill!(Θ, zero(ComplexF64))
        
        # Compute βₙ
        βₙ = beta_n(k, α, n)
        
        # Compute the vector Θ (Fourier coefficients and their conjugates)
        Nmrc.compute_coef!(Θ, fv, dh, facetset, n; period = period)
        
        # Assemble the TBC matrix
        for i in Θ.nzind, j in Θ.nzind
            v = im * βₙ * Θ[i] * conj(Θ[j]) / period
            Ferrite.addindex!(F, v, i, j)
        end
    end
    
    return F
end

struct Nnep{T}
    A₀::SparseMatrixCSC{T, Int}
    A₁::SparseMatrixCSC{T, Int}
    A₂::SparseMatrixCSC{T, Int}
    fv::FacetValues
    dh::DofHandler
    cst::ConstraintHandler
    facetset
    dofsDtN
    N::Int64
    k
    period
end

function (nep::Nnep{T})(α::T) where T
    # Allocate the TBC matrix
    F = allocate_stiff_matrix(nep.dh, nep.cst, nep.dofsDtN)
    
    # Assemble the TBC matrix
    F = assemble_tbc(nep.fv, nep.dh, F, nep.facetset, nep.dofsDtN, nep.N, nep.k, α; period = nep.period)
    
    # Impose the boundary conditions
    apply!(F, nep.cst)
    
    return nep.A₀ + α * nep.A₁ + (α^2) * nep.A₂ - F
end

# The contour integral method for nonlinear eigenvalue problem
function new_cim(ctr::Cim.AbstractContour, nep::Nnep, d::Int, l::Int; n=50, tol=1e-12)
    # Input validation
    d > 0 || throw(ArgumentError("d must be positive"))
    l > 0 || throw(ArgumentError("l must be positive"))

    # Get the quadrature points
    pts = get_quadpts(ctr, n)

    # Preallocate arrays
    A0 = zeros(ComplexF64, d, l)
    A1 = zeros(ComplexF64, d, l)
    Random.seed!(10);
    Vhat = randn(ComplexF64, d, l)

    # Compute A0 and A1 with trapezoid rule
    for j in 1:pts.N - 1
        z = complex(pts.nodes[j, 1], pts.nodes[j, 2])
        z_prime = complex(pts.nodes_prime[j, 1], pts.nodes_prime[j, 2])
        invNEP_Vhat = nep(z) \ Vhat
        A0 .+= invNEP_Vhat * z_prime
        A1 .+= invNEP_Vhat * z * z_prime
    end
    A0 ./= (im * (pts.N - 1))
    A1 ./= (im * (pts.N - 1))

    # Compute the SVD of A0
    (V, Sigma, W) = svd(A0)
    @show Sigma

    # Handle rank deficiency
    if isempty(Sigma)
        @warn "No eigenvalues found!"
        return ComplexF64[]
    end

    # Determine the number of nonzero singular values 
    k = count(Sigma ./ Sigma[1] .> tol)

    # Compute the matrix B 
    Vk = V[:,1:k]
    Sigk = Sigma[1:k]
    Wk = W[:,1:k]

    # Diagonal is more efficient
    B = (Vk' * A1 * Wk) * Diagonal(1 ./ Sigk)

    # Compute the eigenvalues λ and the corresponding eigenvectors s
    λ, s = eigen(B)

    # Check if the eigenvalue lies inside the contour
    # filter!(λ -> Cim.is_inside(λ, ctr), lambda)
    
    # Store the indices of eigenvalues inside the contour
    good = Int[]
    
    for i in eachindex(λ)
        if Cim.is_inside(λ[i], ctr) == true
            push!(good, i)
        end
    end
    

    λ = @view λ[good]
    u = @view s[:, good]
    u = Vk * u

    return λ, u
end

function assemble(cv::CellValues, dh::DofHandler, A::SparseMatrixCSC, n::Function, inc::Incident)
    k = get_wavenumber(inc)
    α = get_alpha(inc)
    
    # Allocate the local stiffness matrix
    n_basefuncs = getnbasefunctions(cv)
    Ae = zeros(ComplexF64, n_basefuncs, n_basefuncs)
    
    # Create an assembler A
    assembler = start_assemble(A)
    
    # Loop over all cells
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for current cell
        reinit!(cv, cell)
        
        # Reset local stiffness matrix to 0.0 + 0.0
        fill!(Ae, 0.0 + 0.0im)
        
        # Get the coordinates of this cell
        coords = getcoordinates(cell)
        
        # Loop over quadrature points
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            
            # Get the coordinates of the quadrature point and 
            # evaluate the refractive index at the quadrature point
            coords_qp = spatial_coordinate(cv, qp, coords)
            ri = n(coords_qp)
            
            # Loop over test shape functions
            for i in 1:n_basefuncs
                v = shape_value(cv, qp, i)
                ∇v = shape_gradient(cv, qp, i)

                # Loop over trial shape functions
                for j in 1:n_basefuncs
                    u = shape_value(cv, qp, j)
                    ∇u = shape_gradient(cv, qp, j)
                    
                    # Assemble the local stiffness matrix
                    Ae[i, j] += (∇u ⋅ ∇v - 2im * α * ∇u[1] * v - ((k^2) * ri - (α^2)) * u * v) * dx
                end
            end
        end
        
        # Assemble Ae into A
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A
end

# TODO: we have implemented functions for computing integrals in Nmrc.jl

# ---------------------------------------------------------
# Define constants and common variables
# ---------------------------------------------------------

# Wavenumber
k = 4.1

# The number of the truncated terms
N = 20

# The height of the periodic cell
h = 2.0

# Period 
p = 2π

# Generate a mesh for the periodic cell
grid = periodic_cell(lc = 0.01, period = p, height = h)

# Interpolation
ip = Lagrange{RefTriangle, 1}()

# Set up CellValues and FacetValues
cv, fv = setup_vals(ip)

# Set up the DofHandler
dh = setup_dh(grid, ip)

# Set up the boundary condition
# Dirichlet boundary condition on "bottom"
# Periodic boundary condition on "left" and "right"
cst = setup_bcs(dh)

# Extract the facetset on the "top" boundary
top = getfacetset(grid, "top")

# Extract the DoFs on the "top" boundary
dofsDtN = dofs_on_dtn(dh, :u, top)

# ---------------------------------------------------------
# Step 1: Compute the exceptional values and eigenfunctions
# ---------------------------------------------------------

# Allocate A₀, A₁, A₂ for the quadratic part
A₀ = allocate_stiff_matrix(dh, cst, dofsDtN)
A₁ = allocate_stiff_matrix(dh, cst, dofsDtN)
A₂ = allocate_stiff_matrix(dh, cst, dofsDtN)

# Assemble A₀, A₁, A₂
A₀ = assemble_A0(cv, dh, A₀, medium, k)
A₁ = assemble_A1(cv, dh, A₁)
A₂ = assemble_A2(cv, dh, A₂)

# Impose the constraint on A₀, A₁, and A₂
apply!(A₀, cst)
apply!(A₁, cst)
apply!(A₂, cst)

# The size of the nonlinear eigenvalue problem
d = size(A₀, 1)

# Construct the nonlinear eigenvalue problem
L = Nnep(A₀, A₁, A₂, fv, dh, cst, top, dofsDtN, N, k, p)

# Define the contour according to the PML result
cir = Cim.circle([0.38, 0.0], 0.01)

# Solve the nonlinear eigenvalue problem by CIM
λ, v = new_cim(cir, L, d, 5; n = 50, tol = 1e-7)

# Get the exceptional value and the corresponding eigenfunction
α = λ[1]
ϕ = v[:, 1]

# Apply the constraint on the eigenfunction
apply!(ϕ, cst)

VTKGridFile("modu_phi", grid) do vtk
    write_solution(vtk, dh, abs.(ϕ))
end;
# 
VTKGridFile("real_phi", grid) do vtk
    write_solution(vtk, dh, real.(ϕ))
end;
# 
# VTKGridFile("imag_phi", grid) do vtk
#     write_solution(vtk, dh, imag.(ϕ))
# end;

# ---------------------------------------------------------
# Step 2: Get the particular solution
# ---------------------------------------------------------

# Compute the incidnet angle in radians based on the exceptional 
# value in Step 1
θ = asin(real(α) / k)

inc = Incident(k, θ)

# Allocate matrices
A = allocate_stiff_matrix(dh, cst, dofsDtN)
F = allocate_stiff_matrix(dh, cst, dofsDtN)
f = zeros(ComplexF64, ndofs(dh))

# Assemble the load vector f
f = assemble_load(fv, dh, top, f, inc, h)

# Assemble the matrix A
A = assemble(cv, dh, A, medium, inc)

F = Nmrc.assemble_tbc(fv, dh, inc, top, F, N, dofsDtN; period = p)

# Add the TBC matrix to A
A = sub_preserve_structure(A, F)

# Impose the boundary
apply!(A, f, cst)

# Solve the linear system
u₀ = A \ f

# Apply the boundary condition again
apply!(u₀, cst)

# VTKGridFile("modu_u0", grid) do vtk
#     write_solution(vtk, dh, abs.(u₀))
# end;
# 
# VTKGridFile("real_u0", grid) do vtk
#     write_solution(vtk, dh, real.(u₀))
# end;
# 
# VTKGridFile("imag_u0", grid) do vtk
#     write_solution(vtk, dh, imag.(u₀))
# end;

# ---------------------------------------------------------
# Step 3: Compute the coefficients in the constraint proposed 
#.        by G. Hu and A. Kirsch
# ---------------------------------------------------------

# Compute constants in the formulation
s = sin(θ)
s² = sin(θ)^2

c1 = integral_∂₁uv̄(cv, dh, ϕ, ϕ)

c2 = real(integral_nuv̄(cv, dh, medium, ϕ, ϕ))

c3 = integral_uv̄(cv, dh, ϕ, ϕ)

c4 = integral_∂₁uv̄(cv, dh, u₀, ϕ)

c5 = integral_nuv̄(cv, dh, medium, u₀, ϕ)

c6 = integral_uv̄(cv, dh, u₀, ϕ)

# Compute the coefficient
coef = (-s * c4 + im * k * c5 - im * k * s² * c6) / (s * c1 - im * k * c2 + im * k * s² * c3)

# ---------------------------------------------------------
# Step 4: Get the numerical solution
# ---------------------------------------------------------

u = u₀ + coef .* ϕ

VTKGridFile("modu_u", grid) do vtk
    write_solution(vtk, dh, abs.(u))
end;

VTKGridFile("real_u", grid) do vtk
    write_solution(vtk, dh, real.(u))
end;

VTKGridFile("imag_u", grid) do vtk
    write_solution(vtk, dh, imag.(u))
end;