# flat.jl
# 
# An flat surface illuminated by an plane wave with Dirichlet boundary condition.

using Nmrc
using Ferrite

## Parameters
# Incident field
k = 1.0
θ = π/3
inc = Incident(k, θ)

# Number of the truncated terms
N = 10
height = 2.0

# Generate the mesh in a periodic cell
grid = periodic_cell(0.1; height=height)

## Set up fevalues(CellValues and FacetValues), DofHandler, and ConstraintHandler
ip = Lagrange{RefTriangle, 1}() 
cv, fv = setup_vals(ip) 
dh = setup_dh(grid, ip)
cst = setup_bcs(dh)

## Allocate the stiffness matrix A, the TBC matrix F and the load vector f 
# Extract dofs on the "top" boundary
top = getfacetset(grid, "top")
dofsDtN = dofs_on_dtn(dh, :u, top)

# Allocate the stiffness matrix and the load vector
A = allocate_stiff_matrix(dh, cst, dofsDtN)
F = allocate_stiff_matrix(dh, cst, dofsDtN)
f = zeros(ComplexF64, ndofs(dh))

# Assemble the matrix A
A = assemble_A(cv, dh, A, inc)

# Assemble the load vector f
f = assemble_load(fv, dh, top, f, inc, height)

# Assemble the TBC matrix
F = assemble_tbc(fv, dh, inc, top, F, N, dofsDtN)

# Add the TBC matrix to A: A - F
A = sub_preserve_structure(A, F)

# Impose the boundary condition
apply!(A, f, cst)

# Solve the linear system
u = A\f 
apply!(u, cst)

# Write the solution to vtk file
VTKGridFile("real_u", grid) do vtk
   write_solution(vtk, dh, real.(u)) 
end;

VTKGridFile("imag_u", grid) do vtk
    write_solution(vtk, dh, imag.(u)) 
end;