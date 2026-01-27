# Sinusoidal perfectly conducting grating
# 
# - A sinusoidal perfectly conducting grating of period 600nm
# and height 180nm. 
# - Incident plane wave with angle 30 degree and wavelength 600nm
# - TM polarization

using Gmsh
using Ferrite, FerriteGmsh
using Nmrc
using SparseArrays

function pec_grid(;lc = 0.1, n = 20)
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    
    # Add the points
    δ = 1.0 / n
    pts = zeros(n + 3)
    
    pts[1] = gmsh.model.geo.addPoint(0, 0, 0, lc)
    pts[n + 1] = gmsh.model.geo.addPoint(1.0, 0, 0, lc)
    
    # Points on the sinusoidal curve
    for i = 2:n
        pts[i] = gmsh.model.geo.addPoint(δ * (i - 1), 0.15 * sin(2π *  δ * (i - 1)), 0.0, lc)
    end
    
    pts[n + 2] = gmsh.model.geo.addPoint(1.0, 1.0, 0.0, lc)
    pts[n + 3] = gmsh.model.geo.addPoint(0.0, 1.0, 0.0, lc)
    
    # Add the lines
    lns = zeros(n + 3)
    
    for i = 1:n + 2
        lns[i] = gmsh.model.geo.addLine(pts[i], pts[i + 1])
    end
    
    lns[n + 3] = gmsh.model.geo.addLine(pts[n + 3], pts[1])
    
    # Create the loop and the surface
    loop = gmsh.model.geo.addCurveLoop(lns)
    surf = gmsh.model.geo.addPlaneSurface([loop])
    
    # Synchronize the model
    gmsh.model.geo.synchronize()
    
    # Create the physical groups
    gmsh.model.addPhysicalGroup(1, lns[1:n], -1, "bottom")
    gmsh.model.addPhysicalGroup(1, [lns[n + 1]], -1, "right")
    gmsh.model.addPhysicalGroup(1, [lns[n + 2]], -1, "top")
    gmsh.model.addPhysicalGroup(1, [lns[n + 3]], -1, "left")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Ω")
    
    # Set periodic boundary 
    transform = [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [lns[n + 1]], [lns[n + 3]], transform)
    
    # Generate a 2D mesh
    gmsh.model.mesh.generate(2)
    
    grid = mktempdir() do dir 
        path = joinpath(dir, "mesh.msh")
        gmsh.write(path)
        togrid(path)
    end
    
    # Finalize the Gmsh library
    gmsh.finalize()
    
    return grid
end

function assemble_rhs(fv::FacetValues, dh::DofHandler, facetset, f, inc::Incident)
    α = get_alpha(inc)
    ia = im * α
    β = get_beta(inc)
    ib = im * β
    
    # Allocate the local load vector fe
    n_basefuncs = getnbasefunctions(fv)
    fe = zeros(ComplexF64, n_basefuncs)
    
    # Loop over all facets on the specific facetset
    for facet in FacetIterator(dh, facetset)
        # Update the fv to the current facet
        reinit!(fv, facet)
        
        # Reset the local vector fe to zero
        fill!(fe, 0.0 + 0.0im)
        
        # Get the coordinates of this facet
        coords = getcoordinates(facet)
        
        # loop over all quadrature points
        for qp in 1:getnquadpoints(fv)
            ds = getdetJdV(fv, qp)
            
            # Get the exterior norm vector for the computational domain
            ν = getnormal(fv, qp)
            
            # Get the coordinates of the quadrature points
            coords_qp = spatial_coordinate(fv, qp, coords)

            g = ib * ν[2] - ia * ν[1] 
            e = exp(-ib * coords_qp[2])

            # Loop over the test shape functions
            for i in 1:n_basefuncs
                v = shape_value(fv, qp, i)
                
                fe[i] += (g * e * v) * ds
            end
        end
        
        assemble!(f, celldofs(facet), fe)
    end
    
    return f
end

function assemble_boundary_term(fv::FacetValues, dh::DofHandler, facetset, B::SparseMatrixCSC, inc::Incident)
    α = get_alpha(inc)
    ia = im * α
    
    # Allocate the local matrix
    n_basefuncs = getnbasefunctions(fv)
    Be = zeros(ComplexF64, n_basefuncs, n_basefuncs)
    
    # Create an assembler B
    assembler = start_assemble(B)
    
    # Loop over all facets
    for facet in FacetIterator(dh, facetset)
        reinit!(fv, facet)
        fill!(Be, 0.0 + 0.0im)
        
        # Loop over quadrature points
        for qp in 1:getnquadpoints(fv)
            ds = getdetJdV(fv, qp)
            ν = getnormal(fv, qp)
            
            # Loop over test functions
            for i in 1:n_basefuncs
                v = shape_value(fv, qp, i)
                for j in 1:n_basefuncs
                    u = shape_value(fv, qp, j)
                    Be[i, j] += (ia * u * v * ν[1]) * ds
                end
            end
        end
        
        assemble!(assembler, celldofs(facet), Be)
    end
    
    return B
end

function plus_preserve_structure(A::SparseMatrixCSC, B::SparseMatrixCSC)
    # TODO: Find a better way. It relates to the functions `allocate_stiff_matrix` and `assemble_tbc`.
    if A.colptr != B.colptr || A.rowval != B.rowval || size(A) != size(B)
        error("Matrices must have the same sparsity structure and dimensions.")
    end

    new_nzval = A.nzval + B.nzval
    return SparseMatrixCSC(size(A,1), size(A,2), A.colptr, A.rowval, new_nzval)
end

# Incident plane wave
k = 2π
θ = π / 6
inc = Incident(k, θ)
α = get_alpha(inc)
β = get_beta(inc)

# Number of the truncation terms
N = 20
height = 1.0

# Period
p = 1.0

# Generate the mesh
grid = pec_grid(lc = 0.01, n = 80)

# Set up CellValues & FacetValues based on the interpolation
ip = Lagrange{RefTriangle, 1}();
cv, fv = setup_vals(ip);

# Set up the DofHandler
# dh = setup_dh(grid, ip);
dh = DofHandler(grid)
add!(dh, :u, ip)
close!(dh)

# For the TM polarization, we only need to impose the 
# periodic boundary conditon on "left" and "right" boundaries
cst = ConstraintHandler(dh)
pfacets = collect_periodic_facets(dh.grid, "right", "left", x -> x + Vec{2}((p, 0.0)))
pbc = PeriodicDirichlet(:u, pfacets)
add!(cst, pbc)
close!(cst)

top = getfacetset(grid, "top")
bottom = getfacetset(grid, "bottom")
dofsDtN = dofs_on_dtn(dh, :u, top)

# Allocate the stiffness matrix and the load vector
A = allocate_stiff_matrix(dh, cst, dofsDtN)
B = allocate_stiff_matrix(dh, cst, dofsDtN)
F = allocate_stiff_matrix(dh, cst, dofsDtN)
f = zeros(ComplexF64, ndofs(dh));

# Assemble the matrix A
A = assemble_A(cv, dh, A, inc)

# Assemble the boundary term
B = assemble_boundary_term(fv, dh, bottom, B, inc)

# Assemble the load vector f
f = assemble_rhs(fv, dh, bottom, f, inc)

# Assemble the TBC matrix F
F = assemble_tbc(fv, dh, inc, top, F, N, dofsDtN; period = p)

# Add the TBC matrix to A: A - F
A = plus_preserve_structure(A, B)
A = sub_preserve_structure(A, F)

# Impose the boundary condition
apply!(A, f, cst)

# Solve the linear system
u = A \ f
apply!(u, cst)

# Compute Rayleigh coefficient
u0 = compute_rayleigh_n(fv, dh, top, u, 0; period = p)
u1 = compute_rayleigh_n(fv, dh, top, u, -1; period = p)

# 
β1 = sqrt(k^2 - (α - 2π)^2)
# 
e0 = abs(u0)^2
e1 = (abs(u1)^2) * β1 / β

ρ = e0 + e1