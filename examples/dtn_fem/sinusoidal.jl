# Sinusoidal perfectly conducting grating
# 
# - A sinusoidal perfectly conducting grating of period 600nm
# and height 180nm. 
# - Incident plane wave with angle 30 degree and wavelength 600nm

using Gmsh
using Ferrite, FerriteGmsh
using Nmrc

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

cst = setup_bcs(dh; period = p) 

top = getfacetset(grid, "top")
dofsDtN = dofs_on_dtn(dh, :u, top)

# Allocate the stiffness matrix and the load vector
A = allocate_stiff_matrix(dh, cst, dofsDtN)
F = allocate_stiff_matrix(dh, cst, dofsDtN)
f = zeros(ComplexF64, ndofs(dh));

# Assemble the matrix A
A = assemble_A(cv, dh, A, inc)

# Assemble the load vector f
f = assemble_load(fv, dh, top, f, inc, height)

# Assemble the TBC matrix F
F = assemble_tbc(fv, dh, inc, top, F, N, dofsDtN; period = p)

# Add the TBC matrix to A: A - F
A = sub_preserve_structure(A, F)

# Impose the boundary condition
apply!(A, f, cst)

# Solve the linear system
u = A \ f
apply!(u, cst)

# Compute Rayleigh coefficient
u0 = compute_rayleigh_n(fv, dh, top, u, 0; period = p) - exp(-im * β * height)
u1 = compute_rayleigh_n(fv, dh, top, u, -1; period = p)
# 
β1 = sqrt(k^2 - (α - 2π)^2)
# 
e0 = abs(u0)^2
e1 = (abs(u1)^2) * β1 / β