# Lamellar perfectly conducting grating
# 
# - A lamellar perfectly conducting grating of period 600nm 
# - The width of the hole t = 300nm
# - The width of the bump b = 300nm
# - Incident plane wave with angle 30 degree and wavelength 600nm
# - TE polarization

using Gmsh
using Ferrite, FerriteGmsh
using Nmrc

function pec_grid(;lc = 0.1)
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    
    # Add the points
    p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, lc)
    p2 = gmsh.model.geo.addPoint(0.5, 0.0, 0.0, lc)
    p3 = gmsh.model.geo.addPoint(0.5, 0.3, 0.0, lc)
    p4 = gmsh.model.geo.addPoint(1.0, 0.3, 0.0, lc)
    p5 = gmsh.model.geo.addPoint(1.0, 1.0, 0.0, lc)
    p6 = gmsh.model.geo.addPoint(0.0, 1.0, 0.0, lc)
    p7 = gmsh.model.geo.addPoint(0.0, 0.3, 0.0, lc)
    
    # Add the lines
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p5)
    l5 = gmsh.model.geo.addLine(p5, p6)
    l6 = gmsh.model.geo.addLine(p6, p7)
    l7 = gmsh.model.geo.addLine(p7, p1)
    
    # Create the loop and the surface
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5, l6, l7])
    surf = gmsh.model.geo.addPlaneSurface([loop])
    
    # Synchronize the model
    gmsh.model.geo.synchronize()
    
    # Create the physical groups
    gmsh.model.addPhysicalGroup(1, [l7, l1, l2, l3], -1, "bottom")
    gmsh.model.addPhysicalGroup(1, [l4], -1, "right")
    gmsh.model.addPhysicalGroup(1, [l5], -1, "top")
    gmsh.model.addPhysicalGroup(1, [l6], -1, "left")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Ω")
    
    # Set periodic boundary 
    transform = [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [l4], [l6], transform)
    
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
grid = pec_grid(lc = 0.01)

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

ρ = e0 + e1