using Gmsh

function write_mesh(d, ĥ, δ; lc=0.5, lp=0.5, filename="mesh")
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)

    # Compute the height of the truncated domain
    b = ĥ + δ

    # Add the points
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(d, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(d, ĥ, 0, lp) 
    p4 = gmsh.model.geo.addPoint(d, b, 0, lp)
    p5 = gmsh.model.geo.addPoint(0, b, 0, lp)
    p6 = gmsh.model.geo.addPoint(0, ĥ, 0, lp)

    # Add the lines 
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p5)
    l5 = gmsh.model.geo.addLine(p5, p6)
    l6 = gmsh.model.geo.addLine(p6, p1)

    # Add the loop and the surface
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5, l6])
    surf = gmsh.model.geo.addPlaneSurface([loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Create the Physical Groups
    gmsh.model.addPhysicalGroup(1, [l1], -1, "bottom")
    gmsh.model.addPhysicalGroup(1, [l2, l3], -1, "right")
    gmsh.model.addPhysicalGroup(1, [l4], -1, "top")
    gmsh.model.addPhysicalGroup(1, [l5, l6], -1, "left")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Ω")

    # Set periodic mesh
    transform = [1, 0, 0, d, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [l2, l3], [l5, l6], transform)

    # Generate a 2D mesh
    gmsh.model.mesh.generate(2)


    vtk_path = joinpath(pwd(), filename * ".vtk")
    msh_path = joinpath(pwd(), filename * ".msh")
    gmsh.write(vtk_path)
    gmsh.write(msh_path)

    # Finalize the Gmsh library
    gmsh.finalize()
    return
end

# d = 2π, ĥ = 1.1, δ = 2.0
# Generate mesh1: lc = 0.01, lp = 0.01
write_mesh(2π, 1.1, 2.0; lc=0.01, lp=0.01, filename="mesh1")

# Generate mesh2: lc = 0.01, lp = 0.008
write_mesh(2π, 1.1, 2.0; lc=0.01, lp=0.008, filename="mesh2")