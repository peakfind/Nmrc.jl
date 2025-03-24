# Basic setup functions to use Ferrite.jl

"""
    setup_grid(;d=2π, ĥ=1.1, δ=0.5, lc=0.5)

Generate a mesh by julia interface for Gmsh and read the `.msh` file by `FerriteGmsh`

# Arguments

- `d`: the period of the periodic layer
- `h₀`: the start of the PML layer
- `δ`: the height of the PML layer
- `lc`: target mesh size close to the point

# Points

```
 4  -------------  3
 |                 |
 1  -------------  2
```

# Lines

```
 .  ----  3  ----  .
 4                 2 
 .  ----  1  ----  .
```

# Misc 

- We can generate a `*.vtk` file by replace the `do end` block with `gmsh.write(*.vtk)`
"""
function setup_grid(;d=2π, ĥ=1.1, δ=0.5, lc=0.5)
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)

    # Compute the height of the truncated domain
    b = ĥ + δ

    # Add the points
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(d, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(d, b, 0, lc) 
    p4 = gmsh.model.geo.addPoint(0, b, 0, lc)

    # Add the lines
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    # Add the loop and the surface
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.geo.addPlaneSurface([loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Create the Physical Groups
    gmsh.model.addPhysicalGroup(1, [l1], -1, "bottom")
    gmsh.model.addPhysicalGroup(1, [l2], -1, "right")
    gmsh.model.addPhysicalGroup(1, [l3], -1, "top")
    gmsh.model.addPhysicalGroup(1, [l4], -1, "left")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Ω")

    # Set periodic mesh
    transform = [1, 0, 0, d, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [l2], [l4], transform)

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

function setup_grid(;d=2π, ĥ=1.1, δ=0.5, lc=0.5, lp=0.5)
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

    path = joinpath(pwd(), "mesh.vtk")
    gmsh.write(path)
end

"""
    setup_fevs(ip)

Define quadrature rules on Reference triangle and setup FE values.
"""
function setup_fevs(ip)
    qr = QuadratureRule{RefTriangle}(2)
    cv = CellValues(qr, ip)
    return cv
end

"""
    setup_dofs(grid::Grid, ip)

Setup degree of freedoms.
"""
function setup_dofs(grid::Grid, ip)
    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)
    return dh
end

"""
    setup_bdcs(dof::DofHandler, d)

Set the periodic boundary condition ("left" and "right") and 
Dirichelt boundary condition ("bottom" and "top").
"""
function setup_bdcs(dof::DofHandler, d)
    cst = ConstraintHandler(dof)

    # Set periodic boundary condition on the "right" and "left"
    pfacets = collect_periodic_facets(dof.grid, "right", "left", x -> x + Vec{2}((d, 0.0)))
    pbc = PeriodicDirichlet(:u, pfacets)
    add!(cst, pbc)

    # Set Dirichlet boundary condition on the "bottom" and "top"
    dfacets = union(getfacetset(dof.grid, "bottom"), getfacetset(dof.grid, "top"))
    dbc = Dirichlet(:u, dfacets, x -> 0)
    add!(cst, dbc)

    close!(cst)
    return cst
end

"""
    allocate_matrices(dof::DofHandler, cst::ConstraintHandler)

Allocate the complex-valued stiffness matrices.
"""
function allocate_matrices(dof::DofHandler, cst::ConstraintHandler)
    sp = init_sparsity_pattern(dof)
    add_cell_entries!(sp, dof)
    add_constraint_entries!(sp, cst)
    M = allocate_matrix(SparseMatrixCSC{ComplexF64, Int}, sp)
    return M
end

function assemble_A2()
    
end

function assemble_A1()
    
end

function assemble_A0()
    
end


"""
    solve_qm(cv::CellValues, dof::DofHandler, cst::ConstraintHandler, A₀::SparseMatrixCSC, A₁::SparseMatrixCSC, A₂::SparseMatrixCSC)

Get the exceptional values (quasi-momentum α) by solving the eigenvalue 
problems generated by FEM.
"""
function solve_qm(cv::CellValues, dof::DofHandler, cst::ConstraintHandler, A₀::SparseMatrixCSC, A₁::SparseMatrixCSC, A₂::SparseMatrixCSC)
    
end
