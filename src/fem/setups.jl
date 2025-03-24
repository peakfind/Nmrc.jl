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
#= function setup_grid(;d=2π, ĥ=1.1, δ=0.5, lc=0.5, vtk=false)
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

    if vtk 
        path = joinpath(pwd(), "mesh.vtk")
        gmsh.write(path)
    end

    # Finalize the Gmsh library
    gmsh.finalize()
    return grid
end
 =#
function setup_grid(;d=2π, ĥ=1.5, δ=2.0, lc=0.5, lp=0.5, vtk=false)
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

    grid = mktempdir() do dir
        path = joinpath(dir, "mesh.msh")
        gmsh.write(path)
        togrid(path) 
    end

    if vtk 
        path = joinpath(pwd(), "mesh.vtk")
        gmsh.write(path)
    end

    # Finalize the Gmsh library
    gmsh.finalize()
    return grid
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

"""
    assemble_A2(cv, dh, A₂, p)

The second order term in the quadratic eigenvalue problem
"""
function assemble_A2(cv::CellValues, dh::DofHandler, A₂::SparseMatrixCSC, p::PML)
    # Allocate the local stiffness matrix
    n_basefuncs = getnbasefunctions(cv)
    Ae = zeros(ComplexF64, n_basefuncs, n_basefuncs)
    # Create an assembler A₂
    assembler = start_assemble(A₂)

    # Loop over all cells 
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        # Reset local stiffness matrix to 0.0 + 0.0im
        fill!(Ae, 0.0 + 0.0im)
        coords = getcoordinates(cell)

        # Loop over quadrature points
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            coords_qp = spatial_coordinate(cv, qp, coords)
            s = coord_transform(coords_qp[2], p)
        
            # Loop over test shape functions
            for i in 1:n_basefuncs
                v = shape_value(cv, qp, i)
                # Loop over trial shape functions
                for j in 1:n_basefuncs
                    u = shape_value(cv, qp, j)
                    # Assemble local stiffness matrix 
                    Ae[i, j] += (s * u * v) * dx
                end
            end
        end
        
        # Assemble Ae into A₂
        assemble!(assembler, celldofs(cell), Ae)
    end

    return A₂
end

"""
    assemble_A1(cv, dh, A₁, p)

The first order term in the quadratic eigenvalue problem
"""
function assemble_A1(cv::CellValues, dh::DofHandler, A₁::SparseMatrixCSC, p::PML)
    # Allocate the local stiffness matrix
    n_basefuncs = getnbasefunctions(cv)
    Ae = zeros(ComplexF64, n_basefuncs, n_basefuncs)

    # Create an assembler
    assembler = start_assemble(A₁)

    # Loop over all cells
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        # Reset local stiffness matrix to 0.0 + 0.0im 
        fill!(Ae, 0.0 + 0.0im)
        coords = getcoordinates(cell)

        # Loop over quadrature points
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            coords_qp = spatial_coordinate(cv, qp, coords)
            s = coord_transform(coords_qp[2], p)

            # Loop over test shape functions 
            for i in 1:n_basefuncs
                v = shape_value(cv, qp, i)
                for j in 1:n_basefuncs 
                    ∇u = shape_gradient(cv, qp, j)
                    # Assemble local stiffness matrix
                    Ae[i, j] += (-2im * s * ∇u[1] * v) * dx
                end
            end
        end

        # Assemble Ae into A₁
        assemble!(assembler, celldofs(cell), Ae)
    end

    return A₁     
end

"""
    assemble_A0(cv, dh, A₀, medium, p, k)

The zero order term in the quadratic eigenvalue problem
"""
function assemble_A0(cv::CellValues, dh::DofHandler, A₀::SparseMatrixCSC, medium::Function, p::PML, k)
    # Allocate the local stiffness matrix
    n_basefuncs = getnbasefunctions(cv)
    Ae = zeros(ComplexF64, n_basefuncs, n_basefuncs)

    # Create an assembler A₀
    assembler = start_assemble(A₀)
    
    # Loop over all cells
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        # Reset local stiffness matrix to 0.0 + 0.0im
        fill!(Ae, 0.0 + 0.0im)
        coords = getcoordinates(cell)
        
        # Loop over quadrature points
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            coords_qp = spatial_coordinate(cv, qp, coords)
            s = coord_transform(coords_qp[2], p)
            n = medium(coords_qp)

            # Loop over test shape functions
            for i in 1:n_basefuncs
                v = shape_value(cv, qp, i)
                ∇v = shape_gradient(cv, qp, i)
                # Loop over trial shape functions
                for j in 1:n_basefuncs
                    u = shape_value(cv, qp, j) 
                    ∇u = shape_gradient(cv, qp, j)
                    # Assemble local stiffness matrix
                    Ae[i, j] += (s*∇u[1]*∇v[1] + ∇u[2]*∇v[2]/s - (k^2)*n*s*u*v) * dx
                end
            end
        end
        
        # Assemble Ae into A₀
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A₀ 
end

function apply_all_bds!(K::Union{SparseMatrixCSC, Symmetric}, ch::ConstraintHandler)
    return apply_all_bds!(K, eltype(K)[], ch, true)
end


function apply_all_bds!(KK::Union{SparseMatrixCSC, Symmetric}, f::AbstractVector, ch::ConstraintHandler, applyzero::Bool = false)
    @assert Ferrite.isclosed(ch)
    sym = isa(KK, Symmetric)
    K = sym ? KK.data : KK
    @assert length(f) == 0 || length(f) == size(K, 1)
    @boundscheck checkbounds(K, ch.prescribed_dofs, ch.prescribed_dofs)
    @boundscheck length(f) == 0 || checkbounds(f, ch.prescribed_dofs)

    m = Ferrite.meandiag(K) # Use the mean of the diagonal here to not ruin things for iterative solver

    # Add inhomogeneities to f: (f - K * ch.inhomogeneities)
    if !applyzero
        @inbounds for i in 1:length(ch.inhomogeneities)
            d = ch.prescribed_dofs[i]
            v = ch.inhomogeneities[i]
            if v != 0
                for j in nzrange(K, d)
                    r = K.rowval[j]
                    sym && r > d && break # don't look below diagonal
                    f[r] -= v * K.nzval[j]
                end
            end
        end
        if sym
            # In the symmetric case, for a constrained dof `d`, we handle the contribution
            # from `K[1:d, d]` in the loop above, but we are still missing the contribution
            # from `K[(d+1):size(K,1), d]`. These values are not stored, but since the
            # matrix is symmetric we can instead use `K[d, (d+1):size(K,1)]`. Looping over
            # rows is slow, so loop over all columns again, and check if the row is a
            # constrained row.
            @inbounds for col in 1:size(K, 2)
                for ri in nzrange(K, col)
                    row = K.rowval[ri]
                    row >= col && break
                    if (i = get(ch.dofmapping, row, 0); i != 0)
                        f[col] -= ch.inhomogeneities[i] * K.nzval[ri]
                    end
                end
            end
        end
    end

    # Condense K (C' * K * C) and f (C' * f)
    Ferrite._condense!(K, f, ch.dofcoefficients, ch.dofmapping, sym)

    # Remove constrained dofs from the matrix
    Ferrite.zero_out_columns!(K, ch.prescribed_dofs)
    Ferrite.zero_out_rows!(K, ch.dofmapping)

    # Add meandiag to constraint dofs
    @inbounds for i in 1:length(ch.inhomogeneities)
        d = ch.prescribed_dofs[i]
        K[d, d] = m

        # Deal with the periodic boundary condition
        dofcoef = Ferrite.coefficients_for_dof(ch.dofmapping, ch.dofcoefficients, d)
        if dofcoef !== nothing
            col = dofcoef[1].first
            K[d, col] = -m
        end

        if length(f) != 0
            vz = applyzero ? zero(eltype(f)) : ch.inhomogeneities[i]
            f[d] = vz * m
        end
    end
    return
end

"""
    solve_qm(cv::CellValues, dof::DofHandler, cst::ConstraintHandler, A₀::SparseMatrixCSC, A₁::SparseMatrixCSC, A₂::SparseMatrixCSC)

Get the exceptional values (quasi-momentum α) by solving the eigenvalue 
problems generated by FEM.
"""
function solve_qm(cv::CellValues, dof::DofHandler, cst::ConstraintHandler, A₀::SparseMatrixCSC, A₁::SparseMatrixCSC, A₂::SparseMatrixCSC)
    
end