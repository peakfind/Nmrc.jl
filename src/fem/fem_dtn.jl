# Use Ferrite.jl to implement DtN-FEM

"""
    periodic_cell(;lc=0.5, period=2π, height=2.0)

Generate a mesh for the periodic cell.

# Arguments

- `lc`: the mesh size near points
- `period`: the period of the periodic cell
- `height`: the height of the periodic cell

# Points

```
 4 ------------ 3
 |              |
 |              |
 1 ------------ 2
```

# Lines

```
 . ---- l3 ---- .
 l4             l2
 . ---- l1 ---- .
```
"""
function periodic_cell(;lc=0.5, period=2π, height=2.0)
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)

    # Add the points
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(period, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(period, height, 0, lc)
    p4 = gmsh.model.geo.addPoint(0, height, 0, lc)

    # Add the lines
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    # Create the loop and the surface
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.geo.addPlaneSurface([loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Create the physical groups (for boundary conditions)
    gmsh.model.addPhysicalGroup(1, [l1], -1, "bottom") 
    gmsh.model.addPhysicalGroup(1, [l2], -1, "right")
    gmsh.model.addPhysicalGroup(1, [l3], -1, "top")
    gmsh.model.addPhysicalGroup(1, [l4], -1, "left")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Ω")

    # Set Periodic boundary condition
    gmsh.model.mesh.setPeriodic(1, [l2], [l4], [1, 0, 0, period, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])

    # Generate a 2D mesh
    gmsh.model.mesh.generate(2)

    # Save the mesh and read the .msh file by using FerriteGmsh
    grid = mktempdir() do dir 
        path = joinpath(dir, "mesh.msh")
        gmsh.write(path)
        togrid(path)
    end

    # Finalize the Gmsh library
    gmsh.finalize()

    return grid
end

"""
    setup_vals(ip)

Set up CellValues and FacetValues by using interpolation `ip`. Here we need 
to define FacetValues because the load vector and the TBC matrix contains 
integral on the boundary.
"""
function setup_vals(ip)
    qr = QuadratureRule{RefTriangle}(2)
    qr_facet = FacetQuadratureRule{RefTriangle}(2)
    cv = CellValues(qr, ip)
    fv = FacetValues(qr_facet, ip)
    
    return cv, fv
end

function setup_dh(grid::Grid, ip)
    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)
    
    return dh
end

"""
    setup_bcs(dh::DofHandler; period=2π)

Set the periodic boundary condition ("left" and "right") and 
Dirichelt boundary condition ("bottom").

# Arguments

- `dh`: DofHandler
- `period`: the period of the periodic cell, see [`periodic_cell`](@ref)
"""
function setup_bcs(dh::DofHandler; period=2π)
    cst = ConstraintHandler(dh)
    
    # Periodic boundary condition
    pfacets = collect_periodic_facets(dh.grid, "right", "left", x -> x + Vec{2}((period, 0.0)))
    pbc = PeriodicDirichlet(:u, pfacets)
    add!(cst, pbc)
    
    # Dirichlet boundary condition
    dbc = Dirichlet(:u, getfacetset(dh.grid, "bottom"), x -> 0)
    add!(cst, dbc)
    
    close!(cst)
    return cst
end

"""
    dofs_on_dtn(dh::DofHandler, field::Symbol, facetset)

Extract the global indices of Dofs associated to the artificial boundary.
"""
function dofs_on_dtn(dh::DofHandler, field::Symbol, facetset)
    # TODO: Maybe find a more natural way to extract the Dofs associated to the artificial boundary
    dtn_ch = ConstraintHandler(dh)
    dbc = Dirichlet(field, facetset, x -> 0)
    add!(dtn_ch, dbc)
    close!(dtn_ch)
    
    return dtn_ch.prescribed_dofs
end

"""
    allocate_stiff_matrix(dofhandler::DofHandler, csthandler::ConstraintHandler, dofs)

Create a sparse pattern for the stiffness matrix. We need to add extra entries 
due to the DtN map by hand.
"""
function allocate_stiff_matrix(dh::DofHandler, cst::ConstraintHandler, dofsDtN)
    sp = init_sparsity_pattern(dh)
    add_cell_entries!(sp, dh)
    # Use add_entry! for DtN term
    for i in dofsDtN, j in dofsDtN
        # if abs(i - j) > 1
            Ferrite.add_entry!(sp, i, j)
        # end
    end

    add_constraint_entries!(sp, cst)
    K = allocate_matrix(SparseMatrixCSC{ComplexF64, Int}, sp)
    return K
end

"""
    assemble_A(cv::CellValues, dh::DofHandler, A::SparseMatrixCSC, inc::Incident)

"""
function assemble_A(cv::CellValues, dh::DofHandler, A::SparseMatrixCSC, inc::Incident)
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
        
        # Reset local stiffness matrix to 0.0 + 0.0im
        fill!(Ae, 0.0 + 0.0im)

        # Loop over quadrature points
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            
            # Loop over test shape functions
            for i in 1:n_basefuncs
                v = shape_value(cv, qp, i)
                ∇v = shape_gradient(cv, qp, i)
                
                # Loop over trial shape functions
                for j in 1:n_basefuncs
                    u = shape_value(cv, qp, j)
                    ∇u = shape_gradient(cv, qp, j)
                    
                    # Assemble the local stiffness matrix
                    Ae[i, j] += (∇u ⋅ ∇v - 2im * α * ∇u[1] * v - (k^2 - α^2) * u * v) * dx
                end
            end
        end
        
        # Assemble Ae into A
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A
end

"""
    assemble_load(fv::FacetValues, dh::DofHandler, facetset, f, inc::Incident, height)

Assemble the load vector due to the incident field.
"""
function assemble_load(fv::FacetValues, dh::DofHandler, facetset, f, inc::Incident, height)
    β = get_beta(inc)
    
    # Allocate the local load vector fe
    n_basefuncs = getnbasefunctions(fv)
    fe = zeros(ComplexF64, n_basefuncs)
    
    # Loop over all facets on the specific facetset
    for facet in FacetIterator(dh, facetset)
        # Update the fv to the current facet
        reinit!(fv, facet)
        
        # Reset the local vector fe to zero
        fill!(fe, 0.0 + 0.0im)
        
        # Loop over quadrature points
        for qp in 1:getnquadpoints(fv)
            ds = getdetJdV(fv, qp)
            
            # right hand side due to the incident field
            # we can also use coords_qp[2] to replace height here
            g = -2im * β * exp(-im * β * height)
            
            # Loop over test functions 
            for i in 1:n_basefuncs
                v = shape_value(fv, qp, i)
                fe[i] += g * v * ds
            end
        end
        
        assemble!(f, celldofs(facet), fe)
    end
    
    return f
end

"""
    assemble_tbc(fv::FacetValues, dh::DofHandler, inc::Incident, F::SparseMatrixCSC, facetset, N, dofsDtN)

Assemble the TBC matrix.
"""
function assemble_tbc(fv::FacetValues, dh::DofHandler, inc::Incident, facetset, F, N, dofsDtN)
    # Allocate the vector Θ 
    Θ = sparsevec(dofsDtN, zeros(ComplexF64, length(dofsDtN)), ndofs(dh))

    # Loop over truncated terms
    for n in -N:N 
        # Reset the vector Θ to zero
        fill!(Θ, 0.0 + 0.0im)

        # Compute βₙ
        βₙ = beta_n(inc, n)

        # Compute the vector Θ (Fourier coefficients and its conjugate) 
        compute_coef!(fv, dh, facetset, Θ, n)
        
        # Assemble the TBC matrix 
        for i in Θ.nzind, j in Θ.nzind
            v = im * βₙ * Θ[i] * conj(Θ[j])/(2π)
            Ferrite.addindex!(F, v, i, j)
        end
    end
    
    # TODO: find a way to avoid compute im/2π in the above iteration
    # F .*= im/(2π)
    
    return F
end

"""
    compute_coef!(fv::FacetValues, dh::DofHandler, θ::SparseVector, facetset, n)

Compute Θⁿ on the `facetset`. Actually the computation of the TBC matrix reduces 
to the computation of the vector Θⁿ.
"""
function compute_coef!(fv::FacetValues, dh::DofHandler, facetset, Θ::SparseVector, n)
    # Allocate the local vector θ
    n_basefuncs = getnbasefunctions(fv)
    θ = zeros(ComplexF64, n_basefuncs)
    
    # Loop over all facets on the specific facetset
    for facet in FacetIterator(dh, facetset)
        # Update the fv to the correct facet
        reinit!(fv, facet)
        
        # Reset the local vector θ to zero
        fill!(θ, 0.0 + 0.0im)
        
        coords = getcoordinates(facet)
        
        # Loop over quadrature points
        for qp in 1:getnquadpoints(fv)
            ds = getdetJdV(fv, qp)

            # Coordinate of the quadrature point
            coords_qp = spatial_coordinate(fv, qp, coords)
            
            # Modes: eⁱⁿˣ
            mode = exp(im * n * coords_qp[1])
            
            for i in 1:n_basefuncs
                ϕ = shape_value(fv, qp, i)
                θ[i] += ϕ * mode * ds
            end
        end
        
        assemble!(Θ, celldofs(facet), θ)
    end
    
    return Θ
end

"""
    sub_preserve_structure(A::SparseMatrixCSC, B::SparseMatrixCSC)

Subtract two sparse matrix without changing the sparse pattern. `A` and `B` 
should have the same sparse pattern.
"""
function sub_preserve_structure(A::SparseMatrixCSC, B::SparseMatrixCSC)
    # TODO: Find a better way. It relates to the functions `allocate_stiff_matrix` and `assemble_tbc`.
    if A.colptr != B.colptr || A.rowval != B.rowval || size(A) != size(B)
        error("Matrices must have the same sparsity structure and dimensions.")
    end

    new_nzval = A.nzval - B.nzval
    return SparseMatrixCSC(size(A,1), size(A,2), A.colptr, A.rowval, new_nzval)
end