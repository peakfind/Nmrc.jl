# This file constains the functions which assemble sparse matrices of quadratic 
# eigenvalue problems from the computation of exceptional values for closed or 
# open waveguides.

"""
    assemble_A0(cv::CellValues, dh::DofHandler, A₀::SparseMatrixCSC, medium::Function, k)

Assemble the zero order term ``\\mathbf{A}_{0}`` in the Nonlinear eigenvalue problem.

``A_{0}(u, v) = \\int \\nabla u \\cdot \\nabla \\bar{v} - k^{2}n(x_{1}, x_{2}) u \\bar{v} d x.``

# Argument

- `cv`: CellValues
- `dh`: DofHandler
- `A₀`: an empty sparse pattern preallocated for A₀
- `medium`: refractive index function which describes the properties of the medium
- `k`: the wavenumber

See also [`assemble_A1`](@ref), [`assemble_A2`](@ref).
"""
function assemble_A0(cv::CellValues, dh::DofHandler, A₀::SparseMatrixCSC, medium::Function, k)
    # Allocate the local matrix
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
                    Ae[i, j] += (∇u ⋅ ∇v - (k^2) * n * u * v) * dx
                end
            end
        end
        
        # Assemble Ae into A₀ 
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A₀
end

"""
    assemble_A1(cv::CellValues, dh::DofHandler, A₁::SparseMatrixCSC)

Assemble the first order term ``\\mathbf{A}_{1}`` in the Nonlinear eigenvalue problem.

``A_{1}(u, v) = -\\int 2i \\partial_{1} u \\bar{v} d x.``

See also [`assemble_A0`](@ref), [`assemble_A2`](@ref).
"""
function assemble_A1(cv::CellValues, dh::DofHandler, A₁::SparseMatrixCSC)
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
        
        # Loop over quadrature points 
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            
            # Loop over test shape functions 
            for i in 1:n_basefuncs 
                v = shape_value(cv, qp, i)

                # Loop over trial shape functions 
                for j in 1:n_basefuncs
                    ∇u = shape_gradient(cv, qp, j)
                    
                    # Assemble local stiffness matrix 
                    Ae[i, j] += (-2im * ∇u[1] * v) * dx
                end
            end
        end
        
        # Assemble Ae into A₁
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A₁
end

"""
    assemble_A2(cv::CellValues, dh::DofHandler, A₂::SparseMatrixCSC)

Assemble the second order term in the Nonlinear eigenvalue problem.

``A_{2}(u, v) = \\int u \\bar{v} d x.``

See also [`assemble_A0`](@ref), [`assemble_A1`](@ref).
"""
function assemble_A2(cv::CellValues, dh::DofHandler, A₂::SparseMatrixCSC)
    # Allocate the local matrix
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
        
        # Loop over quadrature points 
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            
            # Loop over test shape functions 
            for i in 1:n_basefuncs
                v = shape_value(cv, qp, i)
                
                # Loop over trial shape functions 
                for j in 1:n_basefuncs
                    u = shape_value(cv, qp, j)
                    
                    # Assemble local stiffness matrix 
                    Ae[i, j] += (u * v) * dx
                end
            end
        end
        
        # Assemble Ae into A₂ 
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A₂
end 