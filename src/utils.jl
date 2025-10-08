"""
    my_sqrt(z::Complex{T}) where T<:AbstractFloat

Square root for complex numbers which takes the negative imaginary axis as the branch cut.
Note that `sqrt` for complex number in julia takes the negative real axis as the default
branch cut.
"""
function my_sqrt(z::Complex{T}) where T<:AbstractFloat
    # Handle special case of zero input 
    if z == zero(z)
        return zero(z)
    end
    
    # Rotate by -π/2 (multiply by -im)
    z_rot = Complex(imag(z), -real(z))
    
    # Use the default sqrt in julia
    sqrt_rot = sqrt(z_rot)
    
    # Rotate by π/4 (multply by exp(im * π/4))
    e = 1/sqrt(2)
    r = (real(sqrt_rot) - imag(sqrt_rot)) * e
    i = (real(sqrt_rot) + imag(sqrt_rot)) * e
    rst = Complex(r, i)
    
    return rst
end

# Convenience method for other number types
my_sqrt(z::Complex) = my_sqrt(float(z))

"""
    integral_uv(cv::CellValues, dh::DofHandler, u, v)

Compute the integration of the production of function ``u``` and ``v``. Here 
``u`` and ``v`` are vectors indexed by DoFs. Normally, they are results obtained 
by FEM.
"""
function integral_uv(cv::CellValues, dh::DofHandler, u, v)
    # Check if the elements in u and v have the same type
    (eltype(u) == eltype(v)) || throw(ArgumentError("elements in u and v should have the same type!")) 
    
    rst = zero(eltype(u))
    
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        
        # Get the corresponding cell id of the current cell cache
        ci = cellid(cell)
        
        # Get the DoFs of the cell `ci`
        cd = celldofs(dh, ci)

        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            u_qp = function_value(cv, qp, u, cd)
            v_qp = function_value(cv, qp, v, cd)
            rst += u_qp * v_qp * dx
        end
    end
    
    return rst
end

"""
    integral_∂₁uv(cv::CellValues, dh::DofHandler, u, v)

Compute the integration of the product of ``\\partial_{1}u`` (the partial derivative 
of function ``u``) and ``v``. Here ``u`` and ``v`` are vectors indexed by DoFs. Normally, 
they are results obtained by FEM.
"""
function integral_∂₁uv(cv::CellValues, dh::DofHandler, u, v)
    # Check if the elements in u and v have the same type
    (eltype(u) == eltype(v)) || throw(ArgumentError("elements in u and v should have the same type!"))  

    rst = zero(eltype(u))
    
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cv, cell)
        
        # Get the corresponding cell id of the current cell cache
        ci = cellid(cell)
        # Get the DoFs of the cell `ci`
        cd = celldofs(dh, ci)

        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            ∂₁u_qp = function_gradient(cv, qp, u, cd)[1]
            v_qp = function_value(cv, qp, v, cd)
            rst += ∂₁u_qp * v_qp * dx
        end
    end
    
    return rst
end