"""
    csqrt_negimag(z::Complex{T}) where T<:AbstractFloat

Square root for complex numbers which takes the negative imaginary axis as the branch cut.
Note that `sqrt` for complex number in julia takes the negative real axis as the default
branch cut.
"""
function csqrt_negimag(z::Complex{T}) where T<:AbstractFloat
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
csqrt_negimag(z::Complex) = csqrt_negimag(float(z))

"""
    beta_n(k, α, n)

Compute ``\\beta_{n}`` for complex ``\\alpha``. We need this function when 
we construct nonlinear eigenvalue problems. Please note that we take the 
negative imaginary axis as the branch cut for square root. So we use `csqrt_negimag`  
instead of the default square root in Julia. See also [`csqrt_negimag`](@ref).
"""
function beta_n(k, α, n)
    αₙ = α + complex(n)
    βₙ = csqrt_negimag(complex(k)^2 - αₙ^2)
    
    return βₙ
end

"""
    integral_uv(cv::CellValues, dh::DofHandler, u, v)

Compute the integral 
```math
\\int u v dx, 
```
where ``u`` and ``v`` are vectors indexed by degrees of freedom. Normally, they 
are results obtained by FEM.
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

Compute the integral 
```math
\\int \\frac{\\partial u}{\\partial x_{1}} v dx,
```
where ``u`` and ``v`` are vectors indexed by degrees of freedom. Normally, they 
are results obtained by FEM.
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

"""
    integral_uv̄(cv::CellValues, dh::DofHandler, u, v)

Compute the integral
```math
\\int u \\bar{v} dx, 
```
where ``\\bar{v}`` is the conjugate of the function ``v``. See also [`integral_nuv̄`](@ref) 
and [`integral_∂₁uv̄`](@ref).
"""
function integral_uv̄(cv::CellValues, dh::DofHandler, u, v)
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
            rst += u_qp * conj(v_qp) * dx
        end
    end
    
    return rst
end

"""
    integral_nuv̄(cv::CellValues, dh::DofHandler, n::Function, u, v)

Compute the integral
```math
\\int n(x_{1}, x_{2}) u v dx, 
```
where ``n(x_{1}, x_{2})`` is a function defined in our computational domain 
(normally, the refraction index). See also [`integral_uv̄`](@ref) and [`integral_∂₁uv̄`](@ref).
"""
function integral_nuv̄(cv::CellValues, dh::DofHandler, n::Function, u, v)
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
        
        # Get the coordinates of this cell
        coords = getcoordinates(cell)

        # Loop over quadrature points
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            
            # Get the coordinate of the quadrature point 
            coords_qp = spatial_coordinate(cv, qp, coords)

            # Evaluate the function value of n
            ri = n(coords_qp)

            u_qp = function_value(cv, qp, u, cd)
            v_qp = function_value(cv, qp, v, cd)
            rst += ri * u_qp * conj(v_qp) * dx
        end
    end
    
    return rst
end

"""
    integral_∂₁uv̄(cv::CellValues, dh::DofHandler, u, v)

Compute the integral
```math
\\int \\frac{\\partial u}{\\partial x_{1}} \\bar{v} dx, 
```
where ``\\bar{v}`` is the conjugate of the function ``v``. See also [`integral_uv̄`](@ref) 
and [`integral_nuv̄`](@ref).
"""
function integral_∂₁uv̄(cv::CellValues, dh::DofHandler, u, v)
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
            rst += ∂₁u_qp * conj(v_qp) * dx
        end
    end
    
    return rst
end

"""
    compute_rayleigh_n(fv::FacetValues, dh::DofHandler, facetset, u, n)

Compute the ``n``-th order Rayleigh coefficient for the solution ``u``. Be careful that 
we need to compute the Rayleigh coefficients for scattered solutions. 
```math
u_{n} = \\int_{\\text{boundary}} u e^{-inx_{1}} ds.
```
"""
function compute_rayleigh_n(fv::FacetValues, dh::DofHandler, facetset, u, n)
    rst = zero(eltype(u))
    
    for facet in FacetIterator(dh, facetset)
        reinit!(fv, facet)
        cd = celldofs(facet)
        coords = getcoordinates(facet)
        
        for qp in 1:getnquadpoints(fv)
            ds = getdetJdV(fv, qp)
            coords_qp = spatial_coordinate(fv, qp, coords)
            e = exp(-im * n * coords_qp[1])
            u_qp = function_value(fv, qp, u, cd)
            rst += u_qp * e * ds
        end
    end
    
    return rst
end