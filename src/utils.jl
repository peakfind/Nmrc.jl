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