@testset "Complex numbers in quadrants I, II, and IV" begin
    # In this testset, csqrt_negimag() should have the same results as sqrt()
    
    # Complex numbers in quadrants I, II, and IV
    cnumbers = [0.0 + 0.0im, 0.0 + 1.0im, -1.0 + 0.0im, 1.0 - 1.0im]
    
    # Compute the square root by julia's sqrt()
    sqrt_julia = sqrt.(cnumbers)
    
    # Compute the square root by csqrt_negimag()
    sqrt_nmrc = csqrt_negimag.(cnumbers)
    
    # Test the real part
    @test real(sqrt_nmrc) ≈ real(sqrt_julia) atol=1e-15

    # Test the imaginary part
    @test imag(sqrt_nmrc) ≈ imag(sqrt_julia) atol=1e-15
end

@testset "Complex numbers in quadrant III" begin
    # Compute the square root by csqrt_negimag()
    sqrt_nmrc = csqrt_negimag(-√2 - √2im)
    
    # Compute the square root by Euler formula
    r = cos(5π / 8) * √2
    i = sin(5π / 8) * √2
    
    # Test the real part
    @test real(sqrt_nmrc) ≈ r atol=1e-15
    
    # Test the imaginary part
    @test imag(sqrt_nmrc) ≈ i atol=1e-15
end