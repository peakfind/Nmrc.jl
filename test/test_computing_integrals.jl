using Ferrite

function construct_fespace(n)
    # Generate the grid
    grid = generate_grid(Triangle, (n, n))
    
    # Use P1 element
    ip = Lagrange{RefTriangle, 1}()
    
    # Define the quadrature rule
    qr = QuadratureRule{RefTriangle}(2)
    
    # Construct CellValues
    cv = CellValues(qr, ip)
    
    # Construct DofHandler
    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)
    
    return cv, dh
end

cv, dh = construct_fespace(1000)

@testset "Integrate the production of u and v" begin
    u = ones(ndofs(dh))
    v = ones(ndofs(dh))
    rst = integral_uv(cv, dh, u, v)
    
    @test rst ≈ 4.0 atol=1e-9
end

@testset "Integrate the production of ∂₁u and v" begin
    u = ones(ndofs(dh))
    v = ones(ndofs(dh))
    rst = integral_∂₁uv(cv, dh, u, v)
    
    @test rst ≈ 0.0 atol=1e-10
end

@testset "Integrate the production of u and v̄" begin
    u = ones(ComplexF64, ndofs(dh))
    v = ones(ComplexF64, ndofs(dh)) 
    rst = integral_uv̄(cv, dh, u, v)
    
    @test rst ≈ 4.0 + 0.0im atol=1e-9
end

@testset "Integrate the production of n, u and v̄" begin

    function medium(x)
        if x[2] < 0.0
            return 2.0
        else
            return 1.0
        end
    end

    u = ones(ComplexF64, ndofs(dh))
    v = ones(ComplexF64, ndofs(dh)) 
    rst = integral_nuv̄(cv, dh, medium, u, v)
    
    @test rst ≈ 6.0 + 0.0im atol=1e-9
end

@testset "Integrate the production of ∂₁u and v̄" begin
    u = ones(ComplexF64, ndofs(dh))
    v = ones(ComplexF64, ndofs(dh)) 
    rst = integral_∂₁uv̄(cv, dh, u, v)
    
    @test rst ≈ 0.0 + 0.0im atol=1e-10
end