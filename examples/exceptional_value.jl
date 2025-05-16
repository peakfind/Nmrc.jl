# Analytic function for computing exceptional value
# Reference: 
# A. Kirsch & A. Lechleiter, The limiting absorption principle 
# and a radiation condition for the scattering by a periodic layer. SIAM MA.

using CairoMakie
using Roots

"""
    fₒ(ω; k, n, h)

Compute the exceptional values by solving the zeros of fₒ in 
the interval (k, k√n). This function is used in the preprint
"Computation of the exceptional values for an open waveguide"
by A. Kirsch and R. Zhang
"""
function fₒ(ω; k, n, h)
    if ω < k || ω > k*sqrt(n) 
        throw(ArgumentError("ω should be in the interval ((k, k√n))"))
    end

    # Common constants
    c1 = sqrt(ω^2 - k^2)
    c2 = sqrt(n*(k^2) - ω^2)

    return c1*sin(c2*h)/c2 + cos(c2*h)
end

k = 4.1        # wavenumber
n = 3.9        # refractive index
h₀ = 1         # height of the periodic layer
ωₗ = k         # left bound of ω
ωᵣ = k*sqrt(n) # right bound of ω

# Compute the zeros by Roots.jl
fo(ω) = fₒ(ω; k = k, n = n, h = h₀)
ω = range(ωₗ, ωᵣ, length = 200)
fig = Figure()
ax1 = Axis(fig[1, 1], title = L"The function $f_{o}(\omega)$ for a open waveguide")
lines!(ax1, ω, fₒ.(ω; k = k, n = n, h = h₀))
fig

# Show the zeros
@show r₁ = find_zero(fo, (7, 8), Bisection())
@show r₂ = find_zero(fo, (5.5, 6.5), Bisection())

#----------------------------------------------------------
# Another formulation for function f
#----------------------------------------------------------
"""
    f(ω; k, n, h)

Compute the exceptional values by solving the zeros of f in 
the interval (k, k√n)
"""
function f(ω; k, n, h)
    if ω < k || ω > k*sqrt(n) 
        throw(ArgumentError("ω should be in the interval ((k, k√n))"))
    end

    # Common constants
    c1 = sqrt(ω^2 - k^2)
    c2 = sqrt(n*(k^2) - ω^2)

    return c1*sin(c2*h) + c2*cos(c2*h)
end

# ax2 = Axis(fig[1, 2], title = L"The function $f(\omega)$ for a open waveguide")
# lines!(ax2, ω, f.(ω; k = k, n = n, h = h₀))