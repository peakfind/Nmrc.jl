"""
    Pml 

Information on PML layer and coordinate transformation.
"""
struct PML
    ρ::Float64
    χ::ComplexF64
    m::Int64
    start::Float64
    δ::Float64
end

# function PML(ρ, χ, m, start, δ)
#     return PML(ρ, χ, m, start, δ)
# end


"""
    get_width(p::PML)

Get the width of the PML layer.
"""
get_width(p::PML) = p.δ


"""
    coord_transform(x₂, p::PML)

Coordinate transformation of PML method.
"""
function coord_transform(x₂, p::PML)
    if x₂ < p.start
        return 1.0
    else
        return 1.0 + p.ρ * p.χ * ((x₂ - p.start)/p.δ)^p.m
    end 
end