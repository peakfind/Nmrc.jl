"""
    Pml 

Information on PML layer and coordinate transformation.
"""
struct Pml
    ρ::Float64
    χ::ComplexF64
    m::Int64
    start::Float64
    width::Float64
end

function Pml(ρ, χ, m, start, width)
    return Pml(ρ, χ, m, start, width)
end


"""
    get_width(::Pml)

Get the width of the PML layer.
"""
get_width(p::Pml) = p.width


"""
    coord_transform()

Coordinate transformation of PML method.
"""
function coord_transform(x, p::Pml)
    if x < p.start
        return 1.0
    else
        return 1.0 + p.ρ*p.χ*((x - p.start)^p.m)/(p.width^p.m)
    end 
end