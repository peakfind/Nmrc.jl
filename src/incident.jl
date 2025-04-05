"""
    Incident

Some parameters related to the incident field.

# Parameters

- `k`: the wave number
- `θ`: the incident field
- `α`: α = k * sin(θ)
- `β`: β = k * cos(θ)
"""
struct Incident
    k::Float64
    θ::Float64
    α::Float64
    β::Float64
end

"""
    Incident(k::Float64, θ::Float64)

Define the incident field with the wavenumber `k` and the incident field `θ`.
"""
function Incident(k::Float64, θ::Float64)
    α = k * sin(θ)
    β = k * cos(θ)

    return Incident(k, θ, α, β)
end

function get_alpha(inc::Incident)
    return inc.α
end

function get_beta(inc::Incident)
    return inc.β
end

"""
    beta_n(inc::Incident, n)

Compute the `βₙ`.
"""
function beta_n(inc::Incident, n)
    αₙ = inc.α + n
    
    if k > abs(αₙ)
        βₙ = complex(sqrt(k^2 - αₙ^2))
    else
        βₙ = im * sqrt(αₙ^2 - k^2)
    end
    
    return βₙ
end