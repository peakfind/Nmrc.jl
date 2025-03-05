
#----------------------------------------------------------
# Eigenvalue parameter scaling
#----------------------------------------------------------

struct Scaling{T}
    γ::T 
    δ::T
end

get_scaling_factors(s::Scaling) = s.γ, s.δ

"""
    compute_scaling(A₀::M, A₁::M, A₂::M; scaling=:nothing::Symbol) where {T,M<:AbstractMatrix{T}}

TBW
"""
function compute_scaling(A₀::M, A₁::M, A₂::M; scaling=:default::Symbol) where {T,M<:AbstractMatrix{T}}
    scaling == :nothing && return Scaling(T(1), T(1))
    
    # Compute 2-norm of matrices
    n₀ = svdsolve(A₀)[1][1]
    n₁ = svdsolve(A₁)[1][1]
    n₂ = svdsolve(A₂)[1][1]

    τ = n₁/sqrt(n₀*n₂)

    if scaling == :default
        τ >= 10 && return return Scaling(T(1), T(1))

        γ = sqrt(n₀/n₂)
        δ = 2/(n₀ + γ*n₁)
    elseif scaling == :flvs
        γ = sqrt(n₀/n₂)
        δ = 2/(n₀ + γ*n₁)
    elseif scaling == :lts
        γ = n₁/n₂
        δ = 1/max(n₂*γ^2, n₁*γ, n₀)
    elseif scaling == :sts
        γ = n₀/n₁
        δ = 1/max(n₂*γ^2, n₁*γ, n₀)
    else
        return error("Unspported type of scaling!")
    end

    return Scaling(T(γ), T(δ))
end