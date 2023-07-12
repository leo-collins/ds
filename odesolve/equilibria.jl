using Polynomials

par = (a=-1.1, b=-2.0, c=1.0, d=1.0, α=1.0, ϵ=0.01, θ=0.95)

function critical_homophily(par)::Vector{Vector{Float64}}
    result = []
    a, b, c, d, α, ϵ, θ = par
    σₐ = (a*d - 2*b*c + 2*sqrt(b^2*c^2 - a*b*c*d)) / a^2
    # σₐ = ((sqrt(a*d - b*c) - sqrt(-b*c))/a)^2
    σ = θ * σₐ
    D₁ = (σ * a + d) / 4 * σ

    a₁ = σ*ϵ
    b₁ = ϵ^2 - b*σ
    c₁ = ϵ*(d - a*σ)
    d₁ = c*ϵ^2

    sols = real(filter(isreal, roots(Polynomial([d₁, c₁, b₁, a₁]))))

    for sol in sols
        println(ϵ*σ*sol^3 + (ϵ^2 - b*σ)*sol^2 + ϵ*(d - a*σ)*sol + c*ϵ^2)
        push!(result, [ϵ, sol, (a*ϵ + b*sol - ϵ*sol^2) / 2*D₁*ϵ])
    end
    return result
end

function patterned_state(par)::Vector{Vector{Float64}}
    a, b, c, d, α, ϵ, θ = par
    # σₐ = (a*d - 2*b*c + 2*sqrt(b^2*c^2 - a*b*c*d)) / a^2
    σₐ = ((sqrt(a*d - b*c) - sqrt(-b*c))/a)^2
    σ = θ * σₐ

    D₁ = (σ * par.a + par.d) / 4 * σ
    D₂ = σ * D₁

    return [sqrt((((a - 2*D₁)*(d - 2*D₂) - b*c)*(b + d - 2*D₂)) / (2*D₁ - a - c)^2), sqrt(((a - 2*D₁)*(d - 2*D₂) - b*c) / (b + d - 2*D₂)), 1]
end

function attdyn(par, point)
    a, b, c, d, α, ϵ, θ = par
    x, y, A = point

    σₐ = (a * d - 2 * b * c + 2 * sqrt(b^2 * c^2 - a * b * c * d)) / a^2
    # σₐ = ((sqrt(a*d - b*c) - sqrt(-b*c))/a)^2
    σ = θ * σₐ

    D₁ = (σ * a + d) / 4 * σ
    D₂ = σ * D₁
    
    return [
        a * x + b * y - x * y^2 - 2 * D₁ * A * x,
        c * x + d * y + x * y^2 - 2 * D₂ * A * y,
        α * A * (1 - A) * (ϵ - x) * (ϵ + x)
    ]
end