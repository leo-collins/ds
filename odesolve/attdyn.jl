using DifferentialEquations, Plots, Polynomials; plotly()

function opinion!(du, u, p, t)
    a, b, c, d, α, ϵ, θ = p
    x, y, A = u

    σₐ = (a * d - 2 * b * c + 2 * sqrt(b^2 * c^2 - a * b * c * d)) / a^2
    σ = θ * σₐ

    D₁ = (σ * a + d) / 4 * σ
    D₂ = σ * D₁

    du[1] = a * x + b * y - x * y^2 - 2 * D₁ * A * x
    du[2] = c * x + d * y + x * y^2 - 2 * D₂ * A * y
    du[3] = α * A * (1 - A) * (ϵ - x) * (ϵ + x)
end


tspan = (0.0, 100)
p = [-1.1, -2.0, 1.0, 1.0, 1.0, 0.01, 0.3]
u₀ = [0.0, 1.0, 0.1]

prob = ODEProblem(opinion!, u₀, tspan, p)
sol = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-8, abstol=1e-8)
# sol = solve(prob)

scene = plot(sol, idxs=(1, 2, 3), title="u₀=[$(u₀[1]), $(u₀[2]), $(u₀[3])]", xlabel='x', ylabel='y', zlabel='A')