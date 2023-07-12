using BifurcationKit, Plots, LinearAlgebra, Polynomials, Setfield

norminf(x) = norm(x, Inf)

function equilibria(p)::Vector{Vector{Float64}}
    result = []
    a, b, c, d, α, ϵ, θ = p
    # σₐ = (a*d - 2*b*c + 2*sqrt(b^2*c^2 - a*b*c*d)) / a^2
    σₐ = ((sqrt(a*d - b*c) - sqrt(-b*c))/a)^2
    σ = θ * σₐ
    D₁ = (σ * par.a + par.d) / 4 * σ

    # println(σₐ)
    # println(σₛ)

    a₁ = σ
    b₁ = ϵ - b*σ/ϵ
    c₁ = d - a*σ
    d₁ = c*ϵ

    sols = real(filter(isreal, roots(Polynomial([d₁, c₁, b₁, a₁]))))

    for sol in sols
        push!(result, [ϵ, sol, (a*ϵ + b*sol - ϵ*sol^2) / 2*D₁*ϵ])
    end
    return result
end

function opinion!(du, u, p, t)
    a, b, c, d, α, ϵ, θ = p
    x, y, A = u

    # σₐ = (a * d - 2 * b * c + 2 * sqrt(b^2 * c^2 - a * b * c * d)) / a^2
    σₐ = ((sqrt(a*d - b*c) - sqrt(-b*c))/a)^2
    σ = θ * σₐ

    D₁ = (σ * a + d) / 4 * σ
    D₂ = σ * D₁

    du[1] = a * x + b * y - x * y^2 - 2 * D₁ * A * x
    du[2] = c * x + d * y + x * y^2 - 2 * D₂ * A * y
    du[3] = α * A * (1 - A) * (ϵ - x) * (ϵ + x)
    du
end

opinion(u, p) = opinion!(similar(u), u, p, 0)

par = (a=-1.1, b=-2.0, c=1.0, d=1.0, α=1.0, ϵ=0.1, θ=0.95)
# u₀ = equilibria(par)[1]
# u₀ = [0.0, 0.0, 0.0]
# u₀ = [0.1, -0.151173972786538, 0.795717305043761]
println(u₀)

recordFromSolutionOp(x, p) = (X=x[1], Y=x[2], A=x[3])
prob = BifurcationProblem(opinion, u₀, par, (@lens _.θ);
    recordFromSolution=recordFromSolutionOp)

opt_newton = NewtonPar(tol=1e-9, maxIter=2000)

opts_br = ContinuationPar(pMin=0.0, pMax=1.0, dsmin=0.00001, dsmax=0.2, ds=0.001,
    nInversion=8, maxBisectionSteps=30, maxSteps=300, detectBifurcation=3, newtonOptions=opt_newton, nev=10)

br = continuation(prob, PALC(), opts_br; bothside=true)

# diagram = bifurcationdiagram(prob, PALC(), 3, (args...) -> setproperties(opts_br; dsminBisection=1e-18); bothside=true)

scene = plot(br, vars=(:param, :A))

# plot(diagram, vars=(:param, :A))

# automatic branch switching from Hopf point
# opt_po = NewtonPar(tol = 1e-10, maxIter = 100)
# opts_po_cont = ContinuationPar(dsmin = 0.005, dsmax = 0.04, ds = 0.01, pMax = 1.0, maxSteps = 2000, newtonOptions = opt_po, nev = 11, tolStability = 1e-6,
# 	detectBifurcation = 3, dsminBisection = 1e-6, maxBisectionSteps = 15, nInversion = 4)

# br_po = continuation(br, 2, opts_po_cont, PeriodicOrbitTrapProblem(M = 51); δp=0.01, ampfactor=1)
# plot(br_po)

# diagram = bifurcationdiagram(prob, br, 3, (args...) -> setproperties(opts_po_cont))
# plot(diagram)