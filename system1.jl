using BifurcationKit, Plots, LinearAlgebra

norminf(x) = norm(x, Inf)

function opinion!(du, u, p, t)
    a, b, c, d, α, ϵ, θ = p
    x, y, A = u

    σₐ = (a*d - 2*b*c + 2*sqrt(b^2*c^2 - a*b*c*d)) / a^2
    σ = θ*σₐ

    D₁ = (σ*a + d) / 4*σ
    D₂ = σ*D₁

    du[1] = a*x + b*y - x*y^2 -2*D₁*A*x
    du[2] = c*x + d*y + x*y^2 - 2*D₂*A*y
    du[3] = α*A*(1 - A)*(ϵ - x)*(ϵ + x)
    du
end

opinion(u, p) = opinion!(similar(u), u, p, 0)

par = (a=-1.1, b=-2.0, c=1.0, d=1.0, α=1.0, ϵ=0.01, θ=0.5)
u₀ = [0.0, 0.0, 0.0]
        
recordFromSolutionOp(x, p) = (X=x[1], Y=x[2], A=x[3])
prob = BifurcationProblem(opinion, u₀, par, (@lens _.θ);
    recordFromSolution=recordFromSolutionOp)

opt_newton = NewtonPar(tol=1e-9, maxIter=200)

opts_br = ContinuationPar(pMin=0.0, pMax=1.0, dsmin=0.0002, dsmax=0.015, 
    nInversion=6, maxBisectionSteps=25, maxSteps=200, detectBifurcation=3, newtonOptions=opt_newton, nev=10)

# br = continuation(prob, PALC(), opts_br; bothside=true, normC=norminf)

diagram = bifurcationdiagram(prob, PALC(), 2, (args...) -> setproperties(opts_br); bothside=true)

# scene = plot(br, vars=(:param, :A))

plot(diagram, vars=(:param, :A))