using BifurcationKit, Setfield, Plots

F(x, p) = [x[1] * (p.μ - x[1])]

par = (μ = -0.2,)

prob = BifurcationProblem(F, [0.1], par, (@lens _.μ))

opts_br = ContinuationPar(dsmax=0.05, ds=0.01, detectBifurcation=3, nev=2, pMin=-0.5, pMax=1.0)

br = continuation(prob, PALC(), opts_br)