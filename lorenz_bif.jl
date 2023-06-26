using BifurcationKit, Plots, LinearAlgebra

norminf(x) = norm(x, Inf)

function lorenz!(du, u, p, t)
    a, b, F, G = p
    x, y, z = u

    du[1] = -y^2 - z^2 - a * x + a * F
    du[2] = x * y - b * x * z - y + G
    du[3] = b * x * y + x * z - z
    du
end

lorenz(u, p) = lorenz!(similar(u), u, p, 0)

par = (a=0.25, b=4.0, F=0.3, G=0.5)
u0 = [0.0, 0.0, 0.0]

recordFromSolutionLor(x, p) = (X=x[1], Y=x[2], Z=x[3])
prob = BifurcationProblem(lorenz, u0, par, (@lens _.F);
    recordFromSolution=recordFromSolutionLor)

opt_newton = NewtonPar(tol=1e-9, maxIter=200)

opts_br = ContinuationPar(pMin=0.0, pMax=2.0, 
    nInversion=6, 
    maxBisectionSteps=25, 
    nev=4, 
    maxSteps=200,
    newtonOptions=opt_newton, 
    dsmin=0.0002, 
    dsmax=0.015,
    ds=0.01,
    detectBifurcation=3)

println("Starting continuation:")
br = continuation(prob, PALC(), opts_br; 
    bothside=true, verbosity=3)

scene = plot(br, vars=(:X, :Y), ylims=(-1.0, 1.0))

hopf_cont = continuation(br, 4, (@lens _.G), ContinuationPar(opts_br); bothside=true)

plot(hopf_cont)

# br_po = continuation(br, 4, opts_br, PeriodicOrbitTrapProblem(M=51); δp=nothing, ampfactor=1)

# plot(br_po)