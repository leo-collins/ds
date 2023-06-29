using Polynomials

par = (a=-1.1, b=-2.0, c=1.0, d=1.0, α=1.0, ϵ=0.01, θ=0.3)

σₐ = (par.a * par.d - 2 * par.b * par.c + 2 * sqrt(par.b^2 * par.c^2 - par.a * par.b * par.c * par.d)) / par.a^2
σ = par.θ * σₐ

# D₁ = (σ * par.a + par.d) / 4 * σ
# D₂ = σ * D₁

a = σ
b = par.b*σ/par.ϵ
c = par.d - par.a*σ
d = -par.c*par.ϵ

y_pol = Polynomial([d, c, b, a])
Δ = b^2*c^2 - 4*a*c^3 - 4*b^3*d - 27*a^2*d^2 + 18*a*b*c*d
println(Δ)
real(filter(isreal, roots(y_pol)))