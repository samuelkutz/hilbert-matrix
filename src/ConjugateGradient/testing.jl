using Plots, LinearAlgebra

include("../Hilbert.jl")
include("ConjugateGradient.jl") 

n = 100
H = Hilbert(n)
b = ones(Real, n)
max_iter = 10000
tol = 1e-7

x = zeros(Real, n)
r = b - H * x
p = copy(r)

rel_residual_norms = zeros(max_iter)
for iter in 1:max_iter
    global x, r, p

    Hp = H * p

    alpha = dot(p, r) / dot(p, Hp)

    x += alpha * p
    r -= alpha * Hp

    rel_residual_norms[iter] = norm(r) / norm(b)
    if norm(r) < tol
        println("Converged in $iter iterations.")
        break
    end

    beta = -dot(Hp, r) / dot(Hp, p)
    p = r + beta * p
end

println(rel_residual_norms[30])

plot(1:max_iter, rel_residual_norms, xlabel="n de iterações", ylabel="Norma do resíduo relativo", title="Norma do resíduo relativo vs n de iterações", label="Norma relativa do resíduo", c=:magenta, lw=2)
plot!(1:max_iter, tol*ones(max_iter), s=:dash, c=:grey, lw = 2, label="Tolerância")
savefig("./plots/ConjugateGradient/conjugate_gradient_residual_norm.png")