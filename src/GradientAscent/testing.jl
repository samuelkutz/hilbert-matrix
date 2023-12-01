using Plots, LinearAlgebra

include("../Hilbert.jl")
include("GradientAscent.jl") # gradient_ascent()

n = 100
H = Hilbert(n)
b = ones(Real, n)
max_iter = 500
tol = 1e-7

x = zeros(Real, n)
r = b - H * x  

rel_residual_norms = zeros(max_iter)

for iter in 1:max_iter
    global x, r

    alpha = dot(r, r) / dot(r, H * r)
    x += alpha * r
    r -= alpha * H * r

    rel_residual_norms[iter] = norm(r) / norm(b)

    if norm(r) < tol
        println("Converged in $iter iterations.")
        break
    end
end


plot(1:max_iter, rel_residual_norms, xlabel="n de iterações", ylabel="Norma do resíduo relativo", title="Norma do resíduo relativo vs n de iterações", label="Norma relativa do resíduo", c=:magenta, lw=1)
savefig("./src/GradientAscent/gradient_ascent_residual_norm.png")