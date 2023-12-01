using Plots
using LinearAlgebra
include("../Hilbert.jl")
include("Cholesky.jl")

n = 100
rel_residual_norms = []

#  [1.0e-10, 1.29e-10, 1.67e-10, ...]
lambda_values = 10 .^ range(-10, stop = 0, length = 100)

H = Hilbert(n)

b = ones(Real, n)

min_residual = Inf 
best_lambda = NaN

for lambda in lambda_values
    global min_residual, best_lambda

    H_reg = H + lambda * I

    L = chol(H_reg) 
    x = solve_chol(L, b)

    residual_norm = norm(H * x - b) / norm(b)
    push!(rel_residual_norms,  residual_norm)

    if residual_norm < min_residual
        min_residual = residual_norm
        best_lambda = lambda
    end
end

best_lambda = round(best_lambda, digits=4)

println("Best lambda: ", best_lambda)
println("Minimal residual: ", min_residual)

plot(lambda_values, rel_residual_norms, xlabel="Lambda", ylabel="Norma do resíduo relativo",
    title="Norma do resíduo relativo vs Lambda", label="Norma relativa do resíduo", c=:magenta, lw=1)

scatter!([best_lambda], [min_residual], markersize=6, c=:cyan, marker=:diamond, label="Melhor Lambda: $best_lambda")

savefig("./src/Cholesky/residual_norm_vs_lambda.png")