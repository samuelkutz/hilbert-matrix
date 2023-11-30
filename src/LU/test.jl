using Plots
include("../Hilbert.jl")
include("LU.jl")

tol = 1e-7
ns = 1:100
rel_residual_norms = []

for n in ns
    H = zeros(Real, n, n)
    
    for i in 1:n
        for j in 1:n
            H[i, j] = 1 / (i + j- 1)
        end
    end

    b = ones(Real, n)

    L, U = LU(H)    
    x = solve_LU(L, U, b)

    push!(rel_residual_norms,  norm(H * x - b,) / norm(b))
end

plot(xlabel="n", ylabel="Norma do resíduo relativo", title="Norma do resíduo relativo vs n")

plot!(ns, rel_residual_norms, c=:magenta, label="Norma relativa do resíduo")

plot!(ns, tol*ones(100), s=:dash, c=:grey, lw = 2, label="Tolerância")

savefig("./src/LU/residual_norm_plot.png")