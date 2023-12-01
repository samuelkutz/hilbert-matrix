using Plots

include("../Hilbert.jl")
include("SOR.jl")  # Assuming you have a file named "SOR.jl" with the SOR method implementation

n = 100
H = Hilbert(n)
b = ones(Real, n)
rel_residual_norms = []
tol = 1e-7
omega_values = 0.4:0.002:0.5

min_residual = Inf 
best_omega = NaN

for omega in omega_values
    global min_residual, best_omega

    x = SOR(H, b, omega)
    residual_norm = norm(H * x - b) / norm(b)
    push!(rel_residual_norms, residual_norm)

    if residual_norm < min_residual
        min_residual = residual_norm
        best_omega = omega
    end
end

println("Best omega: ", best_omega)
println("Minimal residual: ", min_residual)

plot(omega_values, rel_residual_norms, xlabel="Ômega", ylabel="Norma do resíduo relativo", title="Norma do resíduo relativo vs Ômega", label="Norma relativa do resíduo", c=:magenta, lw=1)
scatter!([best_omega], [min_residual], markersize=6, c=:cyan, marker=:diamond, label="Melhor Ômega: $best_omega")
savefig("./src/SOR/SOR_residual_norm_vs_omega.png")