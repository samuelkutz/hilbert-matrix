using LinearAlgebra

include("./src/Hilbert.jl")
include("./src/LU.jl")

H = Hilbert(20)

show(stdout, "text/plain" , H)

println("\nH is symmetric: ", H' == H)

# grafico do numero de condicionamento de H as n grows

L, U = lu_decomposition(H)


println("\nL: ")
show(stdout, "text/plain" , L)
println("\nU: ")
show(stdout, "text/plain" , U)