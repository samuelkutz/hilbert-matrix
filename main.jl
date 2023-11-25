include("./src/Hilbert.jl")

H = Hilbert(3)

show(stdout, "text/plain" , H)

println("\nH is symmetric: ", H' == H)

# grafico do numero de condicionamento de H as n grows
