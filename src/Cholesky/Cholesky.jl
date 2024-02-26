using LinearAlgebra

function chol(A)
    n, m = size(A)

    if n != m
        error("A matriz precisa ser quadrada")
    end

    L = Matrix{BigFloat}(I, n, n)

    for i in 1:n
        for j in 1:i
            if i == j
                L[i, i] = sqrt(A[i, i] - sum(L[i, 1:i-1].^2))
            else
                L[i, j] = (A[i, j] - sum(L[i, 1:j-1] .* L[j, 1:j-1])) / L[j, j]
            end
        end
    end

    return L
end

function solve_chol(L, b)
    n = size(L, 1)
    y = zeros(BigFloat, n)

    for i in 1:n
        y[i] = (b[i] - sum(L[i, 1:i-1] .* y[1:i-1])) / L[i, i]
    end

    x = zeros(BigFloat, n)

    for i in n:-1:1
        x[i] = (y[i] - sum(L[i+1:end, i] .* x[i+1:end])) / L[i, i]
    end

    return x
end
