using LinearAlgebra

function LU(A::Matrix)
    m, n = size(A)

    if m != n
        error("A matriz precisa ser quadrada.")
    end
    
    L = zeros(m, n)
    U = copy(A)
    
    for k in 1:n
        if U[k, k] == 0
            error("Pivô nulo encontrado, a decomposição não é possível.")
        end
        
        L[k:end, k] = U[k:end, k] / U[k, k]
        
        for i in k+1:n
            L[i, k] = U[i, k] / U[k, k]
            U[i, k:end] -= L[i, k] * U[k, k:end]
        end
    end
    
    for i in 1:m
        L[i, i] = 1.0
    end
    
    return L, U
end

function solve_LU(L, U, b)
    n = size(b, 1)
    y = zeros(Real, n)

    for i in 1:n
        y[i] = b[i] - dot(L[i, 1:i-1], y[1:i-1])
    end

    x = zeros(n)

    for i in n:-1:1
        x[i] = (y[i] - dot(U[i, i+1:end], x[i+1:end])) / U[i, i]
    end

    return x
end