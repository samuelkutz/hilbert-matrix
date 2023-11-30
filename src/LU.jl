export lu_decomposition

using LinearAlgebra

function lu_decomposition(A::Matrix)
    m, n = size(A)

    if m != n
        error("Input matrix must be square.")
    end
    
    L = zeros(m, n)
    U = copy(A)
    
    for k in 1:n
        if U[k, k] == 0
            error("Zero pivot encountered. LU decomposition is not possible.")
        end
        
        # Calculate multipliers and store them in the lower triangular matrix L
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