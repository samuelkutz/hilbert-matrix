using LinearAlgebra

function SOR(A, b, omega, tol=1e-7, max_iter=1000)
    n = length(b)
    x = zeros(n) 

    for iter in 1:max_iter
        x_old = copy(x)
        
        for i in 1:n
            sum1 = dot(A[i, 1:i], x[1:i])  
            sum2 = dot(A[i, i+1:end], x_old[i+1:end])  
            
            x[i] = (omega / A[i, i]) * (b[i] - sum1 - sum2) + (1 - omega) * x_old[i]
        end

        if norm(x - x_old) < tol
            break
        end
    end


    return x
end