using LinearAlgebra

function JOR(A, b, omega, tol=1e-7, max_iter=1000)
    n= length(b)
    x = zeros(n)  

    for k in 1:max_iter
        x_old = copy(x)

        for i in 1:n
            sigma = dot(A[i, :] .* x, .!=(1:n, i))
            x[i] = (omega / A[i, i]) * (b[i] - sigma) + (1 - omega) * x[i]
        end

        if norm(x - x_old) < tol
            break
        end
    end

    return x
end