using LinearAlgebra

function JOR(A, b, omega, max_iter, tol)
    n= length(b)
    x = zeros(n)  # Initial guess

    for k in 1:max_iter
        x_old = copy(x)

        for i in 1:n
            sigma = dot(A[i, :] .* x, .!=(1:n, i))
            x[i] = (omega / A[i, i]) * (b[i] - sigma) + (1 - omega) * x[i]
        end

        # Check for convergence
        if norm(x - x_old) < tol
            return x
        end
    end

    return x
end