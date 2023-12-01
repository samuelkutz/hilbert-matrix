using LinearAlgebra

function gradient_ascent(A, b, tol=1e-7, max_iter=1000)
    n = size(A, 1)
    x = zeros(Real, n)
    r = b - A * x

    for iter in 1:max_iter
        alpha = dot(r, r) / dot(r, A * r)
        x += alpha * r
        r -= alpha * A * r

        if norm(r) < tol
            println("Converged in $iter iterations.")
            break
        end
    end

    return x
end