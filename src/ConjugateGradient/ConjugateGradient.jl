using LinearAlgebra

function conjugate_gradient(A, b, tol=1e-7, max_iter=1000)
    n = size(A, 1)
    x = zeros(Real, n)
    r = b - A * x
    p = copy(r)

    for iter in 1:max_iter
        Ap = A * p

        alpha = dot(p, r) / dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap

        if norm(r) < tol
            println("Converged in $iter iterations.")
            break
        end

        beta = -dot(Ap, r) / dot(Ap, p)
        p = r + beta * p
    end

    return x
end
