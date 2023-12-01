function Hilbert(n::Int64)
    H = zeros(BigFloat, n, n)
    
    for i in 1:n
        for j in 1:n
            H[i, j] = 1 / (i + j- 1)
        end
    end

    return H
end