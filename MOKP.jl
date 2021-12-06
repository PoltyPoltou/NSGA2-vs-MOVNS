struct _bi01IP
    C::Matrix{Int} # objective functions, k=1..2, j=1..n
    A::Matrix{Int} # matrix of constraints, i=1..m, j=1..n
    b::Vector{Int} # right-hand side, i=1..m
end
