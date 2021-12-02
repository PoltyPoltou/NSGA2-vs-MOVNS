struct _bi01IP
    C  :: Matrix{Int} # objective functions, k=1..2, j=1..n
    A  :: Matrix{Int} # matrix of constraints, i=1..m, j=1..n
    b  :: Vector{Int} # right-hand side, i=1..m
end
mutable struct solution
    val_objectif::Vector{Int64}
    sol::Vector{Int64}
    cout::Vector{Int64}
end

C = zeros(2,3)
A = zeros(3,3)
b = zeros(3)
petite_solution= _bi01IP(C,A,b)
# ----------------------------------------------------------------
# Sélectionne les points pour former l'ensemble bornant primal
function kung(E::Array{solution,1}, EPrime::Array{solution,1})
    S = union(E, EPrime)
    sort!(S, by = x -> x.val_objectif[1])
    SN=[] ; push!(SN, S[1]) ; minYFeas = S[1].val_objectif[2]
    for i=2:length(E)
        if S[i].val_objectif[2] < minYFeas
            push!(SN, S[i]) ; minYFeas = S[i].val_objectif[2]
        end
    end
    return SN
end

function verification(prob::_bi01IP, x::solution)
    poids = zeros(size(prob.A,1))
    objectif = [0,0]
    for i in 1:size(prob.A,1)
        for j in 1:size(prob.A,2)
            poids[i] += prob.A[i,j] * x[j]
        end
    end
    for i in 1:size(prob.C,2)
        objectif[1] += prob.C[1,i]
        objectif[2] += prob.C[2,i]
    end
    if poids[1] <= prob.b[1] && poids[2] <= prob.b[2]
        x.val_objectif = objectif
        x.cout = poids
        return true
    else
        return false
    end
end

function no_dominated(x::solution, E::Array{solution,1})
    pos = 0
    while x.val_objectif[1] >= E[pos+1].val_objectif[1]
        pos +=1
    end
    if x.val_objectif[1] == E[pos].val_objectif[1] && x.val_objectif[2] == E[pos].val_objectif[2]
        return false
    end
    if p!=length(E) && x.val_objectif[2] <= E[pos+1].val_objectif[2]
        return true
    end
    if p==length(E) && x.val_objectif[2] <= E[pos].val_objectif[2]
        return true
    end
    if p==0
        return true
    end
end

function improvement(E::Array{solution,1}, EPrime::Array{solution,1})
    for element in EPrime
        if no_dominated(element,E)
            return true
        end
    end
    return false
end

function neighborhood_change(E::Array{solution,1}, EPrime::Array{solution,1}, k::Int)
    if improvemnt(E, EPrime)
        E = kung(E, EPrime)
        return E,1
    end
    return E,k+1
end

function swap(x::solution, k::Int, prob::_bi01IP)
    iter = 0
    xPrime = copy(x)
    while iter < k
        rand1 = rand(1:length(xPrime.sol))
        rand2 = rand(1:length(xPrime.sol))
        xPrime.sol[rand1],xPrime.sol[rand2] = xPrime.sol[rand2],xPrime.sol[rand1]
        if verification(prob, xPrime)
            iter +=1
        else
            xPrime.sol[rand1],xPrime.sol[rand2] = xPrime.sol[rand1],xPrime.sol[rand2]
        end
    end
    return xPrime
end

function replace(x::solution, k::Int, prob::_bi01IP)
    iter = 0
    xPrime = copy(x)
    while iter < k
        random_idx = rand(1:length(xPrime.sol))
        xPrime.sol[random_idx] = (xPrime.sol[random_idx] + 1) % 2
        if verification(prob, xPrime)
            iter +=1
        else
            xPrime.sol[random_idx] = (xPrime.sol[random_idx] + 1) % 2
        end
    end
    return xPrime
end

function shake(x::solution, k::Int, type_shake::Int, prob::_bi01IP)
    if type_shake == 1
        return shake_swap(x,k,prob)
    elseif type_shake == 2
        return shake_replace(x,k,prob)
    elseif type_shake == 3
        xPrime = shake_swap(x,k,prob)
        return shake_replace(xPrime,k,prob)
    end
end

function mo_shake(E::Array{solution,1}, k::Int, type_shake::Int, prob::_bi01IP)
    EPrime = []
    for element in E
        xPrime = shake(element, k, type_shake, prob)
        union!(EPrime,xPrime)
    end
     return EPrime
end

function swap(x::solution, i::Int,j::Int)
    xPrime = copy(x)
    xPrime.sol[i],xPrime.sol[j] = xPrime.sol[j],xPrime.sol[i]
    return xPrime
end

function replace(x::solution, i::Int)
    xPrime = copy(x)
    xPrime.sol[i] = (1 + xPrime.sol[i]) % 2
    return xPrime
end

function voisinage_un_echange(x::solution, prob::_bi01IP)
    N = []
    union!(N,x)
    for i in 1:length(x.sol)
        for j in i+1:length(x.sol)
            xPrime = swap(x,i,j)
            if verification(prob, xPrime)
                union!(N,xPrime)
            end
        end
    end
    return N
end

function voisinage_deux_echange(x::solution, prob::_bi01IP)
    N = []
    N_un_echange = voisinage_un_echange(x,prob)
    for xPrime in N_un_echange
        union!(N,voisinage_un_echange(xPrime,prob))
    end
    return N
end

function replace_neighborhood(x::solution, k::Int, prob::_bi01IP)
    N = []
        for i in 1:length(x.sol)
            xPrime = replace(x,i)
            if verification(prob, xPrime)
                union!(N,xPrime)
            end
        end
    if k != 1
        NPrime = []
        for xPrime in N
            union!(NPrime,replace_neighborhood(xPrime,k-1,prob))
        end
        return NPrime
    else
        return N
    end
end

function swap_neighborhood(x::solution, k::Int, prob::_bi01IP)
    if k == 1
        return voisinage_un_echange(x,prob)
    end
    if k == 2
        return voisinage_deux_echange(x,prob)
    end
end

function VND_i(x::solution, kPrime_max::Int, i::Int)
    k = 1
    E = []
    union!(E,x)
    xPrime = x
    while k < kPrime_max
        N = union(replace_neighborhood(x, k,prob), swap_neighborhood(x,k,prob))
        zPrime = min(x -> x.val_objectif, N)
        for element in N
            if element.val_objectif == zPrime
                xPrime = element
                break
            end
        end
        E = kung(E,N)
        if xPrime.val_objectif[i] < x.val_objectif[i]
            x = xPrime
            k = 1
        else
            k += 1
        end
    end
    return E
end

function VND(E::Array{solution,1}, kPrime_max::Int, r::Int)
    # r = nombre d'objectifs
    S = fill([],r)
    i = 1
    exclusion = setdiff(E, S[i])
    while i <= r
        while length(inter) > 0
            random = rand((1:length(exclusion)))
            xPrime = S[i][random]
            Ei = VND_i(xPrime, kPrime_max, i)
            union!(S[i],E[i])
            union!(S[i],xPrime)
            exclusion = setdiff(E, S[i])
        end
        E,i = neighborhood_change(E, S[i], i)
    end
end

function GVNS(E::Array{solution,1}, k_max::Int, t_max::Int, type_shake::Int, prob::_bi01IP, r=2, kPrime_max=2)
    t = 0.0
    start = time()
    while t < t_max
        k = 1
        while k <= k_max
            EPrime = mo_shake(E, k, type_shake, prob)
            ESecond = VND(EPrime, kPrime_max, r)
            E, k = neighborhood_change(E, ESecond, k)
        end
    t = time() - start
    end
    return E
end

function main(fname::String, E::Array{solution,1}, k_max::Int,  t_max::Int, type_shake::Int,r=2, kPrime_max=2)

    @assert kPrime_max<=2 "Erreur : Au maximum 2 echanges dans la recherche locale"


end
