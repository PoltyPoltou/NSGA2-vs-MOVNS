import Base.copy
import Base.==
import Random.randperm
include("MOKP.jl")
include("tools.jl")



function copy(x::solution)
    return solution(copy(x.val_objectif), copy(x.sol), copy(x.cout))
end

function ==(x::solution, y::solution)
    return x.sol == y.sol
end
# C = [1 2 0 1; 1 1 1 2]
# A = [1 5 3 3; 5 1 1 6; 1 1 9 2]
# b = [5, 10, 10]
# test_pb = _bi01IP(C, A, b)
# ----------------------------------------------------------------
# Sélectionne les points pour former l'ensemble bornant primal
function kung(E::Array{solution,1}, EPrime::Array{solution,1})
    S = union(E, EPrime)
    maxDeltaY2 = maximum(x -> x.val_objectif[2], S) - minimum(x -> x.val_objectif[2], S)
    sort!(S, by = x -> x.val_objectif[1] + x.val_objectif[2] / maxDeltaY2)
    SN::Vector{solution} = []
    push!(SN, S[1])
    minYFeas = S[1].val_objectif
    for i = 2:length(S)
        if S[i].val_objectif[2] < minYFeas[2] && S[i].val_objectif[1] > minYFeas[1]
            push!(SN, S[i])
            minYFeas = S[i].val_objectif
        end
    end
    return SN
end

function verification(prob::_bi01IP, x::solution)
    poids = zeros(size(prob.A, 1))
    objectif = [0, 0]
    for i = 1:size(prob.A, 1)
        for j = 1:size(prob.A, 2)
            poids[i] += prob.A[i, j] * x.sol[j]
        end
    end
    for i = 1:size(prob.C, 2)
        objectif[1] += prob.C[1, i] * x.sol[i]
        objectif[2] += prob.C[2, i] * x.sol[i]
    end
    if sum(poids .<= prob.b) == length(prob.b)
        x.val_objectif = -objectif
        x.cout = poids
        return true
    else
        return false
    end
end

function verification_swap(prob::_bi01IP, x::solution, i::Int, j::Int)
    delta_cout = zeros(size(prob.A, 1))
    for l = 1:size(prob.A, 1)
        delta_cout[l] += prob.A[l, j] * (x.sol[j] - x.sol[i]) + prob.A[l, i] * (x.sol[i] - x.sol[j])
    end

    if sum(x.cout + delta_cout .<= prob.b) == length(prob.b)
        x.val_objectif[1] -= prob.C[1, j] * (x.sol[j] - x.sol[i]) + prob.C[1, i] * (x.sol[i] - x.sol[j])
        x.val_objectif[2] -= prob.C[2, j] * (x.sol[j] - x.sol[i]) + prob.C[2, i] * (x.sol[i] - x.sol[j])
        x.cout += delta_cout
        return true
    else
        return false
    end
end

function no_dominated(x::solution, E::Array{solution,1})
    # E must be sorted lexicographically to work
    pos = 0
    while pos + 1 <= length(E) && x.val_objectif[1] >= E[pos+1].val_objectif[1]
        pos += 1
    end
    if pos == 0
        return true
    elseif pos == length(E)
        return x.val_objectif[2] < E[pos].val_objectif[2]
    else
        x.val_objectif[2] == E[pos+1].val_objectif[2]
        return x.val_objectif[2] < E[pos].val_objectif[2]
    end
end

function improvement(E::Array{solution,1}, EPrime::Array{solution,1})
    sort!(E, by = x -> x.val_objectif[1])
    for element in EPrime
        if no_dominated(element, E)
            return true
        end
    end
    return false
end

function neighborhood_change(E::Array{solution,1}, EPrime::Array{solution,1}, k::Int)
    if improvement(E, EPrime)
        E = kung(E, EPrime)
        return E, 1
    end
    return E, k + 1
end

function swap(x::solution, k::Int, prob::_bi01IP)
    iter = 0
    swap_idx1 = 1
    swap_idx2 = 1
    swap_tuples1 = randperm(length(x.sol))
    swap_tuples2 = randperm(length(x.sol))
    xPrime = copy(x)
    while iter < k
        rand1 = swap_tuples1[swap_idx1]
        rand2 = swap_tuples2[swap_idx2]
        xPrime.sol[rand1], xPrime.sol[rand2] = xPrime.sol[rand2], xPrime.sol[rand1]
        if verification(prob, xPrime)
            iter += 1
        else
            xPrime.sol[rand1], xPrime.sol[rand2] = xPrime.sol[rand1], xPrime.sol[rand2]
        end
        swap_idx1 = swap_idx1 % length(x.sol) + 1
        if swap_idx1 == 1
            swap_idx2 += 1
        end
        if swap_idx2 > length(x.sol)
            return x
        end
    end
    @assert verification(prob, xPrime)
    return xPrime
end

function shake(x::solution, k::Int, type_shake::Int, prob::_bi01IP)
    if type_shake == 1
        return swap(x, k, prob)
    elseif type_shake == 2
        return swap(x, k, prob)
    elseif type_shake == 3
        return swap(x, k, prob)
    elseif type_shake == 4
        return swap(x, k, prob)
    end
end

function mo_shake(E::Array{solution,1}, k::Int, type_shake::Int, prob::_bi01IP)
    EPrime = Vector{solution}(undef, length(E))
    len = 1
    for element in E
        EPrime[len] = shake(element, k, type_shake, prob)
        len += 1
    end
    return EPrime
end

function swap(x::solution, i::Int, j::Int)
    xPrime = copy(x)
    xPrime.sol[i], xPrime.sol[j] = xPrime.sol[j], xPrime.sol[i]
    return xPrime
end

function replace(x::solution, i::Int)
    xPrime = copy(x)
    xPrime.sol[i] = (1 + xPrime.sol[i]) % 2
    return xPrime
end

function voisinage_un_echange(x::solution, prob::_bi01IP)
    N::Vector{solution} = Vector{solution}(undef, length(x.sol) * length(x.sol))
    index = 1
    N[index] = x
    index += 1
    for i = 1:length(x.sol)
        for j = i+1:length(x.sol)
            xPrime = swap(x, i, j)
            if verification_swap(prob, xPrime, i, j) && x != xPrime
                N[index] = xPrime
                index += 1
            end
        end
    end
    return unique(x -> x.sol, N[1:index-1])
end

function replace_neighborhood(x::solution, prob::_bi01IP)
    N = []
    for i = 1:length(x.sol)
        xPrime = replace(x, i)
        if verification(prob, xPrime)
            push!(N, xPrime)
        end
    end
    return N
end

function neighborhood(x::solution, k::Int, prob::_bi01IP)
    if k == 1
        return voisinage_un_echange(x, prob)
    end
    if k == 2
        return replace_neighborhood(x, prob)
    end
    return []
end

function VND_i(x::solution, kPrime_max::Int, i::Int, prob::_bi01IP)
    k = 1
    E::Vector{solution} = [x]
    xPrime = x
    while k <= kPrime_max
        N::Vector{solution} = neighborhood(x, k, prob)
        if length(N) > 0
            zPrime = minimum(x -> x.val_objectif[i], N)
            for element in N
                if element.val_objectif[i] == zPrime
                    xPrime = element
                    break
                end
            end
            E = kung(E, N)
            if xPrime.val_objectif[i] < x.val_objectif[i]
                x = xPrime
                k = 1
            else
                k += 1
            end
        else
            k += 1
        end
    end
    return E
end

function VND(E::Array{solution,1}, kPrime_max::Int, r::Int, prob::_bi01IP)
    # r = nombre d'objectifs
    S::Vector{Vector{solution}} = fill([], r)
    i = 1
    exclusion::Vector{solution} = setdiff(E, S[i])
    while i <= r
        while length(exclusion) > 0
            random = rand((1:length(exclusion)))
            xPrime = exclusion[random]
            Ei = VND_i(xPrime, kPrime_max, i, prob)
            union!(S[i], Ei)
            push!(S[i], xPrime)
            exclusion = setdiff(E, S[i])
        end
        E, i = neighborhood_change(E, S[i], i)
    end
    return E
end

function GVNS(E::Array{solution,1}, k_max::Int, t_max::Int, type_shake::Int, prob::_bi01IP, r = 2, kPrime_max = 2)
    t = 0.0
    loops = 0
    start = time()
    while t < t_max
        k = 1
        while k <= k_max && t < t_max
            EPrime = mo_shake(E, k, type_shake, prob)
            ESecond = VND(EPrime, kPrime_max, r, prob)
            E, k = neighborhood_change(E, ESecond, k)
            t = time() - start
        end
        loops += 1
    end
    #println("time used : ", t, " loops : ", loops)
    return E, t
end
function initPop(nIndiv::Int, prob::_bi01IP)
    population = Vector{solution}(undef, nIndiv)
    len = 1
    for i = 1:nIndiv
        sol::solution = solution([], zeros(size(prob.C, 2)), [])
        rand_indexes = randperm(length(sol.sol))
        idx = 1
        while verification(prob, sol)
            sol.sol[rand_indexes[idx]] = 1
            idx += 1
        end
        population[len] = sol
        len += 1
    end
    return sort(population, by = x -> x.val_objectif)
end

function initPopEpsilon(nIndiv::Int, prob::_bi01IP)
    x1 = solution([], pointsExtremes(prob, 1), [])
    x2 = solution([], pointsExtremes(prob, 2), [])
    pop::Vector{solution} = [x1, x2]
    verification(prob, x1)
    verification(prob, x2)
    delta2 = abs(x2.val_objectif[2] - x1.val_objectif[2])
    delta1 = abs(x2.val_objectif[1] - x1.val_objectif[1])
    step = delta2 / nIndiv / 2
    i = 1
    while i * step < delta2
        x = solution([], solutionEpsilon(prob, 1, -x1.val_objectif[2] + i * step), [])
        @assert verification(prob, x)
        push!(pop, x)
        i += 1
    end
    step = delta1 / nIndiv / 2
    i = 1
    while i * step < delta1
        x = solution([], solutionEpsilon(prob, 2, -x2.val_objectif[1] + i * step), [])
        @assert verification(prob, x)
        push!(pop, x)
        i += 1
    end
    return pop
end

function get_YN_VNS(solutions::Vector{solution})
    return map(x -> -x.val_objectif, solutions)
end
