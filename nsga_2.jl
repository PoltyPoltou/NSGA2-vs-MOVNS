import Random.randperm
import Base.copy

include("MOKP.jl")
include("tools.jl")

mutable struct Individu
    x::Array{Int,1} # vecteur solution
    y::Array{Int,1} # valeurs sur les objectifs
    rang::Int        # rang
    listeDomines::Array{Individu,1}    # liste des individus dominés par celui-ci, vide au début
    nbDomine::Int                    # nombre d'individus qui dominent celui-ci, vide au début
    distanceCrowding::Float64  # valeur de la distance de Crowding
end


function copy(a::Individu)
    return Individu(copy(a.x), copy(a.y), copy(a.rang), copy(a.listeDomines), copy(a.nbDomine), copy(a.distanceCrowding))
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
        x.val_objectif = objectif
        x.cout = poids
        return true
    else
        return false
    end
end

function genererPopulation(nIndividus::Int, prob::_bi01IP)
    pop = Vector{Individu}(undef, nIndividus)
    len = 1

    x1 = Individu(pointsExtremes(prob, 1), [0, 0], 0, [], 0, 0.0)
    x2 = Individu(pointsExtremes(prob, 2), [0, 0], 0, [], 0, 0.0)
    x1 = updateY(x1, prob)
    x2 = updateY(x2, prob)

    pop[len] = x1
    len += 1
    pop[len] = x2
    len += 1

    for i = 1:(nIndividus-2)
        sol = solution([], zeros(size(prob.C, 2)), [])
        rand_indexes = randperm(length(sol.sol))
        idx = 1
        while verification(prob, sol)
            sol.sol[rand_indexes[idx]] = 1
            idx += 1
        end

        indiv = Individu(sol.sol, sol.val_objectif, 0, [], 0, 0.0)
        pop[len] = indiv
        len += 1
    end
    return pop
end

# Prends deux individus et renvoie true si le premier domine le second
function domine(p::Individu, q::Individu)
    if p.y[1] >= q.y[1] && p.y[2] >= q.y[2]
        return true
    end
    return false
end

# Prend un population et renvoie vrai s'il existe un Individu
#avec un rang non assigné
function rangZero(population)
    for p in population
        if p.rang == 0
            return true
        end
    end
    return false
end

# Prend une population et assigne un rang à chaque individus
function triNonDomine(population)
    for p in population
        empty!(p.listeDomines)
        p.nbDomine = 0
        p.rang = 0
    end

    for i = 1:length(population)
        for j = i+1:length(population)
            if domine(population[i], population[j])
                push!(population[i].listeDomines, population[j])
                population[j].nbDomine += 1
            elseif domine(population[j], population[i])
                push!(population[j].listeDomines, population[i])
                population[i].nbDomine += 1
            end
        end
        if population[i].nbDomine == 0
            population[i].rang = 1
        end
    end

    k = 2
    while rangZero(population)
        for p in population
            if p.rang == k - 1
                for q in p.listeDomines
                    q.nbDomine -= 1
                    if q.nbDomine == 0
                        q.rang = k
                    end
                end
            end
        end
        k += 1
    end
end

#Prend une population et renvoie le front maximum
function trouverMaxFront(population)
    max = 0
    for p in population
        if p.rang > max
            max = p.rang
        end
    end
    return max
end

# Prend un frond et assigne la distance de corwding à chaque individu du front
function distanceCrowding(front)
    for i = 1:length(front[1].y)
        sort!(front, lt = (a, b) -> a.y[i] < b.y[i], rev = true)

        max = front[1].y[i]
        min = front[end].y[i]

        for j = 2:(length(front)-1)
            greaterNeighbor = front[j-1].y[i]
            lowerNeighbor = front[j+1].y[i]
            front[j].distanceCrowding += ((greaterNeighbor - lowerNeighbor) / (max - min))
        end
    end
end

#renvoie le gagnant du tournoi entre deux individus
function tournoi(a::Individu, b::Individu)
    if a.rang < b.rang
        return a
    elseif b.rang < a.rang
        return b
    elseif a.distanceCrowding < b.distanceCrowding
        return a
    elseif b.distanceCrowding < a.distanceCrowding
        return b
    end
    return rand((a, b))
end

#teste si deux individus sont égaux
function isequal(a::Individu, b::Individu)
    return a.x == b.x
end

# prends un individu et assigne le y à son vecteur x correspondant
# Si le vecteur ne respecte pas les conditions du probleme alors
# des xi aléatoires passent de 1 à 0 jusqu'à ce que l'individu respecte les conditions
function updateY(a::Individu, prob::_bi01IP)
    a.y[1], a.y[2] = 0, 0
    for i = 1:length(a.x)
        a.y[1] += prob.C[1, i] * a.x[i]
        a.y[2] += prob.C[2, i] * a.x[i]
    end
    poids = zeros(size(prob.A, 1))
    for i = 1:size(prob.A, 1)
        for j = 1:size(prob.A, 2)
            poids[i] += prob.A[i, j] * a.x[j]
        end
    end

    for i = 1:length(prob.b)
        while poids[i] > prob.b[i]
            pos = rand(1:length(a.x))
            if a.x[pos] == 1
                a.x[pos] = 0
                a.y[1] -= prob.C[1, pos]
                a.y[2] -= prob.C[2, pos]
                for j = 1:size(prob.A, 1)
                    poids[j] -= prob.A[j, pos]
                end
            end
        end
    end
    return a
end

# réalise un crossover entre 2 individus
function crossover(a::Individu, b::Individu, probCrossover::Float64, prob::_bi01IP)
    if rand() <= probCrossover
        i = rand(1:length(a.x))
        savea = copy(a.x[i:end])
        a.x[i:end] = b.x[i:end]
        b.x[i:end] = savea

        a = updateY(a, prob)
        b = updateY(b, prob)
    end
    return a, b
end

# réalise une simple mutation en inversant un bit d'indice aléatoire d'un individu a
function mutation(a::Individu, prob::_bi01IP, probMutation::Float64)
    if rand() <= probMutation
        i = rand(1:length(a.x))
        if a.x[i] == 1
            a.x[i] = 0
            a.y[1] -= prob.C[1, i]
            a.y[2] -= prob.C[2, i]
        else
            a.x[i] = 1
            a = updateY(a, prob)
        end
    end
    return a
end

function genererPopulationEnfants(nIndividus::Int, probCrossover::Float64, probMutation::Float64, population, prob::_bi01IP)
    enfants = Vector{Individu}(undef, nIndividus)
    len = 1
    for i = 1:2:nIndividus
        indiv1 = tournoi(population[rand(1:length(population))], population[rand(1:length(population))])
        indiv2 = tournoi(population[rand(1:length(population))], population[rand(1:length(population))])

        a = copy(indiv1)
        b = copy(indiv2)
        a, b = crossover(a, b, probCrossover, prob)

        a = mutation(a, prob, probMutation)
        b = mutation(b, prob, probMutation)


        enfants[len] = a
        len += 1
        enfants[len] = b
        len += 1
    end
    return enfants
end

#Selectionne une population de n individus ou n + epsilon parmis 2xn individus
function selectionIndividus(popTotale, nIndividus::Int)
    if popTotale[nIndividus].rang == 1
        debut = findfirst(a -> a.rang == 2, popTotale)
        fin = findlast(a -> a.rang == 2, popTotale)
        popsecond = view(popTotale, debut:fin)
        distanceCrowding(popsecond)
        sort!(popsecond, by = a -> a.distanceCrowding)
        epsilon = round(Int, length(popsecond) / 5)
        popFinale = union(view(popTotale, 1:debut-1), view(popsecond, 1:epsilon))
    else
        debut = findfirst(a -> a.rang == popTotale[nIndividus].rang, popTotale)
        fin = findlast(a -> a.rang == popTotale[nIndividus].rang, popTotale)
        popsecond = view(popTotale, debut:fin)
        distanceCrowding(popsecond)
        sort!(popsecond, by = a -> a.distanceCrowding)
        popFinale = union(view(popTotale, 1:debut-1), view(popsecond, 1:(nIndividus-debut+1)))
    end
    return popFinale
end

#
function kung(population::Array{Individu,1})
    SN::Vector{Individu} = []
    sort!(population, by = a -> a.y[1])
    push!(SN, population[1])
    minYFeas = population[1].y[2]
    for i = 2:length(population)
        if population[i].y[2] < minYFeas
            push!(SN, population[i])
            minYFeas = population[i].y[2]
        end
    end
    return SN
end

function nsga2(nIndividus::Int, nGenerations::Int, probCrossover::Float64, probMutation::Float64, prob::_bi01IP)
    start = time()
    population = genererPopulation(nIndividus, prob)
    triNonDomine(population)
    t = 0.0
    start = time()
    gen = 0
    while gen < nGenerations && t < 60
        enfants = genererPopulationEnfants(nIndividus, probCrossover, probMutation, population, prob)
        popTotale = union(enfants, population)
        triNonDomine(popTotale)
        sort!(popTotale, by = a -> a.rang)
        population = selectionIndividus(popTotale, nIndividus)
        t = time() - start
        gen += 1
    end
    SN = kung(population)
    println("pts non dominés : ", length(SN))
    println("generations : ", gen)
    println("time : ", time() - start)
    return SN
end

function get_YN_nsga(SN::Vector{Individu})
    return map(x -> x.y, SN)
end
