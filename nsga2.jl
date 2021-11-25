struct_donnee = Vector{Tuple{Tuple{Float,Float},Tuple{Float,Float}}}

function genererPopulation(nIndividus::Int)
end

function classementGoldberg(population::struct_donnee)
end

function calculValeurDensite(elementsVoisins, pointIdeal, PointAntiIdeal)
    result = 0
    for i in 1:2
        result += (elementsVoisins[1][2][i] - elementsVoisins[2][2][i])/(pointIdeal[i] - PointAntiIdeal[i])
    end
    return result
end

function selectionSurpopulation(populationTotale::struct_donnee, populationRangFixe::struct_donnee, nFinal::Int)
    populationTrieeObjectif = sort(populationRangFixe)
    pointIdeal = min(x -> x[2], populationTotale)
    pointAntiIdeal = max(x -> x[2], populationTotale)
    kappa = Dict([(extremaMin, +Inf), (extremaMax, +Inf)])
    for i in 1:length(populationRangFixe)
end


function extrairePopRang(population, ranking::Vector{Int}, rang::Int)
    listeRang::struct_donnee = zeros(length(population))
    idxListe = 1
    for i in 1:length(population)
        if ranking[i] == rang
            listeRang[idxListe] = population[i]
            idxListe += 1
        end
    end
    resize!(listeRang,idxListe-1)
    return listeRang
end

function selectionIndividus(population,ranking, nIndividus::Int)
    nouvelle_pop::struct_donnee = Vector(undef,0)
    rangActuel = 1
    while length(nouvelle_pop) < nIndividus
        popRangFixe = extrairePopRang(population, ranking, rangActuel)
        if length(nouvelle_pop) + length(nouvelle_pop) <= nIndividus
            union!(nouvelle_pop, popRangFixe)
        else
            union!(nouvelle_pop, selectionSurpopulation(popRangFixe, nIndividus - length(nouvelle_pop)))
        end
        rangActuel += 1
    end

end

function binaryTournament(population, ranking)
end
function crossover()
end
function mutation()
end
function genererPopulationEnfants(nIndividus, probCrossover, probMutation, populationAReproduire)
    # tournoi + crossover + mutation et doit renvoyer nIndividus elements
end

function ngsa(nIndividus::Int, nGenerations::Int, probCrossover::Float, probMutation::Float)
    #changer le type quand on a la structure pour
    population::struct_donnee= genererPopulation(nIndividus)
    ranking = classementGoldberg(population)
    for i = 1:nGenerations
        populationEnfants = genererPopulationEnfants(nIndividus, probCrossover, probMutation, population)
        popTotale = union(populationEnfants, population)
        rankingTotal = classementGoldberg(popTotale)
        population = selectionsIndividus(popTotale,rankingTotal)

    end
end
