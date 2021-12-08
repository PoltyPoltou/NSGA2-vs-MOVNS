using JuMP, GLPK, CPLEX
include("MOKP.jl")
function TOOL_read_YN(name)

    f = open(string("YN/", name, ".DAT"), "r")
    lines = readlines(f)
    close(f)

    len = length(lines) - 1

    YN = Vector{Vector{Int}}(undef, len)

    for i = 1:len
        YN[i] = parse.(Int64, split(lines[i+1]))
    end

    return YN
end

function TOOL_write_YN(YN::Vector{Vector{Int}}, name)
    way = string("YN/" * name * ".DAT")
    open(way, "w") do f

        println(f, string(length(YN)))
        for i = 1:length(YN)
            println(f, string(YN[i][1], " ", YN[i][2]))
        end
    end
end

function pointsExtremes(prob::_bi01IP, obj::Int)
    model = Model(with_optimizer(CPLEX.Optimizer))
    set_silent(model)
    @variable(model, x[1:size(prob.C, 2)], Bin)

    @constraint(model, [i = 1:length(prob.b)], sum((prob.A[i, j] * x[j]) for j = 1:size(prob.C, 2)) <= prob.b[i])

    if obj == 1
        @objective(model, Max, sum((prob.C[1, i]) * x[i] for i = 1:size(prob.C, 2)))
    else
        @objective(model, Max, sum((prob.C[2, i]) * x[i] for i = 1:size(prob.C, 2)))
    end
    optimize!(model)
    return round.(value.(x))
end

function solutionEpsilon(prob::_bi01IP, obj::Int, epsilon)
    model = Model(with_optimizer(CPLEX.Optimizer))
    set_silent(model)

    @variable(model, x[1:size(prob.C, 2)], Bin)

    @constraint(model, [i = 1:length(prob.b)], sum((prob.A[i, j] * x[j]) for j = 1:size(prob.C, 2)) <= prob.b[i])

    if obj == 1
        @objective(model, Max, sum((prob.C[1, i]) * x[i] for i = 1:size(prob.C, 2)))
        @constraint(model, sum((prob.C[2, i]) * x[i] for i = 1:size(prob.C, 2)) >= epsilon)
    else
        @objective(model, Max, sum((prob.C[2, i]) * x[i] for i = 1:size(prob.C, 2)))
        @constraint(model, sum((prob.C[1, i]) * x[i] for i = 1:size(prob.C, 2)) >= epsilon)
    end
    optimize!(model)
    return round.(Int, value.(x))
end


function compareHV(YN1::Vector{Vector{Int}}, YN2::Vector{Vector{Int}})
    # must be a maximization problem
    # example of shape of YN = [[5,3],...]
    sort!(by = x -> x[1], YN1)
    sort!(by = x -> x[1], YN2)
    minY1 = minimum(x -> x[1], union(YN1, YN2))
    minY2 = minimum(x -> x[2], union(YN1, YN2))
    hypervolume1::Float64 = 0
    hypervolume2::Float64 = 0
    for j = 1:length(YN1)
        if j == 1
            hypervolume1 += (YN1[1][1] - minY1) * (YN1[1][2] - minY2)
        else
            hypervolume1 += abs(YN1[j][1] - YN1[j-1][1]) * (YN1[j][2] - minY2)
        end
    end
    for j = 1:length(YN2)
        if j == 1
            hypervolume2 += (YN2[1][1] - minY1) * (YN2[1][2] - minY2)
        else
            hypervolume2 += abs(YN2[j][1] - YN2[j-1][1]) * (YN2[j][2] - minY2)
        end
    end
    return hypervolume1, hypervolume2, hypervolume1 / hypervolume2
end

