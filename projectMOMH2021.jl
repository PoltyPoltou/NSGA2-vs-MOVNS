# ==============================================================================
# Xavier Gandibleux - November 2021
#   Implemented in Julia 1.6

# ==============================================================================

println("""\nProject MOMH 2021" --------------------------------\n""")

const verbose = false
const graphic = true

include("parserMomkpZL.jl")
include("parserMomkpPG.jl")
include("vOptMomkp.jl")
include("displayGraphic.jl")
include("VNS.jl")
include("MOKP.jl")

# ==============================================================================
# Datastructure of a generic bi-objective 0/1 IP where all coefficients are integer
#
#   Max  sum{j=1,...,n} C(k,j) x(j)                k=1,...,p
#   st   sum{j=1,...,n} A(i,j) x(j) <= b_{i}       i=1,...,m
#                       x(j) = 0 or 1


# ==============================================================================

# Example ----------------------------------------------------------------------
fname = "instancesZL/knapsack.250.3"

# Read and load an instance of MO-01MKP from the collection of E. Zitzler / M. Laumanns
momkpZL = readInstanceMOMKPformatZL(verbose, fname)

# Reduce the MO-MKP instance to the two first objectives and store it as a generic Bi-01IP
dat1 = _bi01IP(momkpZL.P[1:2, :], momkpZL.W, momkpZL.ω)


# Example ----------------------------------------------------------------------
#fname = "instancesPG/set1/ZL28.DAT"
#fname = "instancesPG/set2/kp28W-Perm.DAT"
#fname = "instancesPG/set3/W7BI-rnd1-1800.DAT"

# Read and load an instance of Bi-01BKP from the collection of O. Perederieieva / X. Gandibleux
#momkpPG = readInstanceMOMKPformatPG(verbose, fname)

# Store it as a generic bi-01IP
#dat2 = _bi01IP(momkpPG.P, momkpPG.W, momkpPG.ω)

# Example ----------------------------------------------------------------------
#solverSelected = GLPK.Optimizer
#YN, XE = vSolveBi01IP(solverSelected, dat1.C, dat1.A, dat1.b)
#YN, XE = vSolveBi01IP(solverSelected, dat2.C, dat2.A, dat2.b)

# ---- Displaying the results (XE and YN)
# for n = 1:length(YN)
#     X = value.(XE, n)
#     print(findall(elt -> elt ≈ 1, X))
#     println("| z = ", YN[n])
# end

# for n = 1:length(YN)
#     println(YN[n])
# end

function TOOL_read_YN(prob::_MOMKP)

    f = open(string("YN/", prob.name, ".DAT"), "r")
    lines = readlines(f)
    close(f)

    len = length(lines) - 1

    YN = Vector{Vector{Int}}(undef, len)

    for i = 1:len
        YN[i] = parse.(Int64, split(lines[i+1]))
    end

    return YN
end
# Example ----------------------------------------------------------------------

function pointsExtremes(prob::_bi01IP, obj::Int)
    model = Model(with_optimizer(GLPK.Optimizer))

    @variable(model, x[1:size(prob.C, 2)], Bin)

    @constraint(model, [i = 1:length(prob.b)], sum((prob.A[i, j] * x[j]) for j = 1:size(prob.C, 2)) <= prob.b[i])

    if obj == 1
        @objective(model, Max, sum((prob.C[1, i]) * x[i] for i = 1:size(prob.C, 2)))
    else
        @objective(model, Max, sum((prob.C[2, i]) * x[i] for i = 1:size(prob.C, 2)))
    end
    optimize!(model)
    return value.(x)
end

pop = initPop(100, dat1)
x1 = solution([], pointsExtremes(dat1, 1), [])
x2 = solution([], pointsExtremes(dat1, 2), [])
verification(dat1, x1)
verification(dat1, x2)
push!(pop, x1)
push!(pop, x2)
sols = GVNS(pop, 2, 30, 1, dat1)
YN_VNS = map(x -> -x.val_objectif, sols)
displayGraphics(fname, TOOL_read_YN(momkpZL), "black")
displayGraphics(fname, YN_VNS, "green")
