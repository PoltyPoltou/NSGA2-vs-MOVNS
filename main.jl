include("parserMomkpZL.jl")
include("parserMomkpPG.jl")
include("displayGraphic.jl")
include("vOptMomkp.jl")
include("VNS.jl")
include("MOKP.jl")
include("tools.jl")
include("nsga_2.jl")
using CPLEX
function main(name::String, flag::Int)
    verbose = false
    nGenerations = 100000
    nIndividus = 100 #size(prob.C)[2]
    probCrossover = 0.9
    probMutation = 0.1
    timeAllowed = 10
    if flag == 1
        fname = "instancesZL/" * name
        momkpZL = readInstanceMOMKPformatZL(verbose, fname)
    else
        fname = "instancesPG/" * name
        momkpZL = readInstanceMOMKPformatPG(verbose, fname)
    end
    prob = _bi01IP(momkpZL.P[1:2, :], momkpZL.W, momkpZL.ω)
    println(name)
    if !isfile("YN/" * name * ".DAT")
        solverSelected = CPLEX.Optimizer
        YN, XE = vSolveBi01IP(solverSelected, prob.C, prob.A, prob.b)
        Yn_int = Vector{Vector{Int}}(undef, length(YN))
        for j = 1:length(YN)
            Yn_int[j] = [0, 0]
            Yn_int[j][1] = round(Int, YN[j][1])
            Yn_int[j][2] = round(Int, YN[j][2])
        end
        TOOL_write_YN(Yn_int, name)
        YN = Yn_int
    else
        YN = TOOL_read_YN(name)
    end
    displayGraphics(fname, YN, "black")
    #VNS
    println("VNS")
    print("temps gen pop : ")
    @time pop = initPopEpsilon(20, prob)
    VNS_sol = GVNS(pop, 2, timeAllowed, 1, prob, 2, 2)
    YN_VNS = get_YN_VNS(VNS_sol)
    # NSGA-II
    println("NSGA")
    nGenerations = 100000
    nIndividus = 100 #size(prob.C)[2]
    probCrossover = 0.9
    probMutation = 0.1
    NSGA_sol = nsga2(nIndividus, nGenerations, probCrossover, probMutation, prob, timeAllowed)
    YN_NSGA = get_YN_nsga(NSGA_sol)


    displayGraphics(fname, YN_VNS, "green")
    displayGraphics(fname, YN_NSGA, "red")
    println("hypervolumes (VNS,exact,ratio) ", compareHV(YN_VNS, YN))
    println("hypervolumes (NSGA-II,exact,ratio) ", compareHV(YN_NSGA, YN))
    println("hypervolumes (VNS,NSGA-II,ratio) ", compareHV(YN_VNS, YN_NSGA))
    println("(|YN_VNS|,|YN_NSGA|,|YN|)", (length(YN_VNS), length(YN_NSGA), length(YN)))

end
main("knapsack.100.4", 1)