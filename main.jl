include("parserMomkpZL.jl")
include("parserMomkpPG.jl")
include("displayGraphic.jl")
include("vOptMomkp.jl")
include("VNS.jl")
include("MOKP.jl")
include("tools.jl")
include("nsga_2.jl")
using CPLEX
function main(name::String, flag::Int, exact::Bool)
    verbose = false
    nGenerations = 10000000
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
    #println("--------------------------------------")
    #println(name)
    if exact
        if !isfile("YN/" * name * ".DAT")
            # we must solve with cplex first
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
        # displayGraphics(fname, YN, "black", L"y \in Y_N")
    end
    #VNS
    #println("VNS")
    #print("temps gen pop : ")
    pop = initPopEpsilon(10, prob)
    VNS_sol, real_time = GVNS(pop, 2, timeAllowed, 1, prob, 2, 2)
    YN_VNS = get_YN_VNS(VNS_sol)
    # displayGraphics(fname, YN_VNS, "green", L"y \in Y_{VNS}")
    # NSGA-II
    #println("NSGA")
    NSGA_sol = nsga2(nIndividus, nGenerations, probCrossover, probMutation, prob, real_time)
    YN_NSGA = get_YN_nsga(NSGA_sol)


    # displayGraphics(fname, YN_NSGA, "red", L"y \in Y_{NSGA}")
    if exact
        #println("hypervolumes (VNS,exact,ratio) ", compareHV(YN_VNS, YN))
        #println("hypervolumes (NSGA-II,exact,ratio) ", compareHV(YN_NSGA, YN))
    else
        #println("hypervolumes (VNS,NSGA-II,ratio) ", compareHV(YN_VNS, YN_NSGA))
    end
    v1 = round(compareHV(YN_VNS, YN)[3], digits = 4)
    v2 = round(compareHV(YN_NSGA, YN)[3], digits = 4)
    real_time = round(digits = 1, real_time)
    #println("(|YN_VNS|,|YN_NSGA|,|YN|)", (length(YN_VNS), length(YN_NSGA), exact ? length(YN) : 0))
    println("\\hline")
    println("$name & $(length(YN)) & $(length(YN_VNS)) & ", v1, " & $(length(YN_NSGA)) & ", v2, " & ", real_time, "\\\\")
    #println("--------------------------------------")
end

#main("knapsack.100.2", 1, true)
#main("knapsack.100.3", 1, true)
#main("knapsack.100.4", 1, true)
#main("knapsack.250.2", 1, true)
#main("knapsack.250.3", 1, true)
#main("knapsack.250.4",1,true)
#main("knapsack.500.2", 1, true)
#main("knapsack.500.3",1,true)
#main("knapsack.500.4",1,true)
#main("knapsack.750.2",1,true)
#main("knapsack.750.3",1,true)
#main("knapsack.750.4",1,true)
#main("set1/ZL28.DAT", 2, true)
#main("set1/ZL100.DAT", 2, true)
#main("set1/ZL105.DAT", 2, true)
#main("set1/ZL250.DAT", 2, true)
#main("set1/ZL500.DAT",2,true)
#main("set1/ZL750.DAT",2,true)
#main("set2/A1.DAT", 2, true)
#main("set2/A2.DAT", 2, true)
#main("set2/A3.DAT", 2, true)
#main("set2/A4.DAT", 2, true)
#main("set2/D1.DAT", 2, true)
#main("set2/D2.DAT", 2, true)
#main("set2/D3.DAT", 2, true)
#main("set2/D4.DAT", 2, true)
#main("set2/kp28.DAT", 2, true)
#main("set2/kp28-2.DAT", 2, true)
main("set2/kp28c1W-c2ZL.DAT", 2, true)
main("set2/kp28cW-WZL.DAT", 2, true)
main("set2/kp28W-Perm.DAT", 2, true)
main("set2/kp28W-ZL.DAT", 2, true)
main("set2/kp28W.DAT", 2, true)
main("set3/W7B1-rnd1-1800.DAT", 2, true)
main("set3/W7B1-rnd1-3000.DAT", 2, true)
main("set3/W7B1-tube1-1800.DAT", 2, true)
main("set3/W7B1-tube1-3000.DAT", 2, true)
main("set3/W7B1-tube1-asyn.DAT", 2, true)
main("set3/W7B1-tube2-1800.DAT", 2, true)
main("set3/Wcollage-tube.DAT", 2, true)
