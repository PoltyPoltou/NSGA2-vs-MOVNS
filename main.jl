include("parserMomkpZL.jl")
include("displayGraphic.jl")
include("VNS.jl")
include("MOKP.jl")
include("tools.jl")
include("nsga_2.jl")



name = "knapsack.100.2"
println(name)
fname = "instancesZL/" * name
momkpZL = readInstanceMOMKPformatZL(false, fname)
dat1 = _bi01IP(momkpZL.P[1:2, :], momkpZL.W, momkpZL.ω)
print("temps gen pop : ")
@time pop = initPopEpsilon(20, dat1)
#VNS
VNS_sol = GVNS(pop, 2, 10, 1, dat1, 2, 2)
YN_VNS = get_YN_VNS(VNS_sol)
# NSGA-II
nGenerations = 100000
nIndividus = 100 #size(prob.C)[2]
probCrossover = 0.9
probMutation = 0.1
NSGA_sol = nsga2(nIndividus, nGenerations, probCrossover, probMutation, dat1)
YN_NSGA = get_YN_nsga(NSGA_sol)


YN = TOOL_read_YN(name)
displayGraphics(fname, YN, "black")
displayGraphics(fname, YN_VNS, "green")
displayGraphics(fname, YN_NSGA, "red")
println("hypervolumes (VNS,exact,ratio) ", compareHV(YN_VNS, YN), " ; (|YN_VNS|,|YN_NSGA|,|YN|)", (length(YN_VNS), length(YN_NSGA), length(YN)))
println("hypervolumes (VNS,NSGA-II,ratio) ", compareHV(YN_VNS, YN_NSGA))