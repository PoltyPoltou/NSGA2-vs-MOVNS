include("parserMomkpZL.jl")
include("parserMomkpPG.jl")
include("vOptMomkp.jl")
include("displayGraphic.jl")
include("VNS.jl")
include("MOKP.jl")
include("tools.jl")
using CPLEX
name = "knapsack.500.2"
fname = "instancesZL/" * name

momkpZL = readInstanceMOMKPformatZL(false, fname)
dat1 = _bi01IP(momkpZL.P[1:2, :], momkpZL.W, momkpZL.ω)
solverSelected = CPLEX.Optimizer
YN, XE = vSolveBi01IP(solverSelected, dat1.C, dat1.A, dat1.b)
Yn_int = Vector{Vector{Int}}(undef, length(YN))
for j = 1:length(YN)
    Yn_int[j] = [0, 0]
    Yn_int[j][1] = round(Int, YN[j][1])
    Yn_int[j][2] = round(Int, YN[j][2])
end
TOOL_write_YN(Yn_int, name)