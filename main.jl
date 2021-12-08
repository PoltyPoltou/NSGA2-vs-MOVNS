include("parserMomkpZL.jl")
include("parserMomkpPG.jl")
include("vOptMomkp.jl")
include("displayGraphic.jl")
include("VNS.jl")
include("MOKP.jl")
include("tools.jl")



name = "knapsack.250.2"
fname = "instancesZL/" * name
momkpZL = readInstanceMOMKPformatZL(false, fname)
dat1 = _bi01IP(momkpZL.P[1:2, :], momkpZL.W, momkpZL.ω)

@time pop = initPopEpsilon(20, dat1)
sols = GVNS(pop, 2, 10, 1, dat1, 2, 2)
YN_VNS = map(x -> -x.val_objectif, sols)
YN = TOOL_read_YN(name)
displayGraphics(fname, YN, "black")
displayGraphics(fname, YN_VNS, "green")
