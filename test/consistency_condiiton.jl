
using Revise, ApproxOperator

elements, nodes = importmsh("./msh/bar.msh")
nₚ = length(nodes[:x])

type = (Node,:Quadratic1D,:□,:CubicSpline)
s = 2.5/20*ones(nₚ)

sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z],n = 3,γ = 5)
elements["Ω"] = ReproducingKernel{type...,:Seg2}(elements["Ω"],sp)
set𝓖!(elements["Ω"],:SegGI5,:∂1)

f = checkConsistency(elements["Ω"],ApproxOperator.get𝝭,ApproxOperator.get𝒑)
