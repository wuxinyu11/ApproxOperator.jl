
using Revise, ApproxOperator

elements, nodes = importmsh("./msh/bar.msh")
nₚ = length(nodes[:x])
nₑ = length(elements["Domain"])

set𝓖!(elements["Domain"],:SegGI2)
set𝓖!(elements["NBC"],:PoiGI1)
set𝓖!(elements["EBC"],:PoiGI1)

union!(elements["NBC"][1],elements["Domain"][nₑ])
