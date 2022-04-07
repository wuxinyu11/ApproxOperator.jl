
using Revise, ApproxOperator

# elements, nodes = importmsh("./msh/bar.msh")
elements, nodes = importmsh("./msh/patchtest.msh")
nₚ = length(nodes[:x])

# type = (Node,:Quadratic2D,:□,:CubicSpline)
# type = (SNode,:Quadratic2D,:□,:QuinticSpline)
type = (SNode,:Cubic2D,:□,:QuinticSpline)
# s = 2.5/20*ones(nₚ)
s = 3.5/20*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z],n = 3,γ = 5)
elements["Ω"] = ReproducingKernel{type...,:Tri3}(elements["Ω"],sp)
set𝓖!(elements["Ω"],:TriRK6,:∂1,:∂x,:∂y,:∂z,:∂x²,:∂x∂y,:∂y²,:∂z²,:∂x∂z,:∂y∂z,:∂x³,:∂x²∂y,:∂x∂y²,:∂y³)
# set𝓖!(elements["Ω"],:TriGI13,:∂1,:∂x,:∂y,:∂z,:∂x²,:∂x∂y,:∂y²,:∂z²,:∂x∂z,:∂y∂z,:∂x³,:∂x²∂y,:∂x∂y²,:∂y³)

# set𝝭!(elements["Ω"])
# set∇𝝭!(elements["Ω"])
# set∇̃𝝭!(elements["Ω"])
# set∇̃²𝝭!(elements["Ω"])
# set∇²𝝭!(elements["Ω"])
set∇³𝝭!(elements["Ω"])
# set∇̃𝝭!(elements["Ω"])
# set∇̃²𝝭!(elements["Ω"])
# f = checkConsistency(elements["Ω"])
# f = checkConsistency(elements["Ω"],ApproxOperator.get∇𝝭,ApproxOperator.get∇𝒑)
# f = checkConsistency(elements["Ω"],ApproxOperator.get∇²𝝭,ApproxOperator.get∇²𝒑)
f = checkConsistency(elements["Ω"],ApproxOperator.get∇³𝝭,ApproxOperator.get∇³𝒑)
