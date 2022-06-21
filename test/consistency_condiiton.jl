
using Revise, ApproxOperator, BenchmarkTools

# elements, nodes = importmsh("./msh/bar.msh")
elements, nodes = importmsh("./msh/patchtest.msh")
nₚ = length(nodes[:x])

# type = (Node,:Quadratic2D,:□,:CubicSpline)
# type = (SNode,:Quadratic2D,:□,:QuinticSpline)
type = (SNode,:Cubic2D,:□,:QuinticSpline)
# s = 2.5/20*ones(nₚ)
s = 3.2/20*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z],n = 2,γ = 5)
elements["Ω"] = ReproducingKernel{type...,:Tri3}(elements["Ω"])
set𝓖!(elements["Ω"],:TriRK6)
set_memory_𝗠!(elements["Ω"],:∂1,:∂x,:∂y,:∂z,:∇̃,:∂x_,:∂y_,:∂z_,:∂x²,:∂x∂y,:∂y²,:∂x∂z,:∂y∂z,:∂z²,:∂x³,:∂x²∂y,:∂x∂y²,:∂y³,:∂x²_,:∂x∂y_,:∂y²_)
set_memory_𝝭!(elements["Ω"],:∂1,:∂x,:∂y,:∂z,:∂x²,:∂x∂y,:∂y²,:∂x∂z,:∂y∂z,:∂z²,:∂x³,:∂x²∂y,:∂x∂y²,:∂y³)
sp(elements["Ω"])
# set𝓖!(elements["Ω"],:TriGI13,:∂1,:∂x,:∂y,:∂z,:∂x²,:∂x∂y,:∂y²,:∂z²,:∂x∂z,:∂y∂z,:∂x³,:∂x²∂y,:∂x∂y²,:∂y³)

# @btime set𝝭!(elements["Ω"])
# @btime set∇𝝭!(elements["Ω"])
# set∇𝝭!(elements["Ω"])
# set∇̃𝝭!(elements["Ω"])
# set∇̃²𝝭!(elements["Ω"])
# @btime set∇²𝝭!(elements["Ω"])
# @btime set∇³𝝭!(elements["Ω"])
set∇̂³𝝭!(elements["Ω"])
# set∇̃𝝭!(elements["Ω"])
# set∇̃²𝝭!(elements["Ω"])
# f = checkConsistency(elements["Ω"])
# f = checkConsistency(elements["Ω"],ApproxOperator.get∇𝝭,ApproxOperator.get∇𝒑)
# f = checkConsistency(elements["Ω"],ApproxOperator.get∇²𝝭,ApproxOperator.get∇²𝒑)
f = checkConsistency(elements["Ω"],ApproxOperator.get∇³𝝭,ApproxOperator.get∇³𝒑)

# err_chol, err, err_x, err_y = ApproxOperator.test_cal𝗠!(elements["Ω"][1],(0.,0.,0.))
# err,err_chol,err_inv,err_I,err1,err2 = ApproxOperator.test_cal𝗠!(elements["Ω"][1],(0.,0.,0.))