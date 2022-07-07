
using Revise, ApproxOperator, BenchmarkTools, YAML

config = YAML.load_file("./fem.yml")
elements = importmsh("./msh/patchtest.msh",config)

nₚ = getnₚ(elements["Ω"])

s = 3.5/20 .* ones(nₚ)
push!(elements["Ω"],:s₁=>s,:s₂=>s,:s₃=>s)
set∇𝝭!(elements["Ω"])
# @btime set𝝭!($elements["Ω"])
# @btime set∇𝝭!(elements["Ω"])
# set∇𝝭!(elements["Ω"])
# set∇̃𝝭!(elements["Ω"])
# set∇̃²𝝭!(elements["Ω"])
# @btime set∇²₂𝝭!(elements["Ω"])
# @btime set∇³𝝭!(elements["Ω"])
# set∇𝝭!(elements["Ω"])
# set∇̂³𝝭!(elements["Ω"])
# set∇²𝝭!(elements["Ω"])
# set∇³𝝭!(elements["Ω"])
# set∇̃𝝭!(elements["Ω"])
# set∇̃²𝝭!(elements["Ω"])
f = check∇𝝭(elements["Ω"])
# f = check∇²𝝭(elements["Ω"])
# f = check∇³𝝭(elements["Ω"])
