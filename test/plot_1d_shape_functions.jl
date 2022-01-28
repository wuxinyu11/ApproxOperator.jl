using Revise, ApproxOperator, BenchmarkTools, Plots

nₚ = 11
nₑ = nₚ-1
x = [1/nₑ*i for i in 0:nₑ]

data = Dict(:x=>x,:y=>zeros(nₚ),:z=>zeros(nₚ))

elements = [Seg2(i,i+1,data) for i in 1:nₑ]

set𝓖!
