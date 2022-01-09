using Revise
using ApproxOperator
using BenchmarkTools

efficiency()


nₚ = 11
nₑ = nₚ - 1
x = [1 / nₑ * i for i = 0:nₑ]
data = Dict(:x => x, :y => zeros(nₚ), :z => zeros(nₚ))
@btime Node(1, data)

a = Node(1, data)
@btime $a.x
@code_warntype a.x
