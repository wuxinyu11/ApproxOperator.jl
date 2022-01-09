using Revise
using ApproxOperator
using BenchmarkTools

efficiency()


nₚ = 11
nₑ = nₚ - 1
x = [1 / nₑ * i for i = 0:nₑ]
data = Dict(:x => x, :y => zeros(nₚ), :z => zeros(nₚ))
paradata = Dict(:ξ => [0.0], :w => [1.0], :σ => [0.0],:b=>[0.0],:t=>[0.0])
# @btime Node(1, data)

𝓒 = [Node(1, data),Node(2, data)]
𝓖 = [Point(1, 1, paradata)]
b = (x, y, z) -> x^2
coefficients = Dict(:k => 1.0)
functions = Dict(:b => b)
ap = Seg2(𝓒,𝓖)
op = Operator(Val(:∇v∇u), coefficients)
k = zeros(nₚ,nₚ)
f = zeros(nₚ)
op1 = Operator(Val(:vt),Dict{Symbol,Float64}())
op1(ap,f)
@btime Operator(Val(:∇v∇u), $coefficients)
@btime op($ap,$k,$f)
# @code_warntype op(ap,k,f,Val(:∇v∇u))
# @btime $a.x
# @btime $ξ.ξ
# @btime getfield($a,:data)[:x][1]
# @code_warntype ξ.ξ

# @btime $a.x
# @btime $data[:x][1]
# @btime $x[1]
