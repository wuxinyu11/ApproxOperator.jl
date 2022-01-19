using Revise
using ApproxOperator
using BenchmarkTools

nₚ = 10
n = 10
x = rand(nₚ)
y = zeros(nₚ)
z = zeros(nₚ)
ξ = rand(n)
data = Dict(:x=>x,:y=>y,:z=>z,:s₁=>rand(nₚ),:s₂=>rand(nₚ),:s₃=>rand(nₚ))
data_ = Dict(:ξ=>ξ)
index = zeros(Int,length(ξ)+1)
𝝭 = Dict(:∂1=>Float64[])
𝓒 = [Node(i,data) for i in 1:nₚ]
𝓖 = [SNode(i,data_,index,𝝭) for i in 1:n]
𝗠 = Dict(:∂1 => SymMat(3))
𝝭ᵉ = Dict(:∂1 => zeros(nₚ))
ap = SegN(𝓒,𝓖,𝗠,𝝭ᵉ,:Linear1D,:□,:CubicSpline)
op = Operator(:𝝭,Any)
op(ap)
