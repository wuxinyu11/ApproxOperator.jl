
using Revise, ApproxOperator, BenchmarkTools

data = Dict([:x=>(:I,rand(10)),:y=>(:J,rand(10))])
# data = Dict([:x=>(1,rand(10)),:y=>(2,rand(10))])

@btime x2 = Node((J=2,),data)
# f(x::Node) = x.x
# @btime a = test($x1,:x)
# x1 = Node((I=2,J=3),data)
@btime begin
    x1 = Node((I=2,J=3),data)
    x1.I
    x1.x
end
# @code_warntype x1.x
# nodes = Node{2}[]
# push!(nodes,x1)
# push!(nodes,x2)
# 𝓒 = [Node((i,1),data) for i in 1:10]
# 𝓖 = [Node((2,),data) for i in 1:10]
# e = Element(𝓒,𝓖)

# @code_warntype test(x1,:x)

# struct Test
#     a::NamedTuple{(:I,:J),NTuple{2,Int}}
#     b::Dict{Symbol,Vector{Float64}}
# end

# a = (I=1,J=2)
# b = Dict([:x=>rand(3)])
# c = Test(a,b)
# @btime $c.b[:x][$c.a.I]