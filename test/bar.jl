
using Revise, ApproxOperator, BenchmarkTools

data = Dict([:x=>(1,rand(10)),:y=>(2,rand(10))])

x1 = Node((2,3),data)
x2 = Node((2,),data)
# nodes = Node{2}[]
# push!(nodes,x1)
# push!(nodes,x2)
𝓒 = [Node((i,1),data) for i in 1:10]
𝓖 = [Node((2,),data) for i in 1:10]
e = Element(𝓒,𝓖)

