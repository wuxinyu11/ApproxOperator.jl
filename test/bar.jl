
using Revise, ApproxOperator, BenchmarkTools

x = rand(100)
y = rand(100)
z = rand(100)
data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
n = Node(2,data)
@btime $n.x
@btime $n.y
@btime $n.I
# @btime ApproxOperator.RV(1,$x)
