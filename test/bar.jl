
using Revise, ApproxOperator, BenchmarkTools

x = rand(100)
y = rand(100)
z = rand(100)
nodes = Node(:x=>x,:y=>y,:z=>z)
# @btime nodes = FENode(:x=>$x,:y=>$y,:z=>$z)
n = nodes[1]
n.x[2]
@btime a = $n.x*1.
# @btime ApproxOperator.RV(1,$x)