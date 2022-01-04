using Revise
using ApproxOperator
using BenchmarkTools

# efficiency_meshfree()

x = [Node(1.0/10*i,0.0,0.0) for i in 0:10]
aps = [Seg2([i,i+1],x) for i in 1:10]
sp = RegularGrid(x)

xₘ = Node(0.5,0.0,0.0)
indices = sp(xₘ)

sp(aps)
