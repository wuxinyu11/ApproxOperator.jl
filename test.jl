using Revise
using ApproxOperator
using BenchmarkTools

efficiency()

#

nₚ = 11
nₑ = nₚ-1

coordinates = [(1/nₑ*i,0.0,0.0) for i in 0:nₑ]
node = Node(1,coordinates)
@btime $node.id
