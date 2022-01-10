# This file contains functions for 1D bar analysis, the whole domain denoted by Ω = (0,L) can be discretized by a set of nodes,

using Revise
using ApproxOperator

nₚ = 11
nₑ = nₚ - 1
x = [1.0/nₑ*i for i in 0:nₑ]
ξ = [-0.5773502691896257645091487805,0.5773502691896257645091487805]
w = [1.0,1.0]

nodedata = Dict(:x=>x,:y=>zeros(nₚ),:z=>zeros(nₚ))
pointdata = Dict(:ξ=>ξ,:w=>w,:b=>zeros(2*nₑ))

ap1 = [Seg2([Node(i,nodedata),Node(i+1,nodedata)],[Point(i,1,pointdata),Point(i,2,pointdata)]) for i in 1:nₑ]
