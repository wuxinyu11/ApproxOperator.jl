# This file contains functions for 1D bar analysis, the whole domain denoted by Ω = (0,L) can be discretized by a set of nodes,

using Revise
using ApproxOperator

nₚ = 11
nₑ = nₚ - 1

# data
data = Dict{Symbol,Vector{Float64}}()
push!(data,:x=>[1.0/nₑ*i for i in 0:nₑ])
push!(data,:y=>zeros(nₚ))
push!(data,:z=>zeros(nₚ))
push!(data,:ξ=>[-0.5773502691896257645091487805,0.5773502691896257645091487805])
push!(data,:w=>[1.0,1.0])


ap1 = [Seg2([Node(i,nodedata),Node(i+1,nodedata)],[Point(i,1,pointdata),Point(i,2,pointdata)]) for i in 1:nₑ]
