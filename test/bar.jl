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
push!(data,:ξ=>Float64[])
push!(data,:w=>Float64[])
for i in 1:nₑ
    push!(data[:ξ],-0.5773502691896257645091487805,0.5773502691896257645091487805)
    push!(data[:w],1.0,1.0)
end
push!(data[:ξ],0.0,0.0)
push!(data[:w],1.0,1.0)
push!(data,:b=>zeros(length(data[:w])))
push!(data,:t=>zeros(length(data[:w])))
push!(data,:g=>zeros(length(data[:w])))

ap = [Seg2([Node(i,data),Node(i+1,data)],[Node(i,data),Node(i+1,data)]) for i in 1:nₑ]
ap1 = Poi1([Node(1,data)],[Node(2*nₑ+1,data)])
apn = Poi1([Node(nₚ,data)],[Node(2*nₑ+2,data)])

coefficients = Dict(:k=>1.0)
op = Operator(:∫∇v∇udΩ,coefficients)
op1 = Operator(:g,coefficients)
opn = Operator(:∫vtdΓ,coefficients)


k = zeros(nₚ,nₚ)
f = zeros(nₚ)
d = zeros(nₚ)
d[nₚ] = 1.0
push!(data,:d=>d)
op(ap,k,f)
op1(apn,k,f)
op1(ap1,k,f)
d .= k\f
