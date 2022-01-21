# This file contains functions for 1D bar analysis, the whole domain denoted by Ω = (0,L) can be discretized by a set of nodes,

using Revise
using ApproxOperator
using BenchmarkTools

ip = Operator(:msh)
ip(:SNode)
ip(:Linear1D,:□,:CubicSpline)
aps,x = ip("./msh/bar.msh")
sp = RegularGrid(x[:x],x[:y],x[:z],n=1,γ=1)
sp(aps["Domain"])
sp(aps["EBC"])
sp(aps["NBC"])
# set𝓖_Ω = Operator(:𝓖)
# set𝓖_Γᵗ = Operator(:𝓖)
# set𝓖_Γᵍ = Operator(:𝓖)
set𝓖_Ω = Operator(:𝓖,:𝝭=>[:∂1,:∂x,:∂y,:∂z])
set𝓖_Γᵗ = Operator(:𝓖,:𝝭=>[:∂1])
set𝓖_Γᵍ = Operator(:𝓖,:𝝭=>[:∂1])
# set𝓖_Ω(aps["Domain"],:SegGI2)
data_Ω = set𝓖_Ω(aps["Domain"],:SegRK3)
data_Γᵗ = set𝓖_Γᵗ(aps["NBC"],:PoiGI1)
data_Γᵍ = set𝓖_Γᵍ(aps["EBC"],:PoiGI1)
nₚ = ip.nₚ
s = 0.15.*ones(nₚ)
push!(x,:s₁=>s,:s₂=>s,:s₃=>s)
cal𝝭_Ω = Operator(:𝝭)
# cal𝝭_Ω = Operator(:𝝭ʳ)
cal𝝭_Γ = Operator(:𝝭)
# cal𝝭 = Operator(:𝝭checkrepeat,Any)
# push!(cal𝝭,:id=>Dict{NTuple{3,Float64},Int}(),:ids=>Int[],:index=>[0])
cal𝝭_Ω(aps["Domain"])
cal𝝭_Γ(aps["NBC"])
cal𝝭_Γ(aps["EBC"])
push!(data_Ω,:n₁=>similar(data_Ω[:w]))
get𝒏(aps["Domain"])
# @btime $cal𝝭_Ω($aps["Domain"])

# coefficients = (:k=>1.0,:α=>1e7)
# op = Operator(:∫∇v∇uvbdΩ,coefficients...)
# op1 = Operator(:∫vgdΓ,coefficients...)
# opn = Operator(:∫vtdΓ,coefficients...)
# r = 3
# prescribe!(aps["NBC"],:t,(x,y,z)->r*x^abs(r-1))
# prescribe!(aps["EBC"],:g,(x,y,z)->x^r)
# prescribe!(aps["Domain"],:b,(x,y,z)->-r*abs(r-1)*x^abs(r-2))
#
# k = zeros(nₚ,nₚ)
# f = zeros(nₚ)
# d = zeros(nₚ)
# push!(ip.nodes,:d=>d)
# op(aps["Domain"],k,f)
# op1(aps["EBC"],k,f)
# opn(aps["NBC"],f)
# d .= k\f
