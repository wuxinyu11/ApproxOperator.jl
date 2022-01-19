# This file contains functions for 1D bar analysis, the whole domain denoted by Ω = (0,L) can be discretized by a set of nodes,

using Revise
using ApproxOperator
using BenchmarkTools

ip = Operator(:msh)
push!(ip,:basisfunction=>:Linear1D,:kerneltype=>:□,:kernelfunction=>:CubicSpline)
push!(ip,:stype=>[:∂1,:∂x,:∂y,:∂z],:spatialpartition=>:RegularGrid,:nᵣ=>1,:γᵣ=>1)
ip.etype[1] = :SegN
ip.etype[15] = :PoiN
ip.qtype[1] = :SegGI5

ip.ntype[1] = :SNode
ip.qtype[1] = :SegRK3
push!(ip,:index=>Int[0],:shapefunctions=>Dict(:∂1=>Float64[],:∂x=>Float64[],:∂y=>Float64[],:∂z=>Float64[]))

aps = ip("./msh/bar.msh")
nₚ = ip.nₚ
s = 0.15.*ones(nₚ)
push!(ip.nodes,:s₁=>s,:s₂=>s,:s₃=>s)

# cal𝝭 = Operator(:∇𝝭,Any)
cal𝝭 = Operator(:𝝭checkrepeat,Any)
push!(cal𝝭,:id=>Dict{NTuple{3,Float64},Int}(),:ids=>Int[],:index=>[0])
cal𝝭(aps["Domain"])
# @btime $cal𝝭($aps["Domain"])

coefficients = Dict(:k=>1.0,:α=>1e7)
op = Operator(:∫∇v∇uvbdΩ,coefficients)
op1 = Operator(:∫vgdΓ,coefficients)
opn = Operator(:∫vtdΓ,coefficients)
r = 3
prescribe!(aps["NBC"],:t,(x,y,z)->r*x^abs(r-1))
prescribe!(aps["EBC"],:g,(x,y,z)->x^r)
prescribe!(aps["Domain"],:b,(x,y,z)->-r*abs(r-1)*x^abs(r-2))

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
d = zeros(nₚ)
push!(ip.nodes,:d=>d)
op(aps["Domain"],k,f)
op1(aps["EBC"],k,f)
opn(aps["NBC"],f)
d .= k\f
