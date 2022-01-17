# This file contains functions for 1D bar analysis, the whole domain denoted by Ω = (0,L) can be discretized by a set of nodes,

using Revise
using ApproxOperator

ip = Operator(:msh)
aps = ip("./msh/bar.msh")
nₚ = ip.nₚ

coefficients = Dict(:k=>1.0)
op = Operator(:∫∇v∇udΩ,coefficients)
op1 = Operator(:g,coefficients)
opn = Operator(:∫vtdΓ,coefficients)
r = 3
prescribe!(aps["NBC"],:t,(x,y,z)->r*x^abs(r-1))
prescribe!(aps["Domain"],:b,(x,y,z)->-r*abs(r-1)*x^abs(r-2))

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
d = zeros(nₚ)
push!(ip.nodes,:d=>d)
op(aps["Domain"],k,f)
op1(aps["EBC"],k,f)
opn(aps["NBC"],f)
d .= k\f
