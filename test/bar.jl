# This file contains functions for 1D bar analysis, the whole domain denoted by Ω = (0,L) can be discretized by a set of nodes,

using Revise
using ApproxOperator

opimport = Operator(:msh,Dict(:nodetype=>:Node,:pointtype=>Dict("EBC"=>:PoiGI1,"NBC"=>:PoiGI1,"Domain"=>:SegGI2)))
aps,datas = op("./msh/bar.msh")
nₚ = datas["nₚ"]
nₑ = datas["nₑ"]

coefficients = Dict(:k=>1.0)
op = Operator(:∫∇v∇udΩ,coefficients)
op1 = Operator(:g,coefficients)
opn = Operator(:∫vtdΓ,coefficients)
r = 1
prescribe = Operator(:prescribe,Dict(:b=>(x,y,z)->-r*abs(r-1)*x^abs(r-2),:t=>(x,y,z)->r*x^abs(r-1),:g=>(x,y,z)->0))
prescribe(aps["Domain"],:b)
prescribe(aps["EBC"],:g)
prescribe(aps["NBC"],:t)
k = zeros(nₚ,nₚ)
f = zeros(nₚ)
op(aps["Domain"],k,f)
op1(aps["EBC"],k,f)
op1(aps["NBC"],k,f)
d .= k\f
