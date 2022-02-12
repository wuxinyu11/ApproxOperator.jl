
using Revise, ApproxOperator

elements, nodes = importmsh("./msh/bar.msh")
nₚ = length(nodes[:x])
nₑ = length(elements["Domain"])

set𝓖!(elements["Domain"],:SegGI2)
set𝓖!(elements["NBC"],:PoiGI1)
set𝓖!(elements["EBC"],:PoiGI1)
elements["EBC"] = Element{:Seg2}(elements["Domain"][1],elements["EBC"][1])

r = 3
prescribe!(elements["Domain"],:u,(x,y,z)->x^r)
prescribe!(elements["Domain"],:∂u∂x,(x,y,z)->r*x^abs(r-1))
prescribe!(elements["Domain"],:b,(x,y,z)->-r*(r-1)*x^abs(r-2))
prescribe!(elements["NBC"],:t,(x,y,z)->r*x^abs(r-1))
prescribe!(elements["EBC"],:g,(x,y,z)->x^r)

coefficient = (:k=>1.0,:α=>1e3)
ops = [Operator(:∫∇v∇uvbdΩ,coefficient...),
       Operator(:∫vtdΓ,coefficient...),
       Operator(:∫∇𝑛vgdΓ,coefficient...),
       Operator(:H₁)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Domain"],k,f)
ops[2](elements["NBC"],f)
ops[3](elements["EBC"],k,f)
d = k\f
push!(nodes,:d=>d)
h1, l2 = ops[4](elements["Domain"])
