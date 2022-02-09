
using Revise, ApproxOperator

elements, nodes = importmsh("./msh/bar.msh")
nₚ = length(nodes[:x])
nₑ = length(elements["Domain"])

type = (SNode,:Quadratic1D,:□,:CubicSpline)
sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z],n = 2,γ = 1)
elements["Domain"] = ReproducingKernel{type...,:Seg2}(elements["Domain"])
elements["NBC"] = ReproducingKernel{type...,:Poi1}(elements["NBC"])
elements["EBC"] = ReproducingKernel{type...,:Poi1}(elements["EBC"])
sp(elements["Domain"])
sp(elements["NBC"])
sp(elements["EBC"])
elements["DomainS"] = ReproducingKernel{type...,:Seg2}(elements["Domain"])
s = 0.25*ones(nₚ)
nodes[:s₁] = s
nodes[:s₂] = s
nodes[:s₃] = s

set𝓖!(elements["Domain"],:SegRK3,:∂1,:∂x,:∂y,:∂z)
set𝓖!(elements["DomainS"],:SegGI2,:∂1,:∂x,:∂y,:∂z)
set𝓖!(elements["NBC"],:PoiGI1,:∂1)
set𝓖!(elements["EBC"],:PoiGI1,:∂1,:∂x,:∂y,:∂z)

elements["EBC"] = ReproducingKernel{type...,:Seg2}(elements["Domain"],elements["EBC"],sharing=true)
elements["EBCD"] = elements["Domain"]∩elements["EBC"]
elements["EBCS"] = elements["DomainS"]∩elements["EBC"]

r = 3
prescribe!(elements["Domain"],:b,(x,y,z)->-r*(r-1)*x^abs(r-2))
prescribe!(elements["NBC"],:t,(x,y,z)->r*x^abs(r-1))
prescribe!(elements["EBC"],:g,(x,y,z)->x^r)

set𝝭!(elements["Domain"])
set∇̃𝝭!(elements["DomainS"],elements["Domain"])
set𝝭!(elements["NBC"])
set∇̃𝝭!(elements["EBCD"])
setg̃!(elements["EBCS"],elements["EBC"])

coefficient = (:k=>1.0,:α=>1e3)
ops = [Operator(:∫∇v∇udΩ,coefficient...),
       Operator(:∫vbdΩ,coefficient...),
       Operator(:∫vtdΓ,coefficient...),
       Operator(:∫∇𝑛vgdΓ,coefficient...),
       Operator(:∫vgdΓ,coefficient...),
       Operator(:H₁)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["DomainS"],k)
ops[2](elements["Domain"],f)
ops[3](elements["NBC"],f)
ops[4](elements["EBC"],k,f)
ops[5](elements["EBCS"],k,f)
d = k\f

push!(nodes,:d=>d)
set𝓖!(elements["Domain"],:SegGI5,:∂1,:∂x,:∂y,:∂z)
prescribe!(elements["Domain"],:u,(x,y,z)->x^r)
prescribe!(elements["Domain"],:∂u∂x,(x,y,z)->r*x^abs(r-1))
set∇𝝭!(elements["Domain"])
h1, l2 = ops[6](elements["Domain"])
