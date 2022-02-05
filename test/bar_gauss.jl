
using Revise, ApproxOperator

elements, nodes = importmsh("./msh/bar.msh")
nₚ = length(nodes[:x])
nₑ = length(elements["Domain"])

sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z],n = 2,γ = 1)
elements["Domain"] = SegN{Node,:Quadratic1D,:□,:CubicSpline}(elements["Domain"])
elements["NBC"] = PoiN{Node,:Quadratic1D,:□,:CubicSpline}(elements["NBC"])
elements["EBC"] = PoiN{Node,:Quadratic1D,:□,:CubicSpline}(elements["EBC"])
sp(elements["Domain"])
sp(elements["NBC"])
s = 0.25*ones(nₚ)
nodes[:s₁] = s
nodes[:s₂] = s
nodes[:s₃] = s

set𝓖!(elements["Domain"],:SegGI5,:∂1,:∂x,:∂y,:∂z)
set𝓖!(elements["NBC"],:PoiGI1,:∂1)
set𝓖!(elements["EBC"],:PoiGI1,:∂1,:∂x,:∂y,:∂z)

elements["EBCL"] = glue(elements["EBC"][1],elements["Domain"][1])

r = 2
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
ops[3](elements["EBCL"],k,f)
d = k\f
push!(nodes,:d=>d)
h1, l2 = ops[4](elements["Domain"])
