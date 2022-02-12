
using Revise, ApproxOperator

elements, nodes = importmsh("./msh/patchtest.msh")
nₚ = length(nodes[:x])
nₑ = length(elements["Ω"])

set𝓖!(elements["Ω"],:TriGI3)
# set𝓖!(elements["Γᵗ₁"],:SegGI2)
# set𝓖!(elements["Γᵗ₂"],:SegGI2)
set𝓖!(elements["Γᵍ"],:SegGI2)
elements["Γᵍ"] = Element{:Tri3}(elements["Ω"],elements["Γᵍ"])

prescribe!(elements["Ω"],:b,(x,y,z)->0.0)
# prescribe!(elements["Γᵗ₁"],:t₁,(x,y,z)->E/(1-ν))
# prescribe!(elements["Γᵗ₁"],:t₂,(x,y,z)->E/(1+ν))
# prescribe!(elements["Γᵗ₂"],:t₁,(x,y,z)->E/(1+ν))
# prescribe!(elements["Γᵗ₂"],:t₂,(x,y,z)->E/(1-ν))
prescribe!(elements["Γᵍ"],:g,(x,y,z)->1.0+x+y)

coefficient = (:k=>1.0,:α=>1e3)
ops = [Operator(:∫∇v∇udΩ,coefficient...),
       Operator(:∫vbdΩ,coefficient...),
       Operator(:∫vtds,coefficient...),
       Operator(:∫∇𝑛vgdΓ,coefficient...),
       Operator(:∫vgdΓ,coefficient...),
       Operator(:H₁,coefficient...)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω"],k)
# ops[3](elements["Γᵗ₁"],f)
# ops[3](elements["Γᵗ₂"],f)
ops[4](elements["Γᵍ"],k,f)
# ops[5](elements["Γᵍ"],k,f)

d = k\f

push!(nodes,:d=>d)
prescribe!(elements["Ω"],:u,(x,y,z)->1.0+x+y)
prescribe!(elements["Ω"],:∂u∂x,(x,y,z)->1.0)
prescribe!(elements["Ω"],:∂u∂y,(x,y,z)->1.0)
h1,l2 = ops[6](elements["Ω"])
