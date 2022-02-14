
using Revise, ApproxOperator

elements, nodes = importmsh("./msh/patchtest.msh")
nₚ = length(nodes[:x])
nₑ = length(elements["Ω"])

elements["Γ̃ᵍ"] = Element{:Seg2}(elements["Γᵍ"])
set𝓖!(elements["Ω"],:TriGI3)
set𝓖!(elements["Γᵍ"],:SegGI2)
set𝓖!(elements["Γ̃ᵍ"],:SegGI2)
elements["Γ̃ᵍ"] = Element{:Tri3}(elements["Ω"],elements["Γ̃ᵍ"])

prescribe!(elements["Ω"],:b,(x,y,z)->0.0)
# prescribe!(elements["Γᵗ₁"],:t₁,(x,y,z)->E/(1-ν))
# prescribe!(elements["Γᵗ₁"],:t₂,(x,y,z)->E/(1+ν))
# prescribe!(elements["Γᵗ₂"],:t₁,(x,y,z)->E/(1+ν))
# prescribe!(elements["Γᵗ₂"],:t₂,(x,y,z)->E/(1-ν))
prescribe!(elements["Γᵍ"],:g,(x,y,z)->1.0+x+y)
prescribe!(elements["Γ̃ᵍ"],:g,(x,y,z)->1.0+x+y)

coefficient = (:k=>1.0,:α=>0e0)
ops = [Operator(:∫∇v∇udΩ,coefficient...),
       Operator(:∫vbdΩ,coefficient...),
       Operator(:∫vtds,coefficient...),
       Operator(:∫∇𝑛vgdΓ,coefficient...),
       Operator(:∫vgdΓ,coefficient...),
       Operator(:H₁,coefficient...)]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](elements["Ω"],k)
ops[4](elements["Γ̃ᵍ"],k,f)
ops[5](elements["Γᵍ"],k,f)

x = nodes[:x]
y = nodes[:y]
test = y'*f

d = k\f

push!(nodes,:d=>d)
prescribe!(elements["Ω"],:u,(x,y,z)->1.0+x+y)
prescribe!(elements["Ω"],:∂u∂x,(x,y,z)->1.0)
prescribe!(elements["Ω"],:∂u∂y,(x,y,z)->1.0)
h1,l2 = ops[6](elements["Ω"])
