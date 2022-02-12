
using Revise, ApproxOperator

elements, nodes = importmsh("./msh/patchtest.msh")
nₚ = length(nodes[:x])
nₑ = length(elements["Ω"])

set𝓖!(elements["Ω"],:TriGI3)
# set𝓖!(elements["Γᵗ₁"],:SegGI2)
# set𝓖!(elements["Γᵗ₂"],:SegGI2)
set𝓖!(elements["Γᵍ"],:SegGI2)
elements["Γᵍ"] = Element{:Tri3}(elements["Ω"],elements["Γᵍ"])

E = 1.0
ν = 0.0
prescribe!(elements["Ω"],:b₁,(x,y,z)->0.0)
prescribe!(elements["Ω"],:b₂,(x,y,z)->0.0)
# prescribe!(elements["Γᵗ₁"],:t₁,(x,y,z)->E/(1-ν))
# prescribe!(elements["Γᵗ₁"],:t₂,(x,y,z)->E/(1+ν))
# prescribe!(elements["Γᵗ₂"],:t₁,(x,y,z)->E/(1+ν))
# prescribe!(elements["Γᵗ₂"],:t₂,(x,y,z)->E/(1-ν))
prescribe!(elements["Γᵍ"],:g₁,(x,y,z)->1.0+x+y)
prescribe!(elements["Γᵍ"],:g₂,(x,y,z)->1.0+x+y)
prescribe!(elements["Γᵍ"],:n₁₁,(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂,(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₂₂,(x,y,z)->1.0)

coefficient = (:E=>E,:ν=>ν,:α=>1e7)
ops = [Operator(:∫∫εᵢⱼσᵢⱼdxdy,coefficient...),
       Operator(:∫∫vᵢbᵢdxdy,coefficient...),
       Operator(:∫vᵢtᵢds,coefficient...),
       Operator(:∫σᵢⱼnⱼgᵢds,coefficient...),
       Operator(:∫vᵢgᵢds,coefficient...),
       Operator(:Hₑ_PlaneStress,coefficient...)]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

ops[1](elements["Ω"],k)
# ops[3](elements["Γᵗ₁"],f)
# ops[3](elements["Γᵗ₂"],f)
ops[4](elements["Γᵍ"],k,f)
# ops[5](elements["Γᵍ"],k,f)

d = k\f

d₁ = zeros(nₚ)
d₂ = zeros(nₚ)
d₁ .= d[1:2:2*nₚ-1]
d₂ .= d[2:2:2*nₚ]
push!(nodes,:d₁=>d₁,:d₂=>d₂)
prescribe!(elements["Ω"],:u,(x,y,z)->1.0+x+y)
prescribe!(elements["Ω"],:v,(x,y,z)->1.0+x+y)
prescribe!(elements["Ω"],:∂u∂x,(x,y,z)->1.0)
prescribe!(elements["Ω"],:∂v∂x,(x,y,z)->1.0)
prescribe!(elements["Ω"],:∂u∂y,(x,y,z)->1.0)
prescribe!(elements["Ω"],:∂v∂y,(x,y,z)->1.0)
h1,l2 = ops[6](elements["Ω"])
