
using Revise,ApproxOperator

elements, nodes = importmsh("./msh/cantilever.msh")
np = length(nodes[:x])

set𝓖!(elements["Ω"],:QuadGI2)
set𝓖!(elements["Γᵗ"],:SegGI2)
set𝓖!(elements["Γᵍ"],:SegGI2)

P = 1000.0
E = 3e6
ν = 0.3
L = 48.0
D = 12.0
I = D^3/12
EI = E*I
t₁(x,y,z) = 0.0
t₂(x,y,z) = P/2/I*(D^2/4-y^2)
g₁(x,y,z) = -P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4))
g₂(x,y,z) = P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2)
prescribe!(elements["Γᵗ"],:t₁,t₁)
prescribe!(elements["Γᵗ"],:t₂,t₂)
prescribe!(elements["Γᵍ"],:g₁,g₁)
prescribe!(elements["Γᵍ"],:g₂,g₂)
prescribe!(elements["Γᵍ"],:n₁₁,(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₂₂,(x,y,z)->1.0)

op_Ω = Operator(:∫∫εᵢⱼσᵢⱼdxdy,:E=>E,:ν=>ν)
op_Γᵗ = Operator(:∫vᵢtᵢds)
op_Γᵍ = Operator(:∫vᵢgᵢds,:α=>1e7*E)

k = zeros(2*np,2*np)
f = zeros(2*np)

op_Ω(elements["Ω"],k)
op_Γᵗ(elements["Γᵗ"],f)
op_Γᵍ(elements["Γᵍ"],k,f)

d = k\f

nodes[:d₁] = d[1:2:2*np-1]
nodes[:d₂] = d[2:2:2*np]

set𝓖!(elements["Ω"],:TriGI3)
prescribe!(elements["Ω"],:u,g₁)
prescribe!(elements["Ω"],:v,g₂)
prescribe!(elements["Ω"],:∂u∂x,(x,y,z)->-P/EI*(L-x)*y)
prescribe!(elements["Ω"],:∂u∂y,(x,y,z)->-P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)))
prescribe!(elements["Ω"],:∂v∂x,(x,y,z)->P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4))
prescribe!(elements["Ω"],:∂v∂y,(x,y,z)->P/EI*(L-x)*y*ν)
op = Operator(:Hₑ_PlaneStress,:E=>E,:ν=>ν)
h1,l2 = op(elements["Ω"])

l2 = log10(l2)
h1 = log10(h1)
