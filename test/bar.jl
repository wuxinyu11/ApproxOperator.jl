# This file contains functions for 1D bar analysis, the whole domain denoted by Ω = (0,L) can be discretized by a set of nodes,

using Revise
using ApproxOperator
using BenchmarkTools

## FEM
# nodes,elements = importdata("./msh/bar.msh")
# set𝓖(elements["Domain"],:SegGI2)
# set𝓖(elements["EBC"],:PoiGI1)
# set𝓖(elements["NBC"],:PoiGI1)

## MESHFREE
# nodes,elements = importdata("./msh/bar.msh",𝒑=:Quadratic1D,𝑠=:□,𝜙=:CubicSpline)
# nₚ = length(nodes[:x])
# sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z],n=2,γ=1)
# sp(elements["Domain"])
# sp(elements["EBC"])
# sp(elements["NBC"])
# s = 0.25*ones(nₚ)
# push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
# set𝓖(elements["Domain"],:SegGI5,:∂1,:∂x,:∂y,:∂z)
# set𝓖(elements["EBC"],:PoiGI1,:∂1)
# set𝓖(elements["NBC"],:PoiGI1,:∂1)

# nodes,elements = importdata("./msh/bar.msh",1=>:SNode,𝒑=:Linear1D,𝑠=:□,𝜙=:CubicSpline)
# nₚ = length(nodes[:x])
# sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z],n=1,γ=1)
# sp(elements["Domain"])
# sp(elements["EBC"])
# sp(elements["NBC"])
# s = 0.15*ones(nₚ)
# push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
# set𝓖(elements["Domain"],:SegGI5,:∂1,:∂x,:∂y,:∂z,ntype=:SNode)
# set𝓖(elements["EBC"],:PoiGI1,:∂1)
# set𝓖(elements["NBC"],:PoiGI1,:∂1)
# set∇𝝭(elements["Domain"])

## RK
nodes, elements = importdata(
  "./msh/bar.msh",
  1 => :SNode,
  𝒑 = :Quadratic1D,
  𝑠 = :□,
  𝜙 = :CubicSpline,
)
nₚ = length(nodes[:x])
sp = RegularGrid(nodes[:x], nodes[:y], nodes[:z], n = 2, γ = 1)
sp(elements["Domain"])
sp(elements["EBC"])
sp(elements["NBC"])
s = 0.25 * ones(nₚ)
push!(nodes, :s₁ => s, :s₂ => s, :s₃ => s)
elements["DomainGS"] = SegN{SNode,:Linear1D,:□,:CubicSpline}(elements["Domain"])
elements["LEBC"] = Poi1(elements["EBC"], renumbering = true)
set𝓖(elements["Domain"], :SegRK3, :∂1, :∂x, :∂y, :∂z, ntype = :SNode)
set𝓖(elements["DomainGS"], :SegGI2, :∂1, :∂x, :∂y, :∂z, ntype = :SNode)
set𝓖(elements["EBC"], :PoiGI1, :∂1)
set𝓖(elements["NBC"], :PoiGI1, :∂1)
set𝓖(elements["LEBC"], :PoiGI1)
prescribe!(elements["Domain"], :n₁, (x, y, z) -> 0.0)
set𝒏(elements["Domain"])
set𝝭(elements["Domain"])
set∇̃𝝭(elements["DomainGS"], elements["Domain"])

coefficients = (:k => 1.0, :α => 1e7)
# op = Operator(:∫∇v∇uvbdΩ,coefficients...)
op = Operator(:∫vbdΩ, coefficients...)
opgs = Operator(:∫∇v∇udΩ, coefficients...)
op1 = Operator(:∫vgdΓ, coefficients...)
opn = Operator(:∫vtdΓ, coefficients...)
opl = Operator(:∫λudΓ, coefficients...)
r = 3
prescribe!(elements["NBC"], :t, (x, y, z) -> r * x^abs(r - 1))
prescribe!(elements["EBC"], :g, (x, y, z) -> x^r)
prescribe!(elements["Domain"], :b, (x, y, z) -> -r * abs(r - 1) * x^abs(r - 2))

k = zeros(nₚ, nₚ)
f = zeros(nₚ)
# d = zeros(nₚ)
d = zeros(nₚ + 1)
q = [0.0]
g = zeros(nₚ, 1)
push!(nodes, :d => d)
# op(elements["Domain"],k,f)
op(elements["Domain"], f)
opgs(elements["DomainGS"], k)
# op1(elements["EBC"],k,f)
opn(elements["NBC"], f)
opl(elements["EBC"], elements["LEBC"], g, q)
# d .= k\f
d .= [k g; g' 0] \ [f; q]

set𝓖(elements["Domain"], :SegGI10, :∂1, :∂x, :∂y, :∂z, ntype = :SNode)
set∇𝝭(elements["Domain"])
prescribe!(elements["Domain"], :u, (x, y, z) -> x^r)
prescribe!(elements["Domain"], :∂u∂x, (x, y, z) -> r * x^abs(r - 1))
prescribe!(elements["Domain"], :∂u∂y, (x, y, z) -> 0.0)
prescribe!(elements["Domain"], :∂u∂z, (x, y, z) -> 0.0)
l2 = Operator(:L₂)
h1 = Operator(:H₁)
l2error = l2(elements["Domain"])
h1error, l2error = h1(elements["Domain"])
