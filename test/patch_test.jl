
include("../mfea.jl")

using .MFEA, BenchmarkTools, Plots

a₁ = 1.
a₂ = 1.

b = (x,y,z)->0.
t₁ = (x,y,z)->2
t₂ = (x,y,z)->3
u = (x,y,z)->1+2x+3y
ū = (x,y,z)->(1+2x+3y,2,3,0.)

aps = import("./msh/patchtest.msh")
nₚ = 9

# ops = [Potential_Ω(b), Potential_Γᵍ_Nitsche(u), H₁Error_scale(ū)]
ops = [Potential_Ω(b), Potential_Γᵗ(t₁), Potential_Γᵗ(t₂), Potential_Γᵍ_penalty(u,1e7), H₁Error_scale(ū)]
# ops = [Potential_Ω(b), EBCDOFS(u), H₁Error_scale(ū)]
k = zeros(nₚ,nₚ)
f = zeros(nₚ)
# assembling
ops[1](aps["Domain"],k,f)
ops[2](aps["Traction1"],f)
ops[3](aps["Traction2"],f)
ops[4](aps["EssentialBC"],k,f)
# ops[2](aps_Ω̄,aps_Γ,k,f)
# solve
d = k\f
# error

H₁Error, L₂Error = ops[5](aps["Domain"],d)
