
include("../mfea.jl")
using .MFEA

P = 1000.0
E = 3e6
ν = 0.3
L = 48.0
D = 12.0
I = D^3/12
EI = E*I

u = (x,y,z) -> (
    -P*y/6/EI*((6*L-3x)*x + (2+ν)*(y^2-D^2/4)),
     P/6/EI*(3*ν*y^2*(L-x) + (4+5*ν)*D^2*x/4 + (3*L-x)*x^2),
    -P/EI*(L-x)*y,
    -P/6/EI*((6*L-3*x)*x + (2+ν)*(3*y^2-D^2/4)),
     P/6/EI*((6*L-3*x)*x - 3*ν*y^2 + (4+5*ν)*D^2/4),
     P/EI*(L-x)*y*ν
)
t = (x,y,z) -> (0.0,P/2/I*(D^2/4-y^2))
b = (x,y,z) -> (0.0,0.0)
n = (x,y,z) -> (1.0,0.0,1.0)

aps = import_msh("./msh/cantilever.msh")
nₚ = length(aps["Domain"][1].nodes)

k = zeros(2nₚ,2nₚ)
f = zeros(2nₚ)

ops = [
    PlaneStress_Ω(b,E,ν),
    PlaneStress_Γᵗ(t),
    PlaneStress_Γᵍ_penalty(u,n,E*1e7),
    HₑError_PlaneStress(u,E,ν)
]

ops[1](aps["Domain"],k,f)
ops[2](aps["Traction"],f)
ops[3](aps["EssentialBC"],k,f)

d = k\f

Hₑ_Error,L₂_Error = ops[4](aps["Domain"],d)

var = Dict{String,Symbol}("disp"=>:uv)
export_VTK(aps["Domain"],d,variables=var)
