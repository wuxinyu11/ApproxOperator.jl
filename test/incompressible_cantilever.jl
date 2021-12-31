
using ApproxOperator

P = 1000.0
Ẽ = 3e6
# ν̃ = 0.3
ν̃ = 0.49999
E = Ẽ/(1-ν̃^2)
ν = ν̃/(1-ν̃)
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

ops = [
    # 1
    PlaneStress_Ω(b,E,ν),
    # 2
    PlaneStrain_Ωᵛ(Ẽ,ν̃),
    # 3
    PlaneStrain_Ωᵈ(b,Ẽ,ν̃),
    # 4
    PlaneStress_Γᵗ(t),
    # 5
    PlaneStress_Γᵍ_penalty(u,n,E*1e7),
    # 6
    HₑError_PlaneStress(u,E,ν),
    # 7
    VTKExport(
                filename="vtk/cantilever_locking.vtk",
                topic="incompressible test for cantilever beam problem",
                pointdata=Dict{String,Symbol}("disp"=>:u_2D),
                celldata=Dict{String,Symbol}("stress"=>:σ_PlaneStrain),
                parameters=Dict{String,Float64}("E"=>Ẽ,"ν"=>ν̃)
            )
]

## locking
k = zeros(2nₚ,2nₚ)
f = zeros(2nₚ)

ops[1](aps["Domain"],k,f)
ops[4](aps["Traction"],f)
ops[5](aps["EssentialBC"],k,f)

d = k\f

Hₑ_Locking,L₂_Locking = ops[6](aps["Domain"],d)

ops[7](aps["Domain"],d)

## locking-free
k = zeros(2nₚ,2nₚ)
f = zeros(2nₚ)

ops[3](aps["Domain"],k,f)
set_integration_rule!(aps["Domain"],:QuadGI1)
ops[2](aps["Domain"],k,f)
ops[4](aps["Traction"],f)
ops[5](aps["EssentialBC"],k,f)

d = k\f

set_integration_rule!(aps["Domain"],:QuadGI2)
Hₑ_Free,L₂_Free = ops[6](aps["Domain"],d)

ops[7].filename = "vtk/cantilever_free.vtk"
ops[7](aps["Domain"],d)
