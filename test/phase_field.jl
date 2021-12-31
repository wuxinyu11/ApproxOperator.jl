
using ApproxOperator, LinearAlgebra

aps = import_msh("./msh/cube_pfm.msh")
nₚ = length(aps["surface"][1].nodes)
set_integration_rule!(aps["surface"],Val(:PhaseField))

E = 2.0e11
ν = 0.3
θ = 0.
u₀ = 5.0e-3
nls = 50
k̄ = 2700
l = 0.01
tol = 1e-12
maxiter = 50

b = (x,y,z)->(0.0,0.0)
n = (x,y,z)->(1.0,0.0,1.0)
ū = (x,y,z)->(0.0,0.0,0.0)

ops = [
    PlaneStress_PhaseField_Ω(b,E,ν),
    NonlinearPlaneStress_Γᵍ_penalty(ū,n,E*1e7),
    SecondOrderPhaseField(k̄,l),
    Update_HistoryField_PlaneStress(E,ν),
    VTKExport(
                filename="figure/movie00.vtk",
                topic="phase field modeling fracture",
                parameters=Dict{String,Float64}("E"=>E,"ν"=>ν)
            )
]

d = zeros(2*nₚ)
k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
kᵥ = zeros(nₚ,nₚ)
fᵥ = zeros(nₚ)
dᵥ = zeros(nₚ)
dᵗ = zeros(nₚ)

ops[5](aps["surface"],d,dᵥ,Val(:PFM_PlaneStress))

for t in 1:nls
    # define essential BC
    u = (x,y,z)->(cos(θ)*u₀/nls*t,sin(θ)*u₀/nls*t)
    op_EBC = NonlinearPlaneStress_Γᵍ_penalty(u,n,E*1e7)
    iter = 0

    # Alternate minimization algorithm
    while iter ≤ maxiter
        dᵗ .= dᵥ
        iter += 1
        # loop of elasticity
        iteru = 0
        while iteru ≤ maxiter
            iteru += 1
            fill!(k,0.0)
            fill!(f,0.0)
            ops[1](aps["surface"],k,f,d,dᵥ)
            ops[2](aps["disp_fix"],k,f,d)
            op_EBC(aps["disp"],k,f,d)
            erru = norm(f)/E/1e7
            print("u: time step = $t, alternate iter = $iter, disp iter = $iteru, erru = $erru\n")
            if erru < tol;break;end
            d .+= k\f
        end

        # update history field
        ops[4](aps["surface"],d)

        # loop of phase field
        iterv = 0
        while iterv ≤ maxiter
            iterv += 1
            fill!(kᵥ,0.0)
            fill!(fᵥ,0.0)
            ops[3](aps["surface"],kᵥ,fᵥ,dᵥ)
            errv = norm(fᵥ)/k̄
            print("v: time step = $t, alternate iter = $iter, damage iter = $iterv, errv = $errv\n")
            if errv < tol;break;end
            dᵥ .+= kᵥ\fᵥ
        end

        err = norm(dᵗ .- dᵥ)/norm(dᵥ)
        print("err = $err\n")
        if err < tol;break;end
    end
    for ap in aps["surface"]
        for qw in ap.qw
            qw.ℋₜ = qw.ℋ
        end
    end
    # export to VTK file
    ops[5].filename = "figure/movie"*string(t,pad=2)*".vtk"
    ops[5](aps["surface"],d,dᵥ,Val(:PFM_PlaneStress))
end
