
include("../mfea.jl")

using .MFEA

orders = ["Linear","Quadratic","Cubic"]
methods = ["Segmental element", "Linear meshfree", "Quadratic meshfree", "Cubic meshfree"]

L = 1.0
nₚ = 21
nₑ = nₚ-1
x = [Node(L/nₑ*i,0.,0.) for i in 0:nₑ]

# u = (x,y,z)-> 1+2x+3x^2

L₂Error = zeros(length(orders),length(methods))
H₁Error = zeros(length(orders),length(methods))

for m in 1:length(methods)
    aps = Vector{AbstractSeg}(undef,nₑ)
    if methods[m] == "Segmental element"
        aps = [Seg2(x,[i,i+1]) for i in 1:nₑ]
    elseif methods[m] == "Linear meshfree"
        bf = Linear1D(:∂1,:∂x,:∂y,:∂z)
        kf = TensorProductKernel(:∂1,:∂x,:∂y,:∂z,ss=[1.5*L/nₑ,1.,1.],nm=nₚ)
        for i in 1:nₑ
            aps[i] = SegM(x,vcat(i,i+1,collect(1:i-1),collect(i+2:nₚ)),bf=bf,kf=kf,qw=Gauss1_4Point)
        end
    elseif methods[m] == "Quadratic meshfree"
        bf = Quadratic1D(:∂1,:∂x,:∂y,:∂z)
        kf = TensorProductKernel(:∂1,:∂x,:∂y,:∂z,ss=[2.5*L/nₑ,1.,1.],nm=nₚ)
        for i in 1:nₑ
            aps[i] = SegM(x,vcat(i,i+1,collect(1:i-1),collect(i+2:nₚ)),bf=bf,kf=kf,qw=Gauss1_4Point)
        end
    elseif methods[m] == "Cubic meshfree"
        bf = Cubic1D(:∂1,:∂x,:∂y,:∂z)
        kf = TensorProductKernel(:∂1,:∂x,:∂y,:∂z,ss=[3.5*L/nₑ,1.,1.],nm=nₚ)
        for i in 1:nₑ
            aps[i] = SegM(x,vcat(i,i+1,collect(1:i-1),collect(i+2:nₚ)),bf=bf,kf=kf,qw=Gauss1_4Point)
        end
    end

    for o in 1:length(orders)
        if orders[o] == "Linear"
            u = (x,y,z)-> (1+2x,2,0,0)
        elseif orders[o] == "Quadratic"
            u = (x,y,z)-> (1+2x+3x^2,2+6x,0,0)
        elseif orders[o] == "Cubic"
            u = (x,y,z)-> (1+2x+3x^2+4x^3,2+6x+12x^2,0,0)
        end
        d = [u(1.0x[i]...)[1] for i in 1:nₚ]
        error_operator = H₁Error_scale(u)
        H₁, L₂ = error_operator(aps,d)
        L₂Error[o,m] = L₂
        H₁Error[o,m] = H₁
    end
end
