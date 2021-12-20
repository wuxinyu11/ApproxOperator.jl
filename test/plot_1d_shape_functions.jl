
include("../mfea.jl")

using .MFEA, Plots

nₚ = 11
nₑ = nₚ-1
x = [Node(1/nₑ*i,0.,0.) for i in 0:nₑ]

# methods = ["Segmental element", "Linear meshfree", "Quadratic meshfree", "Cubic meshfree"]
# methods = ["Segmental element"]
methods = ["Linear meshfree"]
# methods = ["Quadratic meshfree"]
for m in methods
    aps = Vector{AbstractSeg}(undef,nₑ)
## Finite element shape functions
    if m == "Segmental element"
        for i in 1:nₑ
            aps[i] = Seg2(x,[i,i+1])
        end

## Meshfree shpae functions
    elseif m == "Linear meshfree"

        bf = Linear1D(:∂1,:∂x,:∂y,:∂z)
        kf = TensorProductKernel(:∂1,:∂x,:∂y,:∂z,ss=[0.15,1.,1.],nm=nₚ)
        for i in 1:nₑ
            aps[i] = SegM(x,vcat(i,i+1,collect(1:i-1),collect(i+2:nₚ)),bf=bf,kf=kf)
        end
    elseif m == "Quadratic meshfree"

        bf = Quadratic1D(:∂1,:∂x)
        kf = TensorProductKernel(:∂1,:∂x,ss=[0.25,1.,1.],nm=nₚ)
        for i in 1:nₑ
            aps[i] = SegM(x,vcat(i,i+1,collect(1:i-1),collect(i+2:nₚ)),bf=bf,kf=kf)
        end

    elseif m == "Cubic meshfree"

        bf = Cubic1D(:∂1,:∂x)
        kf = TensorProductKernel(:∂1,:∂x,ss=[0.35,1.,1.],nm=nₚ)
        for i in 1:nₑ
            aps[i] = SegM(x,vcat(i,i+1,collect(1:i-1),collect(i+2:nₚ)),bf=bf,kf=kf)
        end
    end

    xᵢ, N, ∂N∂x = export_shape_functions(aps,:∂1,:∂x)
    # xᵢ, N, ∂N∂x, ∂N∂y, ∂N∂z = export_shape_functions(aps,:∂1,:∂x,:∂y,:∂z)

    plt_N = plot(xᵢ,N)
    plt_∂N∂x = plot(xᵢ,∂N∂x)

    display(plt_N)
    display(plt_∂N∂x)
end
