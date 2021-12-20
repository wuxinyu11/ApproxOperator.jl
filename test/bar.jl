# This file contains functions for 1D bar analysis, the whole domain denoted by Ω = (0,L) can be discretized by a set of nodes,

include("../mfea.jl")

using .MFEA, Plots

L = 1. # length of bar
# b = (x,y,z)-> -2.
# t = (x,y,z)-> 2*x
# u = (x,y,z)-> x^2
# ū = (x,y,z)-> (x^2,2*x,0,0)

b = (x,y,z)-> -20*x^3
t = (x,y,z)-> 5*x^4
u = (x,y,z)-> x^5
ū = (x,y,z)-> (x^5,5*x^4,0,0)

ℕ = [5, 10, 20, 40, 80]

methods = ["Segmental element", "Linear meshfree", "Quadratic meshfree", "Cubic meshfree"]
# methods = ["Segmental element","Linear meshfree"]
# methods = ["Linear meshfree"]


h = zeros(length(ℕ))
L₂Error = zeros(length(ℕ),length(methods))
H₁Error = zeros(length(ℕ),length(methods))
L₂slope = zeros(length(ℕ)-1,length(methods))
H₁slope = zeros(length(ℕ)-1,length(methods))
for n in 1:length(ℕ)
    nₑ = ℕ[n] # number of element
    nₚ = nₑ+1 # number of nodes
    x = [Node(L/nₑ*i,0.,0.) for i in 0:nₑ]
    h[n] = L/nₑ
    # ops = [Potential_Ω(b), Potential_Γᵗ(t), Potential_Γᵍ_penalty(u,1e7), H₁Error_scale(ū)]
    # ops = [Potential_Ω(b), Potential_Γᵗ(t), EBCDOFS(u), H₁Error_scale(ū)]
    ops = [Potential_Ω(b), Potential_Γᵗ(t), Potential_Γᵍ_Nitsche(u), H₁Error_scale(ū)]
    for m = 1:length(methods)
        k = zeros(nₚ,nₚ)
        f = zeros(nₚ)
        if methods[m] == "Segmental element"
            aps = [Seg2(x,[i,i+1]) for i in 1:nₑ]
            ape = [Seg2(x,[i,i+1]) for i in 1:nₑ]
            ap1 = Poi1(x,1)
            apnₚ = Poi1(x,nₚ)
        elseif methods[m] == "Linear meshfree"
            rg = RegularGrid(x,n=1,γ=1)
            bf = Linear1D(:∂1,:∂x,:∂y,:∂z)
            kf = TensorProductKernel(:∂1,:∂x,:∂y,:∂z,ss=[1.5*L/nₑ,1.,1.],nm=nₚ)
            aps = Vector{AbstractSeg}(undef,nₑ)
            for i in 1:nₑ
                aps[i] = SegM(x,[i,i+1],bf=bf,kf=kf,qw=:SegGI2,sp=rg)
            end
            ape = Vector{AbstractSeg}(undef,nₑ)
            for i in 1:nₑ
                ape[i] = SegM(x,[i,i+1],bf=bf,kf=kf,qw=:SegGI10,sp=rg)
            end
            ap1 = PoiM(x,1,bf = bf,kf = kf,sp=rg)
            apnₚ = PoiM(x,nₚ,bf = bf,kf = kf,sp=rg)
        elseif methods[m] == "Quadratic meshfree"
            rg = RegularGrid(x,n=2,γ=1)
            bf = Quadratic1D(:∂1,:∂x,:∂y,:∂z)
            kf = TensorProductKernel(:∂1,:∂x,:∂y,:∂z,ss=[2.5*L/nₑ,1.,1.],nm=nₚ)
            aps = Vector{AbstractSeg}(undef,nₑ)
            for i in 1:nₑ
                aps[i] = SegM(x,[i,i+1],bf=bf,kf=kf,qw=:SegGI5,sp=rg)
            end
            ape = aps
            ap1 = PoiM(x,1,bf = bf,kf = kf,sp=rg)
            apnₚ = PoiM(x,nₚ,bf = bf,kf = kf,sp=rg)
        elseif methods[m] == "Cubic meshfree"
            rg = RegularGrid(x,n=3,γ=1)
            bf = Cubic1D(:∂1,:∂x,:∂y,:∂z)
            kf = TensorProductKernel(:∂1,:∂x,:∂y,:∂z,ss=[3.5*L/nₑ,1.,1.],nm=nₚ)
            aps = Vector{AbstractSeg}(undef,nₑ)
            for i in 1:nₑ
                aps[i] = SegM(x,[i,i+1],bf=bf,kf=kf,qw=:SegGI10,sp=rg)
            end
            ape = Vector{AbstractSeg}(undef,nₑ)
            for i in 1:nₑ
                ape[i] = SegM(x,[i,i+1],bf=bf,kf=kf,qw=:SegGI10,sp=rg)
            end
            ap1 = PoiM(x,1,bf = bf,kf = kf,sp=rg)
            apnₚ = PoiM(x,nₚ,bf = bf,kf = kf,sp=rg)
        end
        # assembling
        ops[1](aps,k,f)
        ops[2](apnₚ,f)
        # ops[3](apnₚ,k,f)
        # ops[3](ap1,k,f)
        ops[3](aps[1],ap1,k,f)

        # solve
        d = k\f
        # error
        H₁Error_, L₂Error_ = ops[4](ape,d)
        L₂Error[n,m] = L₂Error_
        H₁Error[n,m] = H₁Error_
    end
end
for i in 1:length(ℕ)-1
    L₂slope[i,:] .= (log10.(L₂Error[i+1,:]) - log10.(L₂Error[i,:]))./(log10.(h[i+1]) - log10.(h[i]))
    H₁slope[i,:] .= (log10.(H₁Error[i+1,:]) - log10.(H₁Error[i,:]))./(log10.(h[i+1]) - log10.(h[i]))
end
L₂Fig = plot(log10.(h),log10.(L₂Error))
H₁Fig = plot(log10.(h),log10.(H₁Error))
display(L₂Fig)
display(H₁Fig)
