
include("../mfea.jl")

using .MFEA, BenchmarkTools, Plots

a₁ = 1.
a₂ = 1.

t1 = (x,y,z)-> -exp(x+y+z)
t2 = (x,y,z)-> exp(x+y+z)
t3 = (x,y,z)-> exp(x+y+z)
t4 = (x,y,z)-> -exp(x+y+z)
b = (x,y,z)-> -2exp(x+y+z)
u = (x,y,z)-> exp(x+y+z)
ū = (x,y,z)-> (exp(x+y+z),exp(x+y+z),exp(x+y+z),0.)

# b = (x,y,z)->0.
# u = (x,y,z)->1+2x+3y
# ū = (x,y,z)->(1+2x+3y,2,3,0.)

ℕ = [5, 10, 20, 40, 80]
# ℕ = [5]
# ℕ = [20]
# ℕ = [40]

# methods = ["Segmental element", "Linear meshfree", "Quadratic meshfree", "Cubic meshfree"]
# methods = ["Segmental element"]
# methods = ["Linear meshfree"]
methods = ["Quadratic meshfree"]

h = zeros(length(ℕ))
L₂Error = zeros(length(ℕ),length(methods))
H₁Error = zeros(length(ℕ),length(methods))
L₂slope = zeros(length(ℕ)-1,length(methods))
H₁slope = zeros(length(ℕ)-1,length(methods))

for n in 1:length(ℕ)
    n₁ = ℕ[n]
    n₂ = ℕ[n]
    nₚ = (n₁+1)*(n₂+1)
    nₑ = 2*n₁*n₂
    h[n] = a₁/n₁
    x = [Node(a₁/n₁*i,a₁/n₂*j,0.) for j in 0:n₂ for i in 0:n₁]
    # x = [Node(a₁/n₁*i + (i≠0 && i≠n₁ && j≠0 && j≠n₂ ? 0.5*a₁/n₁*rand() : 0),a₂/n₂*j + (i≠0 && i≠n₁ && j≠0 && j≠n₂ ? 0.5*a₂/n₂*rand() : 0),0.) for j in 0:n₂ for i in 0:n₁]
    ids_Ω = vcat([[[(n₁+1)*(j-1)+i,(n₁+1)*(j-1)+i+1,(n₁+1)*j+i],
                   [(n₁+1)*j+i+1,(n₁+1)*j+i,(n₁+1)*(j-1)+i+1]] for j in 1:n₂ for i in 1:n₁]...)
    ids_Γ = [[[i,i+1] for i in 1:n₁]...,
             [[(n₁+1)*j,(n₁+1)*(j+1)] for j in 1:n₂]...,
             [[nₚ-k+1,nₚ-k] for k in 1:n₁]...,
             [[(n₁+1)*(n₂-l+1)+1,(n₁+1)*(n₂-l)+1] for l in 1:n₂]...]
    ids_Ω̄ = ids_Ω[vcat([2*i-1 for i in 1:n₁]...,
                       [2*n₁*j for j in 1:n₂]...,
                       [nₑ-2*k+2 for k in 1:n₁]...,
                       [2*n₁*(n₂-l)+1 for l in 1:n₂]...)]
    for m in 1:length(methods)
        if methods[m] == "Segmental element"
            aps = Tri3(x,ids_Ω)
            aps_Γ = Seg2(x,ids_Γ)
            aps_Ω̄ = Tri3(x,ids_Ω̄)
            ape = aps
        elseif methods[m] == "Linear meshfree"
            rg = RegularGrid(x,n=2,γ=1)
            bf = Linear2D(:∂1,:∂x,:∂y,:∂z)
            kf = TensorProductKernel(:∂1,:∂x,:∂y,:∂z,ss=[1.5*a₁/n₁,1.5*a₂/n₂,1.],nm=30)
            aps = TriM(x,ids_Ω,bf=bf,kf=kf,sp=rg,qw=:TriGI7)
            aps_Γ = SegM(x,ids_Γ,bf=bf,kf=kf,sp=rg,qw=:SegGI4)
            aps_Ω̄ = TriM(x,ids_Ω̄,bf=bf,kf=kf)
            ape = TriM(x,ids_Ω,bf=bf,kf=kf,qw=:TriGI16)
        elseif methods[m] == "Quadratic meshfree"
            rg = RegularGrid(x,n=2,γ=1)
            bf = Quadratic2D(:∂1,:∂x,:∂y,:∂z)
            kf = TensorProductKernel(:∂1,:∂x,:∂y,:∂z,ss=[2.5*a₁/n₁,2.5*a₂/n₂,1.],nm=40)
            aps = TriM(x,ids_Ω,bf=bf,kf=kf,sp=rg,qw=:TriGI16)
            aps_Γ = SegM(x,ids_Γ,bf=bf,kf=kf,sp=rg,qw=:SegGI4)
            aps_Ω̄ = TriM(x,ids_Ω̄,bf=bf,kf=kf)
            ape = TriM(x,ids_Ω,bf=bf,kf=kf,qw=:TriGI16)
        end
        # println(length(aps[1].id))
        ops = [Potential_Ω(b), Potential_Γᵍ_Nitsche(u), H₁Error_scale(ū)]
        # ops = [Potential_Ω(b), Potential_Γᵍ_penalty(u,1e7), H₁Error_scale(ū)]
        # ops = [Potential_Ω(b), EBCDOFS(u), H₁Error_scale(ū)]

        k = zeros(nₚ,nₚ)
        f = zeros(nₚ)

        # assembling
        ops[1](aps,k,f)
        # ops[2](aps_Γ,k,f)
        ops[2](aps_Ω̄,aps_Γ,k,f)


        # solve
        d = k\f
        # error
        H₁Error_, L₂Error_ = ops[3](ape,d)
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
