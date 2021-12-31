## DOFs for essential boundary
struct EBCDOFS <: Operator
    g :: Function
end
function (op::EBCDOFS)(ap::Poi1,k::Matrix{Float64},f::Vector{Float64})
    xⱼ = get_coordinates(ap,0.)
    gⱼ = op.g(xⱼ...)
    j = ap.id
    for i in 1:length(f)
        f[i] -= k[i,j]*gⱼ
    end
    k[j,:] .= 0.
    k[:,j] .= 0.
    k[j,j] = 1.
    f[j] = gⱼ
end
## Potential Problem
# k: thermal diffusivity
struct Potential_Ω{F<:Function} <:Operator
    b :: F
    k :: Float64
end
Potential_Ω(b::Function = (x,y,z)->0.) = Potential_Ω(b,1.)
struct Potential_Γᵗ{F<:Function} <:Operator
    t :: F
end
struct Potential_Γᵍ{F<:Function} <:Operator
    g :: F
end
struct Potential_Γᵍ_Lagrange_multiplier{F<:Function} <:Operator
    g :: F
end
struct Potential_Γᵍ_penalty{F<:Function} <:Operator
    g :: F
    α :: Float64
end
Potential_Γᵍ_penalty(g::Function) = Potential_Γᵍ_penalty(g,1e7)
struct Potential_Γᵍ_Nitsche{F<:Function} <:Operator
    g :: F
    k :: Float64
    α :: Float64
end
Potential_Γᵍ_Nitsche(g::Function) = Potential_Γᵍ_Nitsche(g,1.,1e3)
struct Potential_HR_Ω{F<:Function} <:Operator
    b :: F
    k :: Float64
end
struct Potential_HR_Γᵗ{F<:Function} <:Operator
    t :: F
end
struct Potential_HR_Γᵍ{F<:Function} <:Operator
    g :: F
end

function (op::Potential_Ω)(ap::Approximator,k::Matrix{Float64},f::Vector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N,B₁,B₂,B₃ = get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
        Jᵢ = get_jacobe(ap,ξᵢ)
        xᵢ = get_coordinates(ap,ξᵢ)
        W = Jᵢ*wᵢ
        bᵢ = op.b(xᵢ...)
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            for j in 1:get_number_of_indices(ap)
                J = get_global_indice(ap,j)
                k[I,J] += op.k*(B₁[i]*B₁[j] + B₂[i]*B₂[j] + B₃[i]*B₃[j])*W
            end
            f[I] += N[i]*bᵢ*W
        end
    end
end

function (op::Potential_Γᵗ)(ap::Approximator,f::Vector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N = get_shape_functions(ap,ξᵢ,Val(:∂1))
        W = get_jacobe(ap,ξᵢ)*wᵢ
        xᵢ = get_coordinates(ap,ξᵢ)
        tᵢ = op.t(xᵢ...)
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            f[I] = f[I] + N[i]*tᵢ*W
        end
    end
end

function (op::Potential_Γᵍ_penalty)(ap::Approximator,
                                    k ::Matrix{Float64},
                                    f ::Vector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N = get_shape_functions(ap,ξᵢ,Val(:∂1))
        W = get_jacobe(ap,ξᵢ)*wᵢ
        xᵢ = get_coordinates(ap,ξᵢ)
        gᵢ = op.g(xᵢ...)
        α = op.α
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            f[I] += α*N[i]*gᵢ*W
            for j in 1:get_number_of_indices(ap)
                J = get_global_indice(ap,j)
                k[I,J] += α*N[i]*N[j]*W
            end
        end
    end
end

function (op::Potential_Γᵍ_Lagrange_multiplier)(ap1::Approximator,
                                                ap2::Approximator,
                                                g  ::Matrix{Float64},
                                                q  ::Vector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N = get_shape_functions(ap1,ξᵢ,Val(:∂1))
        N̄ = get_shape_functions(ap2,ξᵢ,Val(:∂1))
        W = get_jacobe(ap1,ξᵢ)*wᵢ
        xᵢ = get_coordinates(ap1,ξᵢ)
        gᵢ = op.g(xᵢ...)
        for k in 1:length(ap2.id)
            K = get_global_indice(ap2,k)
            q[K] -= N̄[k]*gᵢ*W
            for i in 1:get_number_of_indices(ap1)
                I = get_global_indice(ap1,i)
                g[I,K] -= N[i]*N̄[k]*W
            end
        end
    end
end

function (op::Potential_Γᵍ_Nitsche)(ap1::Approximator,
                                    ap2::Approximator,
                                    k  ::Matrix{Float64},
                                    f  ::Vector{Float64})
    n₁,n₂,n₃ = get_normal(ap1,ap2)
    for qw in ap.qw
        ηᵢ = qw.ξ
        wᵢ = qw.w
        ξᵢ = get_coordinates(ap1,ap2,ηᵢ)
        N,B₁,B₂,B₃ = get_shape_functions(ap1,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
        W = get_jacobe(ap2,ηᵢ)*wᵢ
        xᵢ = get_coordinates(ap2,ηᵢ)
        gᵢ = op.g(xᵢ...)
        α  = op.α
        for i in 1:get_number_of_indices(ap1)
            I = get_global_indice(ap1,i)
            tᵢ = n₁*B₁[i] + n₂*B₂[i] + n₃*B₃[i]
            f[I] += op.k*(- tᵢ + α*N[i])*gᵢ*W
            for j in 1:get_number_of_indices(ap1)
                J = get_global_indice(ap1,j)
                tⱼ = n₁*B₁[j] + n₂*B₂[j] + n₃*B₃[j]
                k[I,J] += op.k*(- N[i]*tⱼ - tᵢ*N[j] + α*N[i]*N[j])*W
            end
        end
    end
end

## Plane stress
struct Elasticity_Ω{F<:Function} <: Operator
    b::F
    E::Float64
    ν::Float64
end
struct PlaneStress_Ω{F<:Function} <: Operator
    b::F
    E::Float64
    ν::Float64
end
struct NonlinearPlaneStress_C_Ω{F<:Function} <: Operator
    b::F
end
struct PlaneStrain_Ωᵛ <: Operator
    E::Float64
    ν::Float64
end
struct PlaneStrain_Ωᵈ{F<:Function} <: Operator
    b::F
    E::Float64
    ν::Float64
end
struct Elasticity_Γᵗ{F<:Function} <: Operator
    t::F
end
struct PlaneStress_Γᵗ{F<:Function} <: Operator
    t::F
end
struct PlaneStress_Γᵍ_Lagrange_multiplier{F<:Function,N<:Function} <: Operator
    g::F
    n::N
end
struct Elasticity_Γᵍ_penalty{F<:Function,N<:Function} <: Operator
    g::F
    n::N
    α::Float64
end
Elasticity_Γᵍ_penalty(g::Function) = Elasticity_Γᵍ_penalty(g,(x,y,z)->(1.0,0.0,1.0,0.0,0.0,1.0),1e7)

struct NonlinearElasticity_Γᵍ_penalty{F<:Function,N<:Function} <: Operator
    g::F
    n::N
    α::Float64
end
NonlinearElasticity_Γᵍ_penalty(g::Function) = NonlinearElasticity_Γᵍ_penalty(g,(x,y,z)->(1.0,0.0,1.0,0.0,0.0,1.0),1e7)

struct PlaneStress_Γᵍ_penalty{F<:Function,N<:Function} <: Operator
    g::F
    n::N
    α::Float64
end
PlaneStress_Γᵍ_penalty(g::Function) = PlaneStress_Γᵍ_penalty(g,(x,y,z)->(1.0,0.0,1.0),1e7)

struct NonlinearPlaneStress_Γᵍ_penalty{F<:Function,N<:Function} <: Operator
    g::F
    n::N
    α::Float64
end
NonlinearPlaneStress_Γᵍ_penalty(g::Function) = NonlinearPlaneStress_Γᵍ_penalty(g,(x,y,z)->(1.0,0.0,1.0),1e7)

struct PlaneStress_Γᵍ_Nitsche{F<:Function} <: Operator
    g::F
end

function (op::Elasticity_Ω)(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N,B₁,B₂ = get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y))
        Jᵢ = get_jacobe(ap,ξᵢ)
        xᵢ = get_coordinates(ap,ξᵢ)
        W = Jᵢ*wᵢ
        b₁,b₂,b₃ = op.b(xᵢ...)
        E = op.E
        ν = op.ν
        Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-2ν)
        Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-2ν)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            for j in 1:get_number_of_indices(ap)
                J = get_global_indice(ap,j)
                k[3*I-2,3*J-2] += (Cᵢᵢᵢᵢ*B₁[i]*B₁[j] + Cᵢⱼᵢⱼ*(B₂[i]*B₂[j] + B₃[i]*B₃[j]))*W
                k[3*I-2,3*J-1] += (Cᵢᵢⱼⱼ*B₁[i]*B₂[j] + Cᵢⱼᵢⱼ*B₂[i]*B₁[j])*W
                k[3*I-2,3*J]   += (Cᵢᵢⱼⱼ*B₁[i]*B₃[j] + Cᵢⱼᵢⱼ*B₃[i]*B₁[j])*W
                k[3*I-1,3*J-2] += (Cᵢᵢⱼⱼ*B₂[i]*B₁[j] + Cᵢⱼᵢⱼ*B₁[i]*B₂[j])*W
                k[3*I-1,3*J-1] += (Cᵢᵢᵢᵢ*B₂[i]*B₂[j] + Cᵢⱼᵢⱼ*(B₁[i]*B₁[j] + B₃[i]*B₃[j]))*W
                k[3*I-1,3*J]   += (Cᵢᵢⱼⱼ*B₂[i]*B₃[j] + Cᵢⱼᵢⱼ*B₃[i]*B₂[j])*W
                k[3*I,3*J-2]   += (Cᵢᵢⱼⱼ*B₃[i]*B₁[j] + Cᵢⱼᵢⱼ*B₁[i]*B₃[j])*W
                k[3*I,3*J-1]   += (Cᵢᵢⱼⱼ*B₃[i]*B₂[j] + Cᵢⱼᵢⱼ*B₂[i]*B₃[j])*W
                k[3*I,3*J]     += (Cᵢᵢᵢᵢ*B₃[i]*B₃[j] + Cᵢⱼᵢⱼ*(B₁[i]*B₁[j] + B₂[i]*B₂[j]))*W
            end
            f[3*I-2] += N[i]*b₁*W
            f[3*I-1] += N[i]*b₂*W
            f[3*I]   += N[i]*b₃*W
        end
    end
end

function (op::PlaneStrain_Ωᵛ)(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N,B₁,B₂ = get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y))
        Jᵢ = get_jacobe(ap,ξᵢ)
        xᵢ = get_coordinates(ap,ξᵢ)
        W = Jᵢ*wᵢ
        E = op.E
        ν = op.ν
        Cᵛ = E/(1-2*ν)
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            for j in 1:get_number_of_indices(ap)
                J = get_global_indice(ap,j)
                k[2*I-1,2*J-1] += Cᵛ/3*B₁[i]*B₁[j]*W
                k[2*I-1,2*J]   += Cᵛ/3*B₁[i]*B₂[j]*W
                k[2*I,2*J-1]   += Cᵛ/3*B₂[i]*B₁[j]*W
                k[2*I,2*J]     += Cᵛ/3*B₂[i]*B₂[j]*W
            end
        end
    end
end

function (op::PlaneStrain_Ωᵈ)(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N,B₁,B₂ = get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y))
        Jᵢ = get_jacobe(ap,ξᵢ)
        xᵢ = get_coordinates(ap,ξᵢ)
        W = Jᵢ*wᵢ
        b₁,b₂ = op.b(xᵢ...)
        E = op.E
        ν = op.ν
        Cᵈ = E/(1+ν)
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            for j in 1:get_number_of_indices(ap)
                J = get_global_indice(ap,j)
                k[2*I-1,2*J-1] += Cᵈ*( 2/3*B₁[i]*B₁[j]+1/2*B₂[i]*B₂[j])*W
                k[2*I-1,2*J]   += Cᵈ*(-1/3*B₁[i]*B₂[j]+1/2*B₂[i]*B₁[j])*W
                k[2*I,2*J-1]   += Cᵈ*(-1/3*B₂[i]*B₁[j]+1/2*B₁[i]*B₂[j])*W
                k[2*I,2*J]     += Cᵈ*( 2/3*B₂[i]*B₂[j]+1/2*B₁[i]*B₁[j])*W
            end
            f[2*I-1] += N[i]*b₁*W
            f[2*I] += N[i]*b₂*W
        end
    end
end

function (op::PlaneStress_Ω)(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N,B₁,B₂ = get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y))
        Jᵢ = get_jacobe(ap,ξᵢ)
        xᵢ = get_coordinates(ap,ξᵢ)
        W = Jᵢ*wᵢ
        b₁,b₂ = op.b(xᵢ...)
        E = op.E
        ν = op.ν
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            for j in 1:get_number_of_indices(ap)
                J = get_global_indice(ap,j)
                k[2*I-1,2*J-1] += (Cᵢᵢᵢᵢ*B₁[i]*B₁[j] + Cᵢⱼᵢⱼ*B₂[i]*B₂[j])*W
                k[2*I-1,2*J]   += (Cᵢᵢⱼⱼ*B₁[i]*B₂[j] + Cᵢⱼᵢⱼ*B₂[i]*B₁[j])*W
                k[2*I,2*J-1]   += (Cᵢᵢⱼⱼ*B₂[i]*B₁[j] + Cᵢⱼᵢⱼ*B₁[i]*B₂[j])*W
                k[2*I,2*J]     += (Cᵢᵢᵢᵢ*B₂[i]*B₂[j] + Cᵢⱼᵢⱼ*B₁[i]*B₁[j])*W
            end
            f[2*I-1] += N[i]*b₁*W
            f[2*I] += N[i]*b₂*W
        end
    end
end

function (op::NonlinearPlaneStress_C_Ω)(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N,B₁,B₂ = get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y))
        Jᵢ = get_jacobe(ap,ξᵢ)
        xᵢ = get_coordinates(ap,ξᵢ)
        W = Jᵢ*wᵢ
        b₁,b₂ = op.b(xᵢ...)
        C₁₁₁₁ = qw.C[1]
        C₁₁₂₂ = qw.C[2]
        C₂₂₂₂ = qw.C[3]
        C₁₂₁₂ = qw.C[4]
        σ₁₁ = qw.σ[1]
        σ₂₂ = qw.σ[2]
        σ₁₂ = qw.σ[3]
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            for j in 1:get_number_of_indices(ap)
                J = get_global_indice(ap,j)
                k[2*I-1,2*J-1] += (C₁₁₁₁*B₁[i]*B₁[j] + C₁₂₁₂*B₂[i]*B₂[j])*W
                k[2*I-1,2*J]   += (C₁₁₂₂*B₁[i]*B₂[j] + C₁₂₁₂*B₂[i]*B₁[j])*W
                k[2*I,2*J-1]   += (C₁₁₂₂*B₂[i]*B₁[j] + C₁₂₁₂*B₁[i]*B₂[j])*W
                k[2*I,2*J]     += (C₂₂₂₂*B₂[i]*B₂[j] + C₁₂₁₂*B₁[i]*B₁[j])*W
            end
            f[2*I-1] += (N[i]*b₁ - B₁[i]*σ₁₁ - B₂[i]*σ₁₂)*W
            f[2*I] += (N[i]*b₂ - B₁[i]*σ₁₂ - B₂[i]*σ₂₂)*W
        end
    end
end
function (op::PlaneStress_Γᵗ)(ap::Approximator,f::Vector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N = get_shape_functions(ap,ξᵢ,Val(:∂1))
        W = get_jacobe(ap,ξᵢ)*wᵢ
        xᵢ = get_coordinates(ap,ξᵢ)
        t₁,t₂ = op.t(xᵢ...)
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            f[2*I-1] += N[i]*t₁*W
            f[2*I] += N[i]*t₂*W
        end
    end
end

function (op::Elasticity_Γᵗ)(ap::Approximator,f::Vector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N = get_shape_functions(ap,ξᵢ,Val(:∂1))
        W = get_jacobe(ap,ξᵢ)*wᵢ
        xᵢ = get_coordinates(ap,ξᵢ)
        t₁,t₂,t₃ = op.t(xᵢ...)
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            f[3*I-2] += N[i]*t₁*W
            f[3*I-1] += N[i]*t₂*W
            f[3*I]   += N[i]*t₃*W
        end
    end
end

function (op::Elasticity_Γᵍ_penalty)(ap::Approximator,
                                     k ::Matrix{Float64},
                                     f ::Vector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N = get_shape_functions(ap,ξᵢ,Val(:∂1))
        W = get_jacobe(ap,ξᵢ)*wᵢ
        xᵢ = get_coordinates(ap,ξᵢ)
        g₁,g₂,g₃ = op.g(xᵢ...)
        n₁₁,n₁₂,n₂₂,n₁₃,n₂₃,n₃₃ = op.n(xᵢ...)
        α = op.α
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            f[3*I-2] += α*N[i]*g₁*W
            f[3*I-1] += α*N[i]*g₂*W
            f[3*I]   += α*N[i]*g₃*W
            for j in 1:get_number_of_indices(ap)
                J = get_global_indice(ap,j)
                k[3*I-2,3*J-2] += α*n₁₁*N[i]*N[j]*W
                k[3*I-2,3*J-1] += α*n₁₂*N[i]*N[j]*W
                k[3*I-2,3*J]   += α*n₁₃*N[i]*N[j]*W
                k[3*I-1,3*J-2] += α*n₁₂*N[i]*N[j]*W
                k[3*I-1,3*J-1] += α*n₂₂*N[i]*N[j]*W
                k[3*I-1,3*J]   += α*n₂₃*N[i]*N[j]*W
                k[3*I,3*J-2] += α*n₁₃*N[i]*N[j]*W
                k[3*I,3*J-1] += α*n₂₃*N[i]*N[j]*W
                k[3*I,3*J]   += α*n₃₃*N[i]*N[j]*W
            end
        end
    end
end

function (op::PlaneStress_Γᵍ_penalty)(ap::Approximator,
                                     k ::AbstractMatrix{Float64},
                                     f ::AbstractVector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N = get_shape_functions(ap,ξᵢ,Val(:∂1))
        W = get_jacobe(ap,ξᵢ)*wᵢ
        xᵢ = get_coordinates(ap,ξᵢ)
        g₁,g₂ = op.g(xᵢ...)
        n₁₁,n₁₂,n₂₂ = op.n(xᵢ...)
        α = op.α
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            f[2*I-1] += α*N[i]*g₁*W
            f[2*I]   += α*N[i]*g₂*W
            for j in 1:get_number_of_indices(ap)
                J = get_global_indice(ap,j)
                k[2*I-1,2*J-1] += α*n₁₁*N[i]*N[j]*W
                k[2*I-1,2*J]   += α*n₁₂*N[i]*N[j]*W
                k[2*I,2*J-1]   += α*n₁₂*N[i]*N[j]*W
                k[2*I,2*J]     += α*n₂₂*N[i]*N[j]*W
            end
        end
    end
end

function (op::NonlinearPlaneStress_Γᵍ_penalty)(ap::Approximator,
                                               k ::AbstractMatrix{Float64},
                                               f ::AbstractVector{Float64},
                                               d ::AbstractVector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N = get_shape_functions(ap,ξᵢ,Val(:∂1))
        W = get_jacobe(ap,ξᵢ)*wᵢ
        xᵢ = get_coordinates(ap,ξᵢ)
        ḡ₁,ḡ₂ = op.g(xᵢ...)
        n₁₁,n₁₂,n₂₂ = op.n(xᵢ...)
        α = op.α
        g₁ = 0.0
        g₂ = 0.0
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            g₁ += N[i]*d[2*I-1]
            g₂ += N[i]*d[2*I]
        end
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            f[2*I-1] += α*N[i]*(ḡ₁-g₁)*W
            f[2*I]   += α*N[i]*(ḡ₂-g₂)*W
            for j in 1:get_number_of_indices(ap)
                J = get_global_indice(ap,j)
                k[2*I-1,2*J-1] += α*n₁₁*N[i]*N[j]*W
                k[2*I-1,2*J]   += α*n₁₂*N[i]*N[j]*W
                k[2*I,2*J-1]   += α*n₁₂*N[i]*N[j]*W
                k[2*I,2*J]     += α*n₂₂*N[i]*N[j]*W
            end
        end
    end
end

## Phase field model
struct PlaneStress_PhaseField_Ω{F<:Function} <: Operator
    b::F
    E::Float64
    ν::Float64
    η::Float64
end
PlaneStress_PhaseField_Ω(b::Function,E::Float64,ν::Float64) = PlaneStress_PhaseField_Ω(b,E,ν,1e-6)

struct Update_Friction_PhaseField_PlaneStress <: Operator
    E::Float64
    ν::Float64
    c::Float64
    𝜙::Float64
    𝜙ᵣ::Float64
    η::Float64
end
PlaneStress_PhaseField_Ω(b::Function,E::Float64,ν::Float64) = PlaneStress_PhaseField_Ω(b,E,ν,1e-6)

struct SecondOrderPhaseField <: Operator
    k::Float64
    l::Float64
end

struct Update_HistoryField_PlaneStress <: Operator
    E::Float64
    ν::Float64
end

function (op::PlaneStress_PhaseField_Ω)(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64},d::AbstractVector{Float64},dᵥ::AbstractVector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N,B₁,B₂ = get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y))
        Jᵢ = get_jacobe(ap,ξᵢ)
        xᵢ = get_coordinates(ap,ξᵢ)
        W = Jᵢ*wᵢ
        b₁,b₂ = op.b(xᵢ...)
        E = op.E
        ν = op.ν
        v = 0.0
        ε₁₁ = 0.0
        ε₂₂ = 0.0
        ε₁₂ = 0.0
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            v += N[i]*dᵥ[I]
            ε₁₁ += B₁[i]*d[2*I-1]
            ε₂₂ += B₂[i]*d[2*I]
            ε₁₂ += B₂[i]*d[2*I-1] + B₁[i]*d[2*I]
        end
        Cᵢᵢᵢᵢ = (v^2+op.η)*E/(1-ν^2)
        Cᵢᵢⱼⱼ = (v^2+op.η)*E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = (v^2+op.η)*E/2/(1+ν)
        σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        σ₂₂ = Cᵢᵢᵢᵢ*ε₂₂ + Cᵢᵢⱼⱼ*ε₁₁
        σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            for j in 1:get_number_of_indices(ap)
                J = get_global_indice(ap,j)
                k[2*I-1,2*J-1] += (Cᵢᵢᵢᵢ*B₁[i]*B₁[j] + Cᵢⱼᵢⱼ*B₂[i]*B₂[j])*W
                k[2*I-1,2*J]   += (Cᵢᵢⱼⱼ*B₁[i]*B₂[j] + Cᵢⱼᵢⱼ*B₂[i]*B₁[j])*W
                k[2*I,2*J-1]   += (Cᵢᵢⱼⱼ*B₂[i]*B₁[j] + Cᵢⱼᵢⱼ*B₁[i]*B₂[j])*W
                k[2*I,2*J]     += (Cᵢᵢᵢᵢ*B₂[i]*B₂[j] + Cᵢⱼᵢⱼ*B₁[i]*B₁[j])*W
            end
            f[2*I-1] += (N[i]*b₁ - B₁[i]*σ₁₁ - B₂[i]*σ₁₂)*W
            f[2*I] += (N[i]*b₂ - B₁[i]*σ₁₂ - B₂[i]*σ₂₂)*W
        end
    end
end

function (op::Update_Friction_PhaseField_PlaneStress)(ap::Approximator,d::AbstractVector{Float64},dᵥ::AbstractVector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N,B₁,B₂ = get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y))
        Jᵢ = get_jacobe(ap,ξᵢ)
        xᵢ = get_coordinates(ap,ξᵢ)
        W = Jᵢ*wᵢ
        E = op.E
        ν = op.ν
        𝜙 = op.𝜙
        𝜙ᵣ = op.𝜙ᵣ
        c = op.c
        θ = π/4 - 𝜙ᵣ/2
        sinθ = sin(θ)
        cosθ = cos(θ)
        v = 0.0
        ε₁₁ = 0.0
        ε₂₂ = 0.0
        ε₁₂ = 0.0
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            v += N[i]*dᵥ[I]
            ε₁₁ += B₁[i]*d[2*I-1]
            ε₂₂ += B₂[i]*d[2*I]
            ε₁₂ += B₂[i]*d[2*I-1] + B₁[i]*d[2*I]
        end
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        σ₂₂ = Cᵢᵢᵢᵢ*ε₂₂ + Cᵢᵢⱼⱼ*ε₁₁
        σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
        σ₃₃ = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        𝝈 = SMatrix{3,3,Float64}(σ₁₁,σ₁₂,0.0,σ₁₂,σ₂₂,0.0,0.0,0.0,σ₃₃)
        ~,𝒂 = eigen(𝝈)
        𝒂₁ = 𝒂[:,1]
        𝒂₂ = 𝒂[:,2]
        𝒏 = 𝒂₁*sinθ+cross(𝒂₂,𝒂₁)*cosθ
        𝒎 = cross(𝒂₂,𝒏)
        α₁₁ = 𝒏[1]*𝒎[1] + 𝒎[1]*𝒏[1]
        α₁₂ = 𝒏[1]*𝒎[2] + 𝒎[1]*𝒏[2]
        α₂₂ = 𝒏[2]*𝒎[2] + 𝒎[2]*𝒏[2]
        α₃₃ = 𝒏[3]*𝒎[3] + 𝒎[3]*𝒏[3]
        β₁₁ = 𝒏[1]*𝒏[1] + 𝒏[1]*𝒏[1]
        β₁₂ = 𝒏[1]*𝒏[2] + 𝒏[1]*𝒏[2]
        β₂₂ = 𝒏[2]*𝒏[2] + 𝒏[2]*𝒏[2]
        β₃₃ = 𝒏[3]*𝒏[3] + 𝒏[3]*𝒏[3]
        τbulk  = 0.5*(σ₁₁*α₁₁ + 2*σ₁₂*α₁₂ + σ₂₂*α₂₂ + σ₃₃*α₃₃)
        pN = - σ₁₁*β₁₁ + 2*σ₁₂*β₁₂ + σ₂₂*β₂₂ + σ₃₃*β₃₃
        abs(v) < eps()*1e5 ? τY = c + pN*tan(𝜙) : τY = pN*tan(𝜙ᵣ)
        f = abs(τbulk) - τY
        if f < 0
            qw.σ[1] = σ₁₁
            qw.σ[2] = σ₂₂
            qw.σ[3] = σ₁₂
            qw.C[1] = Cᵢᵢᵢᵢ
            qw.C[2] = Cᵢᵢⱼⱼ
            qw.C[3] = Cᵢᵢᵢᵢ
            qw.C[4] = Cᵢⱼᵢⱼ
            if abs(v) < eps()*1e5
                qw.ℋ = (c+pN*tan(𝜙) - pN*tan(𝜙ᵣ))^2/2/Cᵢⱼᵢⱼ
            else
                qw.ℋ = qw.ℋₘ
            end
        else
            qw.σ[1] = σ₁₁ - v^2*(τbulk - qw.τ)*α₁₁
            qw.σ[2] = σ₂₂ - v^2*(τbulk - qw.τ)*α₂₂
            qw.σ[3] = σ₁₂ - v^2*(τbulk - qw.τ)*α₁₂
            qw.C[1] = Cᵢᵢᵢᵢ - v^2*Cᵢⱼᵢⱼ*α₁₁*α₁₁
            qw.C[2] = Cᵢᵢⱼⱼ - v^2*Cᵢⱼᵢⱼ*α₁₁*α₂₂
            qw.C[3] = Cᵢᵢᵢᵢ - v^2*Cᵢⱼᵢⱼ*α₂₂*α₂₂
            qw.C[4] = Cᵢⱼᵢⱼ - v^2*Cᵢⱼᵢⱼ*α₁₂*α₁₂
            Δγ = (ε₁₁*α₁₁ + ε₂₂*α₂₂ + ε₁₂*α₁₂) - qw.τ
            qw.ℋ = qw.ℋₙ + (τbulk - pN*tan(𝜙))*Δγ
        end
    end
end

function (op::SecondOrderPhaseField)(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64},d::AbstractVector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        ℋᵢ = qw.ℋ
        N,B₁,B₂,B₃ = get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
        Jᵢ = get_jacobe(ap,ξᵢ)
        xᵢ = get_coordinates(ap,ξᵢ)
        W = Jᵢ*wᵢ
        v = 0.0
        ∂v∂x = 0.0
        ∂v∂y = 0.0
        ∂v∂z = 0.0
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            v += N[i]*d[I]
            ∂v∂x += B₁[i]*d[I]
            ∂v∂y += B₂[i]*d[I]
            ∂v∂z += B₃[i]*d[I]
        end
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            for j in 1:get_number_of_indices(ap)
                J = get_global_indice(ap,j)
                k[I,J] += ((op.k/2/op.l + 2*ℋᵢ)*(N[i]*N[j]) + op.k*2*op.l*(B₁[i]*B₁[j]+B₂[i]*B₂[j]+B₃[i]*B₃[j]))*W
            end
            f[I] -= (N[i]*(op.k/2/op.l*(v-1.)+2*ℋᵢ*v) + op.k*2*op.l*(B₁[i]*∂v∂x+B₂[i]*∂v∂y+B₃[i]*∂v∂z))*W
        end
    end
end

function (op::Update_HistoryField_PlaneStress)(ap::Approximator,d::AbstractVector{Float64})
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N,B₁,B₂ = get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y))
        Jᵢ = get_jacobe(ap,ξᵢ)
        xᵢ = get_coordinates(ap,ξᵢ)
        W = Jᵢ*wᵢ
        E = op.E
        ν = op.ν
        ε₁₁ = 0.0
        ε₂₂ = 0.0
        ε₁₂ = 0.0
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            ε₁₁ += B₁[i]*d[2*I-1]
            ε₂₂ += B₂[i]*d[2*I]
            ε₁₂ += B₂[i]*d[2*I-1] + B₁[i]*d[2*I]
        end
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        σ₂₂ = Cᵢᵢᵢᵢ*ε₂₂ + Cᵢᵢⱼⱼ*ε₁₁
        σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
        ℋₜ = 0.5*(σ₁₁*ε₁₁ + σ₂₂*ε₂₂ + σ₁₂*ε₁₂)
        qw.ℋ = max(qw.ℋₜ,ℋₜ)
    end
end

## Error measure
struct L₂Error_scale
    ū :: Function
end

struct H₁Error_scale
    ū :: Function
end

struct L₂Error_tensor
    ū :: Function
end

struct L₂Error_2nd_order_tensor
    ū :: Function
end

struct H₁Error_tensor
    ū :: Function
end

struct HₑError
    ū :: Function
    E :: Float64
    ν :: Float64
end

struct HₑError_PlaneStress
    ū :: Function
    E :: Float64
    ν :: Float64
end

function (op::L₂Error_scale)(ap::Approximator,d::Vector{Float64})
    Δu²= 0
    ū² = 0
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N = get_shape_functions(ap,ξᵢ,Val(:∂1))
        W = get_jacobe(ap,ξᵢ)*wᵢ
        xᵢ = get_coordinates(ap,ξᵢ)
        ūᵢ = op.ū(xᵢ...)
        uᵢ = 0
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            uᵢ += N[i]*d[I]
        end
        Δu² += (uᵢ - ūᵢ)^2*W
        ū²  += ūᵢ^2*W
    end
    return Δu², ū²
end

function (op::L₂Error_scale)(aps::Vector{T},d::Vector{Float64}) where T<:Approximator
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δu², ū² = op(ap,d)
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function (op::H₁Error_scale)(ap::Approximator,d::Vector{Float64})
    Δ∇u²= 0
    ∇ū² = 0
    Δu²= 0
    ū² = 0
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N,B₁,B₂,B₃ = get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
        W = get_jacobe(ap,ξᵢ)*wᵢ
        xᵢ = get_coordinates(ap,ξᵢ)
        ūᵢ,∂ūᵢ∂x,∂ūᵢ∂y,∂ūᵢ∂z = op.ū(xᵢ...)
        uᵢ = 0.
        ∂uᵢ∂x = 0.
        ∂uᵢ∂y = 0.
        ∂uᵢ∂z = 0.
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            uᵢ += N[i]*d[I]
            ∂uᵢ∂x += B₁[i]*d[I]
            ∂uᵢ∂y += B₂[i]*d[I]
            ∂uᵢ∂z += B₃[i]*d[I]
        end
        Δ∇u² += ((∂uᵢ∂x - ∂ūᵢ∂x)^2 + (∂uᵢ∂y - ∂ūᵢ∂y)^2 + (∂uᵢ∂z - ∂ūᵢ∂z)^2)*W
        ∇ū² += (∂ūᵢ∂x^2 + ∂ūᵢ∂y^2 + ∂ūᵢ∂z^2)*W
        Δu² += (uᵢ - ūᵢ)^2*W
        ū² += ūᵢ^2*W
    end
    return Δ∇u², ∇ū², Δu², ū²
end

function (op::H₁Error_scale)(aps::Vector{T},d::Vector{Float64}) where T<:Approximator
    H₁Norm_Δu²= 0
    H₁Norm_ū² = 0
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δ∇u², ∇ū², Δu², ū² = op(ap,d)
        H₁Norm_Δu² += Δu² + Δ∇u²
        H₁Norm_ū²  += ū² + ∇ū²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (H₁Norm_Δu²/H₁Norm_ū²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function (op::L₂Error_tensor)(ap::Approximator,d::AbstractVector{Float64})
    Δu²= 0
    ū² = 0
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N = get_shape_functions(ap,ξᵢ,Val(:∂1))
        W = get_jacobe(ap,ξᵢ)*wᵢ
        xᵢ = get_coordinates(ap,ξᵢ)
        ū₁,ū₂,ū₃ = op.ū(xᵢ...)
        u₁ = 0
        u₂ = 0
        u₃ = 0
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            u₁ += N[i]*d[3*I-2]
            u₂ += N[i]*d[3*I-1]
            u₃ += N[i]*d[3*I]
        end
        Δu² += ((u₁ - ū₁)^2+(u₂ - ū₂)^2+(u₃ - ū₃)^2)*W
        ū²  += (ū₁^2+ū₂^2+ū₃^2)*W
    end
    return Δu², ū²
end

function (op::L₂Error_tensor)(aps::Vector{T},d::AbstractVector{Float64}) where T<:Approximator
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δu², ū² = op(ap,d)
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function (op::L₂Error_2nd_order_tensor)(ap::Approximator,d::AbstractVector{Float64})
    Δu²= 0
    ū² = 0
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N = get_shape_functions(ap,ξᵢ,Val(:∂1))
        W = get_jacobe(ap,ξᵢ)*wᵢ
        xᵢ = get_coordinates(ap,ξᵢ)
        ū₁,ū₂ = op.ū(xᵢ...)
        u₁ = 0
        u₂ = 0
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            u₁ += N[i]*d[2*I-1]
            u₂ += N[i]*d[2*I]
        end
        Δu² += ((u₁ - ū₁)^2+(u₂ - ū₂)^2)*W
        ū²  += (ū₁^2+ū₂^2)*W
    end
    return Δu², ū²
end

function (op::L₂Error_2nd_order_tensor)(aps::Vector{T},d::AbstractVector{Float64}) where T<:Approximator
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δu², ū² = op(ap,d)
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function (op::HₑError_PlaneStress)(ap::Approximator,d::Vector{Float64})
    ΔW²= 0
    W̄² = 0
    Δu²= 0
    ū² = 0
    E = op.E
    ν = op.ν
    Cᵢᵢᵢᵢ = E/(1-ν^2)
    Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
    Cᵢⱼᵢⱼ = E/2/(1+ν)
    for qw in ap.qw
        ξᵢ = qw.ξ
        wᵢ = qw.w
        N,B₁,B₂ = get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y))
        W = get_jacobe(ap,ξᵢ)*wᵢ
        xᵢ = get_coordinates(ap,ξᵢ)
        ū₁,ū₂,∂ū₁∂x,∂ū₁∂y,∂ū₂∂x,∂ū₂∂y = op.ū(xᵢ...)
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
        σ̄₁₁ = Cᵢᵢᵢᵢ*ε̄₁₁ + Cᵢᵢⱼⱼ*ε̄₂₂
        σ̄₂₂ = Cᵢᵢᵢᵢ*ε̄₂₂ + Cᵢᵢⱼⱼ*ε̄₁₁
        σ̄₁₂ = Cᵢⱼᵢⱼ*ε̄₁₂
        u₁ = 0.
        u₂ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            u₁ += N[i]*d[2*I-1]
            u₂ += N[i]*d[2*I]
            ε₁₁ += B₁[i]*d[2*I-1]
            ε₂₂ += B₂[i]*d[2*I]
            ε₁₂ += B₂[i]*d[2*I-1] + B₁[i]*d[2*I]
        end
        σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        σ₂₂ = Cᵢᵢᵢᵢ*ε₂₂ + Cᵢᵢⱼⱼ*ε₁₁
        σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
        ΔW² += 0.5*((σ₁₁-σ̄₁₁)*(ε₁₁-ε̄₁₁) + (σ₂₂-σ̄₂₂)*(ε₂₂-ε̄₂₂) + (σ₁₂-σ̄₁₂)*(ε₁₂-ε̄₁₂))*W
        W̄² += 0.5*(σ₁₁*ε₁₁ + σ₂₂*ε₂₂ + σ₁₂*ε₁₂)*W
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*W
        ū² += (ū₁^2 + ū₂^2)*W
    end
    return ΔW², W̄², Δu², ū²
end

function (op::HₑError_PlaneStress)(aps::Vector{T},d::Vector{Float64}) where T<:Approximator
    HₑNorm_ΔW²= 0
    HₑNorm_W̄² = 0
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        ΔW², W̄², Δu², ū² = op(ap,d)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄²  += W̄²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

## General funcitons
for t in subtypes(Operator)
    function (op::t)(aps::Vector{T},f::AbstractVector{Float64}) where T<:Approximator
        for ap in aps
            op(ap,f)
        end
    end
    function (op::t)(aps::Vector{T},k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:Approximator
        for ap in aps
            op(ap,k,f)
        end
    end
    function (op::t)(aps::Vector{T},k::AbstractMatrix{Float64},f::AbstractVector{Float64},d::AbstractVector{Float64}) where T<:Approximator
        for ap in aps
            op(ap,k,f,d)
        end
    end
    function (op::t)(aps::Vector{T},k::AbstractMatrix{Float64},f::AbstractVector{Float64},d::AbstractVector{Float64},p::AbstractVector{Float64}) where T<:Approximator
        for ap in aps
            op(ap,k,f,d,p)
        end
    end
    function (op::t)(aps::Vector{T},d::AbstractVector{Float64},p::AbstractVector{Float64}) where T<:Approximator
        for ap in aps
            op(ap,d,p)
        end
    end
    function (op::t)(aps1::Vector{T} where T<:Approximator,aps2::Vector{U} where U<:Approximator,k::Matrix{Float64},f::Vector{Float64})
        for i in 1:length(aps1)
            op(aps1[i],aps2[i],k,f)
        end
    end
end
