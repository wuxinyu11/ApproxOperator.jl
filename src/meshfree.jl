
struct ReproducingKernel{𝑝,𝑠,𝜙,T}<:AbstractElement{T}
    𝓒::Vector{Node}
    𝓖::Vector{Node}
end

for set𝝭 in (:set𝝭!,:set∇𝝭!,:set∇₁𝝭!,:set∇₂𝝭!,:set∇²𝝭!,:set∇³𝝭!,:set∇̂³𝝭!,:set∇²₂𝝭!)
    @eval begin
        function $set𝝭(aps::Vector{T}) where T<:ReproducingKernel
            for ap in aps
                𝓖 = ap.𝓖
                for 𝒙 in 𝓖
                    $set𝝭(ap,𝒙)
                end
            end
        end
    end
end

for set𝝭 in (:set∇̃𝝭!,:set∇̃²𝝭!,:set∇∇̃²𝝭!)
    @eval begin
        function $set𝝭(gps::Vector{T},aps::Vector{S}) where {T<:ReproducingKernel,S<:ReproducingKernel}
            if length(gps) ≠ length(aps)
                error("Miss match element numbers")
            else
                for i in 1:length(gps)
                    $set𝝭(gps[i],aps[i])
                end
            end
        end
    end
end

for set𝝭 in (:set∇̄𝝭!,:set∇̃₁𝝭!)
    @eval begin
        function $set𝝭(aps::Vector{T}) where T<:ReproducingKernel
            for ap in aps
                $set𝝭(ap)
            end
        end
    end
end


function set∇̄²𝝭!(aps::Vector{T};Γᵍ::Vector{T}=T[],Γᶿ::Vector{T}=T[],Γᴾ::Vector{T}=T[]) where T<:ReproducingKernel
    for ap in aps
        set∇̄²𝝭!(ap,Γᵍ=Γᵍ,Γᶿ=Γᶿ,Γᴾ=Γᴾ)
    end
end

function set∇∇̄²𝝭!(aps::Vector{T};Γᵍ::Vector{T}=T[],Γᶿ::Vector{T}=T[],Γᴾ::Vector{T}=T[]) where T<:ReproducingKernel
    for i in 1:length(aps)
        isempty(Γᵍ) ? a = nothing : a = Γᵍ[i]
        isempty(Γᶿ) ? b = nothing : b = Γᶿ[i]
        set∇∇̄²𝝭!(aps[i],Γᵍ=a,Γᶿ=b,Γᴾ=Γᴾ)
    end
end
