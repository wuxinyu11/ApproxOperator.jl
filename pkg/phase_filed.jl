
## Nodes
# Phase field integration nodes
mutable struct PFNode<:ParametricNode
    ξ::AbstractVector{Float64}
    w::Float64
    ℋ::Float64
    ℋₜ::Float64
end
PFNode(node::ParametricNode) = PFNode(node.ξ,node.w,0.0,0.0)

# Phase field friction integration nodes
mutable struct PFFNode<:ParametricNode
    ξ::AbstractVector{Float64}
    w::Float64
    ℋ::Float64
    ℋₙ::Float64
    ℋₘ::Float64
    τ::Float64
    σ::Vector{Float64}
    C::Vector{Float64}
end
PFFNode(node::ParametricNode) = PFFNode(node.ξ,node.w,0.0,0.0,0.0,0.0,zeros(3),zeros(4))

# actions
function set_integration_rule!(ap::Approximator,::Val{:PhaseField})
    for i in 1:length(ap.qw)
        ap.qw[i] = PFNode(ap.qw[i])
    end
end
function set_integration_rule!(aps::Vector{Approximator},::Val{:PhaseField})
    for ap in aps
        set_integration_rule!(ap,Val(:PhaseField))
    end
end
function set_integration_rule!(ap::Approximator,::Val{:PhaseFieldFriction})
    for i in 1:length(ap.qw)
        ap.qw[i] = PFFNode(ap.qw[i])
    end
end
function set_integration_rule!(aps::Vector{Approximator},::Val{:PhaseFieldFriction})
    for ap in aps
        set_integration_rule!(ap,Val(:PhaseFieldFriction))
    end
end
