
## ----------------- Node ------------------
struct Node
    id::Int
    dp::Dict{Symbol,Vector{Float64}}
end

@inline function getindex(n::Node,s::Symbol)
    return n.dp[s][n.id]
end

## ---------------- DataPool ----------------
# struct DataPool
#     data::Vector{Vector{Float64}}
#     index::Dict{Symbol,Int}
# end
#
# struct Node
#     id::Int
#     dp::DataPool
# end
#
# @inline function getindex(n::Node,s::Symbol)
#     return n.dp.data[n.dp.index[s]][n.id]
# end
