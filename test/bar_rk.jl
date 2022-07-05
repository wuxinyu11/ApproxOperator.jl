using BenchmarkTools

const ref = (I=1,J=2)
struct T1
    index::NTuple{2,Int}
    data::Dict{Symbol,Tuple{Int,Vector{Float64}}}
end

struct T_
    i::Int
    v::Vector{Float64}
end
function test(t::T1,s::Symbol)
    if haskey(ref,s)
        return t.index[ref[s]]
    elseif haskey(t.data,s)
        i,v = t.data[s]
        j = t.index[i]
        return v[j]
    else
        return 0.0
    end
end

data = Dict([:x=>(1,rand(10)),:y=>(2,rand(10))])
t1 = T1((2,3),data)
@btime test($t1,:x)
@btime test($t1,:I)
@btime test($t1,:z)
@code_warntype test(t1,:x)
