
using BenchmarkTools
using TimerOutputs
import Base: getproperty, setindex!

to = TimerOutput()

n = 10000
v = zeros(n)
r = rand(n)

struct Pool{T}
    i::Int
    v::Vector{T}
end
getproperty(p::Pool,f::Symbol) = getfield(p,f)[getfield(p,:i)]

# @timeit to "vector setindex" begin
#     for i in 1:n
#         v[i] = rand()
#     end
# end

p = [Pool(i,v) for i in 1:n]

# for i in 1:n
#     @timeit to "vector setindex" v[i] = rand()
# end

for i in 1:n
    @timeit to "vector getindex" a = v[i]
    @timeit to "pool getindex" a = p[i].v
end

show(to)
