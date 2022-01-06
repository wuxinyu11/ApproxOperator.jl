using BenchmarkTools

x = [(0.1*i) for i in 0:10]
pool = Dict(:x=>x)

node = Node(5,pool)
@btime node = Node(5,pool)

@btime $node[:x]
