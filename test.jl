using Revise
using ApproxOperator
using BenchmarkTools

# efficiency_meshfree()

# x = [Node(1.0/10*i,0.0,0.0) for i in 0:10]
# aps = [Seg2(x,i,i+1) for i in 1:10]
# sp = RegularGrid(x)
#
# xₘ = Node(0.5,0.0,0.0)
# indices = sp(xₘ)
#
# sp(aps)

pool = SparseShapePool(5,Val(:∂1))
@btime SparseShapePool(5,Val(:∂1))
push!(pool.𝝭,rand(10)...)
pool.index[2] = pool.index[1]+10

𝝭 = SparseShape(pool,1,Val(:∂1))
@btime SparseShape(pool,1,Val(:∂1))
a = 𝝭[1]
@btime a = $𝝭[1]
v = Val(:∂1)
# @btime pool(1,v)
@btime a = $pool.𝝭[2]

a = rand(10)
@btime c = $a[1]
