using Revise
using ApproxOperator
using BenchmarkTools

n = 10
a = rand(n)
b = zeros(n)
c = zeros(n)
x = [Node(i,a,b,c) for i in 1:n]
ξ = [Node(i,a,b) for i in 1:5]
@btime Node(1,$a,$b,$c)

ap = Seg2(x,ξ)
t = x[:𝝭,1]
@btime x[:𝝭,1]
