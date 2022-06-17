
using Revise, ApproxOperator, BenchmarkTools

# A = ApproxOperator.SymMat(2)
# A.m .= [2,-2,5]
# A = ApproxOperator.SymMat(3)
# A.m .= [4,12,37,-16,-43,98]
A = ApproxOperator.SymMat(4)
A.m .= [6,12,-8,3,-13,-7,-6,4,1,6]
# @allocated ApproxOperator.cholesky_Gaxpy!(A)
ApproxOperator.cholesky!(A)
# ApproxOperator.inverse!(A)
# @btime ApproxOperator.permute!($A,1,3)
# A.m .= permute!(A.m,[3,2,1])
ApproxOperator.LLᵀ!(A)
# ApproxOperator.UUᵀ!(A)
# ApproxOperator.permute!(A)