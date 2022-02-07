
import Base:*
function *(a::NTuple{N,Float64},b::Float64) where N
    return (a[i]*b for i in 1:N)
end
