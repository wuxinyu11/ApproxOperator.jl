
## PointData
struct Node{N,T}
    i::Int
    data::NTuple{N,Vector{T}}
end
Node(i::Int,datas::AbstractVector...) = Node(i,datas)

getindex(pd::Node,i::Int) = pd.data[i][pd.i]
