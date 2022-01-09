using TimerOutputs

function efficiency()

to = TimerOutput()

nₚ = 11
nₑ = nₚ-1

x = [1/nₑ*i for i in 0:nₑ]

# poolx = DataPool([[1.0/nₑ*i for i in 0:10],zeros(nₚ),zeros(nₚ)],Dict(:x=>1,:y=>2,:z=>3))
# poolξ = DataPool([[0.0],[0.0],[0.0],[1.0]],Dict(:ξ=>1,:η=>2,:γ=>3,:w=>4))

data = Dict(:x=>x,:y=>zeros(nₚ),:z=>zeros(nₚ))
paradata = Dict(:ξ=>[0.0],:w=>[1.0])
physdata = Dict(:σ=>[0.0])
𝓒 = [Node(i,data) for i in 1:2]
𝓖 = [Point(1,paradata,1,physdata)]
@timeit to "node" begin
    ξ = 𝓖[1]
    @timeit to "construct node" node = Node(1,data)
    @timeit to "getindex id" node.I
    @timeit to "getindex node" node.coordinates
    @timeit to "getindex x" node.x
    @timeit to "getindex y" node.y
    @timeit to "getindex z" node.z
    @timeit to "getindex ξ" ξ.ξ
    @timeit to "getindex ξ coordinates" ξ.coordinates
end
@timeit to "cell" begin
    # @timeit to "AbstractPoi" begin
    #     @timeit to "construct Poi1" ap = Poi1(25,x)
    #     @timeit to "get_weight" get_weight(ap,ξ₁)
    #     @timeit to "get_coordinates" get_coordinates(ap,ξ₁)
    #     @timeit to "get_shape_functions Poi1" get_shape_functions(ap,ξ₁,Val(:∂1))
    # end
    # @timeit to "AbstractSeg" begin
    #     ξ = 𝓖[1]
    #     @timeit to "construct Seg2" ap = Seg2(𝓒,𝓖)
    #     @timeit to "get Jacobe" ap.J(0.0)
    #     @timeit to "get coordinates" ap.coordinates(ξ)
    #     @timeit to "get shape functions Seg2 ∂1" ap.𝝭(ξ)
    #     @timeit to "get shape functions Seg2 ∂x" ap.∂𝝭∂x(ξ)
    #     @timeit to "get shape functions Seg2 ∂y" ap.∂𝝭∂y(ξ)
    #     @timeit to "get shape functions Seg2 ∂z" ap.∂𝝭∂z(ξ)
    #     @timeit to "get shape functions Seg2 ∇" ap.∇𝝭(ξ)
        # @timeit to "get_shape_functions Seg2 ∂1" get_shape_functions(ap,ξ,Val(:∂1))
    #     @timeit to "get_shape_functions Seg2 ∂x" get_shape_functions(ap,ξ,Val(:∂x))
    #     @timeit to "get_shape_functions Seg2 ∂y" get_shape_functions(ap,ξ,Val(:∂y))
        # @timeit to "get_shape_functions Seg2 ∂1 ∂x ∂y" get_shape_functions(ap,ξ₁,Val(:∂1),Val(:∂x),Val(:∂y))
        # @timeit to "get_shape_functions Seg2 ∂1 ∂x ∂y ∂z" get_shape_functions(ap,ξ₁,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
    # end
    # 𝒞 = [1,2,12]
    # @timeit to "AbstractTri" begin
    #     @timeit to "construct Tri3" ap = Tri3(𝒞,x)
    #     @timeit to "get_weight" get_weight(ap,ξ₂)
    #     @timeit to "get_coordinates" get_coordinates(ap,ξ₂)
    #     @timeit to "get_shape_functions Tri3 ∂1" get_shape_functions(ap,ξ₂,Val(:∂1))
    #     @timeit to "get_shape_functions Tri3 ∂x" get_shape_functions(ap,ξ₂,Val(:∂x))
    #     @timeit to "get_shape_functions Tri3 ∂y" get_shape_functions(ap,ξ₂,Val(:∂y))
    #     @timeit to "get_shape_functions Tri3 ∂1 ∂x ∂y" get_shape_functions(ap,ξ₂,Val(:∂1),Val(:∂x),Val(:∂y))
    #     @timeit to "get_shape_functions Tri3 ∂1 ∂x ∂y ∂z" get_shape_functions(ap,ξ₂,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
    # end
    # 𝒞 = [1,2,13,12]
    # @timeit to "AbstractQuad" begin
    #     @timeit to "construct Quad" ap = Quad(𝒞,x)
    #     # @timeit to "get_weight" get_weight(ap,ξ₂)
    #     # @timeit to "get_coordinates" get_coordinates(ap,ξ₂)
    #     @timeit to "get_shape_functions Quad ∂1" get_shape_functions(ap,ξ₂,Val(:∂1))
    #     # @timeit to "get_shape_functions Quad ∂x ∂y" get_shape_functions(ap,ξ₂,Val(:∂x),Val(:∂y))
    #     # @timeit to "get_shape_functions Quad ∂x" get_shape_functions(ap,ξ₂,Val(:∂x))
    #     # @timeit to "get_shape_functions Quad ∂y" get_shape_functions(ap,ξ₂,Val(:∂y))
    #     # @timeit to "get_shape_functions Quad ∂1 ∂x ∂y" get_shape_functions(ap,ξ₂,Val(:∂1),Val(:∂x),Val(:∂y))
    #     # @timeit to "get_shape_functions Quad ∂1 ∂x ∂y ∂z" get_shape_functions(ap,ξ₂,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
    # end
end
show(to)

end
