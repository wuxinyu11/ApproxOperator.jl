using TimerOutputs

function efficiency()

to = TimerOutput()

@timeit to "node" begin
    @timeit to "construct node" x = Node(1.0,2.0,3.0)
    @timeit to "getindex node" x₁ = x[1]
end
a₁ = 1.0
a₂ = 1.0
n = 10
x = [Node(a₁/n*i,a₁/n*j,0.) for j in 0:n for i in 0:n]
ξ₁ = 0.0
ξ₂ = (0.0,0.0)
ξ₃ = (0.0,0.0,0.0)
@timeit to "cell" begin
    @timeit to "AbstractPoi" begin
        @timeit to "construct Poi1" ap = Poi1(25,x)
        @timeit to "get_jacobe" get_jacobe(ap,ξ₁)
        @timeit to "get_coordinates" get_coordinates(ap,ξ₁)
        @timeit to "get_shape_functions Poi1" get_shape_functions(ap,ξ₁,Val(:∂1))
    end
    @timeit to "AbstractSeg" begin
        @timeit to "construct Seg2" ap = Seg2([1,2],x)
        @timeit to "get_jacobe" get_jacobe(ap,ξ₁)
        @timeit to "get_coordinates" get_coordinates(ap,ξ₁)
        @timeit to "get_shape_functions Seg2 ∂1" get_shape_functions(ap,ξ₁,Val(:∂1))
        @timeit to "get_shape_functions Seg2 ∂x" get_shape_functions(ap,ξ₁,Val(:∂x))
        @timeit to "get_shape_functions Seg2 ∂y" get_shape_functions(ap,ξ₁,Val(:∂y))
        # @timeit to "get_shape_functions Seg2 ∂1 ∂x ∂y" get_shape_functions(ap,ξ₁,Val(:∂1),Val(:∂x),Val(:∂y))
        # @timeit to "get_shape_functions Seg2 ∂1 ∂x ∂y ∂z" get_shape_functions(ap,ξ₁,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
    end
    @timeit to "AbstractTri" begin
        @timeit to "construct Tri3" ap = Tri3([1,2,12],x)
        @timeit to "get_jacobe" get_jacobe(ap,ξ₂)
        @timeit to "get_coordinates" get_coordinates(ap,ξ₂)
        @timeit to "get_shape_functions Tri3 ∂1" get_shape_functions(ap,ξ₂,Val(:∂1))
        @timeit to "get_shape_functions Tri3 ∂x" get_shape_functions(ap,ξ₂,Val(:∂x))
        @timeit to "get_shape_functions Tri3 ∂y" get_shape_functions(ap,ξ₂,Val(:∂y))
        # @timeit to "get_shape_functions Tri3 ∂1 ∂x ∂y" get_shape_functions(ap,ξ₂,Val(:∂1),Val(:∂x),Val(:∂y))
        # @timeit to "get_shape_functions Tri3 ∂1 ∂x ∂y ∂z" get_shape_functions(ap,ξ₂,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
    end
    @timeit to "AbstractQuad" begin
        @timeit to "construct Quad" ap = Quad([1,2,13,12],x)
        @timeit to "get_jacobe" get_jacobe(ap,ξ₂)
        @timeit to "get_coordinates" get_coordinates(ap,ξ₂)
        @timeit to "get_shape_functions Quad ∂1" get_shape_functions(ap,ξ₂,Val(:∂1))
        # @timeit to "get_shape_functions Quad ∂x" get_shape_functions(ap,ξ₂,Val(:∂x))
        # @timeit to "get_shape_functions Quad ∂y" get_shape_functions(ap,ξ₂,Val(:∂y))
        # @timeit to "get_shape_functions Quad ∂1 ∂x ∂y" get_shape_functions(ap,ξ₂,Val(:∂1),Val(:∂x),Val(:∂y))
        # @timeit to "get_shape_functions Quad ∂1 ∂x ∂y ∂z" get_shape_functions(ap,ξ₂,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
    end
end
show(to)

end
