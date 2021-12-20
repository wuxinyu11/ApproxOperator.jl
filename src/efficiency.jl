using TimerOutputs

function efficiency()

to = TimerOutput()

xᵢ = [0.5,0.5,0.5]
Δx = @SVector [1.,2.,3.]

## ---------------- Operator ----------------
@timeit to "Operator" begin
# @timeit to "Potential Ω" op = Potential_Ω(b)
@timeit to "Potential Ω" op = Potential_Ω((x,y,z)-> -2exp(x+y+z))
@timeit to "force" op.b(xᵢ...)
end # Operator end
## ---------------- 1D --------------------
@timeit to "1D test" begin
ξᵢ = 0.
a = 1.
nₚ = 11
nₑ = nₚ-1

x = [Node(1/nₑ*i,0.,0.) for i in 0:nₑ]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)


@timeit to "SegM" begin
id = [5,6]
@timeit to "Regular grid" rg = RegularGrid(x,n=2,γ=1)
@timeit to "Linear basis function" bf = Linear1D(:∂1,:∂x,:∂y,:∂z)
@timeit to "Kernel function" kf = TensorProductKernel(:∂1,:∂x,:∂y,:∂z,ss=[1.5*a/nₑ,1.,1.],nm=nₚ)
@timeit to "SegM" ap = SegM(x,id,bf=bf,kf=kf,sp=rg,qw=:SegGI2)
# @timeit to "get basis function ∂1" p = get_basis_function(ap,Δx,:∂1)
# @timeit to "get basis function ∂1 debug" p = get_basis_function(ap.bf,Δx,Val(:∂1))
# @timeit to "get basis function ∂x debug" p = get_basis_function(ap.bf,Δx,Val(:∂x))
# @timeit to "get basis function ∂1 ∂x ∂y ∂z" p, px, py, pz = get_basis_function(ap,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
# @timeit to "get kernel function ∂1 ∂x ∂y ∂z" w, wx, wy, wz = get_kernel_function(ap,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
@timeit to "Meshfree shape function ∂1" get_shape_functions(ap,ξᵢ,Val(:∂1))
@timeit to "Meshfree shape function ∂1 ∂x" get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x))
# @code_warntype get_basis_function(ap,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
# @code_warntype get_kernel_function(ap,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
@code_warntype get_shape_functions(ap,ξᵢ,Val(:∂1))
@code_warntype get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
end # SegM end

end # 1D test end
## ---------------- 2D --------------------
@timeit to "2D test" begin
ξᵢ = [0.,0.]
a₁ = 1.
a₂ = 1.

n₁ = 10
n₂ = n₁
nₚ = (n₁+1)*(n₂+1)
nₑ = n₁*n₂

x = [Node(a₁/n₁*i,a₁/n₂*j,0.) for j in 0:n₂ for i in 0:n₁]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

A = SymMat(3,rand(6))
U = SymMat(3,rand(6))

# -------------- Linear algebra --------------
@timeit to "Linear algebra" begin
@timeit to "UUᵀ" UUᵀ!(A)
@timeit to "UᵀAU" UᵀAU!(A,U)
@timeit to "UAUᵀ" UAUᵀ!(A,U)
@timeit to "UUᵀAUUᵀ" UUᵀAUUᵀ!(A,U)
@timeit to "negective" nA = -A

end # linear algebra end

# -------------- Tri3 ----------------
id = [34,35,45]
@timeit to "2D Tri3" begin
@timeit to "Tri3" ap = Tri3(x,id)
@timeit to "get coordinates" get_coordinates(ap,ξᵢ)
@timeit to "get jacobe" get_jacobe(ap,ξᵢ)
@timeit to "get number of indices" get_number_of_indices(ap)
@timeit to "Tri3 shape function ∂1" get_shape_functions(ap,ξᵢ,Val(:∂1))
@timeit to "Tri3 shape function ∂1,∂x" get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y))
@timeit to "Tri3 shape function ∂1,∂x,∂y,∂z" get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
# @code_warntype get_shape_functions(ap,ξᵢ,:∂1)
# @code_warntype get_shape_functions(ap,ξᵢ,:∂1,:∂x,:∂y,:∂z)
@timeit to "get global indices" get_global_indice(ap,1)
@timeit to "Assembly Tri3 potential Ω" op(ap,k,f)

end # 2D Tri3 end

# -------------- TriM ----------------
ξᵢ = [0.,0.,0.]
id = [34,35,45]
@timeit to "2D TriM" begin
@timeit to "Regular grid" rg = RegularGrid(x,n=2,γ=1)
@timeit to "Linear basis function" bf = Linear2D(:∂1,:∂x,:∂y,:∂z)
@timeit to "Kernel function" kf = TensorProductKernel(:∂1,:∂x,:∂y,:∂z,ss=[1.5*a₁/n₁,1.5*a₂/n₂,1.],nm=nₚ)
@timeit to "TriM" ap = TriM(x,id,bf=bf,kf=kf,sp=rg,qw=:TriGI3)
@timeit to "Meshfree shape function ∂1" get_shape_functions(ap,ξᵢ,Val(:∂1))
@timeit to "Meshfree shape function ∂1 ∂x ∂y ∂z" get_shape_functions(ap,ξᵢ,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
# @timeit to "get basis function ∂1" p = get_basis_function(ap,Δx,:∂1)
# @timeit to "get basis function ∂1 debug" p = get_basis_function(ap.bf,Δx,Val(:∂1))
# @timeit to "get basis function ∂x debug" p = get_basis_function(ap.bf,Δx,Val(:∂x))
# @timeit to "get basis function ∂1 ∂x ∂y ∂z" p, px, py, pz = get_basis_function(ap,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
# @timeit to "get kernel function ∂1" w = get_kernel_function(ap,Δx,:∂1)
# @timeit to "get kernel function ∂1 ∂x ∂y ∂z" w, wx, wy, wz = get_kernel_function(ap,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
# @timeit to "calulate moment matrix ∂1" m = cal_moment_matrix!(ap,xᵢ,:∂1)
# @timeit to "calulate moment matrix ∂1 ∂x ∂y ∂z" m, mx, my, mz = cal_moment_matrix!(ap,xᵢ,:∂1,:∂x,:∂y,:∂z)
# @timeit to "Assembly TriM potential Ω" op(ap,k,f)
# @code_warntype get_basis_function(ap,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
# @code_warntype get_kernel_function(ap,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))

end # 2D TriM end

end # 2D end

show(to)

end
