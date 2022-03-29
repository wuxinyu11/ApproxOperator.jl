
using Revise, ApproxOperator, LaTeXStrings, GLMakie

elements, nodes = importmsh("./msh/patchtest.msh")
nₚ = length(nodes[:x])

type = (Node,:Quadratic2D,:□,:CubicSpline)
elements["Ω"] = ReproducingKernel{type...,:Tri3}(elements["Ω"])
s = 2.5*0.25*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝓖!(elements["Ω"],:TriGI3,:∂1)

# n = 12
n = 21
inte = 100
# as = collect(0.0:1/inte:1.0)
as = LinRange(0.0,1.0,inte)

𝗠 = elements["Ω"][1].𝗠
𝝭 = elements["Ω"][1].𝝭
ap = ReproducingKernel{type...,:Node}([Node(i,nodes) for i in 1:nₚ],Node[],𝗠,𝝭)

# x = [x for x in as for y in as]
# y = [y for x in as for y in as]
z = [get𝝭(ap,(x,y,0.0),[n])[1] for x in as, y in as]

# plot(x,y,z,
#     st=:surface,
#     c=:jet1,
#     xlims=(0.0,1.0),
#     ylims=(0.0,1.0),
#     zlims=(-0.1,1.0),
#     xlabel=L"x",
#     ylabel=L"y",
#     zlabel=L"\Psi_I",
#     ztickfontrotation=10.0,
#     framestyle=:box,
#     legend=:none,
#     aspect_ratio=0.8,
#     dpi=300
#     )
f = Figure(fontsize=14)
ax1 = Axis3(f[1,1],
      xlabel="",
      ylabel="",
      zlabel="",
      limits=(0.0,1.0,0.0,1.0,-0.1,1.0),
      aspect=(1,1,0.7)
)
# ax2 = Axis3(f[1,1],limits=(0.0,1.0,0.0,1.0,-0.1,1.0))
surface!(as,as,z,colormap = :jet1)

meshscatter!([nodes[:x][n]], [nodes[:y][n]], [nodes[:z][n]],markersize = 0.05, color=:red)

save("shape1.png",f)
# ax = Axis3(figure[1,1],xlabel="x")
