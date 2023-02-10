
"""
set_memory_𝝭!(ap::T,ss::Symbol...) where T<:AbstractElement
"""
const shape_function = (
    𝝭=(:𝝭,),∇𝝭=(:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z),∇₂𝝭=(:𝝭,:∂𝝭∂x,:∂𝝭∂y),∇̃₂𝝭=(:∂𝝭∂x,:∂𝝭∂y),
    ∇²𝝭=(:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂²𝝭∂x∂z,:∂²𝝭∂y∂z,:∂²𝝭∂z²),
    ∇²₂𝝭=(:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²),∇̃²𝝭=(:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²),
    ∇³𝝭=(:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂³𝝭∂x³,:∂³𝝭∂x²∂y,:∂³𝝭∂x∂y²,:∂³𝝭∂y³),
    ∇∇̃²𝝭=(:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂²𝝭∂x²_,:∂²𝝭∂x∂y_,:∂²𝝭∂y²_,:∂∂²𝝭∂x²∂x,:∂∂²𝝭∂x²∂y,:∂∂²𝝭∂x∂y∂x,:∂∂²𝝭∂x∂y∂y,:∂∂²𝝭∂y²∂x,:∂∂²𝝭∂y²∂y,:∂∂²𝝭∂x²∂x_,:∂∂²𝝭∂x²∂y_,:∂∂²𝝭∂x∂y∂x_,:∂∂²𝝭∂x∂y∂y_,:∂∂²𝝭∂y²∂x_,:∂∂²𝝭∂y²∂y_),
    ∇̃ᵥ₂𝝭=(:𝝭,:∂𝝭∂x,:∂𝝭∂y),
    test=(:𝝭,:∂𝝭∂x,:∂𝝭∂x_)
)
const moment_matrix = (
    𝝭=(:𝗠,),∇𝝭=(:𝗠,:∂𝗠∂x,:∂𝗠∂y,:∂𝗠∂z),∇₂𝝭=(:𝗠,:∂𝗠∂x,:∂𝗠∂y),∇̃₂𝝭=(:∇̃,),
    ∇²𝝭=(:𝗠,:∂𝗠∂x,:∂𝗠∂y,:∂𝗠∂z,:∂²𝗠∂x²,:∂²𝗠∂x∂y,:∂²𝗠∂y²,:∂²𝗠∂x∂z,:∂²𝗠∂y∂z,:∂²𝗠∂z²),
    ∇²₂𝝭=(:𝗠,:∂𝗠∂x,:∂𝗠∂y,:∂²𝗠∂x²,:∂²𝗠∂x∂y,:∂²𝗠∂y²),∇̃²𝝭=(:∇̃²,),
    ∇³𝝭=(:𝗠,:∂𝗠∂x,:∂𝗠∂y,:∂²𝗠∂x²,:∂²𝗠∂x∂y,:∂²𝗠∂y²,:∂³𝗠∂x³,:∂³𝗠∂x²∂y,:∂³𝗠∂x∂y²,:∂³𝗠∂y³),
    ∇∇̃²𝝭=(:𝗠,:∂𝗠∂x,:∂𝗠∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η),
    ∇̃ᵥ₂𝝭=(:𝗠,:∇̃),
    test=(:𝗠,:∂𝗠∂x,:∇̃)
)
function set_memory_𝝭!(aps::Vector{T},ss::Symbol...) where T<:AbstractElement
    n = getnₛ(aps)
    data = getfield(aps[1].𝓖[1],:data)
    for s in ss
        push!(data,s=>(4,zeros(n)))
    end
end

"""
set_memory_𝗠!(aps::Vector{T},ss::Symbol...) where T<:ReproducingKernel
"""
function set_memory_𝗠!(aps::Vector{T},ss::Symbol...) where T<:ReproducingKernel
    data = getfield(aps[1].𝓖[1],:data)
    for s in ss
        if s == :∇̃
            n = get𝑛𝒑₁(aps[1])
        elseif s ∈ (:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
            n = get𝑛𝒑₂(aps[1])
        else
            n = get𝑛𝒑(aps[1])
        end
        m = Int(n*(n+1)/2)
        push!(data,s=>(0,zeros(m)))
    end
end

function set_memory_𝗠!(ap::T,ss::Symbol... = keys(ap[1].𝗠)...) where T<:ReproducingKernel
    n = get𝑛𝒑(ap)
    empty!(ap.𝗠)
    for s in ss
        if s == :∇̃
            n₁ = get𝑛𝒑₁(ap)
            ap.𝗠[s] = SymMat(n₁)
        elseif s ∈ (:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
            n₂ = get𝑛𝒑₂(ap)
            ap.𝗠[s] = SymMat(n₂)
        else
            ap.𝗠[s] = SymMat(n)
        end
    end
end

"""
importmsh
"""
function importmsh(filename::String)
    fid = open(filename,"r")
    readline(fid)
    line = readline(fid)
    v_,f_,d_ = split(line," ")
    version = parse(Float64,v_)
    filetype = parse(Int,f_)
    datasize = parse(Int,d_)
    readline(fid)
    if version == 4.1
        elements,nodes = import_msh_4(fid)
    elseif version == 2.2
        elements,nodes = import_msh_2(fid)
    else
        println("Version does not match!")
    end
    return elements, nodes
end

function import_msh_4(fid::IO) end

function import_msh_2(fid::IO)
    etype = Dict(1=>:Seg2,2=>:Tri3,3=>:Quad,8=>:Seg3,9=>:Tri6,15=>:Poi1)
    nodes = Dict{Symbol,Vector{Float64}}()
    elements = Dict{String,Any}()
    physicalnames = Dict{Int,String}()
    for line in eachline(fid)
        if line == "\$PhysicalNames"
            numPhysicalNames = parse(Int,readline(fid))
            physicalnames=>Dict{Int,String}()
            for i in 1:numPhysicalNames
                line = readline(fid)
                d_,p_,n_ = split(line," ")
                dimension = parse(Int,d_)
                physicalTag = parse(Int,p_)
                name = strip(n_,'\"')
                physicalnames[physicalTag] = name
            end
            readline(fid)
        elseif line == "\$Nodes"
            line = readline(fid)
            nₚ = parse(Int,line)
            x = zeros(nₚ)
            y = zeros(nₚ)
            z = zeros(nₚ)
            for i in 1:nₚ
                line = readline(fid)
                t_,x_,y_,z_ = split(line," ")
                tag = parse(Int,t_)
                x[i] = parse(Float64,x_)
                y[i] = parse(Float64,y_)
                z[i] = parse(Float64,z_)
            end
            nodes[:x] = x
            nodes[:y] = y
            nodes[:z] = z
            nodes = Node(nodes...)
            readline(fid)
        elseif line == "\$Elements"
            line = readline(fid)
            nₑ = parse(Int,line)
            for i in 1:nₑ
                line = readline(fid)
                entries = split(line," ")
                elmN_ = entries[1]
                elmT_ = entries[2]
                numT_ = entries[3]
                phyT_ = entries[4]
                elmE_ = entries[5]
                l_ = entries[6:end]
                elmNumber = parse(Int,elmN_)
                elmType = parse(Int,elmT_)
                numTag = parse(Int,numT_)
                phyTag = parse(Int,phyT_)
                elmEntary = parse(Int,elmE_)
                nodeList = parse.(Int,l_)
                name = physicalnames[phyTag]
                type = etype[elmType]
                haskey(elements,name) ? push!(elements[name],Element{type}([nodes[i] for i in nodeList])) : elements[name]=Element{type}[Element{type}([nodes[i] for i in nodeList])]
            end
        end
    end
    return elements, nodes
end

function importmsh(filename::String,config::Dict{T,Any}) where T<:Any
    elms, nodes = importmsh(filename)
    return generate(elms,nodes,config)
end
function importmsh(file_elements::String,file_nodes::String,config::Dict{T,Any}) where T<:Any
    elms, nodes_ = importmsh(file_elements)
    elms_, nodes = importmsh(file_nodes)
    return generate(elms,nodes,config)
end
function generate(elms::Dict{String,Any},nodes::Vector{Node},config::Dict{T,Any}) where T<:Any
    elements = Dict{String,Any}()
    if haskey(config,"RegularGrid")
        x = getfield(nodes[1],:data)[:x][2]
        y = getfield(nodes[1],:data)[:y][2]
        z = getfield(nodes[1],:data)[:z][2]
        n = config["RegularGrid"]["n"]
        γ = config["RegularGrid"]["γ"]
        sp = RegularGrid(x,y,z,n=n,γ=γ)
        delete!(config,"RegularGrid")
    else
        sp = nothing
    end
    if haskey(config,"BoundaryDofs")
        dofs,ndofs = getboundarydofs(elms["Ω"])
        cfg = config["BoundaryDofs"]
        element_type = eval(Meta.parse(cfg["type"]))
        elements["∂Ω"] = Vector{element_type}(undef,ndofs)
        for (ids,n) in dofs
            elements["∂Ω"][n] = element_type([nodes[i] for i in ids],SNode[])
        end
        integration_type = Meta.parse(cfg["𝓖"]["type"])
        set𝓖!(elements["∂Ω"],integration_type)
        delete!(config,"BoundaryDofs")
    end

    for (name,cfg) in config
         # set𝓖
        element_tag = cfg["𝓒"]["tag"]
        element_type = eval(Meta.parse(cfg["type"]))
        integration_tag = haskey(cfg["𝓖"],"tag") ? cfg["𝓖"]["tag"] : element_tag
        integration_type = Meta.parse(cfg["𝓖"]["type"])
        elements[name] = element_type[]
        if haskey(elms,integration_tag)
            set𝓖!(elms[integration_tag],integration_type)
            if haskey(cfg["𝓖"],"normal") set𝒏!(elms[integration_tag]) end
            if integration_tag ≠ element_tag
                elms[element_tag*"∩"*integration_tag] = unique!(elms[element_tag]∩elms[integration_tag])
                element_tag = element_tag*"∩"*integration_tag
                set𝓖!(elms[element_tag],elms[integration_tag])
            end
            if haskey(cfg["𝓖"],"𝐶")
                𝐶_tag = cfg["𝓖"]["𝐶"]
                set𝐶!(elms[element_tag],elements[𝐶_tag])
            end

            # set 𝓒
            nₑ = length(elms[element_tag])
            if element_type<:Element
                for elm in elms[element_tag]
                    𝓒 = [x for x in elm.𝓒]
                    𝓖 = [ξ for ξ in elm.𝓖]
                    push!(elements[name],element_type(𝓒,𝓖))
                end
            elseif element_type<:ReproducingKernel
                if haskey(cfg["𝓒"],"type")
                    for elm in elms[element_tag]
                        𝓖 = [ξ for ξ in elm.𝓖]
                        push!(elements[name],element_type(Node[],𝓖))
                    end
                    position_type= Meta.parse(cfg["𝓒"]["type"])
                    set𝓖!(elms[element_tag],position_type)
                    for (c,elm) in enumerate(elms[element_tag])
                        𝓒 = [nodes[i] for i in sp(elm.𝓖)]
                        push!(elements[name][c].𝓒,𝓒...)
                    end
                else
                    for elm in elms[element_tag]
                        𝓒 = [nodes[i] for i in sp(elm.𝓒)]
                        𝓖 = [ξ for ξ in elm.𝓖]
                        push!(elements[name],element_type(𝓒,𝓖))
                    end
                end
                s = 0
                for elm in elements[name]
                    𝓖 = elm.𝓖
                    data = getfield(𝓖[1],:data)
                    n = length(elm.𝓒)
                    for (i,ξ) in enumerate(𝓖)
                        g = ξ.𝑔
                        G = ξ.𝐺
                        C = ξ.𝐶
                        𝓖[i] = SNode((g,G,C,s),data)
                        s += n
                    end
                end
            elseif element_type<:TRElement
                data = getfield(nodes[1],:data)
                for elm in elms[element_tag]
                    nodeList = (x.𝐼 for x in elm.𝓒)
                    𝓒 = [GNode((i,dofs[Set(setdiff(nodeList,i))]),data) for i in nodeList]
                    𝓖 = [ξ for ξ in elm.𝓖]
                    push!(elements[name],element_type(𝓒,𝓖))
                end
            end
            if contains(element_tag,"∩")
                element_tag_type = eval(Meta.parse(config[cfg["𝓒"]["tag"]]["type"]))
                element_tag_integration_type = Meta.parse(config[cfg["𝓒"]["tag"]]["𝓖"]["type"])
                set𝓖!(elms[element_tag],element_tag_integration_type)
                elements[element_tag] = element_tag_type[]
                for elm in elms[element_tag]
                    𝓒 = [nodes[i] for i in sp(elm.𝓒)]
                    𝓖 = [ξ for ξ in elm.𝓖]
                    push!(elements[element_tag],element_tag_type(𝓒,𝓖))
                end
                set_memory_𝝭!(elements[element_tag],:𝝭,:∂𝝭∂x,:∂𝝭∂y)
                set_memory_𝗠!(elements[element_tag],:𝗠,:∂𝗠∂x,:∂𝗠∂y)
            end

            # set shape memory
            if haskey(cfg,"𝓖")
                if haskey(cfg["𝓖"],"𝝭") set_memory_𝝭!(elements[name],shape_function[Meta.parse(cfg["𝓖"]["𝝭"])]...) end
                if element_type<:ReproducingKernel set_memory_𝗠!(elements[name],moment_matrix[Meta.parse(cfg["𝓖"]["𝝭"])]...) end
            end
        end
    end
    return elements,nodes
end

function getboundarydofs(elements::Vector{T}) where T<:AbstractElement{:Tri3}
    dofs = Dict{Set{Int},Int}()
    idBoundaries = ((1,2),(2,3),(3,1))
    n = 0
    for elm in elements
        𝓒 = elm.𝓒
        for bc in idBoundaries
            dof = Set(𝓒[i].𝐼 for i in bc)
            if ~haskey(dofs,dof)
                n += 1
                dofs[dof] = n
            end
        end
    end
    return dofs,n
end

function voronoimsh(filename::String)
    elms, nds = importmsh(filename)
    nᵥ = 0
    for (name,elm) in elms
        nᵥ += length(elm)
    end
    data = Dict{Symbol,Tuple{Int,Vector{Float64}}}()
    elements = Dict{String,Any}()
    ind = 0
    xv = Float64[]
    yv = Float64[]
    zv = Float64[]
    push!(data,:x=>(1,xv))
    push!(data,:y=>(1,yv))
    push!(data,:z=>(1,zv))
    nodelist = Dict(node=>Int[] for node in nds)
    nodelist_Γᵗ = Dict{Node,Vector{Int}}()
    nodelist_Γᵍ = Dict{Node,Vector{Int}}()
    normal = Dict{Node,Set{Tuple{Float64,Float64}}}()
    for (name,elm) in elms
        if contains(name,"Ω")
            for el in elm
                ind += 1
                x₁ = el.𝓒[1].x
                y₁ = el.𝓒[1].y
                x₂ = el.𝓒[2].x
                y₂ = el.𝓒[2].y
                x₃ = el.𝓒[3].x
                y₃ = el.𝓒[3].y
                xₘ₁ = 0.5*(x₂+x₃)
                yₘ₁ = 0.5*(y₂+y₃)
                xₘ₂ = 0.5*(x₃+x₁)
                yₘ₂ = 0.5*(y₃+y₁)
                s₁ = (y₃-y₂)/(x₃-x₂)
                s₂ = (y₁-y₃)/(x₁-x₃)
                if x₃ ≈ x₂
                    xc = xₘ₂ + s₂*(yₘ₂-yₘ₁)
                    yc = yₘ₁
                elseif x₁ ≈ x₃
                    xc = xₘ₁ + s₁*(yₘ₁-yₘ₂)
                    yc = yₘ₂
                elseif y₃ ≈ y₂
                    xc = xₘ₁
                    yc = -1/s₂*(xₘ₁-xₘ₂)+yₘ₂
                elseif y₁ ≈ y₃
                    xc = xₘ₂
                    yc = -1/s₁*(xₘ₂-xₘ₁)+yₘ₁
                else
                    xc = (s₁*s₂*(yₘ₁-yₘ₂)+(s₂*xₘ₁-s₁*xₘ₂))/(s₂-s₁)
                    yc = ((s₂*yₘ₂-s₁*yₘ₁)+(xₘ₂-xₘ₁))/(s₂-s₁)
                end
                push!(xv,xc)
                push!(yv,yc)
                push!(zv,0.0)
                for xᵢ in el.𝓒
                    push!(nodelist[xᵢ],ind)
                end
            end
        elseif contains(name,"Γ")
            for el in elm
                ind += 1
                x₁ = el.𝓒[1].x
                y₁ = el.𝓒[1].y
                x₂ = el.𝓒[2].x
                y₂ = el.𝓒[2].y
                L = ((x₁-x₂)^2+(y₁-y₂)^2)^0.5
                normal_ = ((y₂-y₁)/L,(x₁-x₂)/L)
                xc = 0.5*(x₁+x₂)
                yc = 0.5*(y₁+y₂)
                push!(xv,xc)
                push!(yv,yc)
                push!(zv,0.0)
                if name == "Γᵍ"
                    for xᵢ in el.𝓒
                        push!(nodelist[xᵢ],ind)
                        haskey(nodelist_Γᵍ,xᵢ) ? push!(nodelist_Γᵍ[xᵢ],ind) : nodelist_Γᵍ[xᵢ] = [ind]
                        haskey(normal,xᵢ) ? push!(normal[xᵢ],normal_) : normal[xᵢ] = Set([normal_])
                    end
                elseif name == "Γᵗ"
                    for xᵢ in el.𝓒
                        push!(nodelist[xᵢ],ind)
                        haskey(nodelist_Γᵗ,xᵢ) ? push!(nodelist_Γᵗ[xᵢ],ind) : nodelist_Γᵗ[xᵢ] = [ind]
                        haskey(normal,xᵢ) ? push!(normal[xᵢ],normal_) : normal[xᵢ] = Set([normal_])
                    end
                end
            end
        end
    end
    for (xᵢ,normal_) in normal
        if length(normal_) ≠ 1
            ind += 1
            push!(xv,xᵢ.x)
            push!(yv,xᵢ.y)
            push!(zv,0.0)
            push!(nodelist[xᵢ],ind)
            if haskey(nodelist_Γᵍ,xᵢ) push!(nodelist_Γᵍ[xᵢ],ind) end
            if haskey(nodelist_Γᵗ,xᵢ) push!(nodelist_Γᵗ[xᵢ],ind) end
        end
    end
    nodes = [Node(i,data) for i in 1:ind]
    elements["Ω"] = Element{:Vor2}[]
    elements["Γᵗ"] = Element{:Seg2}[]
    elements["Γᵍ"] = Element{:Seg2}[]
    # for (node,list) in nodelist
    for node in nds
        list = nodelist[node]
        # sort
        # cal centroid
        xc = 0.
        yc = 0.
        for i in list
            xc += xv[i]
            yc += yv[i]
        end
        xc = xc/length(list)
        yc = yc/length(list)
        α = [atan(yv[i]-yc,xv[i]-xc) for i in list]
        p = sortperm(α)
        push!(elements["Ω"],Element{:Vor2}([nodes[list[i]] for i in p],SNode[]))
    end
    (xmin,xmax) = extrema(getfield(nds[1],:data)[:x][2])
    (ymin,ymax) = extrema(getfield(nds[1],:data)[:y][2])
    xc = 0.5*(xmin+xmax)
    yc = 0.5*(ymin+ymax)
    for (name,nodelist_) in (("Γᵍ",nodelist_Γᵍ),("Γᵗ",nodelist_Γᵗ))
        for (node,list) in nodelist_
            if length(list) == 2
                x₁ = xv[list[1]]
                x₂ = xv[list[2]]
                y₁ = yv[list[1]]
                y₂ = yv[list[2]]
                xₘ = 0.5*(x₁+x₂)
                yₘ = 0.5*(y₁+y₂)
                n₁ = y₂-y₁
                n₂ = x₁-x₂
                nc₁ = x₁-xc
                nc₂ = y₁-yc
                if n₁*nc₁+n₂*nc₂ > 0
                    push!(elements[name],Element{:Seg2}([nodes[list[1]],nodes[list[2]]],SNode[]))
                else
                    push!(elements[name],Element{:Seg2}([nodes[list[2]],nodes[list[1]]],SNode[]))
                end
            else
                xc_ = 0.
                yc_ = 0.
                for i in list
                    xc_ += xv[i]
                    yc_ += yv[i]
                end
                xc_ = xc_/length(list)
                yc_ = yc_/length(list)
                α = [atan(yv[i]-yc_,xv[i]-xc_) for i in list]
                p = sortperm(α)
                for i in 1:length(p)
                    (I,J) = i ≠ length(p) ? (p[i],p[i+1]) : (p[i],p[1])
                    x₁ = xv[list[I]]
                    x₂ = xv[list[J]]
                    y₁ = yv[list[I]]
                    y₂ = yv[list[J]]
                    xₘ = 0.5*(x₁+x₂)
                    yₘ = 0.5*(y₁+y₂)
                    n₁ = y₂-y₁
                    n₂ = x₁-x₂
                    nc₁ = x₁-xc
                    nc₂ = y₁-yc
                    if n₁*nc₁+n₂*nc₂ > 0
                        push!(elements[name],Element{:Seg2}([nodes[list[I]],nodes[list[J]]],SNode[]))
                    end
                end
            end
        end
    end
    return elements, nodes
end

function voronoimsh(filename::String,config::Dict{T,Any}) where T<:Any
    elms, nds = voronoimsh(filename)
    ~, nodes = importmsh(filename)
    elements = Dict{String,Any}()
    x = getfield(nodes[1],:data)[:x][2]
    y = getfield(nodes[1],:data)[:y][2]
    z = getfield(nodes[1],:data)[:z][2]
    n = config["RegularGrid"]["n"]
    γ = config["RegularGrid"]["γ"]
    sp = RegularGrid(x,y,z,n=n,γ=γ)
    delete!(config,"RegularGrid")

    for (name,cfg) in config
        element_tag = cfg["𝓒"]["tag"]
        element_type = eval(Meta.parse(cfg["type"]))
        integration_tag = haskey(cfg["𝓖"],"tag") ? cfg["𝓖"]["tag"] : element_tag
        integration_type = Meta.parse(cfg["𝓖"]["type"])
        elements[name] = element_type[]
        if haskey(elms,integration_tag)
            if integration_type == :Node
                data_ = getfield(nodes[1],:data)
                data = Dict([:x=>data_[:x],:y=>data_[:y],:z=>data_[:z]])
                G = 0
                for (c,elm) in enumerate(elms[integration_tag])
                    𝓖 = elm.𝓖
                    indices = sp(nodes[c])
                    for g in indices
                        G += 1
                        ξ = SNode((g,G,c,0),data)
                        push!(𝓖,ξ)
                    end
                end
                set𝐴!(elms[integration_tag])
                set𝒙ₘ!(elms[integration_tag])
                setm2!(elms[integration_tag])
            else
                set𝓖!(elms[integration_tag],integration_type)
            end
            if haskey(cfg["𝓖"],"normal") set𝒏!(elms[integration_tag]) end
            if integration_tag ≠ element_tag
                elms[element_tag*"∩"*integration_tag] = unique!(elms[element_tag]∩elms[integration_tag])
                element_tag = element_tag*"∩"*integration_tag
                set𝓖!(elms[element_tag],elms[integration_tag])
            end
            nₑ = length(elms[element_tag])
            if haskey(cfg["𝓒"],"type")
                for elm in elms[element_tag]
                    𝓖 = [ξ for ξ in elm.𝓖]
                    push!(elements[name],element_type(Node[],𝓖))
                end
                position_type= Meta.parse(cfg["𝓒"]["type"])
                if position_type == :Node
                    data_ = getfield(nodes[1],:data)
                    data = Dict([:x=>(2,data_[:x][2]),:y=>(2,data_[:y][2]),:z=>(2,data_[:z][2])])
                    for (c,elm) in enumerate(elms[element_tag])
                        𝓖 = elm.𝓖
                        indices = sp(nodes[c])
                        for (g,G) in enumerate(indices)
                            ξ = SNode((g,G,c,0),data)
                            push!(𝓖,ξ)
                        end
                    end
                else
                    set𝓖!(elms[element_tag],position_type)
                end
                for (c,elm) in enumerate(elms[element_tag])
                    𝓒 = [nodes[i] for i in sp(elm.𝓖)]
                    push!(elements[name][c].𝓒,𝓒...)
                end
            else
                for elm in elms[element_tag]
                    𝓒 = [nodes[i] for i in sp(elm.𝓒)]
                    𝓖 = [ξ for ξ in elm.𝓖]
                    push!(elements[name],element_type(𝓒,𝓖))
                end
            end
            s = 0
            for elm in elements[name]
                𝓖 = elm.𝓖
                data = getfield(𝓖[1],:data)
                n = length(elm.𝓒)
                for (i,ξ) in enumerate(𝓖)
                    g = ξ.𝑔
                    G = ξ.𝐺
                    C = ξ.𝐶
                    𝓖[i] = SNode((g,G,C,s),data)
                    s += n
                end
            end

            # set shape memory
            if haskey(cfg,"𝓖")
                if haskey(cfg["𝓖"],"𝝭") set_memory_𝝭!(elements[name],shape_function[Meta.parse(cfg["𝓖"]["𝝭"])]...) end
                if element_type<:ReproducingKernel set_memory_𝗠!(elements[name],moment_matrix[Meta.parse(cfg["𝓖"]["𝝭"])]...) end
            end
        end
    end
    return elements, nodes
end

