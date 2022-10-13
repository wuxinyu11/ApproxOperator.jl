

"""
setgeometry!(ap::T) where T<:AbstractElement
"""
function setgeometry!(ap::T) where T<:AbstractElement
    𝓖 = ap.𝓖
    for x in 𝓖
        𝒙 = get𝒙(ap,x)
        𝑤 = get𝑤(ap,x)
        x.x = 𝒙[1]
        x.y = 𝒙[2]
        x.z = 𝒙[3]
        x.𝑤 = 𝑤
    end
end

"""
set_memory_𝝭!(ap::T,ss::Symbol...) where T<:AbstractElement
"""
function set_memory_𝝭!(aps::Vector{T},ss::Symbol...) where T<:AbstractElement
    n = sum(length(ap.𝓒)*length(ap.𝓖) for ap in aps)
    for s in ss
        push!(getfield(aps[1].𝓖[1],:data),s=>(3,zeros(n)))
    end
end

"""
set_memory_𝗠!(aps::Vector{T},ss::Symbol... = keys(aps[1].𝗠)...) where T<:ReproducingKernel
"""
function set_memory_𝗠!(aps::Vector{T},ss::Symbol... = keys(aps[1].𝗠)...) where T<:ReproducingKernel
    set_memory_𝗠!(aps[1],ss...)
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

## ---------------- msh ---------------
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
    etype = Dict(1=>:Seg2,2=>:Tri3,3=>:Quad,15=>:Poi1)
    nodes = Dict{Symbol,Vector{Float64}}()
    elements = Dict{String,Vector{Tuple{Symbol,Vector{Int}}}}()
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
                elements[name] = Vector{Tuple{Symbol,Vector{Int}}}()
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
                push!(elements[name],(type,nodeList))
            end
            return elements, nodes
        end
    end
end

function importmsh(filename::String,config::Dict{Any,Any})
    elms, nodes = importmsh(filename)
    elements = Dict{String,Any}()
    if haskey(config,"RegularGrid")
        cfg = config["RegularGrid"]
        sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z];n=cfg["n"],γ=cfg["γ"])
        delete!(config,"RegularGrid")
    else
        sp = nothing
    end
    if haskey(config,"IndependentDofs")
        for (k,v) in config["IndependentDofs"]
            dofs = Set{Int}()
            for (type,nodeList) in elms[v]
                union!(dofs,Set(nodeList))
            end
            elms[k] = [(:Poi1,[dof]) for dof in dofs]
        end
        delete!(config,"IndependentDofs")
    end

    nodes = Node(nodes...)
    for (name,cfg) in config
        Type = eval(Meta.parse(cfg["type"]))
        if Type <: ReproducingKernel
            𝗠 = Dict{Symbol,SymMat}()
            elements[name] = [Type([nodes[i] for i in s[2]],𝗠) for s in elms[cfg["𝓒"]["tag"]]]
        else
            elements[name] = [Type([nodes[i] for i in s[2]]) for s in elms[cfg["𝓒"]["tag"]]]
        end
        sp ≠ nothing ? sp(elements[name]) : nothing
        if haskey(cfg,"𝓖")
            QType = Meta.parse(cfg["𝓖"]["type"])
            if haskey(cfg["𝓖"],"tag")
                elms_𝓖 = [Element{s[1]}([nodes[i] for i in s[2]]) for s in elms[cfg["𝓖"]["tag"]]]
                elements[name] = elements[name]∩elms_𝓖
                set𝓖!(elms_𝓖,QType)
                set𝓖!(elements[name],elms_𝓖)
            else
                set𝓖!(elements[name],QType)
            end
            if haskey(cfg["𝓖"],"𝝭")
                ss = Meta.parse.(cfg["𝓖"]["𝝭"])
                Type<:ReproducingKernel ? set_memory_𝗠!(elements[name],ss...) : nothing
                set_memory_𝝭!(elements[name],ss...)
            end
        end
    end
    return elements, nodes
end

function importmsh(filename1::String,filename2::String,config::Dict{Any,Any})
    elms, nodes_ = importmsh(filename1)
    ~, nodes = importmsh(filename2)
    elements = Dict{String,Any}()
    cfg = config["RegularGrid"]
    sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z];n=cfg["n"],γ=cfg["γ"])
    delete!(config,"RegularGrid")
    nodes = Node(nodes...)
    for (name,cfg) in config
        Type = eval(Meta.parse(cfg["type"]))
        𝗠 = Dict{Symbol,SymMat}()
        QType = Meta.parse(cfg["𝓖"]["type"])
        elms_𝓖 = [Element{s[1]}([nodes_[i] for i in s[2]]) for s in elms[cfg["𝓖"]["tag"]]]
        set𝓖!(elms_𝓖,QType)
        elements[name] = [Type(sp(elm,nodes),𝗠) for elm in elms_𝓖]
        set𝓖!(elements[name],elms_𝓖)

        if haskey(cfg["𝓖"],"𝝭")
            ss = Meta.parse.(cfg["𝓖"]["𝝭"])
            Type<:ReproducingKernel ? set_memory_𝗠!(elements[name],ss...) : nothing
            set_memory_𝝭!(elements[name],ss...)
        end
    end
    return elements, nodes
end

function importmsh(filename::String,::Val{:test})
    elems,nodes = importmsh(filename)
    data = Dict([s=>(2,v) for (s,v) in nodes])
    dofs = getboundarydofs2D(elems["Ω"])
    elements = Dict{String,Any}()
    nodes = Node(nodes...)
    gnodes = GNode[]
    elements["∂Ω"] = Vector{Element{:Seg2}}(undef,length(dofs))
    for (dof,i) in dofs
        elements["∂Ω"][i] = Element{:Seg2}([nodes[j] for j in dof])
    end
    elements["Ω"] = DBelement{:Tri3}[]
    elements["Γ"] = DBelement{:Tri3}[]
    haskey(elems,"Γᵗ") ? elements["Γᵗ"] = DBelement{:Tri3}[] : nothing
    for (type,nodeList) in elems["Ω"]
        𝓒 = [GNode((dofs[Set(setdiff(nodeList,i))],i),data) for i in nodeList]
        union!(gnodes,𝓒)
        push!(elements["Ω"],DBelement{:Tri3}(𝓒))
        push!(elements["Γ"],DBelement{:Tri3}(𝓒))
        haskey(elems,"Γᵗ") ? push!(elements["Γᵗ"],DBelement{:Tri3}(𝓒)) : nothing
    end
    set𝓖!(elements["Ω"],:TriGI13)
    set𝓖_DB!(elements["Γ"],:SegGI2)
    if haskey(elems,"Γᵗ")
        elms_𝓖 = [Element{type}([nodes[i] for i in nodeList])     for (type,nodeList) in elems["Γᵗ"]]
        elements["Γᵗ"] = elements["Γᵗ"]∩elms_𝓖
        set𝓖!(elms_𝓖,:SegGI2)
        set𝓖!(elements["Γᵗ"],elms_𝓖)
    end

    elements["Γᵍ"] = DBelement{:Seg2}[]
    for (type,nodeList) in elems["Γᵍ"]
        𝐼 = dofs[Set(nodeList)]
        𝓒 = [GNode((0,i),data) for i in nodeList]
        push!(𝓒,GNode((𝐼,0),data))
        push!(elements["Γᵍ"],DBelement{:Seg2}(𝓒))
    end
    set𝓖!(elements["Γᵍ"],:SegGI2)

    set_memory_𝝭!(elements["Ω"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
    haskey(elems,"Γᵗ") ? set_memory_𝝭!(elements["Γᵗ"],:𝝭) : nothing
    set_memory_𝝭!(elements["Γᵍ"],:𝝭)
    set_memory_𝝭!(elements["Γ"],:𝝭)
    return elements, gnodes
end

function getboundarydofs2D(elements::Vector{Tuple{Symbol,Vector{Int}}})
    dofs = Dict{Set{Int},Int}()
    idBoundaries = (Tri3=((1,2),(2,3),(3,1)),)
    n = 0
    for (type,nodeList) in elements
        for bc in idBoundaries[type]
            dof = Set(nodeList[i] for i in bc)
            if ~haskey(dofs,dof)
                n += 1
                dofs[dof] = n
            end
        end
    end
    return dofs
end
