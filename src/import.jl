

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
    n₁ = get𝑛𝒑₁(ap)
    n₂ = get𝑛𝒑₂(ap)
    empty!(ap.𝗠)
    for s in ss
        if s == :∇̃
            ap.𝗠[s] = SymMat(n₁)
        elseif s ∈ (:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
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
    elements = Dict{String,Set{Tuple{Symbol,Vector{Int}}}}()
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
                elements[name] = Set{Tuple{Symbol,Vector{Int}}}()
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
                elmN_,elmT_,numT_,phyT_,elmE_,l_... = split(line," ")
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
            nₑ = length(elements[name])
            nᵢ = length(quadraturerule(QType)[:w])
            push!(getfield(elements[name][1].𝓖[1],:data),:x=>(2,zeros(nₑ*nᵢ)),:y=>(2,zeros(nₑ*nᵢ)),:z=>(2,zeros(nₑ*nᵢ)),:𝑤=>(2,zeros(nₑ*nᵢ)))
            setgeometry!.(elements[name])
            if haskey(cfg["𝓖"],"𝝭")
                ss = Meta.parse.(cfg["𝓖"]["𝝭"])
                Type<:ReproducingKernel ? set_memory_𝗠!(elements[name],ss...) : nothing
                set_memory_𝝭!(elements[name],ss...)
            end
        end
    end
    return elements
end
