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
                haskey(elements,name) ? push!(elements[name],Element{type}(nodes,nodeList...)) : elements[name]=Element{type}[Element{type}(nodes,nodeList...)]
            end
            return elements, nodes
        end
    end
end

function importmsh(filename::String,config::Dict{Any,Any})
    elms, nodes = importmsh(filename)
    elements = Dict{String,Any}()
    for (name,cfg) in config
        Type = eval(Meta.parse(cfg["𝓒"]["type"]))
        if haskey(cfg,"𝓖")
            QType = Meta.parse(cfg["𝓖"]["type"])
            if haskey(cfg["𝓖"],"tag")
                elms_𝓖 = elms[cfg["𝓖"]["tag"]]
                elms_𝓒 = elms[cfg["𝓒"]["tag"]]∩elms_𝓖
                set𝓖!(elms_𝓖,QType)
                elems = Type(elms_𝓒,elms_𝓖)
                elements[name] = elems
            else
                elems = Type(elms[cfg["𝓒"]["tag"]])
                set𝓖!(elems,QType)
                elements[name] = elems
            end
            if haskey(cfg["𝓖"],"𝝭")
                ss = cfg["𝓖"]["𝝭"]
                ss = [Meta.parse(s) for s in ss]
                set_storage_𝝭!(elements[name],ss...)
            end
        else
            elements[name] = Type(elms[cfg["𝓒"]["tag"]])
        end
        if haskey(cfg,"𝗠")
            ss = cfg["𝗠"]
            ss = [Meta.parse(s) for s in ss]
            set_memory_𝗠!(elements[name],ss...)
        end
        if haskey(cfg,"𝝭")
            ss = cfg["𝝭"]
            ss = [Meta.parse(s) for s in ss]
            set_memory_𝝭!(elements[name],ss...)
        end
    end
    return elements,nodes
end
