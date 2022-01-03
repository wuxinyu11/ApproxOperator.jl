const etype = Dict{Int,Any}(1=>Seg2, 2=>Tri3, 3=>Quad, 15=>Poi1)

function import_msh(filename::String)
    fid = open(filename,"r")
    readline(fid)
    line = readline(fid)
    v_,f_,d_ = split(line," ")
    version = parse(Float64,v_)
    filetype = parse(Int,f_)
    datasize = parse(Int,d_)
    readline(fid)
    if version == 4.1
        aps = import_msh_4(fid)
    elseif version == 2.2
        aps = import_msh_2(fid)
    else
        println("Version does not match!")
    end
    return aps
end

function import_msh_4(fid::IO)
    aps = Dict{String,Vector{Approximator}}()
    phy = Dict{Int,String}()
    for line in eachline(fid)
        if line == "\$PhysicalNames"
            numPhysicalNames = parse(Int,readline(fid))
            for i in 1:numPhysicalNames
                line = readline(fid)
                d_,p_,name = split(line," ")
                dimension = parse(Int,d_)
                physicalTag = parse(Int,p_)

                phy[physicalTag] = name
            end
            readline(fid)
        elseif line == "\$Nodes"
            line = readline(fid)
            numE_,numN_,minN_,maxN_ = split(line," ")
            numEntityBlocks = parse(Int,numE_)
            numNodes = parse(Int,numN_)
            minNodeTag = parse(Int,minN_)
            maxNodeTag = parse(Int,maxN_)
            x = Vector{PhysicalNode}()
            for i in 1:numEntityBlocks
                line = readline(fid)
                entityD_,entityE_,para_,numN_ = split(line," ")
                # entityD_
            end
        end
    end
    return aps
end

function import_msh_2(fid::IO)
    aps = Dict{String,Vector{Approximator}}()
    nodes = Vector{PhysicalNode}()
    phy = Dict{Int,String}()
    for line in eachline(fid)
        if line == "\$PhysicalNames"
            numPhysicalNames = parse(Int,readline(fid))
            for i in 1:numPhysicalNames
                line = readline(fid)
                d_,p_,n_ = split(line," ")
                dimension = parse(Int,d_)
                physicalTag = parse(Int,p_)
                name = strip(n_,'\"')

                phy[physicalTag] = name
                aps[name] = Vector{Approximator}()
            end
            readline(fid)
        elseif line == "\$Nodes"
            line = readline(fid)
            nₚ = parse(Int,line)
            for i in 1:nₚ
                line = readline(fid)
                t_,x_,y_,z_ = split(line," ")
                tag = parse(Int,t_)
                x = parse(Float64,x_)
                y = parse(Float64,y_)
                z = parse(Float64,z_)
                push!(nodes,Node(x,y,z))
            end
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
                push!(aps[phy[phyTag]],etype[elmType](nodes,nodeList))
            end
        end
    end
    return aps
end
