## Data Pool
# ---------------- msh ---------------
function Operator(type::Val{:msh})
    data = Dict{Symbol,Any}()
    etype = Dict(1=>:Seg2,2=>:Tri3,3=>:Quad,15=>:Poi1)
    ntype = Dict(1=>:Node,2=>:Node,3=>:Node,15=>:Node)
    qtype = Dict(1=>:SegGI2,2=>:TriGI3,3=>:QuadGI2,15=>:PoiGI1)
    push!(data,:etype=>etype)
    push!(data,:ntype=>ntype)
    push!(data,:qtype=>qtype)
    return Operator(type,data)
end

function (op::Operator{:msh})(filename::String)
    fid = open(filename,"r")
    readline(fid)
    line = readline(fid)
    v_,f_,d_ = split(line," ")
    version = parse(Float64,v_)
    filetype = parse(Int,f_)
    datasize = parse(Int,d_)
    readline(fid)
    if version == 4.1
        import_msh_4(fid,op)
    elseif version == 2.2
        import_msh_2(fid,op)
    else
        println("Version does not match!")
    end
    return op.elements
end

function import_msh_4(fid::IO,op::Operator{:msh}) end

function import_msh_2(fid::IO,op::Operator{:msh})
    for line in eachline(fid)
        if line == "\$PhysicalNames"
            numPhysicalNames = parse(Int,readline(fid))
            push!(op,:physicalnames=>Dict{Int,String}())
            push!(op,:parametricnodes=>Dict{String,Dict{Symbol,Vector{Float64}}}())
            push!(op,:elements=>Dict{String,Any}())
            push!(op,:nₑ=>Dict{String,Int}())
            push!(op,:nᵢ=>Dict{String,Int}())
            for i in 1:numPhysicalNames
                line = readline(fid)
                d_,p_,n_ = split(line," ")
                dimension = parse(Int,d_)
                physicalTag = parse(Int,p_)
                name = strip(n_,'\"')

                op.physicalnames[physicalTag] = name
                op.nₑ[name] = 0
                op.nᵢ[name] = 0
                op.parametricnodes[name] = Dict{Symbol,Vector{Float64}}()
            end
            readline(fid)
        elseif line == "\$Nodes"
            line = readline(fid)
            nₚ = parse(Int,line)
            push!(op,:nodes=>Dict(:x=>zeros(nₚ),:y=>zeros(nₚ),:z=>zeros(nₚ)))
            push!(op,:nₚ=>nₚ)
            for i in 1:nₚ
                line = readline(fid)
                t_,x_,y_,z_ = split(line," ")
                tag = parse(Int,t_)
                op.nodes[:x][i] = parse(Float64,x_)
                op.nodes[:y][i] = parse(Float64,y_)
                op.nodes[:z][i] = parse(Float64,z_)
            end
            if haskey(op.data,:spatialpartition)
                if op.spatialpartition == :RegularGrid
                    sp = RegularGrid(op.nodes[:x],op.nodes[:y],op.nodes[:z],n=op.nᵣ,γ=op.γᵣ)
                    op.spatialpartition = sp
                    nₘ = 0
                    for c in sp.cells
                        nₘ = max(length(c),nₘ)
                    end
                    push!(op,:nₘ=>nₘ*op.nᵣ*3)
                end
                n = length(get𝒑(Val(op.basisfunction),(0.0,0.0,0.0)))
                push!(op,:𝗠=>Dict{Symbol,SymMat}())
                push!(op,:𝝭=>Dict{Symbol,Vector{Float64}}())
                for s in op.stype
                    push!(op.data[:𝗠],s=>SymMat(n))
                    push!(op.data[:𝝭],s=>zeros(op.nₘ))
                end
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
                name = op.physicalnames[phyTag]
                quadraturepoints = QuadratureRule[op.qtype[elmType]]
                op.nₑ[name] += 1
                op.nᵢ[name] += length(quadraturepoints)
                haskey(op.elements,name) ? push!(op.elements[name],eval(op.etype[elmType])(op,name,nodeList,quadraturepoints)) : push!(op.elements,name=>eval(op.etype[elmType])[eval(op.etype[elmType])(op,name,nodeList,quadraturepoints)])
            end
        end
    end
end

function Node(op::Operator{:msh},name::String,id::Int)
    data = op.parametricnodes[name]
    return Node(id,data)
end
function SNode(op::Operator{:msh},name::String,id::Int)
    data = op.parametricnodes[name]
    𝝭 = op.shapefunctions
    index = op.index
    return SNode(id,data,index,𝝭)
end
function Poi1(op::Operator{:msh},name::String,id::Vector{Int},quadraturepoints::Tuple)
    𝓒 = [Node(i,op.nodes) for i in id]
    𝓖 = eval(op.ntype[15])[]
    data = op.parametricnodes[name]
    for ξ in quadraturepoints
        haskey(data,:w) ? push!(data[:w],ξ[1]) : push!(data,:w=>[ξ[1]])
        haskey(data,:ξ) ? push!(data[:ξ],ξ[2]) : push!(data,:ξ=>[ξ[2]])
        n = length(data[:w])
        push!(𝓖,eval(op.ntype[15])(op,name,n))
    end
    return Poi1(𝓒,𝓖)
end

function Seg2(op::Operator{:msh},name::String,id::Vector{Int},quadraturepoints::Tuple)
    𝓒 = [Node(i,op.nodes) for i in id]
    𝓖 = eval(op.ntype[1])[]
    data = op.parametricnodes[name]
    for ξ in quadraturepoints
        haskey(data,:w) ? push!(data[:w],ξ[1]) : push!(data,:w=>[ξ[1]])
        haskey(data,:ξ) ? push!(data[:ξ],ξ[2]) : push!(data,:ξ=>[ξ[2]])
        n = length(data[:w])
        push!(𝓖,eval(op.ntype[1])(op,name,n))
    end
    return Seg2(𝓒,𝓖)
end

function Tri3(op::Operator{:msh},name::String,id::Vector{Int},quadraturepoints::Tuple)
    𝓒 = [Node(i,op.nodes) for i in id]
    𝓖 = eval(op.ntype[2])[]
    data = op.parametricnodes[name]
    for ξ in quadraturepoints
        haskey(data,:w) ? push!(data[:w],ξ[1]) : push!(data,:w=>[ξ[1]])
        haskey(data,:ξ) ? push!(data[:ξ],ξ[2]) : push!(data,:ξ=>[ξ[2]])
        haskey(data,:η) ? push!(data[:η],ξ[3]) : push!(data,:η=>[ξ[3]])
        n = length(data[:w])
        push!(𝓖,eval(op.ntype[2])(op,name,n))
    end
    return Tri3(𝓒,𝓖)
end
function PoiN(op::Operator{:msh},name::String,id::Vector{Int},quadraturepoints::Tuple)
    sp = op.spatialpartition
    indices = Set{Int}()
    for i in id
        union!(indices,sp(op.nodes[:x][i],op.nodes[:y][i],op.nodes[:z][i]))
    end
    id = union!(id,collect(indices))
    𝓒 = [Node(i,op.nodes) for i in id]
    𝓖 = eval(op.ntype[15])[]
    data = op.parametricnodes[name]
    for ξ in quadraturepoints
        haskey(data,:w) ? push!(data[:w],ξ[1]) : push!(data,:w=>[ξ[1]])
        haskey(data,:ξ) ? push!(data[:ξ],ξ[2]) : push!(data,:ξ=>[ξ[2]])
        n = length(data[:w])
        push!(𝓖,eval(op.ntype[15])(op,name,n))
    end
    𝗠 = op.𝗠
    𝝭 = op.𝝭
    𝒑 = op.basisfunction
    𝑠 = op.kerneltype
    𝜙 = op.kernelfunction
    return PoiN(𝓒,𝓖,𝗠,𝝭,𝒑,𝑠,𝜙)
end
function SegN(op::Operator{:msh},name::String,id::Vector{Int},quadraturepoints::Tuple)
    sp = op.spatialpartition
    indices = Set{Int}()
    for i in id
        union!(indices,sp(op.nodes[:x][i],op.nodes[:y][i],op.nodes[:z][i]))
    end
    id = union!(id,collect(indices))
    𝓒 = [Node(i,op.nodes) for i in id]
    𝓖 = eval(op.ntype[1])[]
    data = op.parametricnodes[name]
    for ξ in quadraturepoints
        haskey(data,:w) ? push!(data[:w],ξ[1]) : push!(data,:w=>[ξ[1]])
        haskey(data,:ξ) ? push!(data[:ξ],ξ[2]) : push!(data,:ξ=>[ξ[2]])
        if length(ξ) > 2
            haskey(data,:wᵇ) ? push!(data[:wᵇ],ξ[3]) : push!(data,:wᵇ=>[ξ[3]])
        end
        n = length(data[:w])
        push!(𝓖,eval(op.ntype[1])(op,name,n))
    end
    𝗠 = op.𝗠
    𝝭 = op.𝝭
    𝒑 = op.basisfunction
    𝑠 = op.kerneltype
    𝜙 = op.kernelfunction
    return SegN(𝓒,𝓖,𝗠,𝝭,𝒑,𝑠,𝜙)
end
