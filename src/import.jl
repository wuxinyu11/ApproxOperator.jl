## Data Pool
# ---------------- msh ---------------
struct MSHImport{𝒑,𝑠,𝜙}
    etype::Dict{Int,DataType}
    ntype::Dict{Int,DataType}
    mtype::Tuple{Val{𝒑},Val{𝑠},Val{𝜙}}
end

function MSHImport()
    etype = Dict(1=>:Seg2,2=>:Tri3,3=>:Quad,15=>:Poi1)
    ntype = Dict(1=>:Node,2=>:Node,3=>:Node,15=>:Node)
    return MSHImport(etype,ntype,(Val(nothing),Val(nothing),Val(nothing)))
end

function MSHImport(𝒑::Symbol,𝑠::Symbol,𝜙::Symbol)
    etype = Dict(1=>:SegN,2=>:TriN,3=>:QuadN,15=>:PoiN)
    ntype = Dict(1=>:Node,2=>:Node,3=>:Node,15=>:Node)
    return MSHImport(etype,ntype,(Val(𝒑),Val(𝑠),Val(𝜙)))
end

function (op::MSHImport)(::Val{:SNode})
    for key in keys(op.ntype)
        op.ntype[key] = :SNode
    end
end

function (op::MSHImport)(filename::String)
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
    return
end

function import_msh_4(fid::IO,op::Operator{:msh}) end

function import_msh_2(fid::IO,op::Operator{:msh})
    for line in eachline(fid)
        if line == "\$PhysicalNames"
            numPhysicalNames = parse(Int,readline(fid))
            datas = Dict{Symbol,Any}()
            push!(datas,:physicalnames=>Dict{Int,String}())
            push!(datas,:physicalfields=>Dict{Int,String}())
            for i in 1:numPhysicalNames
                line = readline(fid)
                d_,p_,n_ = split(line," ")
                dimension = parse(Int,d_)
                physicalTag = parse(Int,p_)
                name = strip(n_,'\"')
                datas[:physicalnames][physicalTag] = name
                push!(datas[:physicalfields], name=>Dict(:nₑ=>0,:nᵢ=>0,:data=>Dict{Symbol,Vector{Float64}}())
            end
            readline(fid)
        elseif line == "\$Nodes"
            line = readline(fid)
            nₚ = parse(Int,line)
            x=>zeros(nₚ)
            y=>zeros(nₚ)
            z=>zeros(nₚ)
            for i in 1:nₚ
                line = readline(fid)
                t_,x_,y_,z_ = split(line," ")
                tag = parse(Int,t_)
                x[i] = parse(Float64,x_)
                y[i] = parse(Float64,y_)
                z[i] = parse(Float64,z_)
            end
            push!(datas,:nodes=>Dict(:nₚ=>nₚ,:data=>Dict(:x=>x,:y=>y,:z=>z)))
            readline(fid)
        elseif line == "\$Elements"
            line = readline(fid)
            nₑ = parse(Int,line)
            push!(datas,:elements=>Dict{String,Vector{Approximator}}())
            for i in 1:nₑ
                line = readline(fid)
                elmN_,elmT_,numT_,phyT_,elmE_,l_... = split(line," ")
                elmNumber = parse(Int,elmN_)
                elmType = parse(Int,elmT_)
                numTag = parse(Int,numT_)
                phyTag = parse(Int,phyT_)
                elmEntary = parse(Int,elmE_)
                nodeList = parse.(Int,l_)
                name = datas[:physicalnames][phyTag]
                datas[:physicalfields][name][:nₑ] += 1
                elements = datas[:elements]
                etype = op.etype[elmType]
                ntype = op.ntype[elmType]
                physicalfield = datas[:physicalfields][name]
                if ntype == :SNode && haskey(physicalfield,:index)
                    haskey(physicalfield,:index)
                haskey(elements,name) ? push!(elements[name],eval(etype)(physicalfield,nodeList)) : push!(elements,name=>eval(etype)[eval(etype)(physicalfield,nodeList)])
            end
        end
    end
end
#
# function Node(op::Operator{:msh},name::String,id::Int)
#     data = op.parametricnodes[name]
#     return Node(id,data)
# end
# function SNode(op::Operator{:msh},name::String,id::Int)
#     data = op.parametricnodes[name]
#     𝝭 = op.shapefunctions
#     index = op.index
#     return SNode(id,data,index,𝝭)
# end
# function Poi1(op::Operator{:msh},name::String,id::Vector{Int},quadraturepoints::Tuple)
#     𝓒 = [Node(i,op.nodes) for i in id]
#     𝓖 = eval(op.ntype[15])[]
#     data = op.parametricnodes[name]
#     for ξ in quadraturepoints
#         haskey(data,:w) ? push!(data[:w],ξ[1]) : push!(data,:w=>[ξ[1]])
#         haskey(data,:ξ) ? push!(data[:ξ],ξ[2]) : push!(data,:ξ=>[ξ[2]])
#         n = length(data[:w])
#         push!(𝓖,eval(op.ntype[15])(op,name,n))
#     end
#     return Poi1(𝓒,𝓖)
# end
#
# function Seg2(op::Operator{:msh},name::String,id::Vector{Int},quadraturepoints::Tuple)
#     𝓒 = [Node(i,op.nodes) for i in id]
#     𝓖 = eval(op.ntype[1])[]
#     data = op.parametricnodes[name]
#     for ξ in quadraturepoints
#         haskey(data,:w) ? push!(data[:w],ξ[1]) : push!(data,:w=>[ξ[1]])
#         haskey(data,:ξ) ? push!(data[:ξ],ξ[2]) : push!(data,:ξ=>[ξ[2]])
#         n = length(data[:w])
#         push!(𝓖,eval(op.ntype[1])(op,name,n))
#     end
#     return Seg2(𝓒,𝓖)
# end
#
# function Tri3(op::Operator{:msh},name::String,id::Vector{Int},quadraturepoints::Tuple)
#     𝓒 = [Node(i,op.nodes) for i in id]
#     𝓖 = eval(op.ntype[2])[]
#     data = op.parametricnodes[name]
#     for ξ in quadraturepoints
#         haskey(data,:w) ? push!(data[:w],ξ[1]) : push!(data,:w=>[ξ[1]])
#         haskey(data,:ξ) ? push!(data[:ξ],ξ[2]) : push!(data,:ξ=>[ξ[2]])
#         haskey(data,:η) ? push!(data[:η],ξ[3]) : push!(data,:η=>[ξ[3]])
#         n = length(data[:w])
#         push!(𝓖,eval(op.ntype[2])(op,name,n))
#     end
#     return Tri3(𝓒,𝓖)
# end
# function PoiN(op::Operator{:msh},name::String,id::Vector{Int},quadraturepoints::Tuple)
#     sp = op.spatialpartition
#     indices = Set{Int}()
#     for i in id
#         union!(indices,sp(op.nodes[:x][i],op.nodes[:y][i],op.nodes[:z][i]))
#     end
#     id = union!(id,collect(indices))
#     𝓒 = [Node(i,op.nodes) for i in id]
#     𝓖 = eval(op.ntype[15])[]
#     data = op.parametricnodes[name]
#     for ξ in quadraturepoints
#         haskey(data,:w) ? push!(data[:w],ξ[1]) : push!(data,:w=>[ξ[1]])
#         haskey(data,:ξ) ? push!(data[:ξ],ξ[2]) : push!(data,:ξ=>[ξ[2]])
#         n = length(data[:w])
#         push!(𝓖,eval(op.ntype[15])(op,name,n))
#     end
#     𝗠 = op.𝗠
#     𝝭 = op.𝝭
#     𝒑 = op.basisfunction
#     𝑠 = op.kerneltype
#     𝜙 = op.kernelfunction
#     return PoiN(𝓒,𝓖,𝗠,𝝭,𝒑,𝑠,𝜙)
# end
# function SegN(op::Operator{:msh},name::String,id::Vector{Int},quadraturepoints::Tuple)
#     sp = op.spatialpartition
#     indices = Set{Int}()
#     for i in id
#         union!(indices,sp(op.nodes[:x][i],op.nodes[:y][i],op.nodes[:z][i]))
#     end
#     id = union!(id,collect(indices))
#     𝓒 = [Node(i,op.nodes) for i in id]
#     𝓖 = eval(op.ntype[1])[]
#     data = op.parametricnodes[name]
#     for ξ in quadraturepoints
#         haskey(data,:w) ? push!(data[:w],ξ[1]) : push!(data,:w=>[ξ[1]])
#         haskey(data,:ξ) ? push!(data[:ξ],ξ[2]) : push!(data,:ξ=>[ξ[2]])
#         if length(ξ) > 2
#             haskey(data,:wᵇ) ? push!(data[:wᵇ],ξ[3]) : push!(data,:wᵇ=>[ξ[3]])
#         end
#         n = length(data[:w])
#         push!(𝓖,eval(op.ntype[1])(op,name,n))
#     end
#     𝗠 = op.𝗠
#     𝝭 = op.𝝭
#     𝒑 = op.basisfunction
#     𝑠 = op.kerneltype
#     𝜙 = op.kernelfunction
#     return SegN(𝓒,𝓖,𝗠,𝝭,𝒑,𝑠,𝜙)
# end
