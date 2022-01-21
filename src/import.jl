## Data Pool
# ---------------- msh ---------------
function Operator(t::Val{:msh})
    etype = Dict(1=>:Seg2,2=>:Tri3,3=>:Quad,15=>:Poi1)
    ntype = Dict(1=>:Node,2=>:Node,3=>:Node,15=>:Node)
    return Operator(t,Dict{Symbol,Any}(:etype=>etype,:ntype=>ntype))
end

function (op::Operator{:msh})(𝒑::Symbol,𝑠::Symbol,𝜙::Symbol)
    op.etype = Dict(1=>:SegN,2=>:TriN,3=>:QuadN,15=>:PoiN)
    push!(op,:𝒑=>𝒑,:𝑠=>𝑠,:𝜙=>𝜙)
    push!(op,:𝗠=>Dict{Symbol,SymMat}())
    push!(op,:𝝭=>Dict{Symbol,Vector{Float64}}())
end
(op::Operator{:msh})(s::Symbol) = op(Val(s))
function (op::Operator{:msh})(::Val{:SNode})
    for key in keys(op.ntype)
        op.ntype[key] = :SNode
    end
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
    return op.elements, op.nodes
end

function import_msh_4(fid::IO,op::Operator{:msh}) end

function import_msh_2(fid::IO,op::Operator{:msh})
    for line in eachline(fid)
        if line == "\$PhysicalNames"
            numPhysicalNames = parse(Int,readline(fid))
            push!(op,:physicalnames=>Dict{Int,String}())
            push!(op,:nₑ=>Dict{String,Int}())
            for i in 1:numPhysicalNames
                line = readline(fid)
                d_,p_,n_ = split(line," ")
                dimension = parse(Int,d_)
                physicalTag = parse(Int,p_)
                name = strip(n_,'\"')
                op.physicalnames[physicalTag] = name
                op.nₑ[name] = 0
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
            push!(op,:nₚ=>nₚ)
            push!(op,:nodes=>Dict(:x=>x,:y=>y,:z=>z))
            readline(fid)
        elseif line == "\$Elements"
            line = readline(fid)
            nₑ = parse(Int,line)
            # push!(op,:elements=>Dict{String,Vector{Approximator}}())
            push!(op,:elements=>Dict{String,Any}())
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
                op.nₑ[name] += 1
                nodes = op.nodes
                elements = op.elements
                etype = eval(op.etype[elmType])
                ntype = eval(op.ntype[elmType])
                if haskey(op.data,:𝒑)
                    𝒑 = op.𝒑
                    𝑠 = op.𝑠
                    𝜙 = op.𝜙
                    𝗠 = op.𝗠
                    𝝭 = op.𝝭
                    haskey(elements,name) ? push!(elements[name],etype{ntype,𝒑,𝑠,𝜙}(nodeList...,nodes,𝗠,𝝭)) : push!(elements,name=>etype{ntype,𝒑,𝑠,𝜙}[etype{ntype,𝒑,𝑠,𝜙}(nodeList...,nodes,𝗠,𝝭)])
                else
                    haskey(elements,name) ? push!(elements[name],eval(etype)(nodeList...,nodes)) : push!(elements,name=>eval(etype)[eval(etype)(nodeList...,nodes)])
                end
            end
        end
    end
end
