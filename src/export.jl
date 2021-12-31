
function export_shape_functions(ap::AbstractSeg;inte::Int=100)
    # coordinates = Dict{String,Vector{Float64}}()
    # push!(coordinates,"x"=>zeros(inte+1))
    # push!(coordinates,"y"=>zeros(inte+1))
    # push!(coordinates,"z"=>zeros(inte+1))
    # datas = Dict{Int,Vector{Float64}}()
    x = zeros(inte+1)
    max_id = maximum(ap.id)
    # push!(datas,i=>zeros(inte+1) for i in 1:max_id)
    datas = zeros(inte+1,max_id)
    for j in 1:inte+1
        ξⱼ = -1+2/inte*(j-1)
        xⱼ = get_coordinates(ap,ξⱼ)
        # coordinates["x"][j] = xᵢ[1]
        # coordinates["y"][j] = xᵢ[2]
        # coordinates["z"][j] = xᵢ[3]
        x[j] = xⱼ[1]
        N = get_shape_functions(ap,ξⱼ,:∂1)
        for i in 1:get_number_of_indices
            I = get_global_indice(ap,i)
            # datas[I][j] = N[i]
            datas[j,I] = N[i]
        end
    end
    # return coordinates, datas
    return x, datas
end

export_shape_functions(aps::Vector{AbstractSeg},gs::Symbol...;inte::Int=100)= export_shape_functions(aps,inte,Val.(gs)...)

function export_shape_functions(aps::Vector{AbstractSeg},inte::Int,::Val{:∂1})
    nₑ = length(aps)
    x = zeros(nₑ*(inte+1))
    max_id = 0
    for ap in aps
        max_temp = maximum(ap.id)
        max_id = max_id > max_temp ? max_id : max_temp
    end
    N = zeros(nₑ*(inte+1),max_id)
    for n in 1:nₑ
        ap = aps[n]
        for j in 1:inte+1
            ξⱼ = -1+2/inte*(j-1)
            xⱼ = get_coordinates(ap,ξⱼ)
            x[(inte+1)*(n-1)+j] = xⱼ[1]
            N_ = get_shape_functions(ap,ξⱼ,:∂1)
            for i in 1:get_number_of_indices(ap)
                I = get_global_indice(ap,i)
                N[(inte+1)*(n-1)+j,I] = N_[i]
            end
        end
    end
    return x, N
end

function export_shape_functions(aps::Vector{AbstractSeg},inte::Int,::Val{:∂1},::Val{:∂x})
    nₑ = length(aps)
    x = zeros(nₑ*(inte+1))
    max_id = 0
    for ap in aps
        max_temp = maximum(ap.id)
        max_id = max_id > max_temp ? max_id : max_temp
    end
    N = zeros(nₑ*(inte+1),max_id)
    ∂N∂x = zeros(nₑ*(inte+1),max_id)
    for n in 1:nₑ
        ap = aps[n]
        for j in 1:inte+1
            ξⱼ = -1+2/inte*(j-1)
            xⱼ = get_coordinates(ap,ξⱼ)
            x[(inte+1)*(n-1)+j] = xⱼ[1]
            N_, ∂N∂x_ = get_shape_functions(ap,ξⱼ,:∂1,:∂x)
            for i in 1:get_number_of_indices(ap)
                I = get_global_indice(ap,i)
                N[(inte+1)*(n-1)+j,I] = N_[i]
                ∂N∂x[(inte+1)*(n-1)+j,I] = ∂N∂x_[i]
            end
        end
    end
    return x, N, ∂N∂x
end

function export_shape_functions(aps::Vector{AbstractSeg},inte::Int,::Val{:∂1},::Val{:∂x},::Val{:∂y},::Val{:∂z})
    nₑ = length(aps)
    x = zeros(nₑ*(inte+1))
    max_id = 0
    for ap in aps
        max_temp = maximum(ap.id)
        max_id = max_id > max_temp ? max_id : max_temp
    end
    N = zeros(nₑ*(inte+1),max_id)
    ∂N∂x = zeros(nₑ*(inte+1),max_id)
    ∂N∂y = zeros(nₑ*(inte+1),max_id)
    ∂N∂z = zeros(nₑ*(inte+1),max_id)
    for n in 1:nₑ
        ap = aps[n]
        for j in 1:inte+1
            ξⱼ = -1+2/inte*(j-1)
            xⱼ = get_coordinates(ap,ξⱼ)
            x[(inte+1)*(n-1)+j] = xⱼ[1]
            N_, ∂N∂x_, ∂N∂y_, ∂N∂z_ = get_shape_functions(ap,ξⱼ,:∂1,:∂x,:∂y,:∂z)
            for i in 1:get_number_of_indices(ap)
                I = get_global_indice(ap,i)
                N[(inte+1)*(n-1)+j,I] = N_[i]
                ∂N∂x[(inte+1)*(n-1)+j,I] = ∂N∂x_[i]
                ∂N∂y[(inte+1)*(n-1)+j,I] = ∂N∂y_[i]
                ∂N∂z[(inte+1)*(n-1)+j,I] = ∂N∂z_[i]
            end
        end
    end
    return x, N, ∂N∂x, ∂N∂y, ∂N∂z
end

## VTK
const cell_type = Dict{Any,Int}(AbstractPoi=>1,AbstractSeg=>3,AbstractTri=>5,AbstractQuad=>9)

mutable struct VTKExport
    filename::String
    topic::String
    celltype::String
    pointdata::Union{Dict{String,Symbol},Nothing}
    celldata::Union{Dict{String,Symbol},Nothing}
    parameters::Union{Dict{String,Float64},Nothing}
end
VTKExport(;filename::String="default.vtk",topic::String="vtk file",celltype::String="UNSTRUCTURED_GRID",pointdata::Union{Dict{String,Symbol},Nothing}=nothing,celldata::Union{Dict{String,Symbol},Nothing}=nothing,parameters::Union{Dict{String,Float64},Nothing}=nothing) = VTKExport(filename,topic,celltype,pointdata,celldata,parameters)

function (op::VTKExport)(aps::Vector{Approximator},d::AbstractVector{Float64})
    fid = open(op.filename,"w")
    nₑ = length(aps)
    nₚ = length(aps[1].nodes)
    write(fid,"# vtk DataFile Version 2.0\n")
    write(fid,op.topic*"\n")
    write(fid,"ASCII\n")
    write(fid,"DATASET "*op.celltype*"\n")
    # POINTS
    write(fid,"POINTS $nₚ float\n")
    for node in aps[1].nodes
        x = node.x
        y = node.y
        z = node.z
        write(fid,"$x $y $z\n")
    end

    # CELLS
    op(fid,aps,Val(Meta.parse(op.celltype)))

    # POINT_DATA
    if op.pointdata ≠ nothing
        write(fid,"POINT_DATA $nₚ\n")
        for (name,datatype) in op.pointdata
            op(fid,aps,d,name,Val(datatype))
        end
    end
    # CELL_DATA
    if op.celldata ≠ nothing
        write(fid,"CELL_DATA $nₑ\n")
        for (name,datatype) in op.celldata
            op(fid,aps,d,name,Val(datatype))
        end
    end
    close(fid)
end

function (op::VTKExport)(fid::IO,aps::Vector{Approximator},::Val{:UNSTRUCTURED_GRID})
    nₑ = length(aps)
    n = nₑ
    for ap in aps
        n += get_number_of_indices(ap)
    end
    write(fid,"CELLS $nₑ $n\n")
    for ap in aps
        nᵢ = get_number_of_indices(ap)
        str = "$nᵢ"
        for I in ap.id .- 1
            str *= " $I"
        end
        str *= "\n"
        write(fid,str)
    end
    write(fid,"CELL_TYPES $nₑ\n")
    for ap in aps
        ntag = cell_type[supertype(typeof(ap))]
        write(fid,"$ntag\n")
    end
end

function (::VTKExport)(fid::IO,aps::Vector{Approximator},d::AbstractVector{Float64},name::String,::Val{:Scale})
    write(fid,"SCALARS "*name*" float 1\n")
    write(fid,"LOOKUP_TABLE "*name*"\n")
    for u in d
        write(fid,"$u\n")
    end
end

function (::VTKExport)(fid::IO,aps::Vector{Approximator},d::AbstractVector{Float64},name::String,::Val{:u_2D})
    nₚ = length(aps[1].nodes)
    write(fid,"VECTORS "*name*" float\n")
    for i in 1:nₚ
        u = d[2*i-1]
        v = d[2*i]
        w = 0.0
        write(fid,"$u $v $w\n")
    end
end

function (::VTKExport)(fid::IO,aps::Vector{Approximator},d::AbstractVector{Float64},name::String,::Val{:u})
    nₚ = length(aps[1].nodes)
    write(fid,"VECTORS "*name*" float\n")
    for i in 1:nₚ
        u = d[3*i-2]
        v = d[3*i-1]
        w = d[3*i]
        write(fid,"$u $v $w\n")
    end
end

function (op::VTKExport)(fid::IO,aps::Vector{Approximator},d::AbstractVector{Float64},name::String,::Val{:σ_PlaneStress})
    write(fid,"TENSORS "*name*" float\n")
    E = op.parameters["E"]
    ν = op.parameters["ν"]
    Cᵢᵢᵢᵢ = E/(1-ν^2)
    Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
    Cᵢⱼᵢⱼ = E/2/(1+ν)
    for ap in aps
        B₁,B₂ = get_shape_functions(ap,SVector{3,Float64}(0.,0.,0.),Val(:∂x),Val(:∂y))
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            ε₁₁ += B₁[i]*d[2*I-1]
            ε₂₂ += B₂[i]*d[2*I]
            ε₁₂ += B₁[i]*d[2*I] + B₂[i]*d[2*I-1]
        end
        σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        σ₂₂ = Cᵢᵢᵢᵢ*ε₂₂ + Cᵢᵢⱼⱼ*ε₁₁
        σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
        write(fid,"$σ₁₁ $σ₁₂ 0.0\n")
        write(fid,"$σ₁₂ $σ₂₂ 0.0\n")
        write(fid,"0.0 0.0 0.0\n")
        write(fid,"\n")
    end
end

function (op::VTKExport)(fid::IO,aps::Vector{Approximator},d::AbstractVector{Float64},name::String,::Val{:σ_PlaneStrain})
    write(fid,"TENSORS "*name*" float\n")
    E = op.parameters["E"]
    ν = op.parameters["ν"]
    Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-ν^2)
    Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-ν^2)
    Cᵢⱼᵢⱼ = E/2/(1+ν)
    for ap in aps
        B₁,B₂ = get_shape_functions(ap,SVector{3,Float64}(0.,0.,0.),Val(:∂x),Val(:∂y))
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            ε₁₁ += B₁[i]*d[2*I-1]
            ε₂₂ += B₂[i]*d[2*I]
            ε₁₂ += B₁[i]*d[2*I] + B₂[i]*d[2*I-1]
        end
        σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        σ₂₂ = Cᵢᵢᵢᵢ*ε₂₂ + Cᵢᵢⱼⱼ*ε₁₁
        σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
        σ₃₃ = Cᵢᵢⱼⱼ*(ε₁₁+ε₂₂)
        write(fid,"$σ₁₁ $σ₁₂ 0.0\n")
        write(fid,"$σ₁₂ $σ₂₂ 0.0\n")
        write(fid,"0.0 0.0 $σ₃₃\n")
        write(fid,"\n")
    end
end

function (op::VTKExport)(fid::IO,aps::Vector{Approximator},d::AbstractVector{Float64},name::String,::Val{:σ})
    write(fid,"TENSORS "*name*" float\n")
    E = op.parameters["E"]
    ν = op.parameters["ν"]
    Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-ν^2)
    Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-ν^2)
    Cᵢⱼᵢⱼ = E/2/(1+ν)
    for ap in aps
        B₁,B₂,B₃ = get_shape_functions(ap,SVector{3,Float64}(0.,0.,0.),Val(:∂x),Val(:∂y),Val(:∂z))
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₃₃ = 0.
        ε₁₂ = 0.
        ε₁₃ = 0.
        ε₂₃ = 0.
        for i in 1:get_number_of_indices(ap)
            I = get_global_indice(ap,i)
            ε₁₁ += B₁[i]*d[3*I-2]
            ε₂₂ += B₂[i]*d[3*I-1]
            ε₃₃ += B₃[i]*d[3*I]
            ε₁₂ += B₁[i]*d[3*I-1] + B₂[i]*d[3*I-2]
            ε₁₃ += B₁[i]*d[3*I] + B₃[i]*d[3*I-2]
            ε₂₃ += B₂[i]*d[3*I] + B₃[i]*d[3*I-1]
        end
        σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂ + Cᵢᵢⱼⱼ*ε₃₃
        σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢᵢᵢ*ε₂₂ + Cᵢᵢⱼⱼ*ε₃₃
        σ₃₃ = Cᵢᵢⱼⱼ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂ + Cᵢᵢᵢᵢ*ε₃₃
        σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
        σ₁₃ = Cᵢⱼᵢⱼ*ε₁₃
        σ₂₃ = Cᵢⱼᵢⱼ*ε₂₃
        write(fid,"$σ₁₁ $σ₁₂ $σ₁₃\n")
        write(fid,"$σ₁₂ $σ₂₂ $σ₂₃\n")
        write(fid,"$σ₁₃ $σ₂₃ $σ₃₃\n")
        write(fid,"\n")
    end
end
## Phase field export
function (op::VTKExport)(aps::Vector{Approximator},d::AbstractVector{Float64},dᵥ::AbstractVector{Float64},::Val{:PFM_PlaneStress})
    fid = open(op.filename,"w")
    nₑ = length(aps)
    nₚ = length(aps[1].nodes)
    write(fid,"# vtk DataFile Version 2.0\n")
    write(fid,op.topic*"\n")
    write(fid,"ASCII\n")
    write(fid,"DATASET "*op.celltype*"\n")
    # POINTS
    write(fid,"POINTS $nₚ float\n")
    for node in aps[1].nodes
        x = node.x[1]
        y = node.x[2]
        z = node.x[3]
        write(fid,"$x $y $z\n")
    end

    # CELLS
    op(fid,aps,Val(Meta.parse(op.celltype)))

    # POINT_DATA
    write(fid,"POINT_DATA $nₚ\n")
    op(fid,aps,d,"disp",Val(:u_2D))
    op(fid,aps,dᵥ,"damage",Val(:Scale))

    # CELL_DATA
    write(fid,"CELL_DATA $nₑ\n")
    op(fid,aps,d,"stress",Val(:σ_PlaneStress))
    close(fid)
end
