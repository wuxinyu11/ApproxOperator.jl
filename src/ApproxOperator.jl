module ApproxOperator

import Base: +, -, *, getindex, setindex!, fill!
import InteractiveUtils: subtypes
import InteractiveUtils: @code_warntype
# import LinearAlgebra: norm

using StaticArrays
using SparseArrays
using LinearAlgebra
# import LinearAlgebra: inv
# import PyPlot: plot
# abstract type Approximator end
abstract type AbstractNode end
abstract type Approximator end
abstract type BasisFunction end
abstract type KernelFunction end
abstract type ShapeFunction end
abstract type SpatialPartition end
abstract type AbstractFunction end
abstract type AbstractPool end
abstract type AbstractPoi <: Approximator end
abstract type AbstractSeg <: Approximator end
abstract type AbstractTri <: Approximator end
abstract type AbstractQuad <: Approximator end
abstract type Operator end

include("linearalgebra.jl")
include("approximation.jl")
# export Poi1, PoiM, Seg2, SegM, Tri3, TriM
export Linear1D, Linear2D, Quadratic1D, Quadratic2D, Cubic1D, Cubic2D, TensorProductKernel
export RegularGrid
export set_integration_rule!
include("operation.jl")
export EBCDOFS, Potential_Ω, Potential_Γᵍ, Potential_Γᵗ, Potential_Γᵍ_Nitsche, Potential_Γᵍ_penalty, Potential_Γᵍ_Lagrange_multiplier, PlaneStress_Ω, PlaneStrain_Ωᵛ, PlaneStrain_Ωᵈ, PlaneStress_Γᵗ,PlaneStress_Γᵍ_penalty, L₂Error_scale, H₁Error_scale, L₂Error_tensor, L₂Error_2nd_order_tensor, HₑError_PlaneStress
include("import.jl")
export import_msh
include("export.jl")
export export_shape_functions, VTKExport


# debug
export BasisFunction, generate_1d_linear_basis_function
export get_shape_functions, cal_moment_matrix!, Node, get_coordinates, get_basis_function, Linear1D, CubicSplineKernel
export AbstractSeg, Approximator
include("efficiency.jl")
export efficiency
end
