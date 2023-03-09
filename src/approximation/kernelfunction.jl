
function get𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₂
    rz = abs(Δx[3])/x.s₃
    wx = get𝜙ᵣ(ap,rx)
    wy = get𝜙ᵣ(ap,ry)
    wz = get𝜙ᵣ(ap,rz)
    return wx*wy*wz
end

function get∇𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₂
    rz = abs(Δx[3])/x.s₃
    ∂rx = sign(Δx[1])/x.s₁
    ∂ry = sign(Δx[2])/x.s₂
    ∂rz = sign(Δx[3])/x.s₃
    wx = get𝜙ᵣ(ap,rx)
    wy = get𝜙ᵣ(ap,ry)
    wz = get𝜙ᵣ(ap,rz)
    ∂wx = get∂𝜙∂r(ap,rx)*∂rx
    ∂wy = get∂𝜙∂r(ap,ry)*∂ry
    ∂wz = get∂𝜙∂r(ap,rz)*∂rz
    return wx*wy*wz, ∂wx*wy*wz, wx*∂wy*wz, wx*wy*∂wz
end

function get∇₂𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₂
    ∂rx = sign(Δx[1])/x.s₁
    ∂ry = sign(Δx[2])/x.s₂
    wx = get𝜙ᵣ(ap,rx)
    wy = get𝜙ᵣ(ap,ry)
    ∂wx = get∂𝜙∂r(ap,rx)*∂rx
    ∂wy = get∂𝜙∂r(ap,ry)*∂ry
    return wx*wy, ∂wx*wy, wx*∂wy
end

function get∇²𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₂
    rz = abs(Δx[3])/x.s₃
    ∂rx = sign(Δx[1])/x.s₁
    ∂ry = sign(Δx[2])/x.s₂
    ∂rz = sign(Δx[3])/x.s₃
    wx = get𝜙ᵣ(ap,rx)
    wy = get𝜙ᵣ(ap,ry)
    wz = get𝜙ᵣ(ap,rz)
    ∂wx = get∂𝜙∂r(ap,rx)*∂rx
    ∂wy = get∂𝜙∂r(ap,ry)*∂ry
    ∂wz = get∂𝜙∂r(ap,rz)*∂rz
    ∂²wx = get∂𝜙∂r(ap,rx)*∂rx^2
    ∂²wy = get∂𝜙∂r(ap,ry)*∂ry^2
    ∂²wz = get∂𝜙∂r(ap,rz)*∂rz^2
    return wx*wy*wz, ∂wx*wy*wz, wx*∂wy*wz, ∂²wx*wy*wz, ∂wx*∂wy*wz, wx*∂²wy*wz, wx*wy*∂wz, ∂wx*wy*∂wz, wx*∂wy*∂wz, wx*wy*∂²wz
end

function get∇²₂𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₂
    ∂rx = sign(Δx[1])/x.s₁
    ∂ry = sign(Δx[2])/x.s₂
    wx = get𝜙ᵣ(ap,rx)
    wy = get𝜙ᵣ(ap,ry)
    ∂wx = get∂𝜙∂r(ap,rx)*∂rx
    ∂wy = get∂𝜙∂r(ap,ry)*∂ry
    ∂²wx = get∂𝜙∂r(ap,rx)*∂rx^2
    ∂²wy = get∂𝜙∂r(ap,ry)*∂ry^2
    return wx*wy, ∂wx*wy, wx*∂wy, ∂²wx*wy, ∂wx*∂wy, wx*∂²wy
end

function get∇³₂𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₂
    ∂rx = sign(Δx[1])/x.s₁
    ∂ry = sign(Δx[2])/x.s₂
    wx = get𝜙ᵣ(ap,rx)
    wy = get𝜙ᵣ(ap,ry)
    ∂wx = get∂𝜙∂r(ap,rx)*∂rx
    ∂wy = get∂𝜙∂r(ap,ry)*∂ry
    ∂²wx = get∂𝜙∂r(ap,rx)*∂rx^2
    ∂²wy = get∂𝜙∂r(ap,ry)*∂ry^2
    ∂³wx = get∂𝜙∂r(ap,rx)*∂rx^3
    ∂³wy = get∂𝜙∂r(ap,ry)*∂ry^3
    return wx*wy, ∂wx*wy, wx*∂wy, ∂²wx*wy, ∂wx*∂wy, wx*∂²wy, ∂³wx*wy, ∂²wx*∂wy, ∂wx*∂²wy, wx*∂³wy
end

function get∇³𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₂
    rz = abs(Δx[3])/x.s₃
    ∂rx = sign(Δx[1])/x.s₁
    ∂ry = sign(Δx[2])/x.s₂
    ∂rz = sign(Δx[3])/x.s₃
    wx = get𝜙ᵣ(ap,rx)
    wy = get𝜙ᵣ(ap,ry)
    wz = get𝜙ᵣ(ap,rz)
    ∂wx = get∂𝜙∂r(ap,rx)*∂rx
    ∂wy = get∂𝜙∂r(ap,ry)*∂ry
    ∂²wx = get∂𝜙∂r(ap,rx)*∂rx^2
    ∂²wy = get∂𝜙∂r(ap,ry)*∂ry^2
    ∂³wx = get∂𝜙∂r(ap,rx)*∂rx^3
    ∂³wy = get∂𝜙∂r(ap,ry)*∂ry^3
    return wx*wy*wz, ∂wx*wy*wz, wx*∂wy*wz, ∂²wx*wy*wz, ∂wx*∂wy*wz, wx*∂²wy*wz, ∂³wx*wy*wz, ∂²wx*∂wy*wz, ∂wx*∂²wy*wz, wx*∂³wy*wz
end

function get𝜙(ap::ReproducingKernel{𝒑,:○},x::Node,Δx::NTuple{3,Float64}) where 𝒑
    r = (Δx[1]^2+Δx[2]^2+Δx[3]^2)^0.5/x.s
    w = get𝜙ᵣ(ap,r)
    return w
end

function get∇𝜙(ap::ReproducingKernel{𝒑,:○},x::Node,Δx::NTuple{3,Float64}) where 𝒑
    r = (Δx[1]^2+Δx[2]^2+Δx[3]^2)^0.5/x.s
    ∂rx = abs(Δx[1])/x.s/(r+eps())
    ∂ry = abs(Δx[2])/x.s/(r+eps())
    ∂rz = abs(Δx[3])/x.s/(r+eps())
    w = get𝜙ᵣ(ap,r)
    ∂wr = get∂𝜙∂r(ap,r)
    return w, ∂wr*∂rx, ∂wr*∂ry, ∂wr*∂rz
end
function get𝜙ᵣ(::ReproducingKernel{𝒑,𝑠,:CubicSpline},r::Float64) where {𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return 2/3 - 4*r^2 +  4*r^3
    else
        return 4/3 - 4*r + 4*r^2 - 4*r^3/3
    end
end

function get∂𝜙∂r(::ReproducingKernel{𝒑,𝑠,:CubicSpline},r::Float64) where {𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return - 8*r + 12*r^2
    else
        return - 4   + 8*r - 4*r^2
    end
end

function get∂²𝜙∂r²(::ReproducingKernel{𝒑,𝑠,:CubicSpline},r::Float64) where {𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return - 8 + 24*r
    else
        return   8 - 8*r
    end
end

function get𝜙ᵣ(::ReproducingKernel{𝒑,𝑠,:QuinticSpline},r::Float64) where {𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 1/3
        return ((3-3r)^5 - 6(2-3r)^5 + 15(1-3r)^5)/120
    elseif r <= 2/3 && r > 1/3
        return ((3-3r)^5 - 6(2-3r)^5)/120
    else
        return (3-3r)^5/120
    end
end

function get∂𝜙∂r(::ReproducingKernel{𝒑,𝑠,:QuinticSpline},r::Float64) where {𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 1/3
        return -((3-3r)^4 - 6(2-3r)^4 + 15(1-3r)^4)/8
    elseif r <= 2/3 && r > 1/3
        return -((3-3r)^4 - 6(2-3r)^4)/8
    else
        return -(3-3r)^4/8
    end
end

function get∂²𝜙∂r²(::ReproducingKernel{𝒑,𝑠,:QuinticSpline},r::Float64) where {𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 1/3
        return ((3-3r)^3 - 6(2-3r)^3 + 15(1-3r)^3)*1.5
    elseif r <= 2/3 && r > 1/3
        return ((3-3r)^3 - 6(2-3r)^3)*1.5
    else
        return (3-3r)^3*1.5
    end
end

function get∂³𝜙∂r³(::ReproducingKernel{𝒑,𝑠,:QuinticSpline},r::Float64) where {𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 1/3
        return -((3-3r)^2 - 6(2-3r)^2 + 15(1-3r)^2)*13.5
    elseif r <= 2/3 && r > 1/3
        return -((3-3r)^2 - 6(2-3r)^2)*13.5
    else
        return -(3-3r)^2*13.5
    end
end