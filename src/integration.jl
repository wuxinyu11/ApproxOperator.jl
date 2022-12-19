
function set𝓖!(as::Vector{T},ss::Symbol) where T<:AbstractElement
    data_ = quadraturerule(ss)
    n = length(data_[:w])
    G = 0
    s = 0
    data = Dict([s=>(1,v) for (s,v) in data_])
    for (c,a) in enumerate(as)
        empty!(a.𝓖)
        for g in 1:n
            G += 1
            push!(a.𝓖,SNode((g,G,c,s),data))
            s += length(a.𝓒)
        end
    end
    setgeometry!(as)
    if ss ∈ RKQuadratureRule set𝑫!(as) end
end

function set𝓖!(as::Vector{T},bs::Vector{S}) where {T<:AbstractElement,S<:AbstractElement}
    data = getfield(bs[1].𝓖[1],:data)
    s = 0
    nₑ = length(as)
    for c in 1:nₑ
        a = as[c]
        b = bs[c]
        empty!(a.𝓖)
        for ξ_ in b.𝓖
            g = ξ_.𝑔
            G = ξ_.𝐺
            push!(a.𝓖,SNode((g,G,c,s),data))
            s += length(a.𝓒)
        end
    end
end

function set𝓖!(as::Vector{T},bs::Vector{S}) where {T<:AbstractElement{:Seg2},S<:AbstractElement{:Poi1}}
    nₑ = length(bs)
    data = Dict([:ξ=>(1,[-1.0,1.0]),:w=>(1,[1.0,1.0])],:x=>(2,zeros(nₑ)),:y=>(2,zeros(nₑ)),:z=>(2,zeros(nₑ)))
    s = 0
    G = 0
    for (c,a) in enumerate(as)
        empty!(a.𝓖)
        for b in bs
            g = findfirst(x->x.𝐼==b.𝓒[1].𝐼, a.𝓒)
            if i ≠ nothing && i ≤ 2
                G += 1
                ξ = SNode((g,G,c,s),data)
                ξ.x = b.𝓖[1].x
                ξ.y = b.𝓖[1].y
                ξ.z = b.𝓖[1].z
                push!(a.𝓖,)
                s += length(a.𝓒)
            end
        end
    end
end

function set𝓖!(as::Vector{T},bs::Vector{S}) where {T<:AbstractElement{:Tri3},S<:AbstractElement{:Poi1}}
    nₑ = 0
    for b in bs
        for a in as
            g = findfirst(x->x.𝐼==b.𝓒[1].𝐼, a.𝓒)
            g ≠ nothing && g ≤ 3 ? nₑ += 1 : nothing
        end
    end
    data = Dict([:ξ=>(1,[1.0,0.0,0.0]),:η=>(1,[0.0,1.0,0.0]),:w=>(1,[1.0,1.0,1.0])],:x=>(2,zeros(nₑ)),:y=>(2,zeros(nₑ)),:z=>(2,zeros(nₑ)),:Δn₁s₁=>(2,zeros(nₑ)),:Δn₁s₂n₂s₁=>(2,zeros(nₑ)),:Δn₂s₂=>(2,zeros(nₑ)))
    s = 0
    G = 0
    for (c,a) in enumerate(as)
        empty!(a.𝓖)
        x₁ = a.𝓒[1].x
        y₁ = a.𝓒[1].y
        x₂ = a.𝓒[2].x
        y₂ = a.𝓒[2].y
        x₃ = a.𝓒[3].x
        y₃ = a.𝓒[3].y
        n₁₁ = y₃-y₂;n₂₁ = y₁-y₃;n₃₁ = y₂-y₁
        n₁₂ = x₂-x₃;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂
        s₁₁ = -n₁₂;s₂₁ = -n₂₂;s₃₁ = -n₃₂
        s₁₂ =  n₁₁;s₂₂ =  n₂₁;s₃₂ =  n₃₁
        𝐿₁² = n₁₁^2+n₁₂^2
        𝐿₂² = n₂₁^2+n₂₂^2
        𝐿₃² = n₃₁^2+n₃₂^2
        for b in bs
            g = findfirst(x->x.𝐼==b.𝓒[1].𝐼, a.𝓒)
            if g ≠ nothing && g ≤ 3
                G += 1
                ξ = SNode((g,G,c,s),data)
                s += length(a.𝓒)
                if g == 1
                    ξ.x = a.𝓒[1].x
                    ξ.y = a.𝓒[1].y
                    ξ.z = a.𝓒[1].z
                    ξ.Δn₁s₁ = n₂₁*s₂₁/𝐿₂² - n₃₁*s₃₁/𝐿₃²
                    ξ.Δn₁s₂n₂s₁ = n₂₁*s₂₂/𝐿₂² + n₂₂*s₂₁/𝐿₂² - n₃₁*s₃₂/𝐿₃² - n₃₂*s₃₁/𝐿₃²
                    ξ.Δn₂s₂ = n₂₂*s₂₂/𝐿₂² - n₃₂*s₃₂/𝐿₃²
                elseif g == 2
                    ξ.x = a.𝓒[2].x
                    ξ.y = a.𝓒[2].y
                    ξ.z = a.𝓒[2].z
                    ξ.Δn₁s₁ = n₃₁*s₃₁/𝐿₃² - n₁₁*s₁₁/𝐿₁²
                    ξ.Δn₁s₂n₂s₁ = n₃₁*s₃₂/𝐿₃² + n₃₂*s₃₁/𝐿₃² - n₁₁*s₁₂/𝐿₁² - n₁₂*s₁₁/𝐿₁²
                    ξ.Δn₂s₂ = n₃₂*s₃₂/𝐿₃² - n₁₂*s₁₂/𝐿₁²
                else
                    ξ.x = a.𝓒[3].x
                    ξ.y = a.𝓒[3].y
                    ξ.z = a.𝓒[3].z
                    ξ.Δn₁s₁ = n₁₁*s₁₁/𝐿₁² - n₂₁*s₂₁/𝐿₂²
                    ξ.Δn₁s₂n₂s₁ = n₁₁*s₁₂/𝐿₁² + n₁₂*s₁₁/𝐿₁² - n₂₁*s₂₂/𝐿₂² - n₂₂*s₂₁/𝐿₂²
                    ξ.Δn₂s₂ = n₁₂*s₁₂/𝐿₁² - n₂₂*s₂₂/𝐿₂²
                end
                push!(a.𝓖,ξ)
            end
        end
    end
end

function set𝓖!(as::Vector{T},bs::Vector{S}) where {T<:AbstractElement{:Tri3},S<:AbstractElement{:Seg2}}
    nₑ = length(bs)
    nᵢ = length(getfield(bs[1].𝓖[1],:data)[:w][2])
    data = Dict([:ξ=>(2,zeros(nₑ*nᵢ)),:η=>(2,zeros(nₑ*nᵢ)),:w=>(2,zeros(nₑ*nᵢ)),:x=>(2,zeros(nₑ*nᵢ)),:y=>(2,zeros(nₑ*nᵢ)),:z=>(2,zeros(nₑ*nᵢ)),:𝑤=>(2,zeros(nₑ*nᵢ)),:n₁=>(2,zeros(nₑ*nᵢ)),:n₂=>(2,zeros(nₑ*nᵢ)),:𝐴=>(3,zeros(nₑ))])
    G = 0
    s = 0
    for (c,a) in enumerate(as)
        empty!(a.𝓖)
        for b in bs
            i = T<:DiscreteElement ? findfirst(x->x.𝑖==b.𝓒[1].𝐼, a.𝓒) : findfirst(x->x.𝐼==b.𝓒[1].𝐼, a.𝓒)
            j = T<:DiscreteElement ? findfirst(x->x.𝑖==b.𝓒[2].𝐼, a.𝓒) : findfirst(x->x.𝐼==b.𝓒[2].𝐼, a.𝓒)
            if i ≠ nothing && j ≠ nothing && i ≤ 3 && j ≤ 3
                𝐴 = get𝐴(a)
                for ξ_ in b.𝓖
                    G += 1
                    ξ = SNode((ξ_.𝑔,G,c,s),data)
                    s += length(a.𝓒)
                    if i == 1
                        ξ.ξ = (1.0-ξ_.ξ)/2.0
                        ξ.η = 1.0-ξ.ξ
                    elseif i == 2
                        ξ.η = (1.0-ξ_.ξ)/2.0
                        ξ.ξ = 0.0
                    else
                        ξ.ξ = (1.0+ξ_.ξ)/2.0
                        ξ.η = 0.0
                    end
                    ξ.x = ξ_.x
                    ξ.y = ξ_.y
                    ξ.z = ξ_.z
                    ξ.w = 0.5*ξ_.w
                    ξ.𝑤 = 0.5*ξ_.𝑤
                    ξ.n₁ = ξ_.n₁
                    ξ.n₂ = ξ_.n₂
                    ξ.𝐴 = 𝐴
                    push!(a.𝓖,ξ)
                end
            end
        end
    end
end

function set𝓖_DB!(aps::Vector{T},s::Symbol) where T<:AbstractElement
    data_ = quadraturerule(s)
    n = length(data_[:w])
    nₑ = length(aps)
    data = Dict([:w=>(1,zeros(3*n)),:ξ=>(1,zeros(3*n)),:η=>(1,zeros(3*n)),:n₁=>(2,zeros(3*nₑ*n)),:n₂=>(2,zeros(3*nₑ*n)),:x=>(2,zeros(3*nₑ*n)),:y=>(2,zeros(3*nₑ*n)),:z=>(2,zeros(3*nₑ*n)),:𝑤=>(2,zeros(3*nₑ*n))])
    for i in 1:n
        w = data_[:w][i]
        ξ = data_[:ξ][i]
        data[:w][2][3*i-2] = w
        data[:w][2][3*i-1] = w
        data[:w][2][3*i]   = w
        data[:η][2][3*i-2] = (1-ξ)/2
        data[:ξ][2][3*i-1] = (1+ξ)/2
        data[:ξ][2][3*i] = (1-ξ)/2
        data[:η][2][3*i] = (1+ξ)/2
    end
    G = 0
    s = 0
    for ap in aps
        empty!(ap.𝓖)
        for g in 1:3*n
            G += 1
            push!(ap.𝓖,SNode((g,G,s),data))
            s += length(ap.𝓒)
        end
        x₁ = ap.𝓒[1].x
        x₂ = ap.𝓒[2].x
        x₃ = ap.𝓒[3].x
        y₁ = ap.𝓒[1].y
        y₂ = ap.𝓒[2].y
        y₃ = ap.𝓒[3].y
        𝐿₁ = ((x₂-x₃)^2+(y₂-y₃)^2)^0.5
        𝐿₂ = ((x₃-x₁)^2+(y₃-y₁)^2)^0.5
        𝐿₃ = ((x₁-x₂)^2+(y₁-y₂)^2)^0.5
        n₁₁ = (y₃-y₂)/𝐿₁
        n₂₁ = (x₂-x₃)/𝐿₁
        n₁₂ = (y₁-y₃)/𝐿₂
        n₂₂ = (x₃-x₁)/𝐿₂
        n₁₃ = (y₂-y₁)/𝐿₃
        n₂₃ = (x₁-x₂)/𝐿₃
        for i in 1:n
            ξ₁ = ap.𝓖[3*i-2]
            ξ₂ = ap.𝓖[3*i-1]
            ξ₃ = ap.𝓖[3*i]
            ξ₁.n₁ = n₁₁
            ξ₁.n₂ = n₂₁
            ξ₂.n₁ = n₁₂
            ξ₂.n₂ = n₂₂
            ξ₃.n₁ = n₁₃
            ξ₃.n₂ = n₂₃
            ξ₁.x = x₁*ξ₁.ξ + x₂*ξ₁.η + x₃*(1-ξ₁.ξ-ξ₁.η)
            ξ₁.y = y₁*ξ₁.ξ + y₂*ξ₁.η + y₃*(1-ξ₁.ξ-ξ₁.η)
            ξ₂.x = x₁*ξ₂.ξ + x₂*ξ₂.η + x₃*(1-ξ₂.ξ-ξ₂.η)
            ξ₂.y = y₁*ξ₂.ξ + y₂*ξ₂.η + y₃*(1-ξ₂.ξ-ξ₂.η)
            ξ₃.x = x₁*ξ₃.ξ + x₂*ξ₃.η + x₃*(1-ξ₃.ξ-ξ₃.η)
            ξ₃.y = y₁*ξ₃.ξ + y₂*ξ₃.η + y₃*(1-ξ₃.ξ-ξ₃.η)
            ξ₁.𝑤 = 𝐿₁/2*ξ₁.w
            ξ₂.𝑤 = 𝐿₂/2*ξ₂.w
            ξ₃.𝑤 = 𝐿₃/2*ξ₃.w
        end
    end
end
"""
quadraturerule(s::Symbol)
"""
const RKQuadratureRule = (:SegRK2,:SegRK3,:SegRK4,:SegRK5,:TriRK3,:TriRK6,:TriRK13)
function quadraturerule(s::Symbol)
    if s == :PoiGI1
        return Dict([:w=>[1.0],:ξ=>[1.0]])
    elseif s == :SegGI1
        return Dict([:w=>[2.0],:ξ=>[0.0]])
    elseif s == :SegGI2
        return Dict([
            :w=>[1.0,1.0],
            :ξ=>[-0.5773502691896257645091487805,
                  0.5773502691896257645091487805]
        ])
    elseif s == :SegGI3
        return Dict([
            :w=>[0.555555555555555555555555555556,
                 0.88888888888888888888888888889,
                 0.555555555555555555555555555556],
            :ξ=>[-0.774596669241483377035853079957,
                  0.0,
                  0.774596669241483377035853079957]
        ])
    elseif s == :SegGI4
        return Dict([
            :w=>[0.347854845137453857373063949222,
                 0.652145154862546142626936050778,
                 0.652145154862546142626936050778,
                 0.347854845137453857373063949222],
            :ξ=>[-0.861136311594052575223946488893,
                 -0.339981043584856264802665759103,
                  0.339981043584856264802665759103,
                  0.861136311594052575223946488893]
        ])
    elseif s == :SegGI5
        return Dict([
            :w=>[0.23692688505618908751426404072,
                 0.47862867049936646804129151484,
                 0.568888888888888888888888888889,
                 0.47862867049936646804129151484,
                 0.23692688505618908751426404072],
            :ξ=>[-0.906179845938663992797626878299,
                 -0.5384693101056830910363144207,
                  0.0,
                  0.5384693101056830910363144207,
                  0.906179845938663992797626878299]
        ])
    elseif s == :SegGI6
        return Dict([
            :w=>[0.171324492379170345040296142173,
                 0.360761573048138607569833513838,
                 0.46791393457269104738987034399,
                 0.46791393457269104738987034399,
                 0.360761573048138607569833513838,
                 0.171324492379170345040296142173],
            :ξ=>[-0.932469514203152027812301554494,
                 -0.6612093864662645136613995950,
                 -0.238619186083196908630501721681,
                  0.238619186083196908630501721681,
                  0.66120938646626451366139959502,
                  0.932469514203152027812301554494]
        ])
    elseif s == :SegGI7
        return Dict([
            :w=>[0.129484966168869693270611432679,
                 0.27970539148927666790146777142,
                 0.38183005050511894495036977549,
                 0.417959183673469387755102040816,
                 0.381830050505118944950369775489,
                 0.279705391489276667901467771424,
                 0.129484966168869693270611432679],
            :ξ=>[-0.949107912342758524526189684048,
                 -0.741531185599394439863864773281,
                 -0.405845151377397166906606412077,
                  0.0,
                  0.405845151377397166906606412077,
                  0.741531185599394439863864773281,
                  0.949107912342758524526189684048]
        ])
    elseif s == :SegGI8
        return Dict([
            :w=>[0.10122853629037625915253135431,
                 0.22238103445337447054435599443,
                 0.313706645877887287337962201987,
                 0.36268378337836198296515044928,
                 0.362683783378361982965150449277,
                 0.31370664587788728733796220199,
                 0.222381034453374470544355994426,
                 0.10122853629037625915253135431],
            :ξ=>[-0.96028985649753623168356086857,
                 -0.796666477413626739591553936476,
                 -0.525532409916328985817739049189,
                 -0.18343464249564980493947614236,
                  0.18343464249564980493947614236,
                  0.525532409916328985817739049189,
                  0.796666477413626739591553936476,
                  0.96028985649753623168356086857]
        ])
    elseif s == :SegGI9
        return Dict([
            :w=>[0.081274388361574411971892158111,
                 0.180648160694857404058472031243,
                 0.260610696402935462318742869419,
                 0.31234707704000284006863040658,
                 0.330239355001259763164525069287,
                 0.31234707704000284006863040658,
                 0.260610696402935462318742869419,
                 0.180648160694857404058472031243,
                 0.081274388361574411971892158111],
            :ξ=>[-0.968160239507626089835576202904,
                 -0.83603110732663579429942978807,
                 -0.613371432700590397308702039342,
                 -0.32425342340380892903853801464,
                  0.0,
                  0.32425342340380892903853801464,
                  0.613371432700590397308702039342,
                  0.83603110732663579429942978807,
                  0.968160239507626089835576202904]
        ])
    elseif s == :SegGI10
        return Dict([
            :w=>[0.066671344308688137593568809893,
                 0.149451349150580593145776339658,
                 0.219086362515982043995534934228,
                 0.26926671930999635509122692157,
                 0.295524224714752870173892994651,
                 0.295524224714752870173892994651,
                 0.26926671930999635509122692157,
                 0.219086362515982043995534934228,
                 0.149451349150580593145776339658,
                 0.066671344308688137593568809893],
            :ξ=>[-0.973906528517171720077964012085,
                 -0.865063366688984510732096688424,
                 -0.679409568299024406234327365115,
                 -0.433395394129247190799265943166,
                 -0.14887433898163121088482600113,
                  0.14887433898163121088482600113,
                  0.433395394129247190799265943166,
                  0.679409568299024406234327365115,
                  0.865063366688984510732096688424,
                  0.973906528517171720077964012085]
        ])
    elseif s == :SegRG100
        return Dict([
            :w=>ones(100)/100,
            :ξ=>collect(-1.0:2/99:1.0)
        ])
    elseif s == :TriGI1
        return Dict([
            :w=>[1.0],
            :ξ=>[1/3],
            :η=>[1/3]
        ])
    elseif s == :TriGI3
        return Dict([
            :w=>[1/3,1/3,1/3],
            :ξ=>[2/3,1/6,1/6],
            :η=>[1/6,2/3,1/6]
        ])
    elseif s == :TriGI4
        return Dict([
            :w=>[-0.562500000000000,
                  0.520833333333333,
                  0.520833333333333,
                  0.520833333333333],
            :ξ=>[0.333333333333333,
                 0.600000000000000,
                 0.200000000000000,
                 0.200000000000000],
            :η=>[0.333333333333333,
                 0.200000000000000,
                 0.600000000000000,
                 0.200000000000000]
        ])
    elseif s == :TriGI6
        return Dict([
            :w=>[0.223381589678011,
                 0.223381589678011,
                 0.223381589678011,
                 0.109951743655322,
                 0.109951743655322,
                 0.109951743655322],
            :ξ=>[0.108103018168070,
                 0.445948490915965,
                 0.445948490915965,
                 0.816847572980459,
                 0.091576213509771,
                 0.091576213509771],
            :η=>[0.445948490915965,
                 0.108103018168070,
                 0.445948490915965,
                 0.091576213509771,
                 0.816847572980459,
                 0.091576213509771]
        ])
    elseif s == :TriGI7
        return Dict([
            :w=>[0.125939180544800,
                 0.125939180544800,
                 0.125939180544800,
                 0.132394152788500,
                 0.132394152788500,
                 0.132394152788500,
                 0.225000000000000],
            :ξ=>[0.101286507323500,
                 0.797426985353100,
                 0.101286507323500,
                 0.470142064105100,
                 0.470142064105100,
                 0.059715871789800,
                 0.333333333333300],
            :η=>[0.101286507323500,
                 0.101286507323500,
                 0.797426985353100,
                 0.059715871789800,
                 0.470142064105100,
                 0.470142064105100,
                 0.333333333333300]
        ])
    elseif s == :TriGI12
        return Dict([
            :w=>[0.116786275726379,
                 0.116786275726379,
                 0.116786275726379,
                 0.050844906370207,
                 0.050844906370207,
                 0.050844906370207,
                 0.082851075618374,
                 0.082851075618374,
                 0.082851075618374,
                 0.082851075618374,
                 0.082851075618374,
                 0.082851075618374],
            :ξ=>[0.501426509658179,
                 0.249286745170910,
                 0.249286745170910,
                 0.873821971016996,
                 0.063089014491502,
                 0.063089014491502,
                 0.053145049844817,
                 0.053145049844817,
                 0.310352451033784,
                 0.310352451033784,
                 0.636502499121399,
                 0.636502499121399],
            :η=>[0.249286745170910,
                 0.501426509658179,
                 0.249286745170910,
                 0.063089014491502,
                 0.873821971016996,
                 0.063089014491502,
                 0.310352451033784,
                 0.636502499121399,
                 0.053145049844817,
                 0.636502499121399,
                 0.053145049844817,
                 0.310352451033784]
        ])
    elseif s == :TriGI13
        return Dict([
            :w=>[0.053347235608800,
                 0.053347235608800,
                 0.053347235608800,
                 0.077113760890300,
                 0.077113760890300,
                 0.077113760890300,
                 0.077113760890300,
                 0.077113760890300,
                 0.077113760890300,
                 0.175615257433200,
                 0.175615257433200,
                 0.175615257433200,
                 -0.14957004446770],
            :ξ=>[0.065130102902200,
                 0.869739794195600,
                 0.065130102902200,
                 0.312865496004900,
                 0.638444188569800,
                 0.048690315425300,
                 0.638444188569800,
                 0.312865496004900,
                 0.048690315425300,
                 0.260345966079000,
                 0.479308067841900,
                 0.260345966079000,
                 0.333333333333300],
            :η=>[0.065130102902200,
                 0.065130102902200,
                 0.869739794195600,
                 0.048690315425300,
                 0.312865496004900,
                 0.638444188569800,
                 0.048690315425300,
                 0.638444188569800,
                 0.312865496004900,
                 0.260345966079000,
                 0.260345966079000,
                 0.479308067841900,
                 0.333333333333300]
        ])
    elseif s == :TriGI16
        return Dict([
            :w=>[0.144315607677787,
                 0.095091634267285,
                 0.095091634267285,
                 0.095091634267285,
                 0.103217370534718,
                 0.103217370534718,
                 0.103217370534718,
                 0.032458497623198,
                 0.032458497623198,
                 0.032458497623198,
                 0.027230314174435,
                 0.027230314174435,
                 0.027230314174435,
                 0.027230314174435,
                 0.027230314174435,
                 0.027230314174435],
            :ξ=>[0.333333333333333,
                 0.081414823414554,
                 0.459292588292723,
                 0.459292588292723,
                 0.658861384496480,
                 0.170569307751760,
                 0.170569307751760,
                 0.898905543365938,
                 0.050547228317031,
                 0.050547228317031,
                 0.008394777409958,
                 0.008394777409958,
                 0.263112829634638,
                 0.263112829634638,
                 0.728492392955404,
                 0.728492392955404],
            :η=>[0.333333333333333,
                 0.459292588292723,
                 0.081414823414554,
                 0.459292588292723,
                 0.170569307751760,
                 0.658861384496480,
                 0.170569307751760,
                 0.050547228317031,
                 0.898905543365938,
                 0.050547228317031,
                 0.263112829634638,
                 0.728492392955404,
                 0.008394777409958,
                 0.728492392955404,
                 0.008394777409958,
                 0.263112829634638]
        ])
    elseif s == :TriGI25
        return Dict([
            :w=>[0.090817990382754,
                 0.036725957756467,
                 0.036725957756467,
                 0.036725957756467,
                 0.045321059435528,
                 0.045321059435528,
                 0.045321059435528,
                 0.072757916845420,
                 0.072757916845420,
                 0.072757916845420,
                 0.072757916845420,
                 0.072757916845420,
                 0.072757916845420,
                 0.028327242531057,
                 0.028327242531057,
                 0.028327242531057,
                 0.028327242531057,
                 0.028327242531057,
                 0.028327242531057,
                 0.009421666963733,
                 0.009421666963733,
                 0.009421666963733,
                 0.009421666963733,
                 0.009421666963733,
                 0.009421666963733],
            :ξ=>[0.333333333333333,
                 0.028844733232685,
                 0.485577633383657,
                 0.485577633383657,
                 0.781036849029926,
                 0.109481575485037,
                 0.109481575485037,
                 0.141707219414880,
                 0.141707219414880,
                 0.307939838764121,
                 0.307939838764121,
                 0.550352941820999,
                 0.550352941820999,
                 0.025003534762686,
                 0.025003534762686,
                 0.246672560639903,
                 0.246672560639903,
                 0.728323904597411,
                 0.728323904597411,
                 0.009540815400299,
                 0.009540815400299,
                 0.066803251012200,
                 0.066803251012200,
                 0.923655933587500,
                 0.923655933587500],
            :η=>[0.333333333333333,
                 0.485577633383657,
                 0.028844733232685,
                 0.485577633383657,
                 0.109481575485037,
                 0.781036849029926,
                 0.109481575485037,
                 0.307939838764121,
                 0.550352941820999,
                 0.141707219414880,
                 0.550352941820999,
                 0.141707219414880,
                 0.307939838764121,
                 0.246672560639903,
                 0.728323904597411,
                 0.025003534762686,
                 0.728323904597411,
                 0.025003534762686,
                 0.246672560639903,
                 0.066803251012200,
                 0.923655933587500,
                 0.009540815400299,
                 0.923655933587500,
                 0.009540815400299,
                 0.066803251012200]
        ])
    elseif s == :QuadGI1
        return Dict([
            :w=>[2.0],:ξ=>[0.0],:η=>[0.0]
        ])
    elseif s == :QuadGI4
        return Dict([
            :w=>[1.0,1.0,1.0,1.0],
            :ξ=>[-0.5773502691896258,
                  0.5773502691896258,
                  0.5773502691896258,
                 -0.5773502691896258],
            :η=>[-0.5773502691896258,
                 -0.5773502691896258,
                  0.5773502691896258,
                  0.5773502691896258]
        ])
    elseif s == :SegRK2
        return Dict([
            :w=>[1.0,1.0],
            :ξ=>[-1.0,1.0],
            :wᵇ=>[ 1.0,1.0]
        ])
    elseif s == :SegRK3
        return Dict([
            :w=>[1/3,4/3,1/3],
            :ξ=>[-1.0,0.0,1.0],
            :wᵇ=>[ 1.0,0.0,1.0]
        ])
    elseif s == :SegRK4
        return Dict([
            :w=>[1/6,5/6,1/6],
            :ξ=>[-1.0,-0.2^0.5,0.2^0.5,1.0],
            :wᵇ=>[1.0,0.0,0.0,1.0]
        ])
    elseif s == :SegRK5
        return Dict([
            :w=>[1/10,49/90,32/45,49/90,1/10],
            :ξ=>[-1.0,-(3/7)^0.5,0.0,(3/7)^0.5,1.0],
            :wᵇ=>[1.0,0.0,0.0,0.0,1.0]
        ])
    elseif s == :TriRK3
        return Dict([
            :w=>[1/3,1/3,1/3],
            :ξ=>[1.0,0.0,0.0],
            :η=>[0.0,1.0,0.0],
            :wᵇ=>[1/2,1/2,1/2]
        ])
    elseif s == :TriRK6
        return Dict([
            :w=>[0.0,0.0,0.0,1/3,1/3,1/3],
            :ξ=>[1.0,0.0,0.0,0.0,0.5,0.5],
            :η=>[0.0,1.0,0.0,0.5,0.0,0.5],
            :wᵇ=>[1/6,1/6,1/6,2/3,2/3,2/3]
        ])
    elseif s == :TriRK13
        return Dict([
            :w=>[-0.0277777777777778,
                 -0.0277777777777778,
                 -0.0277777777777778,
                  0.0296296296296297,
                  0.0296296296296297,
                  0.0296296296296297,
                  0.0907407407407407,
                  0.0907407407407407,
                  0.0907407407407407,
                  0.0907407407407407,
                  0.0907407407407407,
                  0.0907407407407407,
                  0.4500000000000000],
            :ξ=>[1.0000000000000000,
                 0.0000000000000000,
                 0.0000000000000000,
                 0.0000000000000000,
                 0.5000000000000000,
                 0.5000000000000000,
                 0.0000000000000000,
                 0.0000000000000000,
                 0.1726731646460114,
                 0.8273268353539885,
                 0.8273268353539885,
                 0.1726731646460114,
                 0.3333333333333333],
            :η=>[0.0000000000000000,
                 1.0000000000000000,
                 0.0000000000000000,
                 0.5000000000000000,
                 0.0000000000000000,
                 0.5000000000000000,
                 0.8273268353539885,
                 0.1726731646460114,
                 0.0000000000000000,
                 0.0000000000000000,
                 0.1726731646460114,
                 0.8273268353539885,
                 0.3333333333333333],
            :wᵇ=>[ 1/20,
                   1/20,
                   1/20,
                  16/45,
                  16/45,
                  16/45,
                  49/180,
                  49/180,
                  49/180,
                  49/180,
                  49/180,
                  49/180,
                    0.0]
        ])
    end
end
