
## set𝓖
function set𝓖!(aps::Vector{T},s::Symbol) where T<:AbstractElement
    rule = QuadratureRule[s]
    return set𝓖!(aps,rule)
end
function set𝓖!(aps::Vector{T},s::Symbol,stype::Symbol...) where T<:ReproducingKernel
    rule = QuadratureRule[s]
    isrk = s∈(:SegRK2,:SegRK3,:SegRK4,:SegRK5,:TriRK3,:TriRK6,:TriRK13,:TetRK14,:TetRK27)
    return set𝓖!(aps,rule,stype...;isrk=isrk)
end

function set𝓖!(aps::Vector{Element{T}},𝓖::NTuple{N,NTuple{D,Float64}}) where {T,N,D}
    nₑ = length(aps)
    nᵢ = nₑ*N
    data = Dict(:w=>zeros(nᵢ))
    data[:ξ] = zeros(nᵢ)
    D > 2 ? data[:η] = zeros(nᵢ) : nothing
    D > 3 ? data[:γ] = zeros(nᵢ) : nothing

    n = 0
    for ap in aps
        empty!(ap.𝓖)
        for ξ in 𝓖
            n += 1
            push!(ap.𝓖,Node(n,data))
            set𝓖(n,data,ξ)
        end
    end
end

function set𝓖!(aps::Vector{T},𝓖::NTuple{N,NTuple{D,Float64}},stype::Symbol...;isrk::Bool=false) where {T<:ReproducingKernel{Node},N,D}
    nₑ = length(aps)
    nᵢ = nₑ*N
    data = Dict(:w=>zeros(nᵢ))
    data[:ξ] = zeros(nᵢ)
    if isrk
        data[:wᵇ] = zeros(nᵢ)
        D > 3 ? data[:η] = zeros(nᵢ) : nothing
        D > 4 ? data[:γ] = zeros(nᵢ) : nothing
        nₕ = get𝑛𝒒(aps[1])
        aps[1].𝗠[:∇̃] = SymMat(nₕ)
    else
        D > 2 ? data[:η] = zeros(nᵢ) : nothing
        D > 3 ? data[:γ] = zeros(nᵢ) : nothing
    end

    nₘ = 0
    nₕ = length(get𝒑(aps[1],(0.0,0.0,0.0)))
    for ap in aps
        nₘ = max(nₘ,length(ap.𝓒))
    end
    for s in stype
        s∉(:∂̄x,:∂̄y,∂̄z) ? aps[1].𝗠[s] = SymMat(nₕ) : nothing
        if haskey(aps[1].𝝭,s)
            if nₘ>length(aps[1].𝝭[s])
                aps[1].𝝭[s]=zeros(nₘ)
            end
        else
            aps[1].𝝭[s]=zeros(nₘ)
        end
    end
    n = 0
    for ap in aps
        empty!(ap.𝓖)
        for ξ in 𝓖
            n += 1
            push!(ap.𝓖,Node(n,data))
            set𝓖(n,data,ξ)
        end
    end
end

function set𝓖!(aps::Vector{T},𝓖::NTuple{N,NTuple{D,Float64}},stype::Symbol...;isrk::Bool=false) where {T<:ReproducingKernel{SNode},N,D}
    nₑ = length(aps)
    nᵢ = nₑ*N
    data = Dict(:w=>zeros(nᵢ))
    data[:ξ] = zeros(nᵢ)
    if isrk
        data[:wᵇ] = zeros(nᵢ)
        D > 3 ? data[:η] = zeros(nᵢ) : nothing
        D > 4 ? data[:γ] = zeros(nᵢ) : nothing
        nₕ = get𝑛𝒒(aps[1])
        aps[1].𝗠[:∇̃]=SymMat(nₕ)
    else
        D > 2 ? data[:η] = zeros(nᵢ) : nothing
        D > 3 ? data[:γ] = zeros(nᵢ) : nothing
    end

    n = 0
    nₘ = 0
    nₕ = get𝑛𝒑(aps[1])
    index = zeros(Int,nᵢ)
    𝝭 = Dict{Symbol,Vector{Float64}}()
    for ap in aps
        n += length(ap.𝓒)*N
        nₘ = max(nₘ,length(ap.𝓒))
    end
    for s in stype
        push!(𝝭,s=>zeros(n))
        aps[1].𝗠[s]=SymMat(nₕ)
        if haskey(aps[1].𝝭,s)
            if nₘ>length(aps[1].𝝭[s])
                aps[1].𝝭[s]=zeros(nₘ)
            end
        else
            aps[1].𝝭[s]=zeros(nₘ)
        end
    end
    nₜ = 0
    n = 0
    for ap in aps
        empty!(ap.𝓖)
        for ξ in 𝓖
            n += 1
            index[n] = nₜ
            nₜ += length(ap.𝓒)
            push!(ap.𝓖,SNode(n,data,index,𝝭))
            set𝓖(n,data,ξ,isrk)
        end
    end
end

function set𝓖(i::Int,data::Dict{Symbol,Vector{Float64}},ξ::NTuple{2,Float64},isrk::Bool=false)
    data[:w][i] = ξ[1]
    data[:ξ][i] = ξ[2]
end
function set𝓖(i::Int,data::Dict{Symbol,Vector{Float64}},ξ::NTuple{3,Float64},isrk::Bool=false)
    data[:w][i] = ξ[1]
    data[:ξ][i] = ξ[2]
    isrk ? data[:wᵇ][i] = ξ[3] : data[:η][i] = ξ[3]
end
function set𝓖(i::Int,data::Dict{Symbol,Vector{Float64}},ξ::NTuple{4,Float64},isrk::Bool=false)
    data[:w][i] = ξ[1]
    data[:ξ][i] = ξ[2]
    data[:η][i] = ξ[3]
    isrk ? data[:wᵇ][i] = ξ[4] : data[:γ][i] = ξ[4]
end
function set𝓖(i::Int,data::Dict{Symbol,Vector{Float64}},ξ::NTuple{5,Float64},isrk::Bool=false)
    data[:w][i] = ξ[1]
    data[:ξ][i] = ξ[2]
    data[:η][i] = ξ[3]
    data[:γ][i] = ξ[4]
    data[:wᵇ][i] = ξ[5]
end

## get𝓖
function get𝓖(a::T,b::S) where {T<:AbstractElement{:Seg2},S<:AbstractElement{:Poi1}}
    i = findfirst(x->x.id==b.𝓒[1].id, a.𝓒)
    if i ≠ nothing && i ≤ 2
        for ξ in b.𝓖
            i == 1 ? (ξ.ξ = -1.0;ξ.n₁ = -1.0) : (ξ.ξ = 1.0;ξ.n₁ = 1.0)
        end
        return b.𝓖
    else
        return nothing
    end
end

function get𝓖(a::T,b::S) where {T<:AbstractElement{:Tri3},S<:AbstractElement{:Seg2}}
    i = findfirst(x->x.id==b.𝓒[1].id, a.𝓒)
    j = findfirst(x->x.id==b.𝓒[2].id, a.𝓒)
    if i ≠ nothing && j ≠ nothing && i ≤ 3 && j ≤ 3
        x₁ = a.𝓒[1].x
        y₁ = a.𝓒[1].y
        x₂ = a.𝓒[2].x
        y₂ = a.𝓒[2].y
        x₃ = a.𝓒[3].x
        y₃ = a.𝓒[3].y
        for ξ in b.𝓖
            if i == 1
                ξ.ξ = (1.0-ξ.ξ)/2.0
                ξ.η = 1.0-ξ.ξ
                ξ.n₁ = y₂-y₁
                ξ.n₂ = x₁-x₂
            elseif i == 2
                ξ.η = (1.0-ξ.ξ)/2.0
                ξ.ξ = 0.0
                ξ.n₁ = y₃-y₂
                ξ.n₂ = x₂-x₃
            else
                ξ.ξ = (1.0+ξ.ξ)/2.0
                ξ.η = 0.0
                ξ.n₁ = y₁-y₃
                ξ.n₂ = x₃-x₁
            end
            ξ.w *= 0.5
        end
        return b.𝓖
    else
        return nothing
    end
end

## Quadrature Points
const QuadratureRule = Dict(
:PoiGI1 => ((1.0,1.0),),
:Seg100 => ((1.0,-1.0+0.02*i) for i in 0:100),
:SegGI1 => ((2.0,0.0),),
:SegGI2 =>
(
    (1.0,-0.5773502691896257645091487805),
    (1.0, 0.5773502691896257645091487805)
),
:SegGI3 =>
(
    (0.555555555555555555555555555556,-0.774596669241483377035853079957),
    (0.88888888888888888888888888889,0.0),
    (0.555555555555555555555555555556, 0.774596669241483377035853079957)
),
:SegGI4 =>
(
    (0.347854845137453857373063949222,-0.861136311594052575223946488893),
    (0.652145154862546142626936050778,-0.339981043584856264802665759103),
    (0.652145154862546142626936050778, 0.339981043584856264802665759103),
    (0.347854845137453857373063949222, 0.861136311594052575223946488893)
),
:SegGI5 =>
(
    (0.23692688505618908751426404072,-0.906179845938663992797626878299),
    (0.47862867049936646804129151484,-0.5384693101056830910363144207),
    (0.568888888888888888888888888889, 0.0),
    (0.47862867049936646804129151484, 0.5384693101056830910363144207),
    (0.23692688505618908751426404072, 0.906179845938663992797626878299)
),
:SegGI6 =>
(
    (0.171324492379170345040296142173,-0.932469514203152027812301554494),
    (0.360761573048138607569833513838,-0.6612093864662645136613995950),
    (0.46791393457269104738987034399,-0.238619186083196908630501721681),
    (0.46791393457269104738987034399, 0.238619186083196908630501721681),
    (0.360761573048138607569833513838, 0.66120938646626451366139959502),
    (0.171324492379170345040296142173, 0.932469514203152027812301554494)
),
:SegGI7 =>
(
    (0.129484966168869693270611432679,-0.949107912342758524526189684048),
    (0.27970539148927666790146777142,-0.741531185599394439863864773281),
    (0.38183005050511894495036977549,-0.405845151377397166906606412077),
    (0.417959183673469387755102040816, 0.0),
    (0.381830050505118944950369775489, 0.405845151377397166906606412077),
    (0.279705391489276667901467771424, 0.741531185599394439863864773281),
    (0.129484966168869693270611432679, 0.949107912342758524526189684048)
),
:SegGI8 =>
(
    (0.10122853629037625915253135431,-0.96028985649753623168356086857),
    (0.22238103445337447054435599443,-0.796666477413626739591553936476),
    (0.313706645877887287337962201987,-0.525532409916328985817739049189),
    (0.36268378337836198296515044928,-0.18343464249564980493947614236),
    (0.362683783378361982965150449277, 0.18343464249564980493947614236),
    (0.31370664587788728733796220199, 0.525532409916328985817739049189),
    (0.222381034453374470544355994426, 0.796666477413626739591553936476),
    (0.10122853629037625915253135431, 0.96028985649753623168356086857)
),
:SegGI9 =>
(
    (0.0812743883615744119718921581105,-0.968160239507626089835576202904),
    (0.180648160694857404058472031243,-0.83603110732663579429942978807),
    (0.260610696402935462318742869419,-0.613371432700590397308702039342),
    (0.31234707704000284006863040658,-0.32425342340380892903853801464),
    (0.330239355001259763164525069287, 0.0),
    (0.31234707704000284006863040658, 0.32425342340380892903853801464),
    (0.260610696402935462318742869419, 0.613371432700590397308702039342),
    (0.180648160694857404058472031243, 0.83603110732663579429942978807),
    (0.081274388361574411971892158111, 0.968160239507626089835576202904)
),
:SegGI10 =>
(
    (0.066671344308688137593568809893,-0.973906528517171720077964012085),
    (0.149451349150580593145776339658,-0.865063366688984510732096688424),
    (0.219086362515982043995534934228,-0.679409568299024406234327365115),
    (0.26926671930999635509122692157,-0.433395394129247190799265943166),
    (0.295524224714752870173892994651,-0.14887433898163121088482600113),
    (0.295524224714752870173892994651, 0.14887433898163121088482600113),
    (0.26926671930999635509122692157, 0.433395394129247190799265943166),
    (0.219086362515982043995534934228, 0.679409568299024406234327365115),
    (0.149451349150580593145776339658, 0.865063366688984510732096688424),
    (0.066671344308688137593568809893, 0.973906528517171720077964012085)
),
:TriGI3 =>
(
    (1/3,2/3,1/6),
    (1/3,1/6,2/3),
    (1/3,1/6,1/6)
),
:TriGI4 =>
(
    (-0.562500000000000,0.333333333333333,0.333333333333333),
    (0.520833333333333,0.600000000000000,0.200000000000000),
    (0.520833333333333,0.200000000000000,0.600000000000000),
    (0.520833333333333,0.200000000000000,0.200000000000000)
),
:TriGI6 =>
(
    (0.223381589678011,0.108103018168070,0.445948490915965),
    (0.223381589678011,0.445948490915965,0.108103018168070),
    (0.223381589678011,0.445948490915965,0.445948490915965),
    (0.109951743655322,0.816847572980459,0.091576213509771),
    (0.109951743655322,0.091576213509771,0.816847572980459),
    (0.109951743655322,0.091576213509771,0.091576213509771)
),
:TriGI7 =>
(
    (0.125939180544800,0.101286507323500,0.101286507323500),
    (0.125939180544800,0.797426985353100,0.101286507323500),
    (0.125939180544800,0.101286507323500,0.797426985353100),
    (0.132394152788500,0.470142064105100,0.059715871789800),
    (0.132394152788500,0.470142064105100,0.470142064105100),
    (0.132394152788500,0.059715871789800,0.470142064105100),
    (0.225000000000000,0.333333333333300,0.333333333333300)
),
:TriGI12 =>
(
    (0.116786275726379,0.501426509658179,0.249286745170910),
    (0.116786275726379,0.249286745170910,0.501426509658179),
    (0.116786275726379,0.249286745170910,0.249286745170910),
    (0.050844906370207,0.873821971016996,0.063089014491502),
    (0.050844906370207,0.063089014491502,0.873821971016996),
    (0.050844906370207,0.063089014491502,0.063089014491502),
    (0.082851075618374,0.053145049844817,0.310352451033784),
    (0.082851075618374,0.053145049844817,0.636502499121399),
    (0.082851075618374,0.310352451033784,0.053145049844817),
    (0.082851075618374,0.310352451033784,0.636502499121399),
    (0.082851075618374,0.636502499121399,0.053145049844817),
    (0.082851075618374,0.636502499121399,0.310352451033784)
),
:TriGI13 =>
(
    (0.053347235608800,0.065130102902200,0.065130102902200),
    (0.053347235608800,0.869739794195600,0.065130102902200),
    (0.053347235608800,0.065130102902200,0.869739794195600),
    (0.077113760890300,0.312865496004900,0.048690315425300),
    (0.077113760890300,0.638444188569800,0.312865496004900),
    (0.077113760890300,0.048690315425300,0.638444188569800),
    (0.077113760890300,0.638444188569800,0.048690315425300),
    (0.077113760890300,0.312865496004900,0.638444188569800),
    (0.077113760890300,0.048690315425300,0.312865496004900),
    (0.175615257433200,0.260345966079000,0.260345966079000),
    (0.175615257433200,0.479308067841900,0.260345966079000),
    (0.175615257433200,0.260345966079000,0.479308067841900),
    (-0.14957004446770,0.333333333333300,0.333333333333300)
),
:TriGI16 =>
(
    (0.144315607677787,0.333333333333333,0.333333333333333),
    (0.095091634267285,0.081414823414554,0.459292588292723),
    (0.095091634267285,0.459292588292723,0.081414823414554),
    (0.095091634267285,0.459292588292723,0.459292588292723),
    (0.103217370534718,0.658861384496480,0.170569307751760),
    (0.103217370534718,0.170569307751760,0.658861384496480),
    (0.103217370534718,0.170569307751760,0.170569307751760),
    (0.032458497623198,0.898905543365938,0.050547228317031),
    (0.032458497623198,0.050547228317031,0.898905543365938),
    (0.032458497623198,0.050547228317031,0.050547228317031),
    (0.027230314174435,0.008394777409958,0.263112829634638),
    (0.027230314174435,0.008394777409958,0.728492392955404),
    (0.027230314174435,0.263112829634638,0.008394777409958),
    (0.027230314174435,0.263112829634638,0.728492392955404),
    (0.027230314174435,0.728492392955404,0.008394777409958),
    (0.027230314174435,0.728492392955404,0.263112829634638)
),
:TriGI25 =>
(
    (0.090817990382754,0.333333333333333,0.333333333333333),
    (0.036725957756467,0.028844733232685,0.485577633383657),
    (0.036725957756467,0.485577633383657,0.028844733232685),
    (0.036725957756467,0.485577633383657,0.485577633383657),
    (0.045321059435528,0.781036849029926,0.109481575485037),
    (0.045321059435528,0.109481575485037,0.781036849029926),
    (0.045321059435528,0.109481575485037,0.109481575485037),
    (0.072757916845420,0.141707219414880,0.307939838764121),
    (0.072757916845420,0.141707219414880,0.550352941820999),
    (0.072757916845420,0.307939838764121,0.141707219414880),
    (0.072757916845420,0.307939838764121,0.550352941820999),
    (0.072757916845420,0.550352941820999,0.141707219414880),
    (0.072757916845420,0.550352941820999,0.307939838764121),
    (0.028327242531057,0.025003534762686,0.246672560639903),
    (0.028327242531057,0.025003534762686,0.728323904597411),
    (0.028327242531057,0.246672560639903,0.025003534762686),
    (0.028327242531057,0.246672560639903,0.728323904597411),
    (0.028327242531057,0.728323904597411,0.025003534762686),
    (0.028327242531057,0.728323904597411,0.246672560639903),
    (0.009421666963733,0.009540815400299,0.066803251012200),
    (0.009421666963733,0.009540815400299,0.923655933587500),
    (0.009421666963733,0.066803251012200,0.009540815400299),
    (0.009421666963733,0.066803251012200,0.923655933587500),
    (0.009421666963733,0.923655933587500,0.009540815400299),
    (0.009421666963733,0.923655933587500,0.066803251012200)
),
:QuadGI1 =>
(
    (2.0,0.0,0.0),
),
:QuadGI2 =>
(
    (1.0,-0.5773502691896258,-0.5773502691896258),
    (1.0, 0.5773502691896258,-0.5773502691896258),
    (1.0, 0.5773502691896258, 0.5773502691896258),
    (1.0,-0.5773502691896258, 0.5773502691896258)
),
:SegRK2 =>
(
    (1.0,-1.0,1.0),
    (1.0, 1.0,1.0)
),
:SegRK3 =>
(
    (1/3,-1.0,1.0),
    (4/3, 0.0,0.0),
    (1/3, 1.0,1.0)
),
:SegRK4 =>
(
    (1/6,-1.0,    1.0),
    (5/6,-0.2^0.5,0.0),
    (5/6, 0.2^0.5,0.0),
    (1/6, 1.0,    1.0)
),
:SegRK5 =>
(
    (1/10, -1.0,      1.0),
    (49/90,-(3/7)^0.5,0.0),
    (32/45, 0.0,      0.0),
    (49/90, (3/7)^0.5,0.0),
    (1/10,  1.0,      1.0)
),
:TriRK3 =>
(
    (1/3,1.0,0.0,1/2),
    (1/3,0.0,1.0,1/2),
    (1/3,0.0,0.0,1/2)
),
:TriRK6 =>
(
    (0.0,1.0,0.0,1/6),
    (0.0,0.0,1.0,1/6),
    (0.0,0.0,0.0,1/6),
    (1/3,0.0,0.5,2/3),
    (1/3,0.5,0.0,2/3),
    (1/3,0.5,0.5,2/3)
),
:TriRK13 =>
(
    (-0.0277777777777778,1.0000000000000000,0.0000000000000000, 1/20),
    (-0.0277777777777778,0.0000000000000000,1.0000000000000000, 1/20),
    (-0.0277777777777778,0.0000000000000000,0.0000000000000000, 1/20),
    ( 0.0296296296296297,0.0000000000000000,0.5000000000000000,16/45),
    ( 0.0296296296296297,0.5000000000000000,0.0000000000000000,16/45),
    ( 0.0296296296296297,0.5000000000000000,0.5000000000000000,16/45),
    ( 0.0907407407407407,0.0000000000000000,0.8273268353539885,49/180),
    ( 0.0907407407407407,0.0000000000000000,0.1726731646460116,49/180),
    ( 0.0907407407407407,0.1726731646460116,0.0000000000000000,49/180),
    ( 0.0907407407407407,0.8273268353539885,0.0000000000000000,49/180),
    ( 0.0907407407407407,0.8273268353539885,0.1726731646460116,49/180),
    ( 0.0907407407407407,0.1726731646460116,0.8273268353539885,49/180),
    ( 0.4500000000000000,0.3333333333333333,0.3333333333333333,  0.0)
)
)
