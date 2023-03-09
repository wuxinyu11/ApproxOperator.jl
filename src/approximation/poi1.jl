
@inline get𝒙(ap::T,::Any) where T<:AbstractElement{:Poi1} = (ap.𝓒[1].x,ap.𝓒[1].y,ap.𝓒[1].z)
@inline get𝐽(  ::T,::Any) where T<:AbstractElement{:Poi1} = 1.0
@inline get𝑤(  ::T,::Any) where T<:AbstractElement{:Poi1} = 1.0

function set𝝭!(ap::Element{:Poi1},x::Node)
    𝝭 = x[:𝝭]
    𝝭[1] = 1.0
end