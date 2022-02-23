import Base: +
struct Vertex2D{T<:Real}
    x::T
    y::T
end
Base.show(io::IO, v::Vertex2D) = print(io, "(x = $(v.x),y = $(v.y))")
Base.isless(a::Vertex2D, b::Vertex2D) = Base.isless(a.x, b.x) && Base.isless(a.y, b.y)
Base.isequal(a::Vertex2D, b::Vertex2D) = Base.isequal(a.x, b.x) && Base.isequal(a.y, b.y)
Base.promote_rule(::Type{Vertex2D{T}}, ::Type{Vertex2D{S}}) where {T<:Real,S<:Real} =
    Vertex2D{promote_type(T, S)}
Vertex2D(x::T, y::T) where {T<:Integer} = Vertex2D(convert(Float64, x), convert(Float64, y))
+(a::Vertex2D, b::Vertex2D) = Vertex2D(a.x + b.x, a.y + b.y)

Vertex2D(v::Vector{T}) where {T} =  length(v) == 2 ? Vertex2D(v[1],v[2]) : throw(0Base.DimensionMismatch("only 2-d vectors allowed"))
Base.vec(v::Vertex2D{T} where {T}) = Base.vec([v.x, v.y])