using Plots
include("../src/Types.jl")
@recipe function f(::Type{Vertex2D}, V::Vertex2D, add_marker = false)
    linecolor --> :blue
    seriestype := :scatter
    markershape --> (add_marker ? :circle : :none)
    delete!(plotattributes, :add_marker)
    vec.(V)
end
V = Vertex2D.([rand([-1, -2, 0, -3], 2) for i in 1:5])
plt = plot(; framestyle = :origin,
    aspect_ratio = 1)
# plot!(plt, V)
vec.(V)
plot!(plt, vec.(V), seriestype = :scatter)

