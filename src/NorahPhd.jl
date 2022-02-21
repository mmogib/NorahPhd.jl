module NorahPhd
include("./Example1.jl")
include("./Algorithm1.jl")
import Base: +
using LinearAlgebra
# import SymPy

# import Makie 
# using ImplicitEquations
# using Plots
# using LaTeXStrings
# Write your package code here.
using .Example1

using .Algorithm1
Xbar = algo1(Γ, P1, P2, D2, ϵ = 0.01)
println(Xbar)
end
