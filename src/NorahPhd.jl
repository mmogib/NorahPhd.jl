module NorahPhd
include("./Example1.jl")
include("./Algorithm1.jl")

Xbar, Tbar, V = algo1(Γ, P1, P2, D2, ϵ = 0.05)
print(V)
end
