module NorahPhd
include("./Example1.jl")
include("./Algorithm1.jl")
export primal
c1 = [1; 2]
c = [1; 1]
z1 = [1; 0]
z2 = [0; 1]
C = [c1 c]
Z = [z1 z2]
Xbar, Tbar, V = primal(Γ, P1, P2, D2, C, Z, ϵ = 0.05)
print(V)
end
