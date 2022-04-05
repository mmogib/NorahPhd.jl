using .NorahPhd
using Polyhedra
import CairoMakie as CMakie
include("./Example2.jl")
c1 = [1; 0]
c = [1; 1]
z1 = [1; 0]
z2 = [0; 1]
C = [c1 c]
Z = [z1 z2]
# v0 = [0; 0]
# w, u = D2(v0)
# P1(z1)
# D2(v0)
Xbar, Tbar, V, p = primal(Γ, P1, P2, D2, C, Z, ϵ = 0.01)
# print(p)
m = Polyhedra.Mesh(p)
CMakie.mesh(m, color = :blue)
CMakie.wireframe(m)
