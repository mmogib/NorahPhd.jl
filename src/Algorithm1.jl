include("./Types.jl")
import Polyhedra, CDDLib

"""
    primal(Γ, P1, P2, D2, Cone, Dualcone ; ϵ=0.05)

implements Algorithm1 in 
Löhne, A., Rudloff, B. & Ulus, F. Primal and dual approximation algorithms for convex vector optimization problems. J Glob Optim 60, 713–736 (2014). [Download](https://doi.org/10.1007/s10898-013-0136-0) 
"""
function primal(Γ, P1, P2, D2, Cone, Dualcone; ϵ = 0.05)
    # any linearly independent q vectors in cone C in R^q   
    # c1, c2, ..., cq-1, c where c in intC

    C = Cone
    (_, q) = size(Cone)
    T = C
    Tinv = inv(T)

    # line 1 of Algorithm 1
    # q geberators of vectors of the dual C+ 
    
    Z = Dualcone
    # (x1, _) = P1(Z[:, 1])
    # (x2, _) = P1(Z[:, 2])
    Xbar = [first(P1(Z[:, i])) for i = 1:q]
    # H1 = Polyhedra.HalfSpace(-z1, z1' * x1)
    # H2 = Polyhedra.HalfSpace(-z2, z2' * x2)
    Hs = [Polyhedra.HalfSpace(-Z[:, i], Z[:, i]' * Xbar[i]) for i = 1:q]
    # H-rep of P_0
    # P_hr = Polyhedra.hrep([H1, H2])
    P_hr = Polyhedra.hrep(Hs)
    # P_0 as a polyhedra object
    P_poly = Polyhedra.polyhedron(P_hr)
    # t1 = T' * z1
    # t2 = T' * z2
    Tbar = [Float64.(T' * Z[:, i]) for i = 1:q]
    D = hcat([1 for i = 1:q], [Z[:, i]' * Γ(Xbar[i]) for i = 1:q])
    # stop = 1
    Vknown = Array{Vertex2D{Float64}}(undef, 100)
    j = 1
    go = true
    V = []
    while go
        # global Vknown, P_poly, j, Xbar, Tbar, P_hr
        # V-rep of P_0 
        P_vr = Polyhedra.vrep(P_poly)
        # extract number of vecrtices from P0_vr
        num_v = Polyhedra.npoints(P_vr)
        # extract vecrtices from P0_vr as row
        v = ((Polyhedra.MixedMatVRep(P_vr)).V)'
        vv = [v[:, i] for i = 1:num_v]
        # vs = [Vertex2D(v[:,i]...) for i in 1:num_v]
        # =  vs
        go = false
        for v0 in vv
            if ~(Vertex2D(v0...) ∈ Vknown)
                v0_as_vertic = Vertex2D(v0...)
                Vknown[j] = v0_as_vertic
                push!(V, v0)
                j = j + 1
                # md"__line 7__"
                (xv0, zv0) = P2(v0)
                (wv0, uv0) = D2(v0)
                # md"__line 10__"
                push!(Xbar, xv0)
                tv0 = T'wv0
                push!(Tbar, tv0)
                Dstar_q = (wv0'*v0)[1] + zv0
                D0 = hcat(D, [tv0[1]; Dstar_q])
                # [tv0[1] valtv0], [tv0[1];(wv0'*v0)[1] + zv0]
                #H0 = Polyhedra.HyperPlane(wv)
                # since zv0>ϵ we compute the huperplane in line 13
                if (zv0 > ϵ)
                    H0Mat = [tv0[1] 1] * Tinv
                    H0 = Polyhedra.HalfSpace(-[H0Mat[1], H0Mat[2]], -Dstar_q)
                    P_hr = P_poly.hrep ∩ H0
                    ## we compute P_{k+1} in line 18 and break the loop
                    P_poly = Polyhedra.polyhedron(P_hr)
                    ### we go back to line 7
                    ## we find PV 
                    # P_vr = Polyhedra.vrep(P_poly)
                    # no_v1 = Polyhedra.npoints(P1_vr)
                    # v01 = Polyhedra.MixedMatVRep(P1_vr).V
                    go = true
                end
            end
        end

    end
    return Xbar, Tbar, V
end