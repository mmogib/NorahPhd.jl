module Algorithm1
include("./Types.jl")
import Polyhedra, CDDLib
using .Types
export algo1
function algo1(Γ, P1, P2, D2; ϵ = 0.05)
    c1 = [1; 2]
    c = [1; 1]
    T = [c1 c]
    Tinv = inv(T)

    # line 1 of Algorithm 1
    z1 = [1; 0]
    z2 = [0; 1]
    Z = [z1 z2]
    (x1, _) = P1(z1)
    (x2, _) = P1(z2)
    H1 = Polyhedra.HalfSpace(-z1, z1' * x1)
    H2 = Polyhedra.HalfSpace(-z2, z2' * x2)
    # H-rep of P_0
    P_hr = Polyhedra.hrep([H1, H2])
    # P_0 as a polyhedra object
    P_poly = Polyhedra.polyhedron(P_hr)
    Xbar = [x1 x2]
    t1 = T' * z1
    t2 = T' * z2
    Tbar = [t1 t2]
    D = [1 z1'*Γ(x1); 1 z2'*Γ(x2)]
    # stop = 1
    Vknown = Array{Vertex2D{Float64}}(undef, 100)
    j = 1
    go = true
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
                Vknown[j] = Vertex2D(v0...)
                j = j + 1
                # md"__line 7__"
                (xv0, zv0) = P2(v0)
                (wv0, uv0) = D2(v0)
                # md"__line 10__"
                Xbar = [Xbar xv0]
                tv0 = T'wv0
                Tbar = [Tbar tv0]
                Dstar_q = (wv0'*v0)[1] + zv0
                D0 = [D [tv0[1]; Dstar_q]]
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
    return Xbar
end
end
