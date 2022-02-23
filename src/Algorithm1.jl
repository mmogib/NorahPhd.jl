include("./Types.jl")
import Polyhedra, CDDLib
export primal
"""
    primal(Γ, P1, P2, D2, Cone, Dualcone ; ϵ=0.05)

implements Algorithm1 in 
Löhne, A., Rudloff, B. & Ulus, F. Primal and dual approximation algorithms for convex vector optimization problems. J Glob Optim 60, 713–736 (2014). [Download](https://doi.org/10.1007/s10898-013-0136-0) 
"""
function primal(Γ, P1, P2, D2, Cone, Dualcone; ϵ = 0.05, max_num_vertices = 100)
    # any linearly independent q vectors in cone C in R^q   
    # c1, c2, ..., cq-1, c where c in intC
    C = Cone
    # q geberators of vectors of the dual C+ 
    (_, q) = size(Cone)
    T = C
    Z = Dualcone
    # lines 1 and 3
    Xbar = [first(P1(Z[:, i])) for i = 1:q]
    # line 2
    Hs = [Polyhedra.HalfSpace(-Z[:, i], Z[:, i]' * Xbar[i]) for i = 1:q]
    P_hr = Polyhedra.hrep(Hs)
    P_poly = Polyhedra.polyhedron(P_hr)
    # line 3
    k = 1
    Tinv = inv(T)
    Tbar = [Float64.(T' * Z[:, i]) for i = 1:q]
    # line 4  
    # D = D*(Tbar) = [Tbar₁, Tbar₂, ⋯, Tbar_q,   ]  
    D = hcat([Tbar[i][1:q-1][1] for i in 1:q], [Z[:, i]' * Γ(Xbar[i]) for i = 1:q])
    
    # bags to collect the vertices
    # Vknown = similar(rand(100),Vertex2D{Float64})
    
    # Vknown = Array{Union{Missing,Vertex2D{Float64}}}(missing, max_num_vertices)
    V = []
    # flag to break the loop
    go = true
    while go
        # line 7
        P_vr = Polyhedra.vrep(P_poly)
        # extract number of vecrtices from P_vr
        # that is |Pⱽ|
        num_v = Polyhedra.npoints(P_vr)
        # extract vecrtices from P_vr as row
        v = ((Polyhedra.MixedMatVRep(P_vr)).V)'
        vv = [v[:, i] for i = 1:num_v]
        ## by default stop the loop
        go = false

        # for each vertex found
        for v0 in vv
            # check if this vertex has been visited before? 
            if ~(Vertex2D(v0...) ∈ Vertex2D.(V))
                # if the vertex has not been visited
                # v0_as_vertic = Vertex2D(v0...)
                ## store this vertex as visited
                # Vknown[k] = v0_as_vertic
                ## collect it in the bag
                push!(V, v0)
                # increase the index k
                k = k + 1
                # line 10
                (xv0, zv0) = P2(v0)
                (wv0, uv0) = D2(v0)
                # update Xbar and Tbar line 11
                push!(Xbar, xv0)
                tv0 = T'wv0
                push!(Tbar, tv0)
                Dstar_q = (wv0'*v0)[1] + zv0
                D = hcat(D, [tv0[1]; Dstar_q])
                # since zv0>ϵ we compute the huperplane in line 13
                if (zv0 > ϵ)
                    # H0Mat = [tv0[1:end-1][1] 1] * Tinv
                    H0Mat = (vcat(tv0[1:end-1], 1))' * Tinv
                    # if H = {y∈ Rq | aᵀ*y≤ b}, we construct it as
                    # Polyhedra.HalfSpace(a, b)
                    # if Hⱼ = {y∈ Rq | (zʲ)ᵗy≥ (zʲ)ᵗΓ(xⱼ)}
                    ## This halfspace is constructed as follows 
                    # The HalfSpace of Polyhedra package demands that 
                    #    Hⱼ = {y∈ Rq | ϕ(y,D*(Tᵀzʲ))≥0} = {y∈ Rq | ϕ(y,D*(Tᵀzʲ))≥0}
                    #    Hⱼ = {y∈ Rq | [(Tᵀzʲ)₁ (Tᵀzʲ)₂ ⋯ (Tᵀzʲ)q 1]*Tinv*y -D*(Tᵀzʲ)≥0}
                    #    Hⱼ = {y∈ Rq | [(Tᵀzʲ)₁ (Tᵀzʲ)₂ ⋯ (Tᵀzʲ)q 1]*Tinv*y≥D*(Tᵀzʲ)}
                    #    Hⱼ = {y∈ Rq | H0Mat*y≥Dstar_q}
                    # we write as 
                    #    Hⱼ = {y∈ Rⁿ | - (zⱼ)ᵗy≤ - (zⱼ)ᵗΓ(xⱼ)}
                    # ∴  Hⱼ = {y∈ Rq | -H0Mat*y≤ -Dstar_q}
                    # i.e aᵀ = -H0Mat ⇒ a = -H0Mat'
                    #     b = -Dstar_q
                    H0 = Polyhedra.HalfSpace(-H0Mat', -Dstar_q)
                    P_hr = P_poly.hrep ∩ H0
                    ## we compute P_{k+1} in line 18 and break the loop
                    P_poly = Polyhedra.polyhedron(P_hr)
                    go = true   
                end
            end
        end

    end
    return Xbar, Tbar, V, P_poly
end
    