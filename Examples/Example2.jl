using JuMP
import Ipopt
import NLopt

function P1(z)
    model = Model(Ipopt.Optimizer)
    JuMP.set_silent(model)
    JuMP.@variable(model, x1)
    JuMP.@variable(model, x2)
    JuMP.@NLexpression(model, nl_expr, z[1] * ((x1 - 3)^2 + (x2 - 1)^2) + z[2] * ((x1 - 1)^2 + (x2 - 1)^2))

    JuMP.@NLobjective(model, Min, nl_expr)
    JuMP.@constraint(model, c1, x1 + 2 * x2 <= 2)
    JuMP.@constraint(model, c2, x1 - 2 * x2 <= 2)
    JuMP.@constraint(model, c3, -x1 + 2 * x2 <= 2)
    JuMP.@constraint(model, c4, -x1 - 2 * x2 <= 2)
    JuMP.optimize!(model)
    return value.([x1, x2]), JuMP.objective_value(model)
end
function P2(v)
    model = Model(Ipopt.Optimizer)
    JuMP.set_silent(model)

    JuMP.@variable(model, z)
    JuMP.@variable(model, x1)
    JuMP.@variable(model, x2)
    JuMP.@objective(model, Min, z)
    JuMP.@constraint(model, c1, x1 + 2 * x2 <= 2)
    JuMP.@constraint(model, c2, x1 - 2 * x2 <= 2)
    JuMP.@constraint(model, c3, -x1 + 2 * x2 <= 2)
    JuMP.@constraint(model, c4, -x1 - 2 * x2 <= 2)
    # Gammax = Γ([x1, x2])
    # z[1] * ((x1 - 3)^2 + (x2 - 1)^2) + z[2] * ((x1 - 1)^2 + (x2 - 1)^2)
    JuMP.@NLconstraint(model, c5, (x1 - 3)^2 + (x2 - 1)^2 - z - v[1] <= 0)
    JuMP.@NLconstraint(model, c6, (x1 - 1)^2 + (x2 - 1)^2 - z - v[2] <= 0)
    JuMP.optimize!(model)
    return (value.([x1, x2]), value(z))
end
function dual_objective(u1, u2, u3, u4, w1, w2)
    model2 = Model(Ipopt.Optimizer)
    JuMP.set_silent(model2)
    JuMP.@variable(model2, x1)
    JuMP.@variable(model2, x2)
    JuMP.@NLexpression(model2, my_expr, w1 * (x1 - 3)^2 + (x2 - 1)^2 + w2 * (x1 - 1)^2 + (x2 - 1)^2)
    JuMP.@NLobjective(model2, Min,
        u1 * (x1 + 2 * x2 - 2) +
        u2 * (x1 - 2 * x2 - 2) +
        u3 * (-x1 + 2 * x2 - 2) +
        u4 * (-x1 - 2 * x2 - 2) +
        my_expr
    )
    JuMP.optimize!(model2)
    return JuMP.objective_value(model2)
end
# function D2(v)
#     model = Model(Ipopt.Optimizer)
#     @variable(model, u >= 0)
#     @variable(model, w1 >= 0)
#     @variable(model, w2 >= 0)
#     @NLobjective(
#         model,
#         Max,
#         (-(w1^2 + w2^2) / (4 * u) - u) + w1 * (1 - v[1]) + w2 * (1 - v[2])
#     )
#     @constraint(model, c1, w1 + w2 == 1)
#     optimize!(model)
#     return (value.([w1, w2]), value(u))
# end

# function D2(v)
#     model = Model(Ipopt.Optimizer)
#     @variable(model, u1 >= 0)
#     @variable(model, u2 >= 0)
#     @variable(model, u3 >= 0)
#     @variable(model, u4 >= 0)
#     @variable(model, w1 >= 0)
#     @variable(model, w2 >= 0)

#     @NLexpression(model, my_expr,
#         u3 * (u2 + u4 - u1 - u3 + (u1 + u2 - u3 - u4 - 2w2 - 6w1) / (2w1 + 2w2)) +
#         u1 * (u2 + u4 - u1 - u3 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2)) +
#         u4 * (-4 + u1 + u3 - u2 - u4 + (u1 + u2 - u3 - u4 - 2w2 - 6w1) / (2w1 + 2w2)) +
#         u2 * (-4 + u1 + u3 - u2 - u4 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2)) +
#         w1 * (-3 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2))^2 +
#         (1 / 2) * (-u1 + u2 - u3 + u4)^2 +
#         w2 * (-1 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2))^2)


#     @NLobjective(
#         model,
#         Max,
#         my_expr - w1 * v[1] - w2 * v[2]
#     )
#     @constraint(model, c1, w1 + w2 == 1)
#     optimize!(model)
#     return (value.([w1, w2]), value.([u1, u2, u3, u4]))
# end

function D2(v)
    model = Model(Ipopt.Optimizer)
    @variable(model, u1 >= 0)
    @variable(model, u2 >= 0)
    @variable(model, u3 >= 0)
    @variable(model, u4 >= 0)
    @variable(model, w1 >= 0)
    @variable(model, w2 >= 0)

    @NLexpression(model, my_expr,
        u3 * (u2 + u4 - u1 - u3 + (u1 + u2 - u3 - u4 - 2w2 - 6w1) / (2w1 + 2w2)) +
        u1 * (u2 + u4 - u1 - u3 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2)) +
        u4 * (-4 + u1 + u3 - u2 - u4 + (u1 + u2 - u3 - u4 - 2w2 - 6w1) / (2w1 + 2w2)) +
        u2 * (-4 + u1 + u3 - u2 - u4 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2)) +
        w1 * (-3 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2))^2 +
        (1 / 2) * (-u1 + u2 - u3 + u4)^2 +
        w2 * (-1 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2))^2)


    @NLobjective(
        model,
        Max,
        my_expr - w1 * v[1] - w2 * v[2]
    )
    @constraint(model, c1, w1 + w2 == 1)
    optimize!(model)
    return (value.([w1, w2]), value.([u1, u2, u3, u4]))
end
# function D2para(v)
#     model = Model(Ipopt.Optimizer)
#     @variable(model, u >= 0)
#     @variable(model, w1 >= 0)
#     @variable(model, w2 >= 0)
#     exprrr = ParameterD2()
#     @NLexpression(model,my_expr_obj,exprrr)
#     @NLobjective(
#         model,
#         Max,
#         exprrr - w1*v[1]-w2*v[2]
#     )
#     @constraint(model, c1, w1 + w2 == 1)
#     optimize!(model)
#     return (value.([w1, w2]), value(u))
# end


# function D2(v0::Vector)

#     function ExplicitFunction(uw::Vector, grad::Vector, v::Vector)
#         println("called my own OFunction")
#         if length(grad) > 0
#             grad[1] = 0
#             grad[2] = 0.5 / sqrt(x[2])
#         end
#         w1, w2, u1, u2, u3, u4 = uw
#         return u3 * (u2 + u4 - u1 - u3 + (u1 + u2 - u3 - u4 - 2w2 - 6w1) / (2w1 + 2w2)) +
#                u1 * (u2 + u4 - u1 - u3 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2)) +
#                u4 * (-4 + u1 + u3 - u2 - u4 + (u1 + u2 - u3 - u4 - 2w2 - 6w1) / (2w1 + 2w2)) +
#                u2 * (-4 + u1 + u3 - u2 - u4 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2)) +
#                w1 * (-3 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2))^2 +
#                (1 / 2) * (-u1 + u2 - u3 + u4)^2 +
#                w2 * (-1 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2))^2 - w1 * v[1] - w2 * v[2]
#     end
#     function ImplicitFunction(uw::Vector, grad::Vector, v::Vector)

#         if length(grad) > 0
#             grad[1] = 0
#             grad[2] = 0.5 / sqrt(x[2])
#         end
#         model2 = Model(Ipopt.Optimizer)
#         JuMP.set_silent(model2)
#         JuMP.@variable(model2, x1)
#         JuMP.@variable(model2, x2)
#         JuMP.@NLexpression(model2, my_expr, uw[1] * (x1 - 3)^2 + (x2 - 1)^2 + uw[2] * (x1 - 1)^2 + (x2 - 1)^2)
#         JuMP.@NLobjective(model2, Min,
#             uw[3] * (x1 + 2 * x2 - 2) +
#             uw[4] * (x1 - 2 * x2 - 2) +
#             uw[5] * (-x1 + 2 * x2 - 2) +
#             uw[6] * (-x1 - 2 * x2 - 2) +
#             my_expr
#         )
#         JuMP.optimize!(model2)
#         return JuMP.objective_value(model2) - uw[1] * v[1] - uw[2] * v[2]
#     end


#     function constraint(uw::Vector, _grad::Vector)
#         uw[1] + uw[2] - 1
#     end

#     opt = NLopt.Opt(:LN_COBYLA, 6) # LN_COBYLA LD_MMA
#     opt.lower_bounds = repeat([0], 6)
#     # opt.xtol_rel = 1e-2

#     opt.max_objective = (x, g) -> ImplicitFunction(x, g, v0)
#     NLopt.equality_constraint!(opt, constraint, 1e-8)

#     (minf, minx, ret) = NLopt.optimize(opt, repeat([0.9], 6))
#     numevals = opt.numevals # the number of function evaluations
#     w1, w2, u1, u2, u3, u4 = minx
#     return (value.([w1, w2]), value.([u1, u2, u3, u4]))
# end

function Γ(x)
    [(x[1] - 3)^2 + (x[2] - 1)^2; (x[1] - 1)^2 + (x[2] - 1)^2]
end




