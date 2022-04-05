using JuMP
import NLopt
using Ipopt


function D2(v0::Vector)

    function ExplicitFunction(uw::Vector, grad::Vector, v::Vector)
        println("called my own OFunction")
        if length(grad) > 0
            grad[1] = 0
            grad[2] = 0.5 / sqrt(x[2])
        end
        w1, w2, u1, u2, u3, u4 = uw
        return u3 * (u2 + u4 - u1 - u3 + (u1 + u2 - u3 - u4 - 2w2 - 6w1) / (2w1 + 2w2)) +
               u1 * (u2 + u4 - u1 - u3 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2)) +
               u4 * (-4 + u1 + u3 - u2 - u4 + (u1 + u2 - u3 - u4 - 2w2 - 6w1) / (2w1 + 2w2)) +
               u2 * (-4 + u1 + u3 - u2 - u4 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2)) +
               w1 * (-3 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2))^2 +
               (1 / 2) * (-u1 + u2 - u3 + u4)^2 +
               w2 * (-1 + (u3 + u4 - u1 - u2 + 2w2 + 6w1) / (2w1 + 2w2))^2 - w1 * v[1] - w2 * v[2]
    end
    function ImplicitFunction(uw::Vector, grad::Vector, v::Vector)

        if length(grad) > 0
            grad[1] = 0
            grad[2] = 0.5 / sqrt(x[2])
        end
        model2 = Model(Ipopt.Optimizer)
        JuMP.set_silent(model2)
        JuMP.@variable(model2, x1)
        JuMP.@variable(model2, x2)
        JuMP.@NLexpression(model2, my_expr, uw[1] * (x1 - 3)^2 + (x2 - 1)^2 + uw[2] * (x1 - 1)^2 + (x2 - 1)^2)
        JuMP.@NLobjective(model2, Min,
            uw[3] * (x1 + 2 * x2 - 2) +
            uw[4] * (x1 - 2 * x2 - 2) +
            uw[5] * (-x1 + 2 * x2 - 2) +
            uw[6] * (-x1 - 2 * x2 - 2) +
            my_expr
        )
        JuMP.optimize!(model2)
        return JuMP.objective_value(model2) - uw[1] * v[1] - uw[2] * v[2]
    end


    function constraint(uw::Vector, _grad::Vector)
        uw[1] + uw[2] - 1
    end

    opt = NLopt.Opt(:LN_COBYLA, 6) # LN_COBYLA LD_MMA
    opt.lower_bounds = repeat([0], 6)
    # opt.xtol_rel = 1e-2

    opt.max_objective = (x, g) -> ImplicitFunction(x, g, v0)
    NLopt.equality_constraint!(opt, constraint, 1e-8)

    (minf, minx, ret) = NLopt.optimize(opt, repeat([0.9], 6))
    numevals = opt.numevals # the number of function evaluations
    w1, w2, u1, u2, u3, u4 = minx
    return (value.([w1, w2]), value.([u1, u2, u3, u4]))
end

