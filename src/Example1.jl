module Example1
using JuMP
import Ipopt
export P1, P2, D2, Γ
function P1(z)
    model = Model(Ipopt.Optimizer)
    @variable(model, x1)
    @variable(model, x2)
    @objective(model, Min, z[1] * x1 + z[2] * x2)
    @NLconstraint(model, c1, (x1 - 1)^2 + (x2 - 1)^2 <= 1)
    optimize!(model)
    return value.([x1, x2]), objective_value(model)
end
function P2(v)
    model = Model(Ipopt.Optimizer)
    @variable(model, z)
    @variable(model, x1)
    @variable(model, x2)
    @objective(model, Min, z)
    @NLconstraint(model, c1, (x1 - 1)^2 + (x2 - 1)^2 <= 1)
    @constraint(model, c2, x1 - z - v[1] <= 0)
    @constraint(model, c3, x2 - z - v[2] <= 0)
    optimize!(model)
    return (value.([x1, x2]), value(z))
end
function D2(v)
    model = Model(Ipopt.Optimizer)
    @variable(model, u >= 0)
    @variable(model, w1 >= 0)
    @variable(model, w2 >= 0)
    @NLobjective(
        model,
        Max,
        (-(w1^2 + w2^2) / (4 * u) - u) + w1 * (1 - v[1]) + w2 * (1 - v[2])
    )
    @constraint(model, c1, w1 + w2 == 1)
    optimize!(model)
    return (value.([w1, w2]), value(u))
end
function Γ(x)
    x
end

end
