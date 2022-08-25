using JuMP
using CPLEX
# return true if ||u1(θ)|| > ||u2(θ)|| ∀θ s.t. A' θ ≤ b
function atomic_earlyterminate_qp(A,b,U1,U2;time_limit = 1)
    nth,n = size(U1).-(1,0);
    model = direct_model(CPLEX.Optimizer())
    set_silent(model)

    # Nonconvex QP settings
    set_optimizer_attribute(model, "CPXPARAM_OptimalityTarget", 3)
    # Terminate when we know for sure if the optimal value is 
    # > 0 or < 0
    set_optimizer_attribute(model, "CPXPARAM_MIP_Tolerances_UpperCutoff", 0)
    set_optimizer_attribute(model, "CPX_PARAM_INTSOLLIM", 1)
    set_optimizer_attribute(model, "CPX_PARAM_TILIM", time_limit)

    @variable(model,u1[1:n])
    @variable(model,u2[1:n])
    @variable(model,θ[1:nth])

    @constraint(model, u1.==U1[1:nth,:]'*θ+U1[end,:])
    @constraint(model, u2.==U2[1:nth,:]'*θ+U2[end,:])
    @constraint(model, A'*θ.<=b)

    @objective(model, Min,sum(u1.^2)-sum(u2.^2))

    optimize!(model)
    return termination_status(model)==MOI.INFEASIBLE
end
