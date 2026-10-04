module DisjointReference

using LinearAlgebra, JuMP, Ipopt
using ProximalAlgorithms, ProximalOperators

include(joinpath(@__DIR__, "..", "DisjointProblem", "TreeGeneration.jl"))
include(joinpath(@__DIR__, "..", "DisjointProblem", "CentralSolution.jl"))
include(joinpath(@__DIR__, "..", "DisjointProblem", "hADMM.jl"))

export centralized_reference, final_diagnostics, report_reference

# The standalone hADMM sweep reads this global proximal parameter.
const λₕ = 1.0

function physical_nodes(root)
    nodes = [root]
    if root.children !== nothing
        for child in root.children
            append!(nodes, physical_nodes(child))
        end
    end
    return nodes
end

"""Reuse CentralSolution.jl for primals/objective; solve explicit edge constraints for duals."""
function centralized_reference(root; solve_primal = get_CenVars)
    x_star, J_star = solve_primal(root)
    nodes = physical_nodes(root)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "tol", 1e-10)
    x = Dict(node.ID => @variable(model, [1:node.nV]) for node in nodes)
    constraints = Dict()
    terms = NonlinearExpr[]
    for node in nodes
        local_x = x[node.ID]
        slack = @variable(model, [1:node.nV])
        @constraint(model, slack .>= local_x)
        @constraint(model, slack .>= -local_x)
        push!(terms, node.cost_func.val(local_x; para = node.cost_func.para) +
                     node.cost_func.w * sum(slack))
        if node.children !== nothing
            offset = 1
            for child in node.children
                block = offset:(offset + child.nV - 1)
                constraints[(node.ID, child.ID)] =
                    @constraint(model, local_x[block] .- x[child.ID] .== 0)
                offset += child.nV
            end
        end
    end
    @objective(model, Min, sum(terms))
    optimize!(model)
    is_solved_and_feasible(model; dual = true) ||
        error("Centralized edge-dual model failed: $(termination_status(model))")
    # JuMP equality duals use the opposite sign to the +y'*(parent block-child) Lagrangian.
    y_star = Dict(edge => -dual.(constraint) for (edge, constraint) in constraints)
    return (; x_star, J_star, y_star)
end

"""Report diagnostics for physical disjoint edge state, with unscaled y = rho * stored dual.

The objective is evaluated at each node's own final variable, which may be infeasible.
The signed objective gap is retained. No convergence or monotonicity assertion is made.
"""
function final_diagnostics(root, baseline; rho, node_primal = vect_prime, objective = total_cost)
    nodes = physical_nodes(root)
    primal_error = maximum(norm(node_primal(node) - baseline.x_star[node.ID], Inf) for node in nodes)
    primal_residual = 0.0
    edge_dual_error = 0.0
    for node in nodes
        if node.children !== nothing
            for child in node.children
                primal_residual = max(primal_residual,
                    norm(node.prime[child.ID] - node_primal(child), Inf))
                edge = (node.ID, child.ID)
                edge_rho = rho isa Number ? rho : rho[edge]
                edge_dual_error = max(edge_dual_error,
                    norm(edge_rho * node.dual[child.ID] - baseline.y_star[edge], Inf))
            end
        end
    end
    J = objective(root)
    return (; max_primal_error = primal_error, J, objective_gap = J - baseline.J_star,
            max_primal_residual = primal_residual, edge_dual_error)
end

"""Run a standalone hADMM reporting example on the six-node disjoint topology."""
function report_reference(; io = stdout)
    nodes = [linknode(string(i)) for i in 1:6]
    set_relative!(nodes[1], nodes[2:3])
    set_relative!(nodes[2], nodes[4:5])
    set_relative!(nodes[3], [nodes[6]])
    root = nodes[1]
    assign!(root)
    baseline = centralized_reference(root)
    rho = 1 / λₕ
    println(io, "Standalone disjoint hADMM; rho = ", rho, "; J* = ", baseline.J_star)
    println(io, "K | max |x-x*| | J(x^K) | J(x^K)-J* | max primal residual | edge-dual error")
    for k in 1:2000
        hierarchicalADMM!(root, Float64[])
        if k in (10, 2000)
            metrics = final_diagnostics(root, baseline; rho)
            println(io, k, " | ", join(values(metrics), " | "))
        end
    end
    return baseline
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    DisjointReference.report_reference()
end
