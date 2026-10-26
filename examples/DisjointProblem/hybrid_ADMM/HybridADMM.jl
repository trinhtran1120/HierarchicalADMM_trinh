module HybridADMM

using LinearAlgebra, JuMP, Ipopt
using ProximalAlgorithms, ProximalOperators, Zygote
using DifferentiationInterface: AutoZygote

# Isolate the legacy examples and their global defaults in this module. The hybrid
# drivers below pass their own penalties/tolerances rather than using those globals.
include("../TreeGeneration.jl")
include("../CentralSolution.jl")
include("../hADMM.jl")
include("../fADMM.jl")

export build_problem, SolverConfig, solve_hybrid, centralized_original,
       centralized_reformulation, formulation_residuals, feasible_decisions, solve_oracle

const Node = Tuple{Int,Int}

struct Edge
    parent::Node
    child::Node
    block::UnitRange{Int}
    A::Matrix{Float64}
    B::Matrix{Float64}
    c::Vector{Float64}
end

struct ModuleSpec
    name::Symbol
    vertices::Vector{Node}       # lexicographic (depth, index), root first
    edges::Vector{Edge}
    offsets::Dict{Node,UnitRange{Int}}
    D::Matrix{Float64}
    d::Vector{Float64}
    S::Matrix{Float64}
    blocks::Dict{Node,UnitRange{Int}} # node slices in the module-root vector
end

struct Problem
    root::linknode
    nodes::Dict{Node,linknode}
    edges::Vector{Edge}
    modules::Vector{ModuleSpec}
    boundary::Vector{Edge}
    blocks::Dict{Node,UnitRange{Int}} # node slices in the global-root vector
end

"""Candidate B's immutable tree, with the existing scalar-leaf quartic + L1 costs."""
function build_problem()
    children = Dict(
        (1,1) => [(2,1), (2,2)], (2,1) => [(3,1), (3,2), (3,3)],
        (2,2) => [(3,4), (3,5)], (3,4) => [(4,1), (4,2), (4,3)],
        (3,5) => [(4,4), (4,5)], (4,1) => [(5,1), (5,2)],
        (4,3) => [(5,3)],
    )
    vertices = sort!(unique(vcat(collect(keys(children)), reduce(vcat, values(children)))))
    nodes = Dict(v => linknode(string(i)) for (i, v) in enumerate(vertices))
    for v in sort!(collect(keys(children)))
        set_relative!(nodes[v], [nodes[w] for w in children[v]])
    end
    root = nodes[(1,1)]
    assign!(root)
    # get_postion follows the ordered physical child vectors, never Dict iteration.
    blocks = Dict{Node,UnitRange{Int}}()
    for v in vertices
        start = get_postion(get_path(nodes[v]))
        blocks[v] = start:start+nodes[v].nV-1
    end
    edges = Edge[]
    for v in sort!(collect(keys(children)))
        start = 1
        for w in children[v]
            n = nodes[w].nV
            block = start:start+n-1
            A = Matrix{Float64}(I, nodes[v].nV, nodes[v].nV)[block, :]
            push!(edges, Edge(v, w, block, A, -Matrix{Float64}(I, n, n), zeros(n)))
            start += n
        end
        @assert start == nodes[v].nV + 1
    end
    memberships = [
        [(2,1), (3,1), (3,2), (3,3)],
        [(3,4), (4,1), (4,2), (4,3), (5,1), (5,2), (5,3)],
        [(2,2), (3,5), (4,4), (4,5)],
    ]
    modules = ModuleSpec[]
    for (i, vs) in enumerate(memberships)
        internal = [e for e in edges if e.parent in vs && e.child in vs]
        offsets = Dict{Node,UnitRange{Int}}()
        start = 1
        for v in vs
            offsets[v] = start:start+nodes[v].nV-1
            start += nodes[v].nV
        end
        D = zeros(sum(length(e.c) for e in internal), start-1)
        row = 1
        for e in internal
            rows = row:row+length(e.c)-1
            D[rows, offsets[e.parent]] = e.A
            D[rows, offsets[e.child]] = e.B
            row += length(rows)
        end
        S = zeros(nodes[first(vs)].nV, start-1)
        S[:, offsets[first(vs)]] = Matrix{Float64}(I, size(S, 1), size(S, 1))
        local_blocks = Dict(v => (first(blocks[v])-first(blocks[first(vs)])+1):
                                (last(blocks[v])-first(blocks[first(vs)])+1) for v in vs)
        push!(modules, ModuleSpec(Symbol("M", i), vs, internal, offsets, D,
                                  zeros(size(D, 1)), S, local_blocks))
    end
    boundary = [e for e in edges if !any(e in m.edges for m in modules)]
    @assert length(vertices) == 16 && length(edges) == 15
    @assert length.(memberships) == [4, 7, 4]
    @assert length.(getfield.(modules, :edges)) == [3, 6, 3]
    @assert length(unique(reduce(vcat, memberships))) == 15
    @assert !any((1,1) in vs for vs in memberships)
    @assert [(e.parent, e.child) for e in boundary] ==
            [((1,1),(2,1)), ((1,1),(2,2)), ((2,2),(3,4))]
    return Problem(root, nodes, edges, modules, boundary, blocks)
end

cost(node, x) = node.cost_func.val(x; para=node.cost_func.para) + node.cost_func.w*norm(x, 1)
objective(p, x) = sum(cost(p.nodes[v], x[v]) for v in sort!(collect(keys(p.nodes))))
stack(m, x) = reduce(vcat, (x[v] for v in m.vertices))
maxnorm(x) = norm(x, Inf)

function feasible_decisions(p::Problem, leaf_vector)
    @assert length(leaf_vector) == p.root.nV
    return Dict(v => copy(leaf_vector[b]) for (v, b) in p.blocks)
end

"""Residuals of (9), kept separate to expose boundary-copy mistakes."""
function formulation_residuals(p::Problem, x, z)
    internal = maximum(maxnorm(m.D*stack(m, x) - m.d) for m in p.modules)
    consistency = maximum(maxnorm(m.S*stack(m, x) - z[i]) for (i,m) in enumerate(p.modules))
    boundary_values = Dict((1,1) => x[(1,1)])
    for (i,m) in enumerate(p.modules)
        boundary_values[first(m.vertices)] = z[i]
    end
    boundary = maximum(maxnorm(e.A*boundary_values[e.parent] +
                               e.B*boundary_values[e.child] - e.c) for e in p.boundary)
    physical = maximum(maxnorm(e.A*x[e.parent] + e.B*x[e.child] - e.c) for e in p.edges)
    return (; internal, consistency, boundary, physical)
end

function optimizer_model()
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "tol", 1e-10)
    return model
end

function checked_optimize!(model)
    optimize!(model)
    termination_status(model) in (JuMP.MOI.OPTIMAL, JuMP.MOI.LOCALLY_SOLVED) ||
        error("Central/oracle solve failed: $(termination_status(model))")
end

function model_cost!(model, node, x)
    t = @variable(model, [1:length(x)])
    @constraint(model, t .>= x)
    @constraint(model, t .>= -x)
    return node.cost_func.val(x; para=node.cost_func.para) + node.cost_func.w*sum(t)
end

"""Leaf-only reference using the model builders in CentralSolution.jl, with solver statistics."""
function centralized_original(p::Problem)
    model = optimizer_model()
    variables = Dict()
    setVariables!(p.root, variables, model)
    terms = NonlinearExpr[]
    getCostFunction!(p.root, variables, terms, model)
    @objective(model, Min, sum(terms))
    for v in all_variables(model)
        set_start_value(v, 1.0)
    end
    checked_optimize!(model)
    result = Dict()
    assign_result!(p.root, result, Dict(k => value.(v) for (k,v) in variables))
    x = Dict(v => result[n.ID] for (v,n) in p.nodes)
    return (; x, objective=objective(p,x), solver_objective=objective_value(model),
              iterations=JuMP.MOI.get(model, JuMP.MOI.BarrierIterations()),
              status=termination_status(model), converged=true)
end

"""Explicit stacked model (9), without eliminating any physical node copies."""
function centralized_reformulation(p::Problem)
    model = optimizer_model()
    x = Dict(v => @variable(model, [1:n.nV]) for (v,n) in p.nodes)
    z = [@variable(model, [1:p.nodes[first(m.vertices)].nV]) for m in p.modules]
    for (i,m) in enumerate(p.modules)
        X = stack(m, x)
        @constraint(model, m.D*X .== m.d)
        @constraint(model, m.S*X .== z[i])
    end
    bv = Dict((1,1) => x[(1,1)])
    for (i,m) in enumerate(p.modules)
        bv[first(m.vertices)] = z[i]
    end
    for e in p.boundary
        @constraint(model, e.A*bv[e.parent] + e.B*bv[e.child] .== e.c)
    end
    J = sum(model_cost!(model, p.nodes[v], x[v]) for v in sort!(collect(keys(p.nodes))))
    @objective(model, Min, J)
    checked_optimize!(model)
    return (; x=Dict(v => value.(y) for (v,y) in x), z=[value.(y) for y in z],
              objective=objective_value(model))
end

Base.@kwdef struct SolverConfig
    global_algorithm::Symbol = :fADMM
    module_algorithms::NTuple{3,Symbol} = (:hADMM, :hADMM, :hADMM) # M1,M2,M3
    rho::Float64 = 1000.0
    inner_rho::Float64 = 1000.0
    tol::Float64 = 1e-4
    inner_tol::Float64 = 1e-6
    max_iter::Int = 1000
    inner_max_iter::Int = 1000
    record_messages::Bool = false
end

struct Transmission
    source::Node
    destination::Node
    scalars::Int
    phase::Symbol
    module_name::Symbol
end

function transmit!(p, path, data, phase, name, messages, record)
    lookup = Dict(n.ID => v for (v,n) in p.nodes)
    for (a,b) in zip(path[1:end-1], path[2:end])
        (a.parent === b || b.parent === a) || error("Nonphysical transmission")
        record && push!(messages, Transmission(lookup[a.ID], lookup[b.ID], length(data), phase, name))
    end
    com_cost!(path, data, 1) # includes every forwarding sender
end

mutable struct ModuleState
    x::Dict{Node,Vector{Float64}}
    dual::Dict{Node,Vector{Float64}} # keyed by child; edge or flattened dual
end

function ModuleState(p, m)
    x = Dict(v => ones(p.nodes[v].nV) for v in m.vertices)
    dual = Dict(v => zeros(p.nodes[v].nV) for v in m.vertices[2:end])
    return ModuleState(x, dual)
end

# The legacy hADMM proximal interface handles a uniform quadratic. Split a
# positive diagonal quadratic into that term plus a convex diagonal remainder.
function proximal_node(node, a, b, x0)
    base = minimum(a)
    @assert base > 0
    cf = node.cost_func
    remainder = a .- base
    augmented = CostFunc(
        (x; para) -> cf.val(x; para=cf.para) + dot(x, remainder.*x)/2,
        (x; para) -> cf.grad(x; para=cf.para) + remainder.*x,
        cf.para,
    )
    augmented.w = cf.w
    f, g = Proximal_Iteration(augmented, b/base, 1/base)
    solution, _ = node.solver(f=f, g=g, x0=x0)
    return copy(solution)
end

function module_path(p, m, v)
    path = get_path(p.nodes[v])
    start = findfirst(n -> n === p.nodes[first(m.vertices)], path)
    return path[start:end]
end

function local_flat!(p, m, state, target, cfg, messages)
    r = first(m.vertices)
    beta = cfg.inner_rho
    residual = Inf
    for iteration in 1:cfg.inner_max_iter
        old = deepcopy(state.x)
        a = fill(cfg.rho, p.nodes[r].nV)
        b = cfg.rho*target
        for v in m.vertices[2:end]
            block = m.blocks[v]
            a[block] .+= beta
            b[block] .+= beta*(state.x[v] - state.dual[v])
        end
        state.x[r] = proximal_node(p.nodes[r], a, b, state.x[r])
        p.nodes[r].iteration += 1
        residual = beta*maxnorm(state.x[r] - old[r])
        for v in m.vertices[2:end]
            block = m.blocks[v]
            query = state.x[r][block] + state.dual[v]
            path = module_path(p, m, v)
            transmit!(p, path, query, :inner_query, m.name, messages, cfg.record_messages)
            # Reuse the existing flattened non-root proximal update verbatim.
            state.x[v] = copy(update_leaf_flat(p.nodes[v], query,
                ProximalOperators.NormL1(p.nodes[v].cost_func.w), state.x[v], 1/beta))
            delta = state.x[r][block] - state.x[v]
            state.dual[v] += delta
            residual = max(residual, maxnorm(delta), beta*maxnorm(state.x[v] - old[v]))
            p.nodes[v].iteration += 1
            transmit!(p, reverse(path), state.x[v], :inner_response, m.name, messages, cfg.record_messages)
        end
        residual <= cfg.inner_tol && return (; iterations=iteration, residual, converged=true)
    end
    return (; iterations=cfg.inner_max_iter, residual, converged=false)
end

function local_hierarchical!(p, m, state, target, cfg, messages)
    r = first(m.vertices)
    beta = cfg.inner_rho
    residual = Inf
    # Same preorder primal / return-time edge-dual schedule as hierarchicalADMM!.
    function sweep!(v, incoming)
        node = p.nodes[v]
        old = state.x[v]
        outgoing = [e for e in m.edges if e.parent == v]
        a = fill(v == r ? cfg.rho : beta, node.nV)
        b = v == r ? cfg.rho*target : beta*incoming
        for e in outgoing
            a[e.block] .+= beta
            b[e.block] .+= beta*(state.x[e.child] - state.dual[e.child])
        end
        state.x[v] = proximal_node(node, a, b, old)
        node.iteration += 1
        res = beta*maxnorm(state.x[v] - old)
        for e in outgoing
            query = state.x[v][e.block] + state.dual[e.child]
            path = [node, p.nodes[e.child]]
            transmit!(p, path, query, :inner_query, m.name, messages, cfg.record_messages)
            child_res = sweep!(e.child, query)
            delta = state.x[v][e.block] - state.x[e.child]
            state.dual[e.child] += delta
            # One extra scalar carries the child's subtree stopping residual.
            transmit!(p, reverse(path), vcat(state.x[e.child], child_res),
                      :inner_response, m.name, messages, cfg.record_messages)
            res = max(res, child_res, maxnorm(delta))
        end
        return res
    end
    for iteration in 1:cfg.inner_max_iter
        residual = sweep!(r, nothing)
        residual <= cfg.inner_tol && return (; iterations=iteration, residual, converged=true)
    end
    return (; iterations=cfg.inner_max_iter, residual, converged=false)
end

"""High-accuracy oracle for (12); used to validate distributed module updates."""
function module_oracle(p, m, target, rho)
    model = optimizer_model()
    X = @variable(model, [1:size(m.D, 2)])
    @constraint(model, m.D*X .== m.d)
    J = sum(model_cost!(model, p.nodes[v], X[m.offsets[v]]) for v in m.vertices)
    @objective(model, Min, J + rho/2*sum((m.S*X - target).^2))
    checked_optimize!(model)
    x = Dict(v => value.(X[m.offsets[v]]) for v in m.vertices)
    return (; x, objective=objective_value(model))
end

function global_update(p, roots, responses, u, rho, y)
    # Eliminating G gives z_m = T_m*y. T2 selects inside T3, so M2--M3's
    # original boundary equation is enforced in every global minimization.
    a = zeros(length(y))
    b = zeros(length(y))
    for i in 1:3
        block = p.blocks[roots[i]]
        a[block] .+= rho
        b[block] .+= rho*(responses[i] + u[i])
    end
    return proximal_node(p.root, a, b, y)
end

"""Reference (11)--(13) with centralized, high-accuracy module minimizations.
This diagnostic oracle has no distributed communication statistics.
"""
function solve_oracle(p::Problem; rho=1000.0, tol=1e-4, max_iter=1000)
    rho > 0 && tol > 0 && max_iter > 0 || throw(ArgumentError("Invalid oracle settings"))
    roots = [first(m.vertices) for m in p.modules]
    y = ones(p.root.nV)
    responses = [copy(y[p.blocks[r]]) for r in roots]
    u = [zeros(length(q)) for q in responses]
    history = NamedTuple[]
    x = feasible_decisions(p, y)
    z = copy(responses)
    converged = false
    for iteration in 1:max_iter
        y = global_update(p, roots, responses, u, rho, y)
        z = [copy(y[p.blocks[r]]) for r in roots]
        previous = responses
        x = Dict((1,1) => copy(y))
        for (i,m) in enumerate(p.modules)
            merge!(x, module_oracle(p, m, z[i]-u[i], rho).x)
        end
        responses = [copy(x[r]) for r in roots]
        primal = maximum(maxnorm(responses[i]-z[i]) for i in 1:3)
        dual = rho*maximum(maxnorm(responses[i]-previous[i]) for i in 1:3)
        for i in 1:3
            u[i] += responses[i]-z[i]
        end
        push!(history, (; iteration, objective=objective(p, x), primal, dual))
        if max(primal, dual) <= tol
            converged = true
            break
        end
    end
    return (; x, z, objective=objective(p, x), converged, iterations=length(history), history)
end

"""
Solve (9) with the guide's Z-first scaled ADMM. On this two-level computational
star, global fADMM and hADMM have identical two-block updates (see README.md).
Inner decisions/duals are warm-started; the configured inner budget is inexact.
"""
function solve_hybrid(p::Problem; config=SolverConfig(), reference=nothing, reference_solution=nothing)
    cfg = config
    cfg.global_algorithm in (:fADMM, :hADMM) || throw(ArgumentError("Unsupported global algorithm"))
    all(a in (:fADMM, :hADMM) for a in cfg.module_algorithms) ||
        throw(ArgumentError("Unsupported module algorithm; expected M1,M2,M3 choices"))
    all(isfinite(t) && t > 0 for t in (cfg.rho, cfg.inner_rho, cfg.tol, cfg.inner_tol)) ||
        throw(ArgumentError("Penalties and tolerances must be finite and positive"))
    min(cfg.max_iter, cfg.inner_max_iter) > 0 || throw(ArgumentError("Iteration budgets must be positive"))
    reset!(p.root)
    states = [ModuleState(p, m) for m in p.modules]
    roots = [first(m.vertices) for m in p.modules]
    u = [zeros(p.nodes[r].nV) for r in roots]
    y = ones(p.root.nV)
    z = [copy(y[p.blocks[r]]) for r in roots]
    messages = Transmission[]
    history = NamedTuple[]
    converged = false
    for iteration in 1:cfg.max_iter
        previous = [copy(s.x[r]) for (s,r) in zip(states, roots)]
        y = global_update(p, roots, previous, u, cfg.rho, y)
        p.root.iteration += 1
        z = [copy(y[p.blocks[r]]) for r in roots]
        reports = NamedTuple[]
        for (i,m) in enumerate(p.modules)
            path = get_path(p.nodes[roots[i]])
            target = z[i] - u[i]
            transmit!(p, path, target, :global_query, m.name, messages, cfg.record_messages)
            local_solve! = cfg.module_algorithms[i] == :fADMM ? local_flat! : local_hierarchical!
            report = local_solve!(p, m, states[i], target, cfg, messages)
            push!(reports, report)
            transmit!(p, reverse(path), vcat(states[i].x[roots[i]],
                      report.residual, Float64(report.iterations), Float64(report.converged)),
                      :global_response, m.name, messages, cfg.record_messages)
            u[i] += states[i].x[roots[i]] - z[i]
        end
        x = Dict((1,1) => copy(y))
        for state in states
            merge!(x, state.x)
        end
        primal = maximum(maxnorm(states[i].x[roots[i]] - z[i]) for i in 1:3)
        # Z is the first ADMM block; its dual residual uses changes in S_m X_m.
        dual = cfg.rho*maximum(maxnorm(states[i].x[roots[i]] - previous[i]) for i in 1:3)
        residuals = formulation_residuals(p, x, z)
        J = objective(p, x)
        gap = reference === nothing ? NaN : abs(J-reference)/max(1, abs(reference))
        total, max_node = tt_com_iter(p.root)
        push!(history, (; iteration, objective=J, objective_gap=gap, primal, dual,
                         internal=residuals.internal, physical=residuals.physical,
                         inner_residuals=getfield.(reports, :residual),
                         inner_iterations=getfield.(reports, :iterations),
                         inner_converged=getfield.(reports, :converged),
                         TotalVals=total["com"], MaxNodeVals=max_node["com"],
                         TotalNodeUpdates=sum(n.iteration for n in values(p.nodes)),
                         MaxNodeUpdates=max_node["iter"],
                         solution_metrics(p, x, reference, reference_solution)...))
        if max(primal, dual, residuals.internal) <= cfg.tol && all(r.converged for r in reports)
            converged = true
            break
        end
    end
    x = Dict((1,1) => copy(y))
    for state in states
        merge!(x, state.x)
    end
    # Populate the legacy node dictionaries so total_cost/print_tree still work.
    for (v,node) in p.nodes
        if node.children === nothing
            node.prime[node.ID] = copy(x[v])
        else
            for e in p.edges
                e.parent == v && (node.prime[p.nodes[e.child].ID] = copy(x[v][e.block]))
            end
        end
    end
    total, max_node = tt_com_iter(p.root)
    return (; x, z, objective=objective(p, x), converged, iterations=length(history),
              history, messages, TotalVals=total["com"], MaxNodeVals=max_node["com"],
              node_values=Dict(v => n.com_cost for (v,n) in p.nodes),
              node_updates=Dict(v => n.iteration for (v,n) in p.nodes),
              config=cfg, global_updates_equivalent=true)
end

include("Comparison.jl")

end # module
