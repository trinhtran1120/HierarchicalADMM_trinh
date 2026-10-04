module ScalarMixedADMM

include(joinpath(@__DIR__, "hybrid_mixed_admm.jl"))
using .HybridMixedInfrastructure

const Node = Tuple{Int, Int}
const coefficients = Dict(zip(
    [(d,i) for (d,n) in enumerate((1,2,5,5,3)) for i in 1:n],
    [(1.,1.), (1.1,-.5), (.9,2.), (1.2,0.), (1.3,1.25), (.8,-1.),
     (1.4,.75), (1.,-.25), (1.15,1.5), (.95,-1.25), (1.25,.5),
     (1.05,2.25), (.85,-.75), (1.35,1.), (.75,-1.5), (1.45,.25)]))

mutable struct State
    configuration
    x::Dict{Node,Float64}
    u::Dict{String,Float64}
    local_dual::Dict{Tuple{Node,Node},Float64}
    children::Dict{Node,Vector{Node}}
    tracker::MessageTracker
    events::Vector{NamedTuple}
    iteration::Int
end

function State()
    cfg = m123_configuration()
    inverse = Dict(id => v for (v,id) in cfg.physical_tree.node_ids)
    children = Dict(v => (n.children === nothing ? Node[] :
        [inverse[c.ID] for c in n.children]) for (v,n) in cfg.physical_tree.nodes)
    dual = Dict{Tuple{Node,Node},Float64}()
    for m in values(cfg.modules)
        if m.algorithm == "fadmm"
            for v in setdiff(m.physical_nodes, Set([m.module_root]))
                dual[(m.module_root,v)] = 0.
            end
        else
            for v in m.physical_nodes, c in children[v]
                c in m.physical_nodes && (dual[(v,c)] = 0.)
            end
        end
    end
    State(cfg, Dict(v=>0. for v in keys(coefficients)),
        Dict(id=>0. for id in keys(cfg.modules)), dual, children,
        MessageTracker(), NamedTuple[], 0)
end

function event!(s, kind, module_id, node)
    push!(s.events, (; iteration=s.iteration, kind, module_id, node))
end

function send!(s, origin, destination, kind)
    tree = s.configuration.physical_tree
    route_message!(s.tracker, tree.nodes[origin], tree.nodes[destination];
        iteration=s.iteration, message_type=kind, payload_dimension=1)
end

# Local duals use parent-minus-target orientation and rho_local = 1.
function flat_step!(s, m, query)
    r = m.module_root
    targets = sort(collect(setdiff(m.physical_nodes, Set([r]))))
    a,b = coefficients[r]
    s.x[r] = (a*b + query + sum(s.x[v]-s.local_dual[(r,v)] for v in targets)) /
             (a+1+length(targets))
    event!(s, :LOCAL_ROOT_UPDATE, m.id, r)
    for v in targets
        q = s.x[r]+s.local_dual[(r,v)]
        send!(s,r,v,:local_query)
        a,b = coefficients[v]
        s.x[v] = (a*b+q)/(a+1)
        event!(s,:PRIMAL,m.id,v)
        send!(s,v,r,:local_response)
        s.local_dual[(r,v)] += s.x[r]-s.x[v]
        event!(s,:DUAL,m.id,v)
    end
end

function hierarchical_step!(s, m, v, incoming)
    children = filter(c -> c in m.physical_nodes, s.children[v])
    a,b = coefficients[v]
    s.x[v] = (a*b + incoming + sum((s.x[c]-s.local_dual[(v,c)] for c in children); init=0.)) /
             (a+1+length(children))
    event!(s,:PRIMAL,m.id,v)
    for c in children
        query = s.x[v]+s.local_dual[(v,c)]
        send!(s,v,c,:local_query)
        hierarchical_step!(s,m,c,query)
        send!(s,c,v,:local_response)
        s.local_dual[(v,c)] += s.x[v]-s.x[c]
        event!(s,:DUAL,m.id,c)
    end
end

"""One global iteration, with rho_g = rho_h = rho_local_f = 1."""
function step!(s::State)
    s.iteration += 1
    cfg = s.configuration
    ids = cfg.algorithmic_module_graph["R"]
    a,b = coefficients[cfg.global_root]
    z = (a*b + sum(s.x[cfg.modules[id].module_root]+s.u[id] for id in ids))/(a+length(ids))
    s.x[cfg.global_root] = z
    event!(s,:GLOBAL_ROOT_UPDATE,"R",cfg.global_root)
    for id in ids
        m = cfg.modules[id]
        query = z-s.u[id]
        send!(s,cfg.global_root,m.module_root,:global_query)
        event!(s,:LOCAL_STEP,id,m.module_root)
        if m.algorithm == "fadmm"
            flat_step!(s,m,query)
        else
            hierarchical_step!(s,m,m.module_root,query)
        end
        send!(s,m.module_root,cfg.global_root,:global_response)
    end
    for id in ids
        s.u[id] += s.x[cfg.modules[id].module_root]-z
        event!(s,:GLOBAL_DUAL_UPDATE,id,cfg.global_root)
    end
    return s
end

"""All scalar consensus constraints imply a common value; solve it analytically."""
function scalar_reference()
    x_star = sum(a*b for (a,b) in values(coefficients))/sum(a for (a,b) in values(coefficients))
    J_star = sum(a*(x_star-b)^2/2 for (a,b) in values(coefficients))
    return (; x_star, J_star)
end

function diagnostics(s, reference)
    J = sum(a*(s.x[v]-b)^2/2 for (v,(a,b)) in coefficients)
    residual = maximum(abs(s.x[p]-s.x[c]) for (p,c) in keys(s.local_dual))
    residual = max(residual, maximum(abs(s.x[m.module_root]-s.x[(1,1)]) for m in values(s.configuration.modules)))
    metrics = communication_metrics(s.tracker)
    return (; iteration=s.iteration, J, objective_gap=J-reference.J_star,
        max_primal_error=maximum(abs(x-reference.x_star) for x in values(s.x)),
        max_primal_residual=residual, total_messages=metrics.TotalMessages,
        total_scalars=metrics.TotalScalarValues, relay_messages=metrics.RelayMessages)
end

"""Use the standalone solvers' max(primal residual, dual residual) < tol rule.

Dual residuals measure target-response changes times the corresponding penalty.
The scalar fixture uses unit penalties at both global and local levels.
"""
function run(; iterations=1000, tol=1e-4, io=stdout, output_dir=nothing)
    iterations > 0 || throw(ArgumentError("iterations must be positive"))
    isfinite(tol) && tol >= 0 || throw(ArgumentError("tol must be finite and nonnegative"))
    s = State()
    reference = scalar_reference()
    history = NamedTuple[]
    stop_reason = :iteration_limit
    for _ in 1:iterations
        previous_x = copy(s.x)
        step!(s)
        change = maximum(abs(s.x[v]-previous_x[v]) for v in keys(s.x))
        # Match standalone rho * (target_new - target_old), with rho = 1.
        global_dual_residual = maximum(abs(s.x[m.module_root]-previous_x[m.module_root])
            for m in values(s.configuration.modules))
        local_dual_residual = maximum(abs(s.x[target]-previous_x[target])
            for (_,target) in keys(s.local_dual))
        dual_residual = max(global_dual_residual, local_dual_residual)
        row = merge(diagnostics(s,reference),
            (max_primal_change=change, max_dual_residual=dual_residual,))
        push!(history, row)
        if max(row.max_primal_residual, row.max_dual_residual) < tol
            stop_reason = :tolerance
            break
        end
    end
    println(io,"hybridMixed stopped: ",stop_reason," after ",s.iteration,
            " iterations (limit=",iterations,", tolerance=",tol,")")
    println(io,"hybridMixed scalar problem: reference=",reference)
    isempty(history) || println(io,"Final scalar mixed diagnostics: ",last(history))
    println(io,"Separate scalar reference; objective gaps are not compared with vector-copy experiments. No convergence theorem is asserted.")
    if output_dir !== nothing
        mkpath(output_dir)
        open(joinpath(output_dir,"scalar_hybridMixed.csv"),"w") do file
            if !isempty(history)
                println(file,join(string.(keys(first(history))),','))
                for row in history
                    println(file,join(values(row),','))
                end
            end
        end
    end
    return (; state=s,reference,history,stop_reason)
end

end
