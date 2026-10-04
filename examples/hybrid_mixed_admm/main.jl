import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

include(joinpath(@__DIR__, "scalar_mixed.jl"))

module HybridComparisonExperiment

using Random, Statistics, LinearAlgebra, JuMP, Ipopt, Plots
using Printf
using ProximalAlgorithms, ProximalOperators
import Zygote
using DifferentiationInterface: AutoZygote

include(joinpath(@__DIR__, "..", "DisjointProblem", "TreeGeneration.jl"))
include(joinpath(@__DIR__, "..", "DisjointProblem", "CentralSolution.jl"))
include(joinpath(@__DIR__, "..", "DisjointProblem", "hADMM.jl"))
include(joinpath(@__DIR__, "..", "DisjointProblem", "fADMM.jl"))
include(joinpath(@__DIR__, "..", "hybrid_nested_admm", "hybrid_ADMM.jl"))
include(joinpath(@__DIR__, "reference.jl"))
include(joinpath(@__DIR__, "hybrid_mixed_admm.jl"))

# Fixed M123 physical topology.
const nN = 16
const nD = 5
const mode = 3  # standalone solvers stop on residual tolerance
# Preserve the comparison's existing penalties.
const λₙ = 1e-3
const λₛ = 1.5e-3
const λₕ = 1.2e-3
const max_iter = 1000
const Niter = max_iter
const opt_gap = 0.01/100
const tol = 1e-3

"""Copy the fixed M123 topology into the standalone solvers' linknode type."""
function comparison_tree()
    physical_tree = HybridMixedInfrastructure.m123_tree()
    nodes = Dict(id => linknode(id) for id in values(physical_tree.node_ids))
    for physical_node in values(physical_tree.nodes)
        if physical_node.children !== nothing
            set_relative!(nodes[physical_node.ID], [nodes[child.ID] for child in physical_node.children])
        end
    end
    return nodes[physical_tree.root.ID]
end

function print_metric_header(io)
    println(io, rpad("Metric", 49), " | ", lpad("Min", 14), " | ",
            lpad("Median", 14), " | ", lpad("Max", 14))
    println(io, repeat("-", 100))
end

function print_metric_row(io, label, values)
    @printf(io, "%-49s | %14.7g | %14.7g | %14.7g\n",
            label, minimum(values), median(values), maximum(values))
end

function edge_diagnostics(root, baseline, method)
    nodes = DisjointReference.physical_nodes(root)
    rho = Dict((node.ID, child.ID) =>
        1 / (method == "hybridNested" && node === root ? λₛ : λₕ)
        for node in nodes if node.children !== nothing for child in node.children)
    return DisjointReference.final_diagnostics(root, baseline;
        rho, node_primal = vect_prime, objective = total_cost)
end

# The standalone fADMM coordinator and target copies are supplied by its callback.
function flat_diagnostics(root, state, baseline)
    J = root.cost_func.val(state.z; para = root.cost_func.para) + root.cost_func.w * norm(state.z, 1)
    error = norm(state.z - baseline.x_star[root.ID], Inf)
    residual = 0.0
    for (id, node) in state.targets
        J += node.cost_func.val(state.w[id]; para = node.cost_func.para) + node.cost_func.w * norm(state.w[id], 1)
        error = max(error, norm(state.w[id] - baseline.x_star[id], Inf))
        residual = max(residual, norm(state.blocks[id] - state.w[id], Inf))
    end
    return (; max_primal_error = error, J, objective_gap = J - baseline.J_star,
        max_primal_residual = residual, edge_dual_error = NaN)
end

function run_method(initial_root, baseline, method, iterations)
    root = deepcopy(initial_root)
    reset!(root)
    history = NamedTuple[]
    function record(state)
        diagnostics = method == "fADMM" ? flat_diagnostics(root, state, baseline) :
            edge_diagnostics(root, baseline, method)
        total, maxima = tt_com_iter(root)
        updates = sum(node.iteration for node in DisjointReference.physical_nodes(root))
        inner = method == "hybridNested" ? state.inner_iterations : 0
        push!(history, merge((iteration = state.iteration,), diagnostics,
            (total_scalars = total["com"], max_node_scalars = maxima["com"],
             max_node_updates = maxima["iter"], total_node_updates = updates,
             inner_iterations = inner, stopping_residual = state.stopping_residual)))
    end
    if method == "hADMM"
        hADMM(root, baseline.x_star; tol, max_iter = iterations,
            verbose = true, on_iteration = record)
    elseif method == "fADMM"
        flattenADMM(root; tol, max_iter = iterations,
            verbose = true, on_iteration = record)
    elseif method == "hybridNested"
        hybridADMM(root, baseline.x_star; tol, inner_tol = tol,
            max_iter = iterations, inner_max_iter = max_iter,
            verbose = true, on_iteration = record)
    else
        throw(ArgumentError("Unknown comparison method: $method"))
    end
    return history
end

function print_comparison(results, data_dir)
    methods = ("centralSolution", "fADMM", "hADMM", "hybridNested")
    rows = NamedTuple[]
    for result in results, method in methods
        baseline = result.baseline.J_star
        final = get(result.histories, method, nothing)
        central = method == "centralSolution"
        J = central ? baseline : final === nothing ? missing : last(final).J
        gap = ismissing(J) ? missing : J - baseline
        push!(rows, (
            topology = result.topology, seed = result.seed, method,
            status = central ? "baseline" : final === nothing ? "unavailable" : "completed",
            J_star = baseline, J, objective_gap = gap,
            absolute_relative_gap_percent = ismissing(gap) ? missing : 100 * abs(gap) / abs(baseline),
            total_scalars = final === nothing ? missing : last(final).total_scalars,
            max_node_scalars = final === nothing ? missing : last(final).max_node_scalars,
        ))
    end
    open(joinpath(data_dir, "comparison.csv"), "w") do io
        println(io, join(string.(keys(first(rows))), ','))
        for row in rows
            println(io, join((ismissing(value) ? "" : value for value in values(row)), ','))
        end
    end
    println("\nVector-copy comparison (medians across shared topologies)")
    println("Method | Status | J | Signed J-J* | Absolute relative gap (%) | Total scalars sent | Max scalars sent by a node")
    for method in methods
        method_rows = filter(row -> row.method == method, rows)
        columns = [all(ismissing, getproperty.(method_rows, field)) ? "N/A" :
                   string(median(getproperty.(method_rows, field)))
                   for field in (:J, :objective_gap, :absolute_relative_gap_percent,
                                 :total_scalars, :max_node_scalars)]
        println(join((method, first(method_rows).status, columns...), " | "))
    end
    println("Central solution communication is not modeled; hybridMixed is reported separately for its scalar problem.")
    println("Communication columns count sent scalar values, not physical messages.")
end

function save_plots(results, directory)
    methods = ("hADMM", "fADMM", "hybridNested")
    error_plot = plot(xlabel = "Outer iteration", ylabel = "Mean max |x-x*|")
    gap_plot = plot(xlabel = "Outer iteration", ylabel = "Mean signed J(x)-J*")
    physical_plot = plot(xlabel = "Outer iteration", ylabel = "Mean physical-edge residual")
    consensus_plot = plot(xlabel = "Outer iteration", ylabel = "Mean star consensus residual")
    communication_plot = plot(xlabel = "Mean cumulative scalars sent", ylabel = "Mean signed J(x)-J*")
    for method in methods
        histories = [result.histories[method] for result in results]
        common_length = minimum(length, histories)
        average(field) = mean([getproperty.(history[1:common_length], field) for history in histories])
        k = getproperty.(first(histories)[1:common_length], :iteration)
        plot!(error_plot, k, average(:max_primal_error); label = method)
        plot!(gap_plot, k, average(:objective_gap); label = method)
        plot!(method == "fADMM" ? consensus_plot : physical_plot,
              k, average(:max_primal_residual); label = method)
        plot!(communication_plot, average(:total_scalars), average(:objective_gap); label = method)
    end
    savefig(error_plot, joinpath(directory, "primal_error.pdf"))
    savefig(gap_plot, joinpath(directory, "objective_gap.pdf"))
    savefig(plot(physical_plot, consensus_plot; layout = (2, 1), size = (800, 700)),
            joinpath(directory, "residuals.pdf"))
    savefig(communication_plot, joinpath(directory, "communication.pdf"))
end

"""Run comparisons on the fixed M123 tree; n_topologies repeats that same topology."""
function run_comparison(; n_topologies = 1, iterations = max_iter, seed = 0,
                        output_root = normpath(joinpath(@__DIR__, "..", "..")))
    data_dir = joinpath(output_root, "data", "hybrid-mixed-admm")
    figure_dir = joinpath(output_root, "media", "figs", "hybrid_mixed_admm")
    mkpath(data_dir)
    mkpath(figure_dir)
    println("hybridMixed uses a separate scalar problem and reference, reported separately.")
    println("Disjoint problem on fixed M123 tree: $nN nodes, $nD levels, repeats=$n_topologies, K=$iterations, seed=$seed")
    println("lambda_h=$λₕ, lambda_global=$λₛ; initial primals and scaled duals = 1")
    println("hybridNested: global fADMM over 2 modules rooted at (2,1) and (2,2); local hADMM, inner limit=$max_iter")
    results = NamedTuple[]
    open(joinpath(data_dir, "status.txt"), "w") do io
        println(io, "hybridMixed is reported separately in scalar_hybridMixed.csv against its scalar reference.")
        println(io, "nodes=$nN levels=$nD K=$iterations repeats=$n_topologies seed=$seed lambda_h=$λₕ lambda_global=$λₛ")
        println(io, "Fixed M123 physical topology; comparison penalties are specified explicitly in main.jl.")
        println(io, "fADMM residual is star consensus; its physical edge-dual error is unavailable (NaN).")
        println(io, "Objective uses each method's actual node variables; gaps are signed and unsmoothed.")
        println(io, "Legacy communication counters measure sent scalars only; no message counts are inferred.")
    end
    open(joinpath(data_dir, "topologies.csv"), "w") do topology_io
        println(topology_io, "topology,seed,node_id,parent_id,depth,dimension")
        for topology in 1:n_topologies
            topology_seed = seed + topology - 1
            Random.seed!(topology_seed)
            root = comparison_tree()
            assign!(root)
            nodes = DisjointReference.physical_nodes(root)
            for node in nodes
                parent = node.parent === nothing ? "" : node.parent.ID
                println(topology_io, join((topology, topology_seed, node.ID, parent,
                                          length(get_path(node)), node.nV), ','))
            end
            println("Topology $topology/$n_topologies: $(length(nodes)) nodes; centralized baseline")
            baseline = DisjointReference.centralized_reference(root; solve_primal = get_CenVars)
            println("  Central solution J* = $(baseline.J_star)")
            histories = Dict()
            for method in ("hADMM", "fADMM", "hybridNested")
                println("  Running $method")
                flush(stdout)
                histories[method] = run_method(root, baseline, method, iterations)
                final = last(histories[method])
                relative_gap_percent = 100 * abs(final.objective_gap) / abs(baseline.J_star)
                println("  $method: max|x-x*|=$(final.max_primal_error), J=$(final.J), J-J*=$(final.objective_gap), absolute relative optimality gap=$(relative_gap_percent)%, residual=$(final.max_primal_residual), edge-dual error=$(final.edge_dual_error)")
                open(joinpath(data_dir, "topology_$(topology)_$(method).csv"), "w") do io
                    println(io, join(string.(keys(final)), ','))
                    for row in histories[method]
                        println(io, join(values(row), ','))
                    end
                end
            end
            push!(results, (; topology, seed = topology_seed, baseline, histories))
        end
    end
    open(joinpath(data_dir, "summary.csv"), "w") do io
        fields = keys(last(first(results).histories["hADMM"]))
        println(io, "topology,seed,method,J_star,", join(string.(fields), ','))
        for result in results, method in ("hADMM", "fADMM", "hybridNested")
            println(io, join((result.topology, result.seed, method, result.baseline.J_star,
                             values(last(result.histories[method]))...), ','))
        end
    end
    print_comparison(results, data_dir)
    println("\nFinal metrics across topologies (min / median / max)")
    for method in ("hADMM", "fADMM", "hybridNested")
        println("\n$method")
        print_metric_header(stdout)
        finals = [last(result.histories[method]) for result in results]
        for (field, label) in (
            (:max_node_updates, "Maximum number of iterations in a node"),
            (:max_node_scalars, "Maximum number of scalar values sent by a node"),
            (:total_scalars, "Total scalar values sent in the network"),
            (:total_node_updates, "Total node updates"),
            (:max_primal_error, "Final max |x-x*|"),
            (:J, "Final J(x^K)"),
            (:objective_gap, "Final signed J(x^K)-J*"),
            (:max_primal_residual, method == "fADMM" ?
                "Final star consensus primal residual" : "Final physical-edge primal residual"),
        )
            metric_values = getproperty.(finals, field)
            print_metric_row(stdout, label, metric_values)
        end
        relative_gaps = [100 * abs(last(result.histories[method]).objective_gap) /
                         abs(result.baseline.J_star) for result in results]
        print_metric_row(stdout, "Final absolute relative optimality gap (%)", relative_gaps)
        if method != "fADMM"
            errors = getproperty.(finals, :edge_dual_error)
            print_metric_row(stdout, "Final unscaled edge-dual error", errors)
        end
        if method == "hybridNested"
            inner = [sum(row.inner_iterations for row in result.histories[method]) for result in results]
            print_metric_row(stdout, "Total inner module iterations", inner)
        end
    end
    save_plots(results, figure_dir)
    println("Saved results to $data_dir and plots to $figure_dir")
    return results
end

end

module HybridMixedExample

include(joinpath(@__DIR__, "hybrid_mixed_admm.jl"))
using .HybridMixedInfrastructure

"""Report the fixed M123 configuration and a scalar routing demonstration.

The payload dimension is a communication test assumption, not the undefined
module-level numerical query/response interface. No ADMM updates are executed.
"""
function run_example(; io = stdout)
    configuration = m123_configuration()
    tree = configuration.physical_tree
    indices = sort(collect(keys(tree.nodes)))
    index_by_id = Dict(id => index for (index, id) in tree.node_ids)
    tracker = MessageTracker()

    println(io, "M123 structural example: global ", configuration.global_algorithm)
    println(io, "Physical tree (16 nodes, 15 edges):")
    for index in indices
        node = tree.nodes[index]
        if node.children !== nothing
            for child in node.children
                println(io, "  ", index, " -> ", index_by_id[child.ID])
            end
        end
    end

    println(io, "\nPeer modules; R = ", configuration.global_root, " belongs to no module:")
    for id in configuration.algorithmic_module_graph["R"]
        local_module = configuration.modules[id]
        println(io, "  R -> ", id, "; algorithm=", local_module.algorithm,
                "; root=", local_module.module_root,
                "; members=", sort(collect(local_module.physical_nodes)))
        destination = tree.nodes[local_module.module_root]
        path = physical_path(tree.root, destination)
        println(io, "  Physical route: ", join([index_by_id[node.ID] for node in path], " -> "))
        route_message!(tracker, tree.root, destination;
            iteration = 1, message_type = :query, payload_dimension = 1)
        route_message!(tracker, destination, tree.root;
            iteration = 1, message_type = :response, payload_dimension = 1)
    end

    println(io, "\nRouting demonstration: one scalar query/response per module; no ADMM iteration.")
    println(io, "iteration | logical origin -> destination | physical sender -> receiver | type | dimension | relay")
    for send in tracker.sends
        println(io, send.iteration, " | ", index_by_id[send.logical_origin], " -> ",
                index_by_id[send.logical_destination], " | ", index_by_id[send.physical_sender],
                " -> ", index_by_id[send.physical_receiver], " | ", send.message_type,
                " | ", send.payload_dimension, " | ", send.is_relay)
    end
    metrics = communication_metrics(tracker)
    println(io, "\nTotalMessages=", metrics.TotalMessages,
            "; TotalScalarValues=", metrics.TotalScalarValues,
            "; MaxNodeMessages=", metrics.MaxNodeMessages,
            "; MaxNodeScalarValues=", metrics.MaxNodeScalarValues,
            "; RelayMessages=", metrics.RelayMessages)
    println(io, "Node | sent messages | sent scalars")
    for index in indices
        id = tree.node_ids[index]
        println(io, index, " | ", get(metrics.NodeSentMessages, id, 0),
                " | ", get(metrics.NodeSentScalars, id, 0))
    end
    println(io, "\n(2,2) only relays R/M2 traffic; forwarding triggers no M3 optimization update.")
    println(io, "Local numerical steps and pure-mode regressions: not run.")
    println(io, "The scalar mixed numerical experiment runs separately below.")
    println(io, "See skills/formulation_M123.md Section 10. No mixed correctness or convergence claim.")
    return (; configuration, tracker, metrics)
end

end

"""Print mixed results against their own scalar reference after the comparison report."""
function print_mixed_summary(results; io=stdout)
    rows = map(results) do result
        final = last(result.history)
        state = result.state
        metrics = ScalarMixedADMM.HybridMixedInfrastructure.communication_metrics(state.tracker)
        updates = Dict(node => 0 for node in keys(state.x))
        for event in state.events
            if event.kind in (:GLOBAL_ROOT_UPDATE, :LOCAL_ROOT_UPDATE, :PRIMAL)
                updates[event.node] += 1
            end
        end
        merge(final, (
            max_node_updates=maximum(values(updates)),
            max_node_scalars=metrics.MaxNodeScalarValues,
            total_node_updates=sum(values(updates)),
            x_star=result.reference.x_star, J_star=result.reference.J_star,
            relative_gap=100 * abs(final.objective_gap) / abs(result.reference.J_star)))
    end
    println(io, "\nhybridMixed — scalar final metrics across ", length(results), " run(s) (separate reference problem)")
    println(io, "Stop reasons: ", join(string.(getproperty.(results, :stop_reason)), ", "))
    HybridComparisonExperiment.print_metric_header(io)
    for (field, label) in (
        (:iteration, "Global iterations"),
        (:max_node_updates, "Maximum number of iterations in a node"),
        (:max_node_scalars, "Maximum scalar values sent by a node"),
        (:total_scalars, "Total scalar values sent in the network"),
        (:total_messages, "Total physical messages"),
        (:relay_messages, "Relay transmissions"),
        (:total_node_updates, "Total node updates"),
        (:x_star, "Scalar reference x*"),
        (:J_star, "Scalar reference J*"),
        (:max_primal_error, "Final max |x-x*|"),
        (:J, "Final J(x^K)"),
        (:objective_gap, "Final signed J(x^K)-J*"),
        (:max_primal_residual, "Final consensus primal residual"),
        (:max_dual_residual, "Final dual residual"),
        (:max_primal_change, "Final maximum primal change"),
        (:relative_gap, "Final absolute relative optimality gap (%)"),
    )
        HybridComparisonExperiment.print_metric_row(io, label, getproperty.(rows, field))
    end
    println(io, "History: data/hybrid-mixed-admm/scalar_hybridMixed.csv")
    println(io, "Scalar objective and gap are not comparable to the vector-copy results above.")
end

# Execute both from the command line and from the IDE's include-based runner.
HybridMixedExample.run_example()
mixed_result = ScalarMixedADMM.run(; iterations=HybridComparisonExperiment.max_iter,
    tol=HybridComparisonExperiment.tol,
    output_dir=joinpath(@__DIR__, "..", "..", "data", "hybrid-mixed-admm"))
HybridComparisonExperiment.run_comparison()
print_mixed_summary([mixed_result])
nothing
