using Printf

export solve_hierarchical, compare_solvers, save_comparison, print_comparison

function solution_metrics(p, x, reference, reference_solution)
    feasible_objective = objective(p, feasible_decisions(p, x[(1,1)]))
    feasible_gap = reference === nothing ? NaN : abs(feasible_objective-reference)/max(1,abs(reference))
    solution_error = NaN
    relative_solution_error = NaN
    if reference_solution !== nothing
        vertices = sort!(collect(keys(p.nodes)))
        delta = reduce(vcat, (x[v]-reference_solution[v] for v in vertices))
        ref = reduce(vcat, (reference_solution[v] for v in vertices))
        solution_error = maxnorm(delta)
        relative_solution_error = norm(delta)/max(1,norm(ref))
    end
    return (; feasible_objective, feasible_gap, solution_error, relative_solution_error)
end

"""Standalone physical-tree hADMM using the existing hierarchicalADMM! sweep."""
function solve_hierarchical(p::Problem; rho=1000.0, tol=1e-4, max_iter=1000,
                            reference=nothing, reference_solution=nothing)
    all(isfinite(t) && t > 0 for t in (rho,tol)) && max_iter > 0 ||
        throw(ArgumentError("Invalid standalone hADMM settings"))
    reset!(p.root)
    for n in values(p.nodes), u in values(n.dual)
        fill!(u, 0.0)
    end
    history = NamedTuple[]
    converged = false
    x = Dict(v => copy(vect_prime(n)) for (v,n) in p.nodes)
    for iteration in 1:max_iter
        previous = x
        residuals = Float64[]
        hierarchicalADMM!(p.root, residuals; λ=1/rho)
        x = Dict(v => copy(vect_prime(n)) for (v,n) in p.nodes)
        primal = maximum(maxnorm(e.A*x[e.parent]+e.B*x[e.child]-e.c) for e in p.edges)
        dual = rho*maximum(maxnorm(x[e.child]-previous[e.child]) for e in p.edges)
        # The legacy sweep counts decision vectors only. Add the scalar subtree
        # stopping report on each return edge, matching local hybrid hADMM.
        for e in p.edges
            com_cost!(p.nodes[e.child], 0.0, 1)
        end
        total, max_node = tt_com_iter(p.root)
        J = objective(p,x)
        gap = reference === nothing ? NaN : abs(J-reference)/max(1,abs(reference))
        push!(history, (; iteration, objective=J, objective_gap=gap, primal, dual,
                         physical=primal, TotalVals=total["com"], MaxNodeVals=max_node["com"],
                         TotalNodeUpdates=sum(n.iteration for n in values(p.nodes)),
                         MaxNodeUpdates=max_node["iter"],
                         solution_metrics(p,x,reference,reference_solution)...))
        if max(primal,dual) <= tol
            converged = true
            break
        end
    end
    h = last(history)
    return (; x, objective=objective(p,x), converged, iterations=length(history), history,
              TotalVals=h.TotalVals, MaxNodeVals=h.MaxNodeVals,
              node_values=Dict(v => n.com_cost for (v,n) in p.nodes),
              node_updates=Dict(v => n.iteration for (v,n) in p.nodes))
end

function comparison_row(p, name, result, reference, seconds, kind)
    gap = abs(result.objective-reference.objective)/max(1,abs(reference.objective))
    metrics = solution_metrics(p,result.x,reference.objective,reference.x)
    central = kind == :Ipopt
    physical = maximum(maxnorm(e.A*result.x[e.parent]+e.B*result.x[e.child]-e.c) for e in p.edges)
    return (; method=name, converged=result.converged, objective=result.objective,
              objective_gap=gap, metrics..., physical,
              primal=central ? missing : last(result.history).primal,
              dual=central ? missing : last(result.history).dual,
              iterations=result.iterations, iteration_kind=kind,
              TotalNodeUpdates=central ? missing : sum(values(result.node_updates)),
              MaxNodeUpdates=central ? missing : maximum(values(result.node_updates)),
              inner_sweeps=kind == :outer_ADMM ? sum(sum(h.inner_iterations) for h in result.history) : missing,
              TotalVals=central ? 0 : result.TotalVals,
              MaxNodeVals=central ? 0 : result.MaxNodeVals, seconds)
end

"""Compare the central reference, physical-tree hADMM, and requested hybrids.
Warm-up runs are excluded from wall time; times include model/solver setup and
diagnostics, and execute sequentially. The first config sets baseline rho/tol.
"""
function compare_solvers(p::Problem; configs=[SolverConfig(module_algorithms=(:hADMM,:fADMM,:hADMM))], warmup=true)
    isempty(configs) && throw(ArgumentError("At least one hybrid config is required"))
    cfg = first(configs)
    all(c.rho == cfg.rho && c.tol == cfg.tol && c.max_iter == cfg.max_iter for c in configs) ||
        throw(ArgumentError("Compared hybrids must share outer rho, tolerance, and budget"))
    warmup && centralized_original(p)
    central_time = @elapsed reference = centralized_original(p)
    reform = centralized_reformulation(p)
    isapprox(reference.objective,reform.objective;rtol=1e-8) || error("Central formulations disagree")
    run_h = () -> solve_hierarchical(p; rho=cfg.rho, tol=cfg.tol, max_iter=cfg.max_iter,
                                     reference=reference.objective, reference_solution=reference.x)
    warmup && run_h()
    h_time = @elapsed hierarchical = run_h()
    rows = NamedTuple[comparison_row(p,"Centralized",reference,reference,central_time,:Ipopt),
            comparison_row(p,"hADMM",hierarchical,reference,h_time,:tree_sweep)]
    runs = NamedTuple[(name="hADMM",result=hierarchical)]
    for c in configs
        name = "Hybrid[$(c.global_algorithm);$(join(c.module_algorithms,'/'))]"
        run_hybrid = () -> solve_hybrid(p; config=c, reference=reference.objective, reference_solution=reference.x)
        warmup && run_hybrid()
        seconds = @elapsed hybrid = run_hybrid()
        push!(rows,comparison_row(p,name,hybrid,reference,seconds,:outer_ADMM))
        push!(runs,(;name,result=hybrid))
    end
    return (; problem=p, reference, rows, runs, warmed_up=warmup)
end

function print_comparison(report)
    display_count(x) = ismissing(x) ? "N/A" : string(x)
    println("Method                                       Conv  Iters      RelGap     FeasGap      SolInf    Physical  TotalVals MaxNodeVals  NodeUpdates    MaxNode     Seconds")
    for r in report.rows
        @printf("%-44s %-5s %5d %11.3e %11.3e %11.3e %11.3e %10d %11d %12s %10s %11.4f\n",
                r.method,r.converged,r.iterations,r.objective_gap,r.feasible_gap,r.solution_error,
                r.physical,r.TotalVals,r.MaxNodeVals,display_count(r.TotalNodeUpdates),display_count(r.MaxNodeUpdates),r.seconds)
    end
    println("Iters: central=Ipopt iterations; hADMM=tree sweeps; hybrid=outer iterations.")
    println("NodeUpdates counts node proximal solves, including hybrid inner sweeps; central=N/A.")
    println("Communication counts optimization payloads and stopping reports; central data is already resident.")
    println("Wall time includes setup/diagnostics; warm-up excluded: ", report.warmed_up)
end

function write_csv(path, rows)
    isempty(rows) && return
    columns = propertynames(first(rows))
    cell(x) = ismissing(x) ? "" : "\""*replace(string(x),'"'=>"\"\"")*"\""
    open(path,"w") do io
        println(io,join(string.(columns),','))
        for row in rows
            println(io,join((cell(getproperty(row,c)) for c in columns),','))
        end
    end
end

"""Save summary, common convergence histories, and all physical node solutions."""
function save_comparison(report, output_dir)
    mkpath(output_dir)
    write_csv(joinpath(output_dir,"summary.csv"),report.rows)
    histories = [(; method=run.name, iteration=h.iteration, objective=h.objective,
                   objective_gap=h.objective_gap, feasible_gap=h.feasible_gap,
                   solution_error=h.solution_error, relative_solution_error=h.relative_solution_error,
                   primal=h.primal, dual=h.dual, physical=h.physical,
                   TotalVals=h.TotalVals, MaxNodeVals=h.MaxNodeVals,
                   TotalNodeUpdates=h.TotalNodeUpdates, MaxNodeUpdates=h.MaxNodeUpdates)
                 for run in report.runs for h in run.result.history]
    write_csv(joinpath(output_dir,"history.csv"),histories)
    solutions = NamedTuple[]
    all_runs = vcat([(name="Centralized",result=report.reference)],report.runs)
    for run in all_runs, v in sort!(collect(keys(report.problem.nodes))), j in eachindex(run.result.x[v])
        value = run.result.x[v][j]
        reference = report.reference.x[v][j]
        push!(solutions,(; method=run.name, depth=v[1], node=v[2], coordinate=j,
                           value, reference, absolute_error=abs(value-reference)))
    end
    write_csv(joinpath(output_dir,"solutions.csv"),solutions)
    return output_dir
end
