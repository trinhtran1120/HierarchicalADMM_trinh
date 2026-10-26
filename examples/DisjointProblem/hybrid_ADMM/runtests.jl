using Test, LinearAlgebra

include("HybridADMM.jl")
using .HybridADMM
const H = HybridADMM

@testset "Candidate B hybrid ADMM" begin
    p = build_problem()
    @testset "Physical tree, module partition and operators" begin
        @test length(p.nodes) == 16
        @test length(p.edges) == 15
        @test count(n.parent === nothing for n in values(p.nodes)) == 1
        @test all(count(e.child == v for e in p.edges) == 1 for v in keys(p.nodes) if v != (1,1))
        @test length.(getfield.(p.modules, :vertices)) == [4,7,4]
        @test length.(getfield.(p.modules, :edges)) == [3,6,3]
        @test [(e.parent,e.child) for e in p.boundary] ==
              [((1,1),(2,1)), ((1,1),(2,2)), ((2,2),(3,4))]
        @test sort(reduce(vcat, getfield.(p.modules, :vertices))) == sort(setdiff(collect(keys(p.nodes)), [(1,1)]))
        for e in p.edges
            @test size(e.A) == (p.nodes[e.child].nV,p.nodes[e.parent].nV)
            @test e.B == -Matrix{Float64}(I, length(e.c), length(e.c))
            probe = collect(1.0:p.nodes[e.parent].nV)
            @test e.A*probe == probe[e.block]
        end
        x = feasible_decisions(p, collect(1.0:p.root.nV))
        z = [copy(x[first(m.vertices)]) for m in p.modules]
        for (i,m) in enumerate(p.modules)
            @test m.D*H.stack(m,x) == m.d
            @test m.S*H.stack(m,x) == z[i]
        end
        @test all(==(0), formulation_residuals(p,x,z))
        z[2][1] += 1
        res = formulation_residuals(p,x,z)
        @test res.consistency == 1
        @test res.boundary == 1
        @test res.internal == 0
        # M3's M2 block is free internally, while its M3-child block is constrained.
        m = p.modules[3]
        x = feasible_decisions(p, ones(p.root.nV))
        x[(2,2)][1] += 1
        @test H.maxnorm(m.D*H.stack(m,x)) == 0
        x[(2,2)][end] += 1
        @test H.maxnorm(m.D*H.stack(m,x)) == 1
    end

    reference = centralized_original(p)
    @testset "Independent centralized formulations and exact-block oracle" begin
        reform = centralized_reformulation(p)
        @test isapprox(reference.objective,reform.objective;rtol=1e-8)
        @test maximum(formulation_residuals(p,reform.x,reform.z)) < 1e-8
        @test maximum(H.maxnorm(reference.x[v]-reform.x[v]) for v in keys(p.nodes)) < 1e-4
        oracle = solve_oracle(p)
        @test oracle.converged
        @test abs(oracle.objective-reference.objective)/abs(reference.objective) < 1e-7
        @test maximum(formulation_residuals(p,oracle.x,oracle.z)) < 1e-4
    end

    @testset "Both inner algorithms solve the assigned penalized problem" begin
        cfg = SolverConfig()
        for m in p.modules
            target = fill(0.5, p.nodes[first(m.vertices)].nV)
            oracle = H.module_oracle(p,m,target,cfg.rho)
            for solver! in (H.local_flat!,H.local_hierarchical!)
                state = H.ModuleState(p,m)
                report = solver!(p,m,state,target,cfg,H.Transmission[])
                @test report.converged
                @test report.residual <= cfg.inner_tol
                @test maximum(H.maxnorm(state.x[v]-oracle.x[v]) for v in m.vertices) < 1e-5
                @test H.maxnorm(m.D*H.stack(m,state.x)) < 1e-6
                J = sum(H.cost(p.nodes[v],state.x[v]) for v in m.vertices) +
                    cfg.rho/2*sum(abs2,state.x[first(m.vertices)]-target)
                @test isapprox(J,oracle.objective;rtol=1e-8)
            end
        end
    end

    @testset "All 16 solver assignments" begin
        results = Dict()
        for g in (:fADMM,:hADMM), a in (:fADMM,:hADMM),
            b in (:fADMM,:hADMM), c in (:fADMM,:hADMM)
            cfg = SolverConfig(global_algorithm=g,module_algorithms=(a,b,c))
            s = solve_hybrid(p;config=cfg,reference=reference.objective)
            h = last(s.history)
            @test s.converged
            @test h.objective_gap < 1e-7
            @test max(h.primal,h.dual,h.internal,h.physical) <= cfg.tol
            @test all(h.inner_converged)
            @test maximum(h.inner_residuals) <= cfg.inner_tol
            @test maximum(H.maxnorm(s.x[v]-reference.x[v]) for v in keys(p.nodes)) < 1e-4
            @test isapprox(H.total_cost(p.root),s.objective;rtol=1e-12)
            @test s.TotalVals == sum(values(s.node_values))
            @test s.MaxNodeVals == maximum(values(s.node_values))
            results[(g,a,b,c)] = s
        end
        # The two global labels intentionally implement the same star updates.
        for a in (:fADMM,:hADMM), b in (:fADMM,:hADMM), c in (:fADMM,:hADMM)
            f, h = results[(:fADMM,a,b,c)],results[(:hADMM,a,b,c)]
            @test f.objective == h.objective
            @test f.TotalVals == h.TotalVals
            @test f.iterations == h.iterations
        end
        @test results[(:fADMM,:hADMM,:hADMM,:hADMM)].TotalVals !=
              results[(:fADMM,:fADMM,:fADMM,:fADMM)].TotalVals
    end

    @testset "Physical relay accounting and inner-budget reporting" begin
        cfg = SolverConfig(module_algorithms=(:hADMM,:fADMM,:hADMM),
                           max_iter=1,inner_max_iter=1,record_messages=true)
        s = solve_hybrid(p;config=cfg)
        @test !s.converged
        @test !all(last(s.history).inner_converged)
        @test s.TotalVals == sum(t.scalars for t in s.messages)
        @test s.MaxNodeVals == maximum(sum(t.scalars for t in s.messages if t.source == v)
                                      for v in keys(p.nodes))
        physical = Set((e.parent,e.child) for e in p.edges)
        @test all((t.source,t.destination) in physical || (t.destination,t.source) in physical
                  for t in s.messages)
        query = filter(t -> t.module_name == :M2 && t.phase == :global_query,s.messages)
        reply = filter(t -> t.module_name == :M2 && t.phase == :global_response,s.messages)
        @test [(t.source,t.destination) for t in query] == [((1,1),(2,2)),((2,2),(3,4))]
        @test [(t.source,t.destination) for t in reply] == [((3,4),(2,2)),((2,2),(1,1))]
        n = p.nodes[(3,4)].nV
        @test [t.scalars for t in query] == [n,n]
        @test [t.scalars for t in reply] == [n+3,n+3]
        @test any(t.source == (4,1) && t.destination == (5,1) && t.phase == :inner_query for t in s.messages)
        # M3's coordinator relays M2 traffic without extra optimization updates.
        @test p.nodes[(2,2)].iteration == 1
        for v in keys(p.nodes)
            @test s.node_values[v] == sum(t.scalars for t in s.messages if t.source == v; init=0)
        end
        @test_throws ErrorException H.transmit!(p,[p.root,p.nodes[(3,4)]],[1.0],:global_query,:M2,H.Transmission[],false)
        @test_throws ArgumentError solve_hybrid(p;config=SolverConfig(global_algorithm=:unknown))
        @test_throws ArgumentError solve_hybrid(p;config=SolverConfig(module_algorithms=(:hADMM,:unknown,:fADMM)))
        @test_throws ArgumentError solve_hybrid(p;config=SolverConfig(rho=0.0))
        @test_throws ArgumentError solve_hybrid(p;config=SolverConfig(inner_max_iter=0))
    end

    @testset "Central, standalone hADMM, and hybrid comparison" begin
        cfg = SolverConfig(module_algorithms=(:hADMM,:fADMM,:hADMM))
        report = compare_solvers(p;configs=[cfg],warmup=false)
        central, hierarchical, hybrid = report.rows
        @test length(report.rows) == 3
        @test central.objective == H.objective(p,report.reference.x)
        @test central.iterations > 0
        @test central.iteration_kind == :Ipopt
        @test central.objective_gap == central.solution_error == 0
        @test central.feasible_gap < 1e-14
        @test central.TotalVals == central.MaxNodeVals == 0
        @test ismissing(central.TotalNodeUpdates)
        for row in (hierarchical,hybrid)
            @test row.converged
            @test row.objective_gap < 1e-7
            @test row.feasible_gap < 1e-7
            @test row.solution_error < 1e-4
            @test row.relative_solution_error < 1e-5
            @test max(row.primal,row.dual,row.physical) <= cfg.tol
            @test row.TotalVals > row.MaxNodeVals > 0
            @test row.seconds > 0
        end
        @test hierarchical.TotalNodeUpdates == length(p.nodes)*hierarchical.iterations
        @test hierarchical.MaxNodeUpdates == hierarchical.iterations
        @test hierarchical.TotalVals == hierarchical.iterations*(2*sum(length(e.c) for e in p.edges)+length(p.edges))
        @test hybrid.TotalNodeUpdates == hybrid.iterations +
              sum(sum(h.inner_iterations .* length.(getfield.(p.modules,:vertices))) for h in report.runs[2].result.history)
        @test hybrid.inner_sweeps == sum(sum(h.inner_iterations) for h in report.runs[2].result.history)
        @test all(h.solution_error >= 0 for run in report.runs for h in run.result.history)
        @test last(report.runs[1].result.history).TotalVals == hierarchical.TotalVals
        # A nondefault penalty exercises propagation through every legacy sweep.
        alternate = solve_hierarchical(p;rho=500.0,reference=report.reference.objective)
        @test alternate.converged
        @test abs(alternate.objective-report.reference.objective)/report.reference.objective < 1e-7
        @test alternate.iterations != hierarchical.iterations
        @test_throws ArgumentError solve_hierarchical(p;rho=0.0)
        @test_throws ArgumentError compare_solvers(p;configs=SolverConfig[])
        @test_throws ArgumentError compare_solvers(p;configs=[cfg,SolverConfig(rho=500.0)])
        mktempdir() do dir
            save_comparison(report,dir)
            @test length(readlines(joinpath(dir,"summary.csv"))) == 4
            @test length(readlines(joinpath(dir,"history.csv"))) == 1+hierarchical.iterations+hybrid.iterations
            @test length(readlines(joinpath(dir,"solutions.csv"))) == 1+3*sum(n.nV for n in values(p.nodes))
            @test occursin("relative_solution_error",first(readlines(joinpath(dir,"summary.csv"))))
        end
    end
end
