include("HybridADMM.jl")
using .HybridADMM

include("PlotComparison.jl")

function main(; all_combinations="--all" in ARGS,
                output_dir=joinpath(@__DIR__,"results"))
    configs = if all_combinations
        [SolverConfig(global_algorithm=g, module_algorithms=(m1,m2,m3))
         for g in (:fADMM,:hADMM), m1 in (:fADMM,:hADMM),
             m2 in (:fADMM,:hADMM), m3 in (:fADMM,:hADMM)]
    else
        [SolverConfig(global_algorithm=:fADMM, module_algorithms=(:hADMM,:fADMM,:hADMM))]
    end
    report = compare_solvers(build_problem(); configs)
    print_comparison(report)
    save_comparison(report,output_dir)
    plot_comparison(report,output_dir)
    println("Centralized, hADMM, fADMM, nested hybrid, and Candidate B results saved to ", output_dir)
    println("Global fADMM/hADMM use equivalent updates on this two-level star.")
    return report
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
