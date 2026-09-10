import Pkg
Pkg.activate(".")
Pkg.instantiate()

using Pkg, Plots, Graphs, GraphRecipes, NPZ, Statistics
using LinearAlgebra, JuMP, Ipopt
using Zygote
using ProximalAlgorithms
using ProximalOperators
using ProximalCore
using DifferentiationInterface: AutoZygote

include("../DisjointProblem/TreeGeneration.jl")
include("../DisjointProblem/CentralSolution.jl")
include("../DisjointProblem/hADMM.jl")
include("../DisjointProblem/fADMM.jl")
include("hybrid_ADMM.jl")

const mode = 2  # 1: optimality gap, 2: run in specific iterations, 3: stopping criteria
const nN = 20
const nD = 5

include("../DisjointProblem/Parameters.jl")

global countID
global opt_value

methods = ["hADMM", "fADMM", "hybridADMM"]

node_iter = Dict(method => Int64[] for method in methods)
tt_com = Dict(method => Int64[] for method in methods)
max_com = Dict(method => Int64[] for method in methods)
inner_iter = Int64[]

nTestTopo = 10

fontsize = 16
figPrime = plot(framestyle = :box, guidefont = font(16), tickfontsize = fontsize, xlabel = "Number of iterations", yticks = [1, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6])
figRes = plot(framestyle = :box, guidefont = font(16), tickfontsize = fontsize, xlabel = "Number of iterations", yticks = [1e2, 1e1, 1., 0.1, 1e-2, 1e-3, 1e-4])
figJ = plot(framestyle = :box, guidefont = font(16), tickfontsize = fontsize, xlabel = "Number of iterations", yticks = [1, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6])

avg_gap_hADMM = zeros(Niter)
avg_gap_fADMM = zeros(Niter)
avg_gap_hybrid = zeros(Niter)

mkpath(joinpath("media", "figs", "hybrid_admm"))
mkpath(joinpath("data", "hybrid-admm"))

for tp in 1:nTestTopo
    println("Solving topology $tp")

    global countID = 0
    global opt_value = 0

    root = linknode(string(countID += 1))
    topo_gen!(root, nN - 1, nD - 1)
    assign!(root)

    dict_result, opt_value = get_CenVars(root)

    reset!(root)
    traj_err_hADMM, traj_res_hADMM, traj_opt_hADMM = hADMM(root, dict_result)
    total, max_num = tt_com_iter(root)
    push!(node_iter["hADMM"], max_num["iter"])
    push!(max_com["hADMM"], max_num["com"])
    push!(tt_com["hADMM"], total["com"])

    if mode == 2
        avg_gap_hADMM .+= abs.(traj_opt_hADMM .- opt_value)/maximum(abs.(traj_opt_hADMM .- opt_value))/nTestTopo
        avg_gap_hADMM .= monotone_func(avg_gap_hADMM)
    end

    plot!(figPrime, 1:length(traj_err_hADMM), traj_err_hADMM, yscale = :log10, grid = true, label = "")
    plot!(figRes, 1:length(traj_res_hADMM), traj_res_hADMM, yscale = :log10, grid = true, label = "")
    plot!(figJ, 1:length(traj_opt_hADMM), abs.(traj_opt_hADMM .- opt_value)/maximum(abs.(traj_opt_hADMM[1] .- opt_value)), yscale = :log10, grid = true, label = "")

    reset!(root)
    traj_opt_fADMM = flattenADMM(root)
    total, max_num = tt_com_iter(root)
    push!(node_iter["fADMM"], max_num["iter"])
    push!(max_com["fADMM"], max_num["com"])
    push!(tt_com["fADMM"], total["com"])

    if mode == 2
        avg_gap_fADMM .+= abs.(traj_opt_fADMM .- opt_value)/maximum(abs.(traj_opt_fADMM[1] .- opt_value))/nTestTopo
        avg_gap_fADMM .= monotone_func(avg_gap_fADMM)
    end

    reset!(root)
    traj_err_hybrid, traj_res_hybrid, traj_opt_hybrid, traj_inner_hybrid = hybridADMM(root, dict_result)
    total, max_num = tt_com_iter(root)
    push!(node_iter["hybridADMM"], max_num["iter"])
    push!(max_com["hybridADMM"], max_num["com"])
    push!(tt_com["hybridADMM"], total["com"])
    append!(inner_iter, traj_inner_hybrid)

    if mode == 2
        avg_gap_hybrid .+= abs.(traj_opt_hybrid .- opt_value)/maximum(abs.(traj_opt_hybrid[1] .- opt_value))/nTestTopo
        avg_gap_hybrid .= monotone_func(avg_gap_hybrid)
    end

    plot!(figPrime, 1:length(traj_err_hybrid), traj_err_hybrid, yscale = :log10, grid = true, label = "")
    plot!(figRes, 1:length(traj_res_hybrid), traj_res_hybrid, yscale = :log10, grid = true, label = "")
    plot!(figJ, 1:length(traj_opt_hybrid), abs.(traj_opt_hybrid .- opt_value)/maximum(abs.(traj_opt_hybrid[1] .- opt_value)), yscale = :log10, grid = true, label = "")
end

savefig(figPrime, joinpath("media", "figs", "hybrid_admm", string("HY-Prime-Conver-D=", nD, "-N=", nN, ".pdf")))
savefig(figRes, joinpath("media", "figs", "hybrid_admm", string("HY-Res-Conver-D=", nD, "-N=", nN, ".pdf")))
savefig(figJ, joinpath("media", "figs", "hybrid_admm", string("HY-Cost-Conver-D=", nD, "-N=", nN, ".pdf")))

medstep = Dict()
for (key, value) in node_iter
    medstep[key] = round(median(value))
    min_value, _ = findmin(value)
    max_value, _ = findmax(value)
    println("$key Maximum number of iterations in a node (min) avg. (max): ($min_value) $(medstep[key]) ($max_value)")
end

println()
for (key, value) in max_com
    med = round(median(value))
    min_value, _ = findmin(value)
    max_value, _ = findmax(value)
    println("$key Maximum number of scalar variables sent by a node (min) avg. (max): ($min_value) $med ($max_value)")
end

println()
figComCost = plot(framestyle = :box, guidefont = fontsize, legendfontsize = fontsize, tickfontsize = fontsize, xlabel = "TotalVals", ylabel = "Optimality Gap", yscale = :log10, ylims = [1e-7, 2])
for (key, value) in tt_com
    med = round(median(value))
    min_value, _ = findmin(value)
    max_value, _ = findmax(value)
    println("$key Total number of scalar variables sent in network (min) avg. (max): ($min_value) $med ($max_value)")
end

if mode == 2
    plot!(figComCost, (0:(Niter - 1))*mean(tt_com["hADMM"])/(Niter - 1), avg_gap_hADMM, label = "hADMM", linewidth = 4, seriestype = :steppost)
    plot!(figComCost, (0:(Niter - 1))*mean(tt_com["fADMM"])/(Niter - 1), avg_gap_fADMM, label = "fADMM", linewidth = 4, seriestype = :steppost)
    plot!(figComCost, (0:(Niter - 1))*mean(tt_com["hybridADMM"])/(Niter - 1), avg_gap_hybrid, label = "hybridADMM", linewidth = 4, seriestype = :steppost)
    savefig(figComCost, joinpath("media", "figs", "hybrid_admm", string("HY-Cost-Com-ADMM-D=", nD, "-N=", nN, ".pdf")))
end

npzwrite(joinpath("data", "hybrid-admm", string("Max-iter-D=", nD, "-N=", nN, ".npz")), node_iter)
npzwrite(joinpath("data", "hybrid-admm", string("Max-com-D=", nD, "-N=", nN, ".npz")), max_com)
npzwrite(joinpath("data", "hybrid-admm", string("Tot-com-D=", nD, "-N=", nN, ".npz")), tt_com)
npzwrite(joinpath("data", "hybrid-admm", string("Inner-iter-D=", nD, "-N=", nN, ".npz")), Dict("hybridADMM" => inner_iter))
