using Plots

function plot_comparison(report, output_dir)
    mkpath(output_dir)
    short_name(name) = startswith(name,"Hybrid[") ?
        replace(name,"Hybrid["=>"Hybrid ","fADMM"=>"F","hADMM"=>"H",";"=>":","]"=>"") : name
    gap_iter = plot(xlabel="Outer iteration / tree sweep",ylabel="Feasible relative objective gap",yscale=:log10)
    gap_com = plot(xlabel="Transmitted scalars (all physical hops)",ylabel="Feasible relative objective gap",yscale=:log10)
    solution = plot(xlabel="Outer iteration / tree sweep",ylabel="Maximum solution error",yscale=:log10)
    physical = plot(xlabel="Outer iteration / tree sweep",ylabel="Physical constraint residual",yscale=:log10)
    for run in report.runs
        h = run.result.history
        positive(values) = max.(values,1e-16)
        label = short_name(run.name)
        plot!(gap_iter,getfield.(h,:iteration),positive(getfield.(h,:feasible_gap));label)
        plot!(gap_com,getfield.(h,:TotalVals),positive(getfield.(h,:feasible_gap));label)
        plot!(solution,getfield.(h,:iteration),positive(getfield.(h,:solution_error));label)
        plot!(physical,getfield.(h,:iteration),positive(getfield.(h,:physical));label)
    end
    # The central reference is a single solution; it has no ADMM trajectory.
    hline!(gap_iter,[1e-16];label="Central reference (zero gap)",linestyle=:dash,color=:black)
    for panel in (gap_iter,gap_com,solution,physical)
        plot!(panel;legend=:bottomleft,legendfontsize=7)
    end
    plot!(gap_com;legend=:topright)
    many_runs = length(report.runs) > 4
    if many_runs
        plot!(gap_iter;legend=:outerbottom,legend_columns=2,xlabel="")
        for panel in (gap_com,solution,physical)
            plot!(panel;legend=false)
        end
    end
    plot!(gap_iter;title="Hybrid global:M1/M2/M3 (F=fADMM, H=hADMM)",titlefontsize=10)
    convergence = plot(gap_iter,gap_com,solution,physical;layout=(2,2),
                       size=many_runs ? (1800,1400) : (1400,950),margin=7Plots.mm)
    savefig(convergence,joinpath(output_dir,"convergence.pdf"))

    names = short_name.(getfield.(report.rows,:method))
    positions = collect(eachindex(names))
    communications = bar(positions,getfield.(report.rows,:TotalVals);ylabel="Total scalar transmissions",legend=false,
                         title="Hybrid global:M1/M2/M3 (F=fADMM, H=hADMM)",titlefontsize=10)
    peak = bar(positions,getfield.(report.rows,:MaxNodeVals);ylabel="Peak node scalar transmissions",legend=false)
    updates = bar(positions,[ismissing(r.TotalNodeUpdates) ? NaN : r.TotalNodeUpdates for r in report.rows];
                  ylabel="Node proximal solves (central N/A)",legend=false)
    time = bar(positions,getfield.(report.rows,:seconds);ylabel="Wall time after warm-up (s)",legend=false)
    resources = plot(communications,peak,updates,time;layout=(2,2),size=(1400,1000),
                     xticks=(positions,names),xrotation=70,bottom_margin=15Plots.mm,
                     left_margin=12Plots.mm,top_margin=5Plots.mm,tickfontsize=7)
    savefig(resources,joinpath(output_dir,"resources.pdf"))
    return (; convergence, resources)
end
