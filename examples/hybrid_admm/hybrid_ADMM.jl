function direct_child_indices(root::linknode)
    indices = Dict()
    start_index = 1

    for child in root.children
        indices[child.ID] = start_index:(start_index + child.nV - 1)
        start_index += child.nV
    end

    return indices
end


function update_root_hybrid!(
    root::linknode,
    module_prime::Dict,
    module_dual::Dict,
    index::Dict,
    query::Dict,
    g,
    x0::Vector{Float64};
    λ = λₛ,
)
    para = root.cost_func.para
    func = root.cost_func.val

    for child in root.children
        query[child.ID] = module_prime[child.ID] - module_dual[child.ID]
    end

    f = ProximalAlgorithms.AutoDifferentiable(
        v -> func(v; para = para) + (1/(2λ))*sum(
            dot(v[index[child.ID]] - query[child.ID], v[index[child.ID]] - query[child.ID])
            for child in root.children
        ),
        AutoZygote(),
    )

    solution, _ = root.solver(f = f, g = g, x0 = x0)
    copyto!(x0, solution)

    for child in root.children
        root.prime[child.ID] = solution[index[child.ID]]
    end
end


function update_module_root_hybrid!(
    node::linknode,
    top_query::Vector{Float64};
    λ_global = λₛ,
    λ_inner = λₕ,
)
    para = node.cost_func.para
    func = node.cost_func.val
    g = ProximalOperators.NormL1(node.cost_func.w)
    x0 = vect_prime(node)

    if node.children === nothing
        f = ProximalAlgorithms.AutoDifferentiable(
            v -> func(v; para = para) + (1/(2λ_global))*dot(v - top_query, v - top_query),
            AutoZygote(),
        )

        solution, _ = node.solver(f = f, g = g, x0 = x0)
        node.prime[node.ID] = solution
        return
    end

    child_query = Dict(
        child.ID => vect_prime(child) - node.dual[child.ID]
        for child in node.children
    )

    index = Dict()
    start_index = 1
    for child in node.children
        index[child.ID] = start_index:(start_index + child.nV - 1)
        start_index += child.nV
    end

    f = ProximalAlgorithms.AutoDifferentiable(
        v -> func(v; para = para)
            + (1/(2λ_global))*dot(v - top_query, v - top_query)
            + (1/(2λ_inner))*sum(
                dot(v[index[child.ID]] - child_query[child.ID], v[index[child.ID]] - child_query[child.ID])
                for child in node.children
            ),
        AutoZygote(),
    )

    solution, _ = node.solver(f = f, g = g, x0 = x0)

    for child in node.children
        node.prime[child.ID] = solution[index[child.ID]]
    end
end


function hierarchical_sweep_inside_module!(
    node::linknode,
    ter::Vector{Float64};
    λ_inner = λₕ,
)
    node.iteration += 1
    proxAlg!(node; λ = λ_inner)

    if node.children !== nothing
        for child in node.children
            child_prime_old = vect_prime(child)
            hierarchical_sweep_inside_module!(child, ter; λ_inner = λ_inner)

            child_prime = vect_prime(child)
            prime_res = node.prime[child.ID] - child_prime
            dual_res = (child_prime - child_prime_old)/λ_inner
            node.dual[child.ID] += prime_res

            com_cost!(child, child_prime, 1)
            push!(ter, max(norm(prime_res, Inf), norm(dual_res, Inf)))
        end
    end
end


function solve_module_hADMM!(
    module_root::linknode,
    top_query::Vector{Float64};
    λ_global = λₛ,
    λ_inner = λₕ,
    tol = tol,
    max_iter = max_iter,
    verbose = false,
)
    final_residual = 0.0

    for iteration in 1:max_iter
        ter = Float64[]
        module_root.iteration += 1

        com_cost!(module_root.parent, top_query, 1)
        update_module_root_hybrid!(
            module_root,
            top_query;
            λ_global = λ_global,
            λ_inner = λ_inner,
        )

        if module_root.children !== nothing
            for child in module_root.children
                child_prime_old = vect_prime(child)
                hierarchical_sweep_inside_module!(child, ter; λ_inner = λ_inner)

                child_prime = vect_prime(child)
                prime_res = module_root.prime[child.ID] - child_prime
                dual_res = (child_prime - child_prime_old)/λ_inner
                module_root.dual[child.ID] += prime_res

                com_cost!(child, child_prime, 1)
                push!(ter, max(norm(prime_res, Inf), norm(dual_res, Inf)))
            end
        end

        final_residual = isempty(ter) ? 0.0 : maximum(ter)
        if final_residual < tol
            return final_residual, iteration
        end
    end

    if verbose
        println("Module $(module_root.ID) reached inner max_iter with residual $final_residual")
    end

    return final_residual, max_iter
end


function hybridADMM(
    root::linknode,
    dict_result::Dict;
    tol = tol,
    λ_global = λₛ,
    λ_inner = λₕ,
    inner_tol = tol,
    max_iter = max_iter,
    inner_max_iter = max_iter,
    verbose = false,
)
    if root.children === nothing
        error("hybridADMM requires a root with at least one child module")
    end

    traj_err = Float64[]
    traj_res = Float64[]
    opt_value_hybrid = Float64[]
    inner_iterations = Int64[]

    index = direct_child_indices(root)
    root_g = ProximalOperators.NormL1(root.cost_func.w)
    root_x0 = vect_prime(root)
    query = Dict(child.ID => similar(root.prime[child.ID]) for child in root.children)
    module_prime = Dict(child.ID => copy(vect_prime(child)) for child in root.children)
    module_prime_old = deepcopy(module_prime)

    for iteration in 1:max_iter
        root.iteration += 1

        update_root_hybrid!(
            root,
            module_prime,
            root.dual,
            index,
            query,
            root_g,
            root_x0;
            λ = λ_global,
        )

        max_residual = 0.0
        inner_count_this_iter = 0

        for child in root.children
            top_query = root.prime[child.ID] + root.dual[child.ID]
            old_module_prime = copy(module_prime_old[child.ID])

            inner_residual, inner_count = solve_module_hADMM!(
                child,
                top_query;
                λ_global = λ_global,
                λ_inner = λ_inner,
                tol = inner_tol,
                max_iter = inner_max_iter,
                verbose = verbose,
            )

            module_prime[child.ID] = copy(vect_prime(child))
            prime_res = root.prime[child.ID] - module_prime[child.ID]
            dual_res = (module_prime[child.ID] - old_module_prime)/λ_global
            root.dual[child.ID] += prime_res

            com_cost!(child, module_prime[child.ID], 1)

            max_residual = max(max_residual, norm(prime_res, Inf), norm(dual_res, Inf), inner_residual)
            inner_count_this_iter += inner_count
        end

        module_prime_old = deepcopy(module_prime)

        current_cost = total_cost(root)
        push!(opt_value_hybrid, current_cost)
        push!(traj_res, max_residual)
        push!(inner_iterations, inner_count_this_iter)

        err = Float64[]
        get_err!(root, dict_result, err)
        push!(traj_err, sum(err))

        if mode == 1
            if abs((opt_value - current_cost)/opt_value) < opt_gap
                if verbose println("hybridADMM converged after $iteration outer iterations") end
                return traj_err, traj_res, opt_value_hybrid, inner_iterations
            end
        elseif mode == 2
            if iteration >= Niter
                if verbose println("hybridADMM stopped after $iteration outer iterations") end
                return traj_err, traj_res, opt_value_hybrid, inner_iterations
            end
        else
            if max_residual < tol
                if verbose println("hybridADMM converged after $iteration outer iterations") end
                return traj_err, traj_res, opt_value_hybrid, inner_iterations
            end
        end
    end

    return traj_err, traj_res, opt_value_hybrid, inner_iterations
end
