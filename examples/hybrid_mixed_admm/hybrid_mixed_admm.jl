module HybridMixedInfrastructure

import ProximalAlgorithms
include(joinpath(@__DIR__, "..", "..", "src", "linknodes.jl"))

export HybridModule, Dependency, MessageTracker, PhysicalSend,
       m123_tree, m123_configuration, physical_path, route_message!, communication_metrics

const NodeIndex = Tuple{Int, Int}

struct HybridModule
    id::String
    physical_nodes::Set{NodeIndex}
    module_root::NodeIndex
    algorithm::Union{String, Nothing}

    function HybridModule(; id, physical_nodes, module_root, algorithm)
        nodes = Set{NodeIndex}(physical_nodes)
        module_root in nodes || throw(ArgumentError("module_root must belong to physical_nodes"))
        algorithm in ("hadmm", "fadmm", nothing) ||
            throw(ArgumentError("algorithm must be hadmm, fadmm, or nothing"))
        # nothing denotes a Stage A routing fixture only; no solver is executed here.
        new(id, nodes, module_root, algorithm)
    end
end

struct Dependency
    producer_module::String
    consumer_module::String
    quantity::String
    iterate_index::String
    consumer_update::String

    function Dependency(; producer_module, consumer_module, quantity, iterate_index, consumer_update)
        isempty(strip(quantity)) && throw(ArgumentError("dependency quantity must be specified"))
        isempty(strip(consumer_update)) && throw(ArgumentError("consumer_update must be specified"))
        iterate_index in ("k", "k+1") || throw(ArgumentError("iterate_index must be k or k+1"))
        new(producer_module, consumer_module, quantity, iterate_index, consumer_update)
    end
end

"""Build the exact 16-node M123 physical tree with canonical (d,i) -> string IDs.

Parent links describe communication routes, not cross-module optimization dependencies.
"""
function m123_tree()
    indices = [(d, i) for (d, count) in enumerate((1, 2, 5, 5, 3)) for i in 1:count]
    node_ids = Dict(index => string(i) for (i, index) in enumerate(indices))
    nodes = Dict(index => linknode(node_ids[index]) for index in indices)
    set_relative!(nodes[(1, 1)], [nodes[(2, 1)], nodes[(2, 2)]])
    set_relative!(nodes[(2, 1)], [nodes[(3, 1)], nodes[(3, 2)], nodes[(3, 3)]])
    set_relative!(nodes[(2, 2)], [nodes[(3, 4)], nodes[(3, 5)]])
    set_relative!(nodes[(3, 4)], [nodes[(4, 1)], nodes[(4, 2)], nodes[(4, 3)]])
    set_relative!(nodes[(3, 5)], [nodes[(4, 4)], nodes[(4, 5)]])
    set_relative!(nodes[(4, 1)], [nodes[(5, 1)], nodes[(5, 2)]])
    set_relative!(nodes[(4, 3)], [nodes[(5, 3)]])
    return (; root = nodes[(1, 1)], nodes, node_ids)
end

"""Return the fixed physical tree and independent global fADMM peer-module graph.

This is structural configuration only. In particular, (2,2) relays traffic to M2;
it is not M2's algorithmic parent. The q2 -> M2 root update -> xi2 numerical
interface remains undefined in formulation_M123.md Section 10.
"""
function m123_configuration()
    physical_tree = m123_tree()
    modules = Dict(
        "M1" => HybridModule(id = "M1", module_root = (2, 1), algorithm = "fadmm",
            physical_nodes = [(2, 1), (3, 1), (3, 2), (3, 3)]),
        "M2" => HybridModule(id = "M2", module_root = (3, 4), algorithm = "hadmm",
            physical_nodes = [(3, 4), (4, 1), (4, 2), (4, 3), (5, 1), (5, 2), (5, 3)]),
        "M3" => HybridModule(id = "M3", module_root = (2, 2), algorithm = "fadmm",
            physical_nodes = [(2, 2), (3, 5), (4, 4), (4, 5)]),
    )
    algorithmic_module_graph = Dict("R" => ["M1", "M2", "M3"])
    return (; physical_tree, modules, algorithmic_module_graph,
            global_root = (1, 1), global_algorithm = "fadmm")
end

"""Return the physical path, including endpoints, using parent links and a common ancestor."""
function physical_path(origin::linknode, destination::linknode)
    ancestors = linknode[]
    pointer = origin
    while pointer !== nothing
        push!(ancestors, pointer)
        pointer = pointer.parent
    end

    destination_branch = linknode[]
    pointer = destination
    while pointer !== nothing
        common = findfirst(node -> node === pointer, ancestors)
        if common !== nothing
            return vcat(ancestors[1:common], reverse(destination_branch))
        end
        push!(destination_branch, pointer)
        pointer = pointer.parent
    end
    throw(ArgumentError("origin and destination must belong to the same physical tree"))
end

struct PhysicalSend
    iteration::Int
    logical_origin::String
    logical_destination::String
    physical_sender::String
    physical_receiver::String
    message_type::Symbol
    payload_dimension::Int
    is_relay::Bool
end

struct MessageTracker
    sends::Vector{PhysicalSend}
end

MessageTracker() = MessageTracker(PhysicalSend[])

"""Record a logical message as separate physical sends; no optimization or legacy accounting."""
function route_message!(tracker::MessageTracker, origin::linknode, destination::linknode;
                        iteration::Int, message_type::Symbol, payload_dimension::Int)
    payload_dimension >= 0 || throw(ArgumentError("payload_dimension must be nonnegative"))
    path = physical_path(origin, destination)
    for i in 1:(length(path) - 1)
        sender, receiver = path[i], path[i + 1]
        push!(tracker.sends, PhysicalSend(
            iteration, origin.ID, destination.ID, sender.ID, receiver.ID,
            message_type, payload_dimension, sender !== origin && sender !== destination,
        ))
    end
    return path
end

"""Compute aggregate and per-sender metrics over recorded physical sends.

Per-node dictionaries use physical string IDs; an absent sender has count zero.
"""
function communication_metrics(tracker::MessageTracker)
    sent_messages = Dict{String, Int}()
    sent_scalars = Dict{String, Int}()
    total_scalars = 0
    relay_messages = 0
    for send in tracker.sends
        sender = send.physical_sender
        sent_messages[sender] = get(sent_messages, sender, 0) + 1
        sent_scalars[sender] = get(sent_scalars, sender, 0) + send.payload_dimension
        total_scalars += send.payload_dimension
        relay_messages += send.is_relay
    end
    return (
        TotalMessages = length(tracker.sends),
        TotalScalarValues = total_scalars,
        NodeSentMessages = sent_messages,
        MaxNodeMessages = maximum(values(sent_messages); init = 0),
        NodeSentScalars = sent_scalars,
        MaxNodeScalarValues = maximum(values(sent_scalars); init = 0),
        RelayMessages = relay_messages,
    )
end

end
