# kidney matching fallback method

"""
To dos:
- add function to create random edge lists - need to select pairs without replacement
- add MIP for selecting two-cycles
- add comparison between MIP and greedy algorithm
- add other distributions of weights
"""


# read in packages
using CSV, LightGraphs, DataFrames

# include("match.jl")

# read in list of edges
edge_list = DataFrame(CSV.read("fake_edge_list.csv",header=false))

# run functions for fallback method
edge_weights = generate_random_weights(edge_list)
graph = generate_graph(edge_weights)
cycles = find_2cycles(edge_weights)
cycles_chosen, total_weight = greedy_2cycle_match(cycles)

# functions
"""
find_2cycles() identifies cycles of length two using an edge-weight matrix as input;
the function returns a dataframe with a list of cycles and their corresponding weight;
cycles returned are sorted in descending order by weight
"""
function find_2cycles(weights::Array{Float64,2}; graph::DiGraph=DiGraph())
    if nv(graph) == 0
      graph = generate_graph(weights)
    end
    cycle_vector = simplecycles_limited_length(graph,2)
    cycles = DataFrame(Vert1 = Int64[], Vert2 = Int64[], Weight = Float64[])
    for i in 1:length(cycle_vector)
      cycle_val = weights[cycle_vector[i][1],cycle_vector[i][2]]+weights[cycle_vector[i][2],cycle_vector[i][1]]
      push!(cycles,(cycle_vector[i][1],cycle_vector[i][2],cycle_val))
    end
    sort!(cycles,[:Weight], rev=true)
    return cycles
end

"""
greedy_2cycle_match() uses a greedy algorithm to select two-cycles,
taking as input a list of cycles with weights and returning a list
of selected cycles and the total value of their weights
"""
function greedy_2cycle_match(cycles::DataFrame)
    cycles_chosen = DataFrame(Vert1 = Int64[], Vert2 = Int64[], Weight = Float64[])
    edges_used = Int64[]
    for i in 1:nrow(cycles)
        if (!(cycles[i,:Vert1] in edges_used) & !(cycles[i,:Vert2] in edges_used))
            append!(cycles_chosen,cycles[i,:])
            push!(edges_used,cycles[i,:Vert1])
            push!(edges_used,cycles[i,:Vert2])
        end
    end
    total_weight = sum(cycles_chosen[:,:Weight])
    return cycles_chosen, total_weight
end

"""
generate_random_weights() creates random weights using a list of edges as input;
weights are real values sampled uniformly between 0 and 100;
cycles of length one are automatically removed to avoid future errors
"""
function generate_random_weights(edge_list::DataFrame)
    num_vertices = max(maximum(edge_list[:,1]),maximum(edge_list[:,2]))
    weights = zeros(Float64,num_vertices,num_vertices)
    for i in 1:nrow(edge_list)
        weights[edge_list[i,1],edge_list[i,2]] = rand()*100
    end
    for i in 1:size(edge_weights,1)
        weights[i,i] = 0
    end
    return weights
end
