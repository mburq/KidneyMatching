"""
basic_MIP
Returns a JuMP model with the basic matching MIP. On its own, solves the problem with
unbounded chains and cycles.
"""
function basic_MIP(weights::Array{Float64,2}, graph::DiGraph, ndds::Array{Int64, 1})
    m = Model(solver=GLPKSolverMIP())
    n = nv(graph)
    @variable(m, 1 >= x[collect(edges(graph))] >= 0, Int)
    @variable(m, 1 >= in[1:n] >= 0, Int)
    @variable(m, 1 >= out[1:n] >= 0, Int)
    @objective(m, Max, sum(x[Edge((i, j))] * weights[i, j] for i in 1:n, j in outneighbors(graph, i)))

    @constraint(m, c1[k=1:n], sum(x[Edge((i,k))] for i=inneighbors(graph, k)) == in[k])
    @constraint(m, c2[k=1:n], sum(x[Edge((k,j))] for j=outneighbors(graph, k)) == out[k])
    @constraint(m, c3[k=setdiff(1:n, ndds)], out[k] <= in[k])
    @constraint(m, c4[k=ndds], in[k] == 0)
    return (m, x, in)
end

"""
chain_cycle_match
Returns a maximum-weight matching in the form of three elements:
- An array of matched vertices
- An array of matched edges
- The corresponding match value.

Iteratively re-solves because callbacks are not available in GLPK.jl.
"""
function chain_cycle_match(weights::Array{Float64,2};
                           graph::DiGraph=DiGraph(),
                           verbose::Int64 = 0,
                           max_cycle_length::Int64 = 2,
                           max_chain_length::Float64 = Inf64,
                           solver::String = "GLPK",
                           ndds::Array{Int64, 1} = zeros(Int64, 0))

  if nv(graph) == 0
      graph = generate_graph(weights)
  end
  m, x, in = basic_MIP(weights, graph, ndds)
  n = nv(graph)
  status = solve(m)
  x_val = getvalue(x)
  found, length, cycle = check_cycle_length(x_val, n, max_cycle_length, graph)
  while found
      @constraint(m, sum(x[e] for e in cycle) <= length - 1)
      status = solve(m)
      x_val = getvalue(x)
      found, length, cycle = check_cycle_length(x_val, n, max_cycle_length, graph)
  end
  match_edges = getvalue(x)
  match_vertices = getvalue(in)
  value = getobjectivevalue(m)

  return  (match_vertices, match_edges, value)
end

"""
Checks that a matching does not contain cycles of length greater than `max_cycle_length`.
Returns such a cycle if one is found as an array of edges.
"""
function check_cycle_length(x_val, n::Int64, max_cycle_length::Int64, graph::DiGraph)
    flag = zeros(Int64, n)
    for i in 1:n
        if flag[i] == 0
            current_vertex = i
            cycle_length = 0
            cycle = []
            while true
                next_vertex = find_recipient(current_vertex, x_val, graph)
                if next_vertex == 0 # no recipient, end of a chain
                    break
                end
                push!(cycle, Edge((current_vertex, next_vertex)))
                flag[next_vertex] = 1
                cycle_length += 1
                if next_vertex == i
                    if cycle_length > max_cycle_length
                        return (true, cycle_length, cycle)
                    else
                        break # found the cycle
                    end
                end
                current_vertex = next_vertex
            end
        end
    end
    return (false, 0, [])
end

function find_recipient(i::Int64, x_val, graph::DiGraph)
    recipient = 0
    for j in outneighbors(graph, i)
        if x_val[Edge((i,j))] == 1
            recipient = j
            break
        end
    end
    return recipient
end

"""
generate_graph
Creates a LightGraph directed graph from an adjacency edge-weight matrix.
"""
function generate_graph(weights::Array{Float64,2})
    @assert size(weights)[1] == size(weights)[2]
    n = size(weights)[1]
    graph = DiGraph(n)
    for i in 1:n
        for j in 1:n
            if weights[i,j] > 0
                add_edge!(graph, i,j)
            end
        end
    end
    return graph
end
