"""
kidney_match
Returns a maximum-weight matching in the form of three elements:
- An array of matched vertices
- An array of matched edges
- The corresponding match value.
"""
function kidney_match(weights::Array{Float64,2};
                      graph::DiGraph=DiGraph(),
                      verbose::Int64 = 0,
                      max_cycle_length::Int64 = 2,
                      chain_max_length::Int64 = 0,
                      solver::String = "GLPK",
                      method::String = "MIP_callbacks", # Can be "MIP_callbacks" or "hungarian"
                      ndds::Array{Int64, 1} = zeros(Int64, 0))

  if nv(graph) == 0
      graph = generate_graph(weights)
  end
  m = init_model(solver)
  n = nv(graph)
  @variable(m, 1 >= x[collect(edges(graph))] >= 0, Int)
  @variable(m, 1 >= in[1:n] >= 0, Int)
  @variable(m, 1 >= out[1:n] >= 0, Int)
  @objective(m, Max, sum(x[Edge((i, j))] * weights[i, j] for i in 1:n, j in outneighbors(graph, i)))

  @constraint(m, c1[k=1:n], sum(x[Edge((i,k))] for i=inneighbors(graph, k)) == in[k])
  @constraint(m, c2[k=1:n], sum(x[Edge((k,j))] for j=outneighbors(graph, k)) == out[k])
  @constraint(m, c3[k=setdiff(1:n, ndds)], out[k] <= in[k])
  @constraint(m, c4[k=ndds], in[k] == 0)
  if method == "MIP_callbacks"
      status = solve(m)
      x_val = getvalue(x)
      found, length, cycle = check_cycle_length(x_val, n, max_cycle_length, graph)
      while found
          @constraint(m, sum(x[e] for e in cycle) <= length - 1)
          status = solve(m)
          x_val = getvalue(x)
          found, length, cycle = check_cycle_length(x_val, n, max_cycle_length, graph)
      end
  elseif method == "LP_2_cycles"
      two_cycles_constraint(m, n, graph, x)
      status = solve(m)
  elseif method == "LP_3_cycles"
      two_three_cycles_constraint(m, n, graph, x)
      status = solve(m)
  end
  match_edges = getvalue(x)
  match_vertices = getvalue(in)
  value = getobjectivevalue(m)

  return  (match_vertices, match_edges, value)
end

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

function init_model(solver::String)
  if solver == "gurobi"
    error("No gurobi at this time") # m = Model(solver=GurobiSolver(OutputFlag = 0, PreCrush=1))
  elseif solver == "mosek"
    m = Model(solver=MosekSolver(LOG = 0))
  elseif solver == "GLPK"
    m = Model(solver=GLPKSolverMIP())
  else
    error("Solver $(solver) is not supported")
  end
  return m
end

function two_three_cycles_constraint(m::JuMP.Model, n::Int64, graph::DiGraph,
          x#::JuMP.JuMPArray{JuMP.Variable,1,Tuple{Array{LightGraphs.SimpleGraphs.SimpleEdge{Int64},1}}}
          ) # TODO: Breaking change in JuMP syntax, will fix later

  # If there is a path of length 2, it has to be closed as a 3-cycle
  # 2-cycles are not impacted because i is not in ingoing[i]

  @constraint(m,
              c5[i=1:n,
                 j=outneighbors(graph, i),
                 k=intersect(outneighbors(graph, j), inneighbors(graph, i))],
              x[Edge((k,i))] >= x[Edge((i,j))] + x[Edge((j,k))] - 1
             )
 @constraint(m,
             c6[i=1:n, j=outneighbors(graph, i), k=setdiff(outneighbors(graph, j), inneighbors(graph, i))],
             Int(k == i) >= x[Edge((i,j))] + x[Edge((j,k))] - 1
            )
  # if k == i, then no constraint.
end

function two_cycles_constraint(m::JuMP.Model, n::Int64, graph::DiGraph,
  x#::JuMP.JuMPArray{JuMP.Variable,1,Tuple{Array{LightGraphs.SimpleGraphs.SimpleEdge{Int64},1}}}
  )
  # For 2-cycles, edges.list only contains pairs for which both ways exist
  # this enables faster model creation.
  @constraint(m,
              c5[i=1:n, j=intersect(outneighbors(graph, i), inneighbors(graph, i))],
              x[Edge((i,j))] == x[Edge((j,i))])
  @constraint(m,
              c6[i=1:n, j=setdiff(outneighbors(graph, i), inneighbors(graph, i))],
              x[Edge((i,j))] == 0)
  # TODO: Check that the following constraint is redundant
  # @constraint(m,
  #             c6[i=1:n, j=setdiff(inneighbors(graph, i), outneighbors(graph, i))],
  #             x[Edge((j,i))] == 0)
end

"""
Checks that a matching does not contain cycles of length greater than `max_cycle_length`.
Returns such a cycle if one is found as an array of edges.
"""
function check_cycle_length(x_val, n, max_cycle_length, graph)
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

function find_recipient(i, x_val, graph)
    recipient = 0
    for j in outneighbors(graph, i)
        if x_val[Edge((i,j))] == 1
            recipient = j
            break
        end
    end
    return recipient
end
