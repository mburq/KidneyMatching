using JuMP, LightGraphs, GLPK, GLPKMathProgInterface
using Test, KidneyMatching

g1 = DiGraph(4)
add_edge!(g1, 1,2)
add_edge!(g1, 2,1)
add_edge!(g1, 1,3)
add_edge!(g1, 3,4)
add_edge!(g1, 4,1)
# l = Log(1, 1, false)
w1 = ones(Float64, 4, 4)

v, e, val = kidney_match(g1, w1, techno = "cycles", max_cycle_length = 2, solver = "GLPK")
@test v == [1,1,0,0]
@test val == 2

v, e, val = kidney_match(g1, w1, techno = "cycles", max_cycle_length = 3, solver = "GLPK")
@test v == [1,0,1,1]
@test val == 3

v, e, val = kidney_match(g1, w1; techno = "chains", chain_max_length = 3,
             max_cycle_length = 0, solver = "GLPK", ndds= 3 * ones(Int64, 1))
@test v == [1,1,0,1]
@test val == 3

# v, e, val = kidney_match(g1, w1, techno = "chains", chain_max_length = 3,
#              max_cycle_length = 0, solver = "gurobi")
# @test v == [0,0,0,0]
# @test val == 0

v, e, val = kidney_match(g1, w1, techno = "chains", chain_max_length = 3,
             max_cycle_length = 0, solver = "GLPK")
@test v == [0,0,0,0]
@test val == 0

g2 = DiGraph(4)
add_edge!(g2, 1,2)
add_edge!(g2, 2,1)
add_edge!(g2, 1,3)
add_edge!(g2, 1,4)

v, e, val = kidney_match(g2, w1, techno = "cycles", max_cycle_length = 3, solver = "GLPK")
@test v == [1,1,0,0]
@test val == 2
