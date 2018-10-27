@testset "LP_2_cycles formulation" begin

v, e, val = cycle_match(w1, graph = g1, method="LP_2_cycles", max_cycle_length = 2)
@test (v == [1,1,0,0]) & (val == 2)
end;

@testset "LP_3_cycles formulation" begin

v, e, val = cycle_match(w1, graph = g1, max_cycle_length = 3)
@test (v == [1,0,1,1]) & (val == 3)
end;
