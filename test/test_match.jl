g1 = DiGraph(4)
add_edge!(g1, 1,2)
add_edge!(g1, 2,1)
add_edge!(g1, 1,3)
add_edge!(g1, 3,4)
add_edge!(g1, 4,1)
w1 = ones(Float64, 4, 4)

@testset "chain_cycle_match" begin
    @testset "MIP_callbacks" begin

    v, e, val = chain_cycle_match(w1, graph = g1, max_cycle_length = 2)
    @test (v == [1,1,0,0]) & (val == 2)

    v, e, val = chain_cycle_match(w1, graph = g1, max_cycle_length = 3)
    @test (v == [1,0,1,1]) & (val == 3)

    v, e, val = chain_cycle_match(w1, graph = g1;
                             max_cycle_length = 0,
                             ndds= 3 * ones(Int64, 1))
    @test (v == [1,1,0,1]) & (val == 3)

    v, e, val = chain_cycle_match(w1, graph = g1, max_cycle_length = 0)
    @test (v == [0,0,0,0]) & (val == 0)

    g2 = DiGraph(4)
    add_edge!(g2, 1,2)
    add_edge!(g2, 2,1)
    add_edge!(g2, 1,3)
    add_edge!(g2, 1,4)

    v, e, val = chain_cycle_match(w1, graph = g2, max_cycle_length = 3)
    @test (v == [1,1,0,0]) & (val == 2)
    end;

    @testset "matrix input" begin
    w2 = zeros(Float64, 4, 4)
    w2[1,2] = 1
    w2[2,1] = 1
    w2[1,3] = 1
    w2[3,4] = 1
    w2[4,1] = 1
    v, e, val = chain_cycle_match(w2, max_cycle_length = 3)
    @test (v == [1,0,1,1]) & (val == 3)
    end;

end; # test chain_cycle_match
