using SyntheticGrids
using Base.Test
using JSON

@testset "SyntheticGrids" begin
    const SEED = 666
    const GENPATH = joinpath(dirname(@__FILE__), "..", "data", "GenData.json")
    const CENSUSPATH = joinpath(dirname(@__FILE__), "..", "data", "Census_data.dat")
    const GENCOORDPATH = joinpath(dirname(@__FILE__), "..", "data", "Generator_coord.dat")
    const GENDATAPATH = joinpath(dirname(@__FILE__), "..", "data", "Generator_data.dat")
    const DUMMYPATH = joinpath(dirname(@__FILE__), "..", "data", "dummy.json")
    grid = Grid(SEED)
    place_loads_from_zips!(grid; latlim = (38, 40), longlim = (-89, -88))
    @test length(grid.buses) == 130
    place_gens_from_data!(grid; latlim = (38, 40), longlim = (-89, -88))
    @test length(grid.buses) == 147
    connect!(grid)
    @test sum(grid.bus_conn) == 358
    count = SyntheticGrids.count_bus_type(grid)
    @test count == (130, 17)
    cluster!(
        grid,
        round(Int, 0.3 * count[1]),
        round(Int, 0.05 * count[2]),
        round(Int, 0.5 * count[2])
    )
    @test length(grid.substations) == 47
    @test grid.substations[1].population == 6193
    @test grid.substations[47].generation == 16.1
    create_lines!(grid)
    @test grid.trans_lines[1].capacity == 4900
    @test length(grid.trans_lines) == 179
    pgrid = to_pandapower(grid)
    @test length(pgrid[:trafo]) == 36
    @test test_connectivity(grid.bus_conn, false)
    @test test_connectivity(grid.sub_conn, false)
    @test total_links(grid.bus_conn) == 179
    @test total_links(grid.sub_conn) == 69
    @test mean_node_deg(grid.bus_conn) ≈ 2.435374149659
    @test mean_shortest_path(sub_connectivity(grid) .> 0) ≈ 4.107741059302851
    @test mean_shortest_path(adjacency(grid)) ≈ 6.764218612615115
    @test mean_shortest_path(adjacency(grid), distance(buses(grid))) ≈ 133.1803034218734
    @test mean_shortest_path(
        (sub_connectivity(grid) .> 0),
        distance(substations(grid))
        ) ≈ 113.21283328057127
    @test robustness_line(grid.bus_conn, 10) > 1
    @test robustness_node(grid.bus_conn, 10) > 1
    @test robustness_line(grid.sub_conn, 10) > 1
    @test robustness_node(grid.sub_conn, 10) > 1
    SyntheticGrids.prepare_gen_data(GENDATAPATH, GENCOORDPATH, DUMMYPATH)
    genjson = JSON.parsefile(GENPATH)
    dummyjson = JSON.parsefile(DUMMYPATH)
    @test genjson == dummyjson
    grid1 = Grid(SEED)
    place_loads_from_zips!(
        grid1;
        latlim = (38, 40),
        longlim = (-89, -88),
    )
    place_gens_from_data!(
        grid1;
        latlim = (38, 40),
        longlim = (-89, -88),
    )
    connect!(grid1)
    count = SyntheticGrids.count_bus_type(grid1)
    cluster!(
        grid1,
        round(Int, 0.3 * count[1]),
        round(Int, 0.05 * count[2]),
        round(Int, 0.5 * count[2])
    )
    create_lines!(grid1)
    @test grid == grid1
    rm(DUMMYPATH)
    grid = Grid(true)
    @test size(grid.bus_conn) == (137, 137)
    @test size(grid.sub_conn) == (43, 43)
    add_gen!(grid, (11., 11.), 100, [11], ["ki"], reconnect = false)
    add_load!(grid, (19.,12.), 100, 100, 100, reconnect = true)
    add_substation!(
        grid,
        (11., 18.),
        [11],
        1000,
        6000,
        500,
        Set(),
        [grid.buses[139], grid.buses[138]]
    )
    @test size(grid.bus_conn) == (139, 139)
    @test size(grid.sub_conn) == (44, 44)
end
