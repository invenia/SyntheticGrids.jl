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
    const DUMMYPPC = joinpath(dirname(@__FILE__), "..", "data", "dummy.p")
    const GRIDPATH = joinpath(dirname(@__FILE__), "..", "data", "testgrid.json")

    grid = Grid(SEED)
    grid1 = Grid(SEED)

    @testset "Basic grid creation" begin
        place_loads_from_zips!(grid; latlim = (38, 40), longlim = (-89, -88))
        @test length(buses(grid)) == 130
        place_gens_from_data!(grid; latlim = (38, 40), longlim = (-89, -88))
        @test length(buses(grid)) == 147
        connect!(grid)
        @test sum(adjacency(grid)) == 358
        count = SyntheticGrids.count_bus_type(grid)
        @test count == (130, 17)
        create_lines!(grid)
        @test grid.trans_lines[1].capacity == 4900
        @test length(trans_lines(grid)) == 179
        @test test_connectivity(grid.bus_conn, false)
        @test total_links(grid.bus_conn) == 179
        @test mean_node_deg(grid.bus_conn) ≈ 2.435374149659
        @test mean_shortest_path(adjacency(grid)) ≈ 6.764218612615115
        @test mean_shortest_path(adjacency(grid), distance(buses(grid))) ≈ 133.1803034218734
        @test robustness_line(grid.bus_conn, 10) > 1
        @test robustness_node(grid.bus_conn, 10) > 1
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
        create_lines!(grid1)
        @test grid == grid1

        grid2 = Grid()
        function simple_func(t::Tuple)
            return t[1] > -t[2]/2
        end
        place_loads_from_zips!(grid2, simple_func)
        @test length(buses(grid2)) == 8171
        place_gens_from_data!(grid2, simple_func)
        @test length(buses(grid2)) == 9957

        grid3 = Grid(true)
        @test size(grid3.bus_conn) == (137, 137)
        @test size(grid3.sub_conn) == (43, 43)
        add_gen!(grid3, (11., 11.), 100, [11], ["ki"], reconnect = false)
        add_load!(grid3, (19.,12.), 100, 100, 100, reconnect = true)
        add_substation!(
            grid3,
            (11., 18.),
            [11],
            1000,
            6000,
            500,
            Set(),
            [grid3.buses[139], grid3.buses[138]]
        )
        @test size(grid3.bus_conn) == (139, 139)
        @test size(grid3.sub_conn) == (44, 44)
    end

    @testset "Clustering grids" begin
        count = SyntheticGrids.count_bus_type(grid)
        cluster!(
            grid,
            round(Int, 0.3 * count[1]),
            round(Int, 0.05 * count[2]),
            round(Int, 0.5 * count[2])
        )
        @test length(substations(grid)) == 47
        @test grid.substations[1].population == 6193
        @test grid.substations[47].generation == 16.1
        @test test_connectivity(grid.sub_conn, false)
        @test total_links(grid.sub_conn) == 69
        @test mean_shortest_path(sub_connectivity(grid) .> 0) ≈ 4.107741059302851
        @test mean_shortest_path(
            (sub_connectivity(grid) .> 0),
            distance(substations(grid))
            ) ≈ 113.21283328057127
        @test robustness_line(grid.sub_conn, 10) > 1
        @test robustness_node(grid.sub_conn, 10) > 1
        count = SyntheticGrids.count_bus_type(grid1)
        cluster!(
            grid1,
            round(Int, 0.3 * count[1]),
            round(Int, 0.05 * count[2]),
            round(Int, 0.5 * count[2])
        )
        @test (grid == grid1)
    end

    @testset "Input and Output" begin
        @testset "Pandapower interface" begin
            pgrid = to_pandapower(grid, DUMMYPPC)
            @test length(pgrid[:trafo]) == 36
            pgrid = SyntheticGrids.load_pp_grid(DUMMYPPC)
            @test length(pgrid[:trafo]) == 36
        end

        @testset "Generator data" begin
            SyntheticGrids.prepare_gen_data(GENDATAPATH, GENCOORDPATH, DUMMYPATH)
            genjson = JSON.parsefile(GENPATH)
            dummyjson = JSON.parsefile(DUMMYPATH)
            @test genjson == dummyjson
        end

        @testset "Save and Load Grid" begin
            save(grid, GRIDPATH)
            lgrid = load_grid(GRIDPATH)
            @test adjacency(grid) == adjacency(lgrid)
            @test sub_connectivity(grid) == sub_connectivity(lgrid)
            lgrid = 0 # trying to solve issues with Windows locking permissions to GRIDPATH
        end

        rm(DUMMYPATH)
        rm(DUMMYPPC)
        rm(GRIDPATH)
    end
end
