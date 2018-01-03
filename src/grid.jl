@auto_hash_equals mutable struct Grid
    seed::Integer # to allow reproducibility
    buses::Vector{Bus}
    trans_lines::Vector{TransLine}
    substations::Vector{Substation}
    bus_conn::AbstractArray{Bool,2} # bus connectivity matrix
    sub_conn::AbstractArray{Int,2} # this is slightly more than a substation connectivity
                # matrix, as it counts how many connections are there between substations.

    function Grid(
        seed::Integer,
        bus = [],
        line = [],
        sub = [],
        conm = falses(0, 0),
        subm = zeros(Int, 0, 0),
    )
         return new(seed, bus, line, sub, conm, subm)
    end
end

buses(g::Grid) = g.buses
trans_lines(g::Grid) = g.trans_lines
substations(g::Grid) = g.substations
adjacency(g::Grid) = g.bus_conn
sub_connectivity(g::Grid) = g.sub_conn

"""
    add_load!(grid::Grid, args...; reconnect = false)

Add LoadBus to `grid` by calling the `LoadBus(coords, load, volt, pop, connected_to =
Set(), connections = Set())` method with `args...` as arguments. If `reconnect = true`, all
connections will be remade.
"""
function add_load!(grid::Grid, args...; reconnect = false)
    bus = LoadBus(length(grid.buses)+1, args...)
    add_bus!(grid, bus, reconnect = reconnect)
end

"""
    add_gen!(grid::Grid, args...; reconnect = false)

Add Genbus to `grid` by calling the `GenBus(coords, gen, volt, tech, connected_to =
Set(), connections = Set(), pfactor = -1, summgen = -1, wintgen = -1, gens = [])` method
with `args...` as arguments. If `reconnect = true`, all connections will be remade.
"""
function add_gen!(grid::Grid, args...; reconnect = false)
    bus = GenBus(length(grid.buses)+1, args...)
    add_bus!(grid, bus, reconnect = reconnect)
end

function add_bus!(grid::Grid, bus::Union{LoadBus, GenBus}; reconnect = false)
    # The ids are tricky in that they have to correspond to the actual position of the bus
    # in the vector, but, at the same time, it is convenient to have all buses of the same
    # type grouped and ordered.
    loads = [b for b in grid.buses if isa(b, LoadBus)]
    gens = [b for b in grid.buses if isa(b, GenBus)]
    load_ids = [b.id for b in loads]
    gen_ids = [b.id for b in gens]
    ls = sortperm(load_ids)
    gs = sortperm(gen_ids)
    loads = loads[ls]
    gens = gens[gs]
    ml_id = isempty(loads) ? 0 : loads[end].id
    mg_id = isempty(gens) ? 0 : gens[end].id
    first = mg_id > ml_id ? LoadBus : GenBus
    second_empty = isempty(loads) || isempty(gens)
    if !second_empty && isa(bus, first)
        if first == LoadBus
            # Add new LoadBus
            push!(
                loads,
                LoadBus(
                    ml_id + 1,
                    bus.coords,
                    bus.load,
                    bus.voltage,
                    bus.population,
                    bus.connected_to,
                    bus.connections
                )
            )
            # Update IDs of all GenBuses
            new_gens = GenBus[]
            for g in gens
                push!(
                    new_gens,
                    GenBus(
                        g.id + 1,
                        g.coords,
                        g.generation,
                        g.voltage,
                        g.tech_type,
                        g.connected_to,
                        g.connections,
                        g.pfactor,
                        g.summgen,
                        g.wintgen,
                        g.gens
                    )
                )
            end
            # Update grid
            grid.buses = vcat(loads, new_gens)
        else
            # Add new GenBus
            push!(
                gens,
                GenBus(
                    mg_id + 1,
                    bus.coords,
                    bus.generation,
                    bus.voltage,
                    bus.tech_type,
                    bus.connected_to,
                    bus.connections,
                    bus.pfactor,
                    bus.summgen,
                    bus.wintgen,
                    bus.gens
                )
            )
            # Update the IDs of all LoadBuses
            new_loads = LoadBus[]
            for l in loads
                push!(
                    new_loads,
                    LoadBus(
                        l.id + 1,
                        l.coords,
                        l.load,
                        l.voltage,
                        l.population,
                        l.connected_to,
                        l.connections
                    )
                )
            end
            # Update grid
            grid.buses = vcat(gens, new_loads)
        end
    elseif bus.id != length(grid.buses) + 1
        if isa(bus, LoadBus)
            push!(
                grid.buses,
                LoadBus(
                    length(grid.buses) + 1,
                    bus.coords,
                    bus.load,
                    bus.voltage,
                    bus.population,
                    bus.connected_to,
                    bus.connections
                )
            )
        elseif isa(bus, GenBus)
            push!(
                buses(grid),
                GenBus(
                    length(grid.buses) + 1,
                    bus.coords,
                    bus.generation,
                    bus.voltage,
                    bus.tech_type,
                    bus.connected_to,
                    bus.connections,
                    bus.pfactor,
                    bus.summgen,
                    bus.wintgen,
                    bus.gens
                )
            )
        end
    else
        push!(grid.buses, bus)
    end
    if reconnect
        # Break previous connections
        grid.bus_conn = falses(length(grid.buses), length(grid.buses))
        for b in grid.buses
            empty!(b.connected_to)
            empty!(b.connections)
        end
        # Redo connections including the new bus
        connect!(grid)
    elseif size(grid.bus_conn)[1] > 0 # Grid has already been connected
        # Expand connectivity matrix to support the new bus
        bottom_row = zeros(length(grid.buses) - 1)
        right_col = zeros(length(grid.buses))
        grid.bus_conn = vcat(grid.bus_conn, bottom_row')
        grid.bus_conn = hcat(grid.bus_conn, right_col)
        # Check if there are connections to be made
        for b in bus.connected_to
            if b.id > 0
                grid.bus_conn[b.id, bus.id] = 1
                grid.bus_conn[bus.id, b.id] = 1
            end
        end
    end
end

"""
    add_substation!(grid::Grid, args...; reconnect = false)

Add Substation to `grid` by calling the `Substation(coords, volts, load, gen, pop,
con = Set(), group = [])` method with `args...` as arguments.
"""
function add_substation!(grid::Grid, args...)
    sub = Substation(length(grid.substations)+1, args...)
    push!(grid.substations, sub)
    # Break previous connections
    grid.sub_conn = zeros(length(grid.substations), length(grid.substations))
    for s in grid.substations
        s.connected_to = Set()
    end
    # Reconnect
    connect_subs!(grid)
end

Grid() = Grid(rand(0:round(Int, typemax(Int) / 2)))

function mean(v::Vector{LatLon})
    xmean = mean([c.lat for c in v])
    ymean = mean([c.lon for c in v])
    return LatLon(xmean, ymean)
end

"""
    twst!(grid, k)

Execute Tunable Weight Spannning Tree algorithm for the grid using weight 'k'.

REFERENCE: Soltan, Saleh, and Gil Zussman. "Generation of synthetic spatially
embedded power grid networks." arXiv:1508.04447 [cs.SY], Aug. 2015.
"""
function twst!(grid::Grid, k)
    node_locations = [b.coords for b in buses(grid)]
    n = length(buses(grid))
    ind = collect(1:n)
    permut = Vector{Int}(n)
    unorm_prob = [haversine(bus.coords, mean(node_locations))^(-k) for bus in buses(grid)]
    srand(grid.seed + 15)
    rand_range = 0:PREC:(1-PREC)
    for step in 1:n # sets node permutation
        probs = unorm_prob / sum(unorm_prob)
        rd = rand(rand_range)
        cumul = 0.0
        bus_ind = -1
        for i in 1:length(probs) # selects bus
            cumul += probs[i]
            if cumul >= rd
                bus_ind = i
                break
            end
        end
        permut[step] = ind[bus_ind]
        deleteat!(unorm_prob, bus_ind)
        deleteat!(ind, bus_ind)
    end

    bus_conn = zeros(Bool, n, n)
    for i in 2:length(permut) # builds connections
        distances = [
            haversine(node_locations[permut[i]],
            node_locations[permut[j]]) for j in 1:(i-1)
            ]
        ind = findmin(distances)[2] # nearest neighbor with j < i
        bus_conn[permut[i], permut[ind]] = true
        bus_conn[permut[ind], permut[i]] = true
    end

    grid.bus_conn = sparse(bus_conn)
end

"""
    reinforce!(grid, m, a, b, g, t, N)

Execute the Reinforcement procedure to increase robustness of a grid.

# Arguments:

* m: total number of connections
* a, b, g, t: paramaters of the model (see REFERENCE)
* N: number of nearest neighbors in the average distance computation

REFERENCE: Soltan, Saleh, and Gil Zussman. "Generation of synthetic spatially
embedded power grid networks." arXiv:1508.04447 [cs.SY], Aug. 2015.
"""
function reinforce!(grid::Grid, m, a, b, g, t, N)
    n = length(buses(grid))
    if N > (n - 1)
        throw(ArgumentError("ERROR in Reinforcement procedure (reinforce!()):",
              "Number of nearest neighbors has to be less than the number or nodes"))
    end
    if m < 0
        m = round(Int, 1.22 * n) # using Overbye's approximation here
    end
    av_dist = sizehint!(Real[], length(buses(grid)))
    for bus in buses(grid) # computes averge distances to neighbors
        dists = [haversine(bus, b2) for b2 in buses(grid)]
        sort!(dists)
        deleteat!(dists, 1) # minimum distance will always be the self-distance = 0
        push!(av_dist, mean(dists[1:N]))
    end
    unorm_probs = Real[]
    inds = Int[]
    degrees = sum(grid.bus_conn, 1) # this gives the node degrees
    degrees = float(reshape(degrees, length(degrees))) # Just so we can zip it later. Using
    # float here in order to avoid problems when doing Int^(-t).
    if n >= LARGE_GRID_N # only sample nodes with degree < 3
        unorm_probs = [av_dist[i]^(-a) for i in 1:n if degrees[i] < 3]
        inds = [i for i in 1:n if degrees[i] < 3]
    else # sample all nodes
        unorm_probs = [degrees[i]^(-t)*av_dist[i]^(-a) for i in 1:n]
        inds = collect(1:n)
    end
    srand(grid.seed + 479) # Just so we don't keep returning to the same point.
    rand_range = 0:PREC:(1-PREC)
    count = 1
    while count <= (m - n + 1)
        probs = unorm_probs / sum(unorm_probs)
        rd = rand(rand_range)
        cumul = 0.0
        bus_ind = -1
        for i in 1:length(probs) #selects bus
            cumul += probs[i]
            if cumul >= rd
                bus_ind = i # selected bus is inds[bus_ind]
                break
            end
        end
        b_dists = [haversine(buses(grid)[inds[bus_ind]], bus) for bus in buses(grid)]
        b_probs = [dist^(-b) * deg^g for (dist, deg) in zip(b_dists, degrees)]
        b_probs = [bp == Inf ? 0 : bp for bp in b_probs] # If buses overlap we may
                      # have multiple Infs (besides the Inf from the self-distance).
                      # Ideally buses should never overlap!
        b_probs = b_probs / sum(b_probs)
        rd = rand(rand_range)
        cumul = 0.0
        b_ind = -1
        for i in 1:length(b_probs) # draws node to which connect
            cumul += b_probs[i]
            if cumul >= rd
                b_ind = i
                break
            end
        end
        if b_ind == inds[bus_ind] # can't connect a bus to itself
            if b_ind == n
                b_ind -= 1
            else
                b_ind += 1
            end
        end
        if grid.bus_conn[b_ind, inds[bus_ind]] == false # not yet connected
            grid.bus_conn[b_ind, inds[bus_ind]] = true
            grid.bus_conn[inds[bus_ind], b_ind] = true
            deleteat!(unorm_probs, bus_ind)
            deleteat!(inds, bus_ind)
            count += 1
        end
    end
end

"""
    connect_buses!(grid)

Create reference links between buses that are connected according to the
connectivity matrix of the grid.
"""
function connect_buses!(grid::Grid)
    for j in 1:(length(buses(grid)) - 1), i in (j + 1):length(buses(grid))
        if grid.bus_conn[i, j] == true
            push!(buses(grid)[i].connected_to, buses(grid)[j])
            push!(buses(grid)[j].connected_to, buses(grid)[i])
        end
    end
end

"""
    connect!(grid::Grid; k=2.5, m=-1, a=1, b=3.2, g=2.5, t=2, N=10)

Connect buses in an electric grid.

# Arguments:
* k: weight of the spanning tree (see REFERENCE)
* m: total number of connections (computed from the # of nodes if default)
* a, b, g, t: paramaters of the model (see REFERENCE)
* N: number of nearest neighbors in the average distance computation

REFERENCE: Soltan, Saleh, and Gil Zussman. "Generation of synthetic spatially
embedded power grid networks." arXiv:1508.04447 [cs.SY], Aug. 2015.
"""
function connect!(
    grid::Grid;
    k=2.5,
    m=-1,
    a=1,
    b=3.2,
    g=2.5,
    t=2.0,
    N=10,
)
    twst!(grid, k)
    reinforce!(grid, m, a, b, g, t, N)
    connect_buses!(grid)
end


"""
    create_lines!(grid::Grid; impedfunc = linear_imped, capfunc = volt_cap)

Create transmission lines for a synthetic grid. Line impedancies are determined
by impedfunc and line capacities are calculated via capfunc.
This function will only work after the grid has been connected through the
connect!() function.
"""
function create_lines!(grid::Grid; impedfunc = linear_imped, capfunc = volt_cap)
    srand(grid.seed + 16) # Just so we don't keep returning to the same point.
    for j in 1:(length(buses(grid)) - 1), i in (j + 1):length(buses(grid))
        if grid.bus_conn[i, j] == true
            a, b = buses(grid)[i], buses(grid)[j]
            push!(grid.trans_lines, TransLine(a, b; impedfunc=impedfunc, capfunc=capfunc))
            push!(buses(grid)[i].connections, grid.trans_lines[end])
            push!(buses(grid)[j].connections, grid.trans_lines[end])
        end
    end
end

function merge!(grid::Grid, sub1::Substation, sub2::Substation)
    if sub1.load > 0 || sub2.load > 0
        tot_pop = (sub1.population + sub2.population)
        lat = (sub1.population * sub1.coords.lat + sub2.population
            * sub2.coords.lat) / tot_pop
        lon = (sub1.population * sub1.coords.lon + sub2.population
            * sub2.coords.lon) / tot_pop
        sub1.coords = LatLon(lat, lon)
    else
        tot_gen = sub1.generation + sub2.generation
        lat = (sub1.generation * sub1.coords.lat + sub2.generation
            * sub2.coords.lat) / tot_gen
        lon = (sub1.generation * sub1.coords.lon + sub2.generation
            * sub2.coords.lon) / tot_gen
        sub1.coords = LatLon(lat, lon)
    end
    for v in sub2.voltages
        if !(v in sub1.voltages)
            push!(sub1.voltages, v)
        end
    end
    sub1.load += sub2.load
    sub1.generation += sub2.generation
    sub1.population += sub2.population
    sub1.grouping = vcat(sub1.grouping, sub2.grouping)
    deleteat!(substations(grid), findfirst(substations(grid), sub2))
end

function merge!(grid::Grid, i1::Int, i2::Int)
    sub1 = substations(grid)[i1]
    sub2 = substations(grid)[i2]
    if sub1.load > 0 || sub2.load > 0
        tot_pop = (sub1.population + sub2.population)
        lat = (sub1.population * sub1.coords.lat + sub2.population
            * sub2.coords.lat) / tot_pop
        lon = (sub1.population * sub1.coords.lon + sub2.population
            * sub2.coords.lon) / tot_pop
        sub1.coords = LatLon(lat, lon)
    else
        tot_gen = sub1.generation + sub2.generation
        lat = (sub1.generation * sub1.coords.lat + sub2.generation
            * sub2.coords.lat) / tot_gen
        lon = (sub1.generation * sub1.coords.lon + sub2.generation
            * sub2.coords.lon) / tot_gen
        sub1.coords = LatLon(lat, lon)
    end
    for v in sub2.voltages
        if !(v in sub1.voltages)
            push!(sub1.voltages, v)
        end
    end
    sub1.load += sub2.load
    sub1.generation += sub2.generation
    sub1.population += sub2.population
    sub1.grouping = vcat(sub1.grouping, sub2.grouping)
    deleteat!(substations(grid), i2)
end

"""
    connect_subs!(grid)

Determine which substations are connected to each other through their buses.
"""
function connect_subs!(grid::Grid)
    n = length(substations(grid))
    grid.sub_conn = zeros(Int, n, n)
    aux = copy(grid.bus_conn)
    aux = convert(Matrix{Int}, aux)
    for s in 1:n
        for j in 1:length(buses(grid)), b in substations(grid)[s].grouping
            if grid.bus_conn[b.id, j] == true
                aux[b.id, j] = s
            end
        end
    end

    for i in 1:(length(buses(grid)) - 1)
        for j in (i + 1):length(buses(grid))
            if aux[i, j] != aux[j, i]
                grid.sub_conn[aux[i, j], aux[j, i]] += 1 # This approach
                grid.sub_conn[aux[j, i], aux[i, j]] += 1 # not only gives
                # which substation is connected to which, but also gives how many
                # connections are there between them, i.e., gives a measure of
                # strength for the connection.
                push!(
                    substations(grid)[aux[i, j]].connected_to,
                    substations(grid)[aux[j, i]]
                )
                push!(
                    substations(grid)[aux[j, i]].connected_to,
                    substations(grid)[aux[i, j]]
                )
            end
        end
    end
    grid.sub_conn = sparse(grid.sub_conn)
end

function Grid(
        sub::Bool;
        lat = (33,35),
        long = (-95,-93),
    )
    grid = Grid()
    place_loads_from_zips!(grid; latlim = lat, longlim = long)
    place_gens_from_data!(grid; latlim = lat, longlim = long)
    connect!(grid)
    count = count_bus_type(grid)
    if sub # This is just one very crude approach to having substations automatically
           # created. Good for testing.
        cluster!(
            grid,
            round(Int, 0.3 * count[1]),
            round(Int, 0.05 * count[2]),
            round(Int, 0.5 * count[2])
        )
    end
    create_lines!(grid)
    return grid
end
