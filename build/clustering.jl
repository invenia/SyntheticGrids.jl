"""
    cluster_loads!(grid,nload)

Execute aglomerative clustering of all load buses into 'nload' substations
using geographic distances as metric.
"""
function cluster_loads!(grid::Grid, nload)
    loads = [bus for bus in buses(grid) if isa(bus, LoadBus)]
    for b in loads # creates 1 substation for each load bus
        push!(
            substations(grid),
            Substation(
                length(substations(grid)) + 1,
                b.coords,
                [b.voltage],
                b.load,
                0,
                b.population
            )
        )
        substations(grid)[end].grouping = [b]
    end
    load_dist = distance(loads)
    while length(substations(grid)) > nload
        sub_dist = load_sub_dist(substations(grid), loads, load_dist)
        n = length(substations(grid))
        ok = false
        const POPLIMIT = 200000 # population limit for a substation (=> 400MW load)
        i = j = -1
        while true
            m = indmin(sub_dist)
            if sub_dist[m] == Inf # Can't cluster anymore without going over POPLIMIT
                warn("Can't cluster anymore without going over population limit!")
                warn("Stopping clustering process before target number is reached")
                nload = Inf
                break
            end
            i, j = ind2sub(size(sub_dist), m)
            if (substations(grid)[i].population + substations(grid)[j].population) < POPLIMIT
                break
            else # this keeps loads from clustering beyond the population limit
                sub_dist[i, j] = Inf
                sub_dist[j, i] = Inf
            end
        end
        if nload == Inf
            break
        end
        merge!(grid, substations(grid)[i], substations(grid)[j])
    end
end

"""
    fill_gen!(sub, gens, MW)

Group generators ('gens') into substations ('subs') by geographic proximity
until generation threshold 'MW' is reached.
"""
function fill_gen!(sub, gens, MW)
    dists = [haversine(sub, g) for g in gens]
    while sub.generation < MW
        m = indmin(dists)
        sub.generation += gens[m].generation
        sub.voltages = [v for v in Set(vcat(sub.voltages, gens[m].voltage))]
        push!(sub.grouping, gens[m])
        deleteat!(dists, m)
        deleteat!(gens, m)
    end
end

"""
    cluster_load_gen!(grid, nboth)

Create 'nboth' substations with both load and generation, chosen randomly.
"""
function cluster_load_gen!(grid::Grid, nboth)
    gens = [bus for bus in buses(grid) if isa(bus, GenBus)]
    nb = 0
    srand(grid.seed + 66) # Just so we don't keep returning to the same point.
    s = rand(1:length(substations(grid)))
    rrange = e^(-1):PREC:e
    while nb < nboth
        while substations(grid)[s].generation > 0 # ensures we are drawing a
            s = rand(1:length(substations(grid))) # load-only substation
        end
        MW = rand(rrange) * substations(grid)[s].load
        fill_gen!(substations(grid)[s], gens, MW)
        nb += 1
    end
end

"""
    cluster_gens!(grid, ngen)

Execute aglomerative clustering of generators into 'ngen' substations.
"""
function cluster_gens!(grid::Grid, ngen)
    gen_both = GenBus[]
    for sub in substations(grid), bus in sub.grouping
        if isa(bus, GenBus)
            push!(gen_both, bus)
        end
    end
    gens = filter(buses(grid)) do bus
        isa(bus, GenBus) && !(bus in gen_both)
    end
    g_subs = Array{Substation}(length(gens))
    st = length(substations(grid))
    for (i, g) in enumerate(gens) # creates substations for the remaining generators
        substation = Substation(
            length(substations(grid)) + 1,
            g.coords, g.voltage,
            0,
            g.generation,
            0
        )
        substation.grouping = [g]
        push!(substations(grid), substation)
        g_subs[i] = substation
    end
    while length(g_subs) > ngen
        dists = distance(g_subs)
        m = indmin(dists)
        n = length(g_subs)
        ii, jj = ind2sub(size(dists), m)
        merge!(grid, substations(grid)[ii + st], substations(grid)[jj + st])
        deleteat!(g_subs, jj)
    end
end

"""
    cluster!(grid::Grid, nloads, nboth, ngens)

Cluster all nodes of a grid into 'nloads' load substations, 'ngens' generation
substations and 'nboth' substations with both load and generation.

REFERENCE: Birchfield, Adam B., et al. "Grid Structural Characteristics as
Validation Criteria for Synthetic Networks."
IEEE Transactions on Power Systems (2016).
"""
function cluster!(grid::Grid, nloads, nboth, ngens)
    nloads > 0 ? cluster_loads!(grid, nloads) : nothing
    nboth > 0 ? cluster_load_gen!(grid, nboth) : nothing
    ngens > 0 ? cluster_gens!(grid, ngens) : nothing
    (nloads + nboth + ngens) > 0 ? connect_subs!(grid) : nothing
    # The processes above should mess substation ids, so we need to correct them
    for i in length(substations(grid))
        substations(grid)[i].id = i
    end
end
