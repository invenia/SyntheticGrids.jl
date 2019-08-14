"""
    cluster_loads!(grid,nload)

Execute agglomerative clustering of all load buses into 'nload' substations
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
    # Now let's avoid computing some things several times
    load_dist = distance(loads)
    all_ids = [l.id for l in loads]
    id2ind = Dict{Int, Int}()
    for id in all_ids
        id2ind[id] = findfirst(==(id), all_ids)
    end
    m = length(loads)
    sum_pops = fill(0, (m, m))
    for s1 in 1:(m - 1), s2 in (s1 + 1):m
        sum_pops[s2, s1] = loads[s1].population + loads[s2].population
    end
    sum_pops = Symmetric(sum_pops, :L)
    sub_dist = Matrix(load_sub_dist(substations(grid), load_dist, id2ind, sum_pops))
    POPLIMIT = 200000 # population limit for a substation (=> 400MW load)

    while length(substations(grid)) > nload
        n = length(substations(grid))
        ok = false
        i = j = -1
        while true
            m = argmin(sub_dist)
            if sub_dist[m] == Inf # Can't cluster anymore without going over POPLIMIT
                warn("Can't cluster anymore without going over population limit!")
                warn("Stopping clustering process before target number is reached")
                nload = Inf
                break
            end
            i, j = Tuple(m)
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
        merge!(grid, i, j)
        # Update substation distances
        mask = fill(true, (n, n))
        mask[j, :] .= false # Eliminate rows and columns corresponding to the
        mask[:, j] .= false # deleted substation
        sub_dist = reshape(sub_dist[mask], (n - 1, n - 1)) # Remove old dists
        i = i > j ? i - 1 : i # Update i to the new length
        update_sub_dist!(sub_dist, i, substations(grid), load_dist, id2ind, sum_pops)
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
        m = argmin(dists)
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
    seed!(grid.seed + 66) # Just so we don't keep returning to the same point.
    s = rand(1:length(substations(grid)))
    e = MathConstants.e
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
    g_subs = Array{Substation}(undef, length(gens))
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
    dists = distance(g_subs)
    while length(g_subs) > ngen
        m = argmin(dists)
        n = length(g_subs)
        ii, jj = Tuple(m)
        merge!(grid, ii + st, jj + st)
        deleteat!(g_subs, jj)
        # Update substation distances
        mask = fill(true, (n, n))
        mask[jj, :] .= false # Eliminate rows and columns corresponding to the
        mask[:, jj] .= false # deleted substation
        dists = reshape(dists[mask], (n - 1, n - 1)) # Remove old dists
        ii = ii > jj ? ii - 1 : ii # Update i to the new length
        update_distance!(dists, g_subs, ii)
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
