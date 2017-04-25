"""
    bus_dist(b1::Bus, b2::Bus)

OBSOLETE. Compute geometrical distance between two buses. Obsolete now that
buses are located by latitude and longitude.
"""
function bus_dist(b1::Bus, b2::Bus)
    return sqrt((b1.coords[1] - b2.coords[1])^2 + (b1.coords[2] - b2.coords[2])^2)
end


"""
    haversine{T<:Real}(c1::Tuple{T,T}, c2::Tuple{T,T})

Compute distance between two pairs of latitude and longitude using the
haversine formula.
"""
function haversine{T<:Real}(c1::Tuple{T,T}, c2::Tuple{T,T})
    const AVG_RADIUS = 6369.783 # Earth radius at 38.8 latitude (in km).
    d = 2 * AVG_RADIUS * asin(sqrt(sind((c2[1] - c1[1]) / 2)^2
        + cosd(c2[1]) * cosd(c1[1]) * sind((c2[2] - c1[2]) / 2)^2))
    return d
end

function haversine(c1::LatLon, c2::LatLon)
    return haversine((c1.lat, c1.lon), (c2.lat, c2.lon))
end

"""
    haversine(b1::Bus, b2::Bus)

Compute distance between two buses from their latitude and longitude using the
haversine formula.
"""
function haversine(b1::Bus, b2::Bus)
    return haversine(b1.coords, b2.coords)
end

function mean{T <: Real}(v::Vector{LatLon{T}})
    xmean = mean([c.lat for c in v])
    ymean = mean([c.lon for c in v])
    return LatLon(xmean, ymean)
end

"""
    subs_dist(inds1, inds2, buses, b_dist)

Compute distance between two substations using equation 1 from the REFERENCE.

REFERENCE: Birchfield, Adam B., et al. "Grid Structural Characteristics as
Validation Criteria for Synthetic Networks."
IEEE Transactions on Power Systems (2016).
"""
function subs_dist{T <: Real}(
    inds1::Vector{Int},
    inds2::Vector{Int},
    buses::Vector{LoadBus},
    b_dist::Matrix{T},
)
    dist = 0
    totw = 0
    for i1 in inds1, i2 in inds2
        totw += (buses[i1].population + buses[i2].population)
    end
    for i1 in inds1, i2 in inds2
        dist += ((buses[i1].population + buses[i2].population)
                * b_dist[i2,i1] / totw)
    end
    return dist
end

function load_sub_dist(subs, buses, b_dist)
    b_inds = sizehint!([], length(subs)) # indexes of buses in each substations
    for sub in subs
        temp = sizehint!(Int[],length(sub.grouping))
        for load in sub.grouping
          push!(temp,load.id)
        end
        push!(b_inds, temp)
    end
    n = length(subs)
    dist_mat = fill(Inf, (n, n))
    for s1 in 1:(n - 1), s2 in (s1 + 1):n
        dist_mat[s1, s2] = subs_dist(b_inds[s1], b_inds[s2], buses, b_dist)
        dist_mat[s2, s1] = dist_mat[s1, s2]
    end
    return dist_mat
end

function distance{T<:Bus}(buses::Vector{T})
    n = length(buses)
    dist_mat = fill(Inf, (n, n))
    for b1 in 1:(n - 1), b2 in (b1 + 1):n
        dist_mat[b1, b2] = haversine(buses[b1], buses[b2])
        dist_mat[b2, b1] = dist_mat[b1, b2]
    end
    return dist_mat
end
