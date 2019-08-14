"""
    bus_dist(b1::Bus, b2::Bus)

OBSOLETE. Compute geometrical distance between two buses. Obsolete now that
buses are located by latitude and longitude.
"""
function bus_dist(b1::Bus, b2::Bus)
    return sqrt((b1.coords[1] - b2.coords[1])^2 + (b1.coords[2] - b2.coords[2])^2)
end


"""
    haversine(c1::Tuple{T,T}, c2::Tuple{T,T}) where T <: Real

Compute distance between two pairs of latitude and longitude using the
haversine formula.
"""
function haversine(c1::Tuple{T,T}, c2::Tuple{T,T}) where T <: Real
    AVG_RADIUS = 6369.783 # Earth radius at 38.8 latitude (in km).
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

function Statistics.mean(v::Vector{<:LatLon})
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
function subs_dist(
    inds1::Vector{Int},
    inds2::Vector{Int},
    b_dist::Matrix{T},
    id2ind::Dict{Int, Int},
    sum_pops::Symmetric{Int}
) where T <: Real
    dist = 0
    totw = 0
    for i1 in inds1
        j1 = id2ind[i1]
        for i2 in inds2
            j2 = id2ind[i2]
            dist += sum_pops[j2, j1] * b_dist[j2, j1]
            totw += sum_pops[j2, j1]
        end
    end
    return dist/totw
end


function load_sub_dist(
    subs::Vector{Substation},
    b_dist::Matrix{Float64},
    id2ind::Dict{Int, Int},
    sum_pops::Symmetric{Int}
)
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
    for s1 in 1:(n - 1)
        for s2 in (s1 + 1):n
            dist_mat[s2, s1] = subs_dist(b_inds[s1], b_inds[s2], b_dist, id2ind, sum_pops)
        end
    end
    return Symmetric(dist_mat, :L)
end

function update_sub_dist!(
    dist_mat::Matrix{Float64},
    i::Int,
    subs::Vector{Substation},
    b_dist::Matrix{Float64},
    id2ind::Dict{Int, Int},
    sum_pops::Symmetric{Int}
)
    b_inds = sizehint!([], length(subs)) # indexes of buses in each substations
    for sub in subs
        temp = sizehint!(Int[],length(sub.grouping))
        for load in sub.grouping
          push!(temp,load.id)
        end
        push!(b_inds, temp)
    end

    n = length(subs)
    ran = [c for c in collect(1:n) if c != i]
    for s2 in ran
        dist_mat[s2, i] = subs_dist(b_inds[i], b_inds[s2], b_dist, id2ind, sum_pops)
        dist_mat[i, s2] = dist_mat[s2, i]
    end
end

function distance(buses::Vector{T}) where T <: Bus
    n = length(buses)
    dist_mat = fill(Inf, (n, n))
    for b1 in 1:(n - 1), b2 in (b1 + 1):n
        dist_mat[b1, b2] = haversine(buses[b1], buses[b2])
        dist_mat[b2, b1] = dist_mat[b1, b2]
    end
    return dist_mat
end

function update_distance!(dist_mat::Matrix{Float64}, buses::Vector{T}, i::Int) where T <: Bus
    n = length(buses)
    ran = [c for c in collect(1:n) if c != i]
    for b in ran
        dist_mat[b, i] = haversine(buses[b], buses[i])
        dist_mat[i, b] = dist_mat[b, i]
    end
    return dist_mat
end
