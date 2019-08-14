nload(buses::Vector{<:Bus}) = length([b for b in buses if isa(b, LoadBus)])

ngen(buses::Vector{<:Bus}) = length([b for b in buses if isa(b, GenBus)])

count_bus_type(grid::Grid) = (nload(buses(grid)), ngen(buses(grid)))

#######################
### CHECK FUNCTIONS ###
#######################
# Everything here uses connectivty matrices as input instead of grids because
# all functions are supposed to work both for buses and for substations,
# effectively giving us two representations of the network (fine and coarse)
"""
    laplacian(con_mat::AbstractMatrix{<:Real})

Return the Laplacian matrix for a graph defined by adjacency matrix `con_mat`.
"""
function laplacian(conmat::AbstractMatrix{Bool})
    lap = zeros(Int, size(conmat))
    degree = sum(conmat, dims=1) # this gives the node degrees
    n = length(degree)
    for j in 1:(n - 1), i in (j + 1):n
        if conmat[i,j] > 0
            lap[i,j] = -1
        end
        lap[j,j] = degree[j]
    end
    lap[n,n] = degree[n]
    return Symmetric(lap, :L)
end

function laplacian(con_mat::AbstractMatrix{T}) where T <: Real
    return laplacian(con_mat .> zero(T)) # treating the substation graph as a simple graph.
end

# Shorter code but worse performance.
# function laplacian(con_mat::AbstractMatrix{<:Bool})
#     return sum(con_mat, 1) .* eye(Int, size(con_mat)[1]) - con_mat
# end

"""
    total_links(con_mat::AbstractMatrix{Bool})

Return the total number of links in a system with connectivity matrix `con_mat`.
"""
function total_links(con_mat::AbstractMatrix{Bool})
    return (1 / 2) * sum(con_mat)
end

"""
    total_links(con_mat::AbstractMatrix{Integer})

Return the total number of links in a system with connectivity matrix `con_mat`. The graph
will be treated as a simple graph.
"""
function total_links(con_mat::AbstractMatrix{T}) where T <: Real
    return total_links(con_mat .> zero(T)) # treating the substation graph as a simple graph.
end

"""
    mean_node_deg(con_mat::AbstractMatrix)

Return the average nodal degree of a system with connectivity matrix `con_mat`.
"""
function mean_node_deg(con_mat::AbstractArray{Bool, 2})
    return mean(sum(con_mat, dims=1))
end

function mean_node_deg(con_mat::AbstractArray{Int, 2})
    return mean_node_deg(con_mat .> 0) # treating the substation graph as a simple graph.
end

"""
    cluster_coeff(con_mat::AbstractMatrix)

Return the clustering coefficient of a system with connectivity matrix `con_mat`.
"""
function cluster_coeff(conmat::AbstractMatrix{Bool})
    degree = sum(conmat, 1) # this gives the node degrees
    clust_coeff = Vector{Float64}(length(degree))
    for i in 1:length(degree)
        maxedges = degree[i] * (degree[i] - 1) / 2
        cumul = 0
        for j in 1:length(degree)
            if conmat[i,j]
                cumul += dot(conmat[i,:],conmat[:,j])
            end
        end
        push!(clust_coeff, maxedges > 0 ? cumul/maxedges : 0)
    end
    return mean(clust_coeff)
end

function cluster_coeff(con_mat::AbstractMatrix{Int})
    return cluster_coeff(con_mat .> 0) # treating the substation graph as a simple graph.
end

"""
    cluster_coeff_degw(con_mat::AbstractMatrix{<:Real})

Return the degree-weighted clustering coefficient of a system with connectivity
matrix `con_mat`.
"""
function cluster_coeff_degw(conmat::AbstractMatrix{Bool})
    degree = sum(conmat, 1) # this gives the node degrees
    clust_coeff = Vector{Float64}(length(degree))
    for i in 1:length(degree)
        maxedges = degree[i] * (degree[i] - 1) / 2
        cumul = 0
        for j in 1:length(degree)
            if conmat[i, j]
                cumul += dot(conmat[i,:], conmat[:,j])
            end
        end
        push!(clust_coeff, maxedges > 0 ? degree[i] * cumul / maxedges : 0)
    end
    return mean(clust_coeff) / mean(degree)
end

function cluster_coeff_degw(con_mat::AbstractMatrix{Int})
    return cluster_coeff(con_mat .> 0) # treating the substation graph as a simple graph.
end

"""
    test_connectivity(con_mat::AbstractMatrix{<:Integer}, verb=true)

Return `true` if the system with connectivity matrix `con_mat`is connected and `false`
otherwise. If `verb=true` the result will be printed onscreen together with the Fiedler
eigenvalue.
"""
function test_connectivity(con_mat::AbstractMatrix{<:Integer}, verb=true)
    ZERO_TOL = 1e-13 # tolerance for zero due to numerical noise
    lap = laplacian(con_mat)
    evals = eigvals(lap)
    sort!(evals)
    if length(evals) > 1 && evals[2] > ZERO_TOL
        fiedler = evals[2]
        if verb
            println("Fiedler eigenvalue: $fiedler. System is connected.")
        end
        return true
    else
        if verb
            println("Fiedler eigenvalue equals zero. Network is NOT connected!")
        end
        return false
    end
end

"""
    dijkstra(con_mat::AbstractMatrix{<:Real}, v1)

Compute minimum distances to each node in the network with connectivity matrix
con_mat, starting from node v1, using Dijkstra's algorithm.
"""
function dijkstra(con_mat::AbstractMatrix{<:Real}, v1)
    n = size(con_mat, 1)
    c = v1
    visited = [v1]
    dist = [i == v1 ? 0 : Inf for i in 1:n]
    r_dist = [Inf for d in dist]
    while length(visited) < n
        neigh = [i for i in 1:n if con_mat[c, i] > 0 && !(i in visited)]
        for i in neigh
            if dist[i] > (dist[c] + con_mat[c, i])
                dist[i] = (dist[c] + con_mat[c, i])
                r_dist[i] = dist[i]
            end
        end
        c = indmin(r_dist)
        push!(visited, c)
        r_dist[c] = Inf
    end
    return dist
end

"""
    floyd_warshall(weights::AbstractMatrix{<:Real})

Compute minimum distances to each node in the network using Floyd-Warshall's algorithm. Here
we assume that the weigths matrix is symmetric.
"""
function floyd_warshall(weights::AbstractMatrix{<:Real})
    s = size(weights)
    dist = zeros(Float64, s)
    for j in 2:s[1]
        for i in 1:(j - 1)
            if weights[i, j] > 0
                dist[i, j] = weights[i, j] # we are assuming that weights is symmetric as it
                dist[j, i] = weights[i, j] # should be in our case.
            else
                dist[i, j] = Inf
                dist[j, i] = Inf
            end
        end
    end

    for k in 1:s[1]
        for i in 1:s[1]
            for j in 1:s[1]
                if dist[i, j] > dist[i, k] + dist[k, j]
                    dist[i, j] = dist[i, k] + dist[k, j]
                end
            end
        end
    end
    return dist
end

"""
    mean_shortest_path(conmat::AbstractMatrix{<:Real})

Compute the mean shortest distance between nodes for a network represented by `conmat`.
The `conmat` matrix should specify the distance between each connected node and be equal to
zero whenever two nodes are not connected. Using the adjacency matrix will return the mean
shortest path in hops. `conmat` is assumed to be symmetric. Floyd-Warshall's algorithm is
used to compute the distances.
"""
function mean_shortest_path(conmat::AbstractMatrix{<:Real})
    n = size(conmat, 1)
    return sum(floyd_warshall(conmat))/(n * n)
end

"""
    mean_shortest_path(conmat::AbstractMatrix{<:Integer}, distmat::AbstractMatrix{Float64})

Compute the mean shortest distance between nodes for a network with connectivity matrix
`conmat` and distances between nodes given by `distmat`. Distances are computed using
Floyd-Warshall's algorithm.
"""
function mean_shortest_path(
        conmat::AbstractMatrix{<:Integer},
        distmat::AbstractMatrix{Float64}
    )
    return mean_shortest_path(conmat.*distmat)
end

# """
#     mean_shortest_path(conmat::AbstractMatrix{<:Real})
#
# Compute the mean shortest distance between nodes for a network represented by `conmat`.
# The `conmat` matrix should specify the distance between each connected node and be equal to
# zero whenever two nodes are not connected. Distances are computed using Dijkstra's algorithm.
# """
# function mean_shortest_path(conmat::AbstractMatrix{<:Real})
#     n = size(conmat, 1)
#     mdist = [mean(dijkstra(conmat, i)) for i in 1:n]
#     return mean(mdist)
# end

# """
#     mean_shortest_path(conmat::AbstractMatrix{<:Real}, distmat::AbstractMatrix{<:Real})
#
# Compute the mean shortest distance between nodes for a network with connectivity matrix
# `conmat` and distances between nodes given by `distmat`. Distances are computed using
# Dijkstra's algorithm.
# """
# function mean_shortest_path(
#         conmat::AbstractMatrix{<:Real},
#         distmat::AbstractMatrix{<:Real}
#     )
#     return mean_shortest_path(conmat.*distmat)
# end

"""
    robustness_line(con_mat::AbstractMatrix{<:Integer}, n=1000)

Return the average number of transmission lines that have to be randomly removed from a
system with connectivity matrix `con_mat` before it becomes disconnected. Average is
computed over n iterations.
"""
function robustness_line(con_mat::AbstractMatrix{Bool}, n=1000)
    s = size(con_mat, 1)
    conns = []
    for i in 1:(s - 1), j in (i + 1):s
        con_mat[j, i] > 0 ? push!(conns, (i, j)) : nothing
    end
    total = 0
    for rep in 1:n
        connscopy = copy(conns)
        conmat = copy(con_mat)
        step = 0
        while test_connectivity(conmat, false)
            k = rand(1:length(connscopy))
            ij = connscopy[k]
            conmat[ij[1],ij[2]] = false
            conmat[ij[2],ij[1]] = false
            deleteat!(connscopy, k)
            step += 1
        end
        total += step
    end
    return total / n
end

function robustness_line(con_mat::AbstractMatrix{Int}, n=1000)
    return robustness_line(con_mat .> 0 , n)
end

"""
    robustness_node(con_mat::AbstractMatrix{<:Integer}, n=1000)

Return the average number of nodes that have to be randomly removed from a system with
connectivity matrix `con_mat` before it becomes disconnected. Average is computed over n
iterations.
"""
function robustness_node(con_mat::AbstractMatrix{Bool}, n=1000)
    s = size(con_mat, 1)
    total = 0
    for rep in 1:n
        len = s
        conmat = copy(con_mat)
        step = 0
        while test_connectivity(conmat, false)
            k = rand(1:len)
            mask = ones(Bool, len)
            mask[k] = false
            conmat = conmat[mask, mask]
            step += 1
            len -= 1
        end
        total += step
    end
    return total / n
end

function robustness_node(con_mat::AbstractMatrix{Int}, n=1000)
    return robustness_node(con_mat .> 0, n)
end
