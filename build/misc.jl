nload{T<:Bus}(buses::Vector{T}) = length([b for b in buses if isa(b, LoadBus)])

ngen{T<:Bus}(buses::Vector{T}) = length([b for b in buses if isa(b, GenBus)])

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
function laplacian(con_mat::AbstractMatrix{<:Real})
    lap = zeros(Int, size(con_mat))
    conmat = con_mat .> 0 # adjusting for the case of substation connectivty matrix
    degree = sum(conmat, 1) # this gives the node degrees
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

"""
    total_links(con_mat::AbstractMatrix{<:Real})

Return the total number of links in a system with connectivity matrix `con_mat`.
"""
function total_links(con_mat::AbstractMatrix{<:Real})
    return (1 / 2) * trace(laplacian(con_mat))
end

"""
    mean_node_deg(con_mat::AbstractMatrix{<:Real})

Return the average nodal degree of a system with connectivity matrix `con_mat`.
"""
function mean_node_deg{T}(con_mat::AbstractArray{T,2})
    return (1 / size(con_mat, 1)) * trace(laplacian(con_mat))
end

"""
    cluster_coeff(con_mat::AbstractMatrix{<:Real})

Return the clustering coefficient of a system with connectivity matrix `con_mat`.
"""
function cluster_coeff(con_mat::AbstractMatrix{<:Real})
    conmat = con_mat .> 0 # adjusting for the case of substation connectivty matrix
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

"""
    cluster_coeff_degw(con_mat::AbstractMatrix{<:Real})

Return the degree-weighted clustering coefficient of a system with connectivity
matrix `con_mat`.
"""
function cluster_coeff_degw(con_mat::AbstractMatrix{<:Real})
    conmat = con_mat .> 0 # adjusting for the case of substation connectivty matrix
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

"""
    test_connectivity(con_mat::AbstractMatrix{<:Real}, verb=true)

Return `true` if the system with connectivity matrix `con_mat`is connected and `false`
otherwise. If `verb=true` the result will be printed onscreen together with the Fiedler
eigenvalue.
"""
function test_connectivity(con_mat::AbstractMatrix{<:Real}, verb=true)
    lap = laplacian(con_mat)
    evals = eigvals(lap)
    unique_vals = [v for v in Set(evals)]
    sort!(unique_vals)
    if length(unique_vals) > 1 && unique_vals[2] > 0
        fiedler = unique_vals[2]
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

function net_adj_mat(con_mat::AbstractMatrix{<:Real})
    lap = laplacian(con_mat)
    return Diagonal(lap) - lap
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
    mean_shortest_path(conmat::AbstractMatrix{<:Real})

Compute the mean shortest distance between nodes for a network represented by `conmat`.
The `conmat` matrix should specify the distance between each connected node and be equal to
zero whenever two nodes are not connected. Distances are computed using Dijkstra's algorithm.
"""
function mean_shortest_path(conmat::AbstractMatrix{<:Real})
    n = size(conmat, 1)
    mdist = [mean(dijkstra(conmat, i)) for i in 1:n]
    return mean(mdist)
end

"""
    mean_shortest_path(conmat::AbstractMatrix{<:Real}, distmat::AbstractMatrix{<:Real})

Compute the mean shortest distance between nodes for a network with connectivity matrix
`conmat` and distances between nodes given by `distmat`. Distances are computed using
Dijkstra's algorithm.
"""
function mean_shortest_path(
        conmat::AbstractMatrix{<:Real},
        distmat::AbstractMatrix{<:Real}
    )
    return mean_shortest_path(conmat.*distmat)
end

"""
    robustness_line(con_mat::AbstractMatrix{<:Real}, n=1000)

Return the average number of transmission lines that have to be randomly removed from a
system with connectivity matrix `con_mat` before it becomes disconnected. Average is
computed over n iterations.
"""
function robustness_line(con_mat::AbstractMatrix{<:Real}, n=1000)
    s = size(con_mat, 1)
    conns = []
    for i in 1:(s - 1), j in (i + 1):s
        con_mat[j, i] > 0 ? push!(conns, (i, j)) : nothing
    end
    total = 0
    for rep in 1:n
        connscopy = copy(conns)
        conmat = con_mat .> 0 # adjusting for the case of substation connectivty matrix
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

"""
    robustness_node(con_mat::AbstractMatrix{<:Real}, n=1000)

Return the average number of nodes that have to be randomly removed from a system with
connectivity matrix `con_mat` before it becomes disconnected. Average is computed over n
iterations.
"""
function robustness_node(con_mat::AbstractMatrix{<:Real}, n=1000)
    s = size(con_mat, 1)
    total = 0
    for rep in 1:n
        nodes = collect(1:s)
        conmat = con_mat .> 0 # adjusting for the case of substation connectivty matrix
        step = 0
        while test_connectivity(conmat, false)
            k = rand(1:length(nodes))
            i = nodes[k]
            for j in 1:s
                conmat[i,j] = false
                conmat[j,i] = false
            end
            deleteat!(nodes, k)
            step += 1
        end
        total += step
    end
    return total / n
end
