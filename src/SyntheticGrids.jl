module SyntheticGrids

# This needs to be set before `using PyCall` otherwise the numba library won't see it
# at load time ? :/
ENV["NUMBA_DISABLE_JIT"] = 1

using PyCall
using JSON
using JLD
using CSV
using AutoHashEquals
using Geodesy

const pp = PyNULL()

# https://github.com/JuliaPy/PyCall.jl#using-pycall-from-julia-modules
function __init__()
    copy!(pp, pyimport("pandapower"))
end

export
###### Types ######
Bus,
LoadBus,
GenBus,
Generator,
TransLine,
Substation,
Grid,
###### Main Functions ######
connect!,
cluster!,
create_lines!,
place_loads_from_zips!,
place_gens_from_data!,
add_load!,
add_gen!,
add_substation!,
###### Export Functions ######
to_pandapower,
###### Check Functions ######
total_links,
mean_node_deg,
cluster_coeff,
test_connectivity,
mean_shortest_path,
robustness_line,
robustness_node,
###### Getters ######
distance,
coords,
voltage,
technology,
capacity,
power_factor,
minimum_load,
summer_capacity,
winter_capacity,
status,
generators,
trans_lines,
population,
neighbors,
buses,
substations,
adjacency,
sub_connectivity,
end_points,
impedance

import Base.==
import Base.show
import Base.hash
import Base.mean
import Base.merge!
global const PREC = 1.0e-6 # Precision for random number generation
global const LARGE_GRID_N = 2000 # defines what is considered a large grid
abstract type Bus end

include("transline.jl")
include("generator.jl")
include("bus.jl")
include("distances.jl")
include("substation.jl")
include("grid.jl")
include("clustering.jl")
include("io.jl")
include("misc.jl")

end
