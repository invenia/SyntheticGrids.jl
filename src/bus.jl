struct LoadBus <: Bus
    id::Int # for internal use only, don't mess with this.
    coords::LatLon # latitude and longitude
    load::Real #in MW
    voltage::Real #in kV
    population::Integer
    connected_to::Set{Bus}
    connections::Set{TransLine}

    function LoadBus(id, coords, load, volt, pop, connected_to = Set(), connections = Set())
        return new(id, coords, load, volt, pop, connected_to, connections)
    end
end

function LoadBus(
        id,
        coords::Tuple{T,T} where T <: Real,
        load,
        volt,
        pop,
        connected_to = Set(),
        connections = Set()
    )
    return LoadBus(
        id,
        LatLon(coords[1], coords[2]),
        load,
        volt,
        pop,
        connected_to,
        connections
    )
end

function LoadBus(coords::LatLon, load, volt, pop, connected_to = Set(), connections = Set())
    return LoadBus(-1, coords, load, volt, pop, connected_to, connections)
end

function LoadBus(
        coords::Tuple{T,T} where T <: Real,
        load,
        volt,
        pop,
        connected_to = Set(),
        connections = Set()
    )
    return LoadBus(
        -1,
        LatLon(coords[1], coords[2]),
        load,
        volt,
        pop,
        connected_to,
        connections
    )
end

function show(io::IO, bus::LoadBus)
    println(io, "LoadBus(")
    println(io, "\tid=$(bus.id),")
    println(io, "\tcoords=$(bus.coords),")
    println(io, "\tload=$(bus.load)")
    println(io, "\tvoltage=$(bus.voltage),")
    println(io, "\tpopulation=$(bus.population), ")
    if length(bus.connected_to) > 0
        println(io, "\tconnected_to=Set{Bus}(...)")
    else
        println(io, "\tconnected_to=Set{Bus}()")
    end
    if length(bus.connections) > 0
        println(io, "\tconnections=Set{TransLine}(...)")
    else
        println(io, "\tconnections=Set{Transline}()")
    end
    print(io, ")")
end

function hash(a::LoadBus, h::UInt)
    h = hash(a.id, hash(a.coords, hash(a.load, hash(a.voltage, hash(a.population, h)))))
end

struct GenBus <: Bus # this actually represents a power plant
    id::Int # for internal use only, don't mess with this.
    coords::LatLon # latitude and longitude
    generation::Real # MW
    voltage::Vector{Real} # kV
    tech_type::Vector{AbstractString} # for all generators in a plant
    connected_to::Set{Bus}
    connections::Set{TransLine}
    pfactor::Real
    summgen::Real # MW
    wintgen::Real # MW
    gens::Vector{Generator}

    function GenBus(
        id,
        coords,
        gen,
        volt,
        tech,
        connected_to = Set(),
        connections = Set(),
        pfactor = -1,
        summgen = -1,
        wintgen = -1,
        gens = []
    )
        return new(
            id,
            coords,
            gen,
            volt,
            tech,
            connected_to,
            connections,
            pfactor,
            summgen,
            wintgen,
            gens
        )
    end
end

function GenBus(
    id,
    coords::Tuple{T,T} where T <: Real,
    gen,
    volt,
    tech,
    connected_to = Set(),
    connections = Set(),
    pfactor = -1,
    summgen = -1,
    wintgen = -1,
    gens = []
)
    return GenBus(
        id,
        LatLon(coords[1], coords[2]),
        gen,
        volt,
        tech,
        connected_to,
        connections,
        pfactor,
        summgen,
        wintgen,
        gens
    )
end

function GenBus(
    coords::LatLon,
    gen,
    volt,
    tech,
    connected_to = Set(),
    connections = Set(),
    pfactor = -1,
    summgen = -1,
    wintgen = -1,
    gens = []
)
    return GenBus(
        -1,
        coords,
        gen,
        volt,
        tech,
        connected_to = Set(),
        connections = Set(),
        pfactor,
        summgen,
        wintgen,
        gens
    )
end

function GenBus(
    coords::Tuple{T,T}  where T <: Real,
    gen,
    volt,
    tech,
    connected_to = Set(),
    connections = Set(),
    pfactor = -1,
    summgen = -1,
    wintgen = -1,
    gens = []
)
    return GenBus{T <: Real}(
        -1,
        coords::Tuple{T,T},
        gen,
        volt,
        tech,
        connected_to,
        connections,
        pfactor,
        summgem,
        wintgen,
        gens = []
    )
end

function GenBus(
    id,
    plant::AbstractArray{Generator},
    on::AbstractArray{Generator},
    pfactorfunc
)
    return GenBus(
        id,
        plant[1].coords,
        sum([g.cap for g in on]),
        plant[1].volt,
        [t for t in Set([g.tech for g in plant])],
        Set(),
        Set(),
        pfactorfunc(plant),
        sum([g.scap for g in on]),
        sum([g.wcap for g in on]),
        plant,
    )
end

function GenBus(id, plant::AbstractArray{Generator}, pfactorfunc)
    on = filter(g -> g.status in ["OP", "SP"], plant)
    return GenBus(id, plant, on, pfactorfunc)
end

function hash(a::GenBus, h::UInt)
    h = hash(a.id, h)
    h = hash(a.gens, h)
    h = hash(a.wintgen, h)
    h = hash(a.summgen, h)
    h = hash(a.pfactor, h)
    h = hash(a.tech_type, h)
    h = hash(a.voltage, h)
    h = hash(a.generation, h)
    return hash(a.coords, h)
end

function show(io::IO, bus::GenBus)
    println(io, "GenBus(")
    println(io, "\tid=$(bus.id)")
    println(io, "\tcoords=$(bus.coords),")
    println(io, "\tgeneration=$(bus.generation)")
    println(io, "\tvoltage=$(bus.voltage),")
    println(io, "\ttech_type=$(bus.tech_type), ")
    if length(bus.connected_to) > 0
        println(io, "\tconnected_to=Set{Bus}(...)")
    else
        println(io, "\tconnected_to=Set{Bus}()")
    end
    if length(bus.connections) > 0
        println(io, "\tconnections=Set{TransLine}(...)")
    else
        println(io, "\tconnections=Set{Transline}()")
    end
    println(io, "\tpfactor=$(bus.pfactor)")
    println(io, "\tsummgen=$(bus.summgen)")
    println(io, "\twintgen=$(bus.wintgen)")
    println(io, "\tgens=$(bus.gens)")
    print(io, ")")
end


function ==(a::Union{LoadBus, GenBus}, b::Union{LoadBus, GenBus})
    typeof(a) != typeof(b) && return false
    compared = Set{Tuple{UInt, UInt}}()

    bus_equal(a::GenBus, b::LoadBus) = false
    bus_equal(a::LoadBus, b::GenBus) = false

    function bus_equal(a::LoadBus, b::LoadBus)
        a.id == b.id &&
        a.coords == b.coords &&
        a.load == b.load &&
        a.voltage == b.voltage &&
        a.population == b.population &&
        hash(a.connections) == hash(b.connections) &&
        length(a.connected_to) == length(b.connected_to) || return false

        pair = (objectid(a), objectid(b))
        if !(pair in compared)
            push!(compared, pair)

            all(bus_equal(a, b) for (a, b) in zip(a.connected_to, b.connected_to)) ||
            return false
        end
        return true
    end

    function bus_equal(a::GenBus, b::GenBus)
        a.id == b.id &&
        a.coords == b.coords &&
        a.generation == b.generation &&
        a.voltage == b.voltage &&
        a.tech_type == b.tech_type &&
        hash(a.connections) == hash(b.connections) &&
        length(a.connected_to) == length(b.connected_to) &&
        a.pfactor == b.pfactor &&
        a.summgen == b.summgen &&
        a.wintgen == b.wintgen &&
        a.gens == b.gens || return false

        pair = (objectid(a), objectid(b))
        if !(pair in compared)
            push!(compared, pair)

            all(bus_equal(a, b) for (a, b) in zip(a.connected_to, b.connected_to)) ||
            return false
        end
        return true
    end

    return bus_equal(a, b)
end

capacity(b::GenBus) = b.generation
technology(b::GenBus) = b.tech_type
power_factor(b::GenBus) = b.pfactor
summer_capacity(b::GenBus) = b.summgen
winter_capacity(b::GenBus) = b.wintgen
generators(b::GenBus) = b.gens
trans_lines(b::Union{LoadBus, GenBus}) = b.connections
coords(b::Bus) = b.coords
load(b::LoadBus) = b.load
voltage(b::Bus) = b.voltage
population(b::LoadBus) = b.population
neighbors(b::Bus) = b.connected_to
