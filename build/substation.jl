mutable struct Substation <: Bus
    id::Int # for internal use only, don't mess with this.
    coords::LatLon # latitude and longitude
    voltages::Vector{Real} # kV
    load::Real # MW
    generation::Real # MW
    population::Integer
    connected_to::Set{Substation}
    grouping::Vector{Bus}

    function Substation(id, coords, volts, load, gen, pop, con = Set(), group = [])
        return new(id, coords, volts, load, gen, pop, con, group)
    end
end

function Substation{T<:Real}(
    id,
    coords::Tuple{T,T},
    volts,
    load,
    gen,
    pop,
    con=Set(),
    group=[]
)
    return Substation(
        id,
        LatLon(coords[1], coords[2]),
        volts,
        load,
        gen,
        pop,
        con,
        group
    )
end

function Substation(coords::LatLon, volts, load, gen, pop, con = Set(), group = [])
    return Substation(-1, coords, volts, load, gen, pop, con, group)
end

function Substation{T<:Real}(
    coords::Tuple{T,T},
    volts,
    load,
    gen,
    pop,
    con=Set(),
    group=[]
)
    return Substation(
        -1,
        LatLon(coords[1], coords[2]),
        volts,
        load,
        gen,
        pop,
        con,
        group
    )
end

function show(io::IO, bus::Substation)
    println(io, "Substation(")
    println(io, "\tid=$(bus.id)")
    println(io, "\tcoords=$(bus.coords),")
    println(io, "\tvoltages=$(bus.voltages),")
    println(io, "\tload=$(bus.load), ")
    println(io, "\tgeneration=$(bus.generation), ")
    println(io, "\tpopulation=$(bus.population), ")
    if length(bus.connected_to) > 0
        println(io, "\tconnected_to=Set{Substation}(...)")
    else
        println(io, "\tconnected_to=Set{Substation}()")
    end
    println(io, "\tgrouping=$(bus.grouping)")
    print(io, ")")
end

function hash(a::Substation, h::UInt)
    h = hash(a.id, h)
    h = hash(a.grouping, h)
    h = hash(a.generation, h)
    h = hash(a.population, h)
    h = hash(a.voltages, h)
    h = hash(a.load, h)
    return hash(a.coords, h)
end

function ==(a::Substation, b::Substation)
    compared = Set{Tuple{UInt, UInt}}()

    function subs_equal(a::Substation, b::Substation)
        a.id == b.id &&
        a.coords == b.coords &&
        a.load== b.load &&
        a.voltages == b.voltages &&
        a.population == b.population &&
        length(a.connected_to) == length(b.connected_to) &&
        a.grouping == b.grouping &&
        a.generation == b.generation || return false

        pair = (object_id(a), object_id(b))
        if !(pair in compared)
            push!(compared, pair)

            all(subs_equal(a, b) for (a, b) in zip(a.connected_to, b.connected_to)) ||
            return false
        end
        return true
    end

    return subs_equal(a, b)
end

voltage(s::Substation) = s.voltages
capacity(s::Substation) = s.generation
population(s::Substation) = s.population
buses(s::Substation) = s.grouping
