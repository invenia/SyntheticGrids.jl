@auto_hash_equals struct Generator
    coords::LatLon #latitude and longitude
    volt::Vector{Real} #voltages (may be multiplie values,
                                # our source gives no breakdown) # kV
    tech::AbstractString #technology type
    cap::Real #nominal capacity # MW
    pfactor::Real #nominal power factor
    minload::Real #minimum load # MW
    scap::Real #summer nominal capacity # MW
    wcap::Real #winter nominal capacity # MW
    shut2loadtime::AbstractString #time from cold shutdown to full load
    status::AbstractString #"OP" - Operational; "SB" - Standby;
                           # "OA" and "OS" - Out of service
end

function Generator{T <: Real}(
        coords::Tuple{T,T},
        volt,
        tech,
        cap,
        pfactor,
        minload,
        scap,
        wcap,
        shut2loadtime,
        status
    )
    return Generator(
        LatLon(coords[1], coords[2]),
        volt,
        tech,
        cap,
        pfactor,
        minload,
        scap,
        wcap,
        shut2loadtime,
        status
    )
end

coords(g::Generator) = g.coords
voltage(g::Generator) = g.volt
technology(g::Generator) = g.tech
capacity(g::Generator) = g.cap
power_factor(g::Generator) = g.pfactor
minimum_load(g::Generator) = g.minload
summer_capacity(g::Generator) = g.scap
winter_capacity(g::Generator) = g.wcap
status(g::Generator) = g.status
