const CENSUSPATH = joinpath(dirname(@__FILE__), "..", "data", "Census_data.dat")
const GENCOORDPATH = joinpath(dirname(@__FILE__), "..", "data", "Generator_coord.dat")
const GENDATAPATH = joinpath(dirname(@__FILE__), "..", "data", "Generator_data.dat")
const GENJSONPATH = joinpath(dirname(@__FILE__), "..", "data", "GenData.json")

function zipcode_builder(datapath = CENSUSPATH)
    df = CSV.read(datapath, delim='\t')
    zipcodes = Vector(undef, size(df, 1))
    for i in 1:size(df, 1)
        zip = [df[i, 1], df[i, 2], df[i, 8], df[i, 9]]
        zipcodes[i] = zip
    end
    return zipcodes
end

function linear_load(bus::LoadBus)
    return 0.002 * bus.population # Overbye's first approxmation
end

function rand_volt(bus::Bus)
    volts = [100, 200, 350, 500]
    return volts[rand(1:length(volts))]
end

"""
    place_loads_from_zips!(
        grid::Grid;
        latlim = (-Inf,Inf),
        longlim = (-Inf,Inf),
        datapath = CENSUSPATH,
        loadfunc = linearload,
        voltfunc = randvolt
    )

Create load buses based on zipcode information.

# Arguments
* grid: The synthetic grid instance that will be altered.
* latlim: Tuple of the form (minimum latitude, maximum latitude).
* longlim: Tuple of the form (minimum longitude, maximum longitude).
* datapath: Path to the text file containing census information.
* loadfunc: Function for computing bus load based on some rule.
* voltfunc: Function for computing bus voltage based on some rule.

REFERENCE: https://www.census.gov/geo/maps-data/data/gazetteer2010.html
"""
function place_loads_from_zips!(
    grid::Grid;
    latlim = (-Inf, Inf),
    longlim = (-Inf, Inf),
    datapath = CENSUSPATH,
    loadfunc = linear_load,
    voltfunc = rand_volt,
)
    allzips = zipcode_builder(datapath)
    fzips = [
        z for z in allzips if z[3]>latlim[1] &&
        z[3]<latlim[2] &&
        z[4]>longlim[1] &&
        z[4]<longlim[2]
    ]

    seed!(grid.seed + 23) # Just so we don't keep returning to the same point.
    for zip in fzips
        dummy = LoadBus(-1, (zip[3],zip[4]), -1, -1, zip[2])
        push!(
            buses(grid),
            LoadBus(
                length(buses(grid)) + 1,
                (zip[3],zip[4]),
                loadfunc(dummy), # This is pretty convoluted for rigid values, but makes it
                voltfunc(dummy), # extremely easy to extend in the future for more
                zip[2]           # ellaborate functions.
            )
        )
    end
end

"""
    place_loads_from_zips!(
        grid::Grid,
        geolim::Function;
        datapath = CENSUSPATH,
        loadfunc = linearload,
        voltfunc = randvolt
    )

Create load buses based on zipcode information.

# Arguments
* grid: The synthetic grid instance that will be altered.
* geolim: Function that receives a latitude-longitude pair and returns `true` or `false`.
Specifies the region of interest.
* datapath: Path to the text file containing census information.
* loadfunc: Function for computing bus load based on some rule.
* voltfunc: Function for computing bus voltage based on some rule.

REFERENCE: https://www.census.gov/geo/maps-data/data/gazetteer2010.html
"""
function place_loads_from_zips!(
    grid::Grid,
    geolim::Function;
    datapath = CENSUSPATH,
    loadfunc = linear_load,
    voltfunc = rand_volt,
)
    allzips = zipcode_builder(datapath)
    fzips = [z for z in allzips if geolim((z[3], z[4]))]

    seed!(grid.seed + 23) # Just so we don't keep returning to the same point.
    for zip in fzips
        dummy = LoadBus(-1, (zip[3],zip[4]), -1, -1, zip[2])
        push!(
            buses(grid),
            LoadBus(
                length(buses(grid)) + 1,
                (zip[3],zip[4]),
                loadfunc(dummy), # This is pretty convoluted for rigid values, but makes it
                voltfunc(dummy), # extremely easy to extend in the future for more
                zip[2]           # ellaborate functions.
            )
        )
    end
end

function get_plant_data(coordpath=GENCOORDPATH)
    column_names = ["Grid Voltage (kV)", "Grid Voltage 2 (kV)", "Grid Voltage 3 (kV)"]
    column_types = Dict(Pair.(column_names, Float64))
    df = CSV.read(coordpath, types=column_types)
    pcodes = sizehint!(Int[], size(df, 1))
    pcoords = sizehint!([], size(df, 1))
    pvolts = sizehint!([], size(df, 1))
    for i in 1:size(df, 1)
        push!(pcodes, df[i, Symbol("Plant Code")])
        push!(
            pcoords,
            (df[i, Symbol("Latitude")], df[i, Symbol("Longitude")])
        )
        vs = Real[]

        for s in Symbol.(column_names)
            if !ismissing(df[i, s])
                push!(vs, df[i, s])
            end
        end
        push!(pvolts, vs)
    end
    return pcodes, pcoords, pvolts
end

function get_gen_data(
    datapath = GENDATAPATH,
    coordpath = GENCOORDPATH,
)
    # Generator data is separated in two spreadsheets.
    # We have to parse both and then combine.

    #Grabbing plant information
    pcodes, pcoords, pvolts = get_plant_data(coordpath)
    #Grabbing generator information
    keys2keep = [
        "Plant Code",
        "Technology",
        "Nameplate Capacity (MW)",
        "Nameplate Power Factor",
        "Minimum Load (MW)",
        "Summer Capacity (MW)",
        "Winter Capacity (MW)",
        "Time from Cold Shutdown to Full Load",
        "Status"
    ]
    df = CSV.read(datapath)
    tempdata = sizehint!([], size(df, 1))
    for i in 1:size(df, 1)
        gdata = []
        for k in keys2keep[2:end]
            push!(gdata, df[i, Symbol(k)])
        end
        for j in 1:length(pcodes)
            if pcodes[j] == df[i, Symbol(keys2keep[1])]
                gdata = vcat([pcoords[j], pvolts[j]], gdata)
                break
            end
        end
        push!(tempdata, gdata)
    end
    plants = []
    for c in pcoords
        p = [g for g in tempdata if g[1] == c]
        length(p) > 0 ? push!(plants, p) : nothing #lots of plants have no
        # corresponding generators. This means our survey data is not so good.
    end
    return plants
end

function format_gen!(gen::Dict)
    gen["cap"] = gen["cap"]
    gen["pfactor"] = !ismissing(gen["pfactor"]) ? gen["pfactor"] : 1
    gen["minload"] = !ismissing(gen["minload"]) ? gen["minload"] : 0
    gen["scap"] = !ismissing(gen["scap"]) ? gen["scap"] : gen["cap"]
    gen["wcap"] = !ismissing(gen["wcap"]) ? gen["wcap"] : gen["cap"]
    gen["shut2loadtime"] = !ismissing(gen["shut2loadtime"]) ? gen["shut2loadtime"] : " "
    gen["status"] = !ismissing(gen["status"]) ? gen["status"] : " "
end

"""
    prepare_gen_data(indata = "Generator_data.dat",
                    incoord = "Generator_coord.dat",
                    outfile = "GenData.json")

Prepare json (outfile) with generator data gathered from the survey spreadsheets
(indata and incoord). Output is organized as an array of plants. Each plant is
an array of generators. Each generator is a dictionary.
"""
function prepare_gen_data(
    indata = GENDATAPATH,
    incoord = GENCOORDPATH,
    outfile = GENJSONPATH,
)
    plants = get_gen_data(indata, incoord)
    fplants = sizehint!([], length(plants))
    _keys = ["coords", "volt", "tech", "cap", "pfactor", "minload", "scap",
            "wcap", "shut2loadtime", "status"]
    for p in plants
        plant_dict = [Dict(k => v for (k,v) in zip(_keys,g)) for g in p]

        for gen in plant_dict
            format_gen!(gen)
        end

        push!(fplants, plant_dict)
    end

    open(outfile, "w") do f
        JSON.print(f, fplants)
    end
end

function format_gen!(gen::Generator)
    gen.cap = gen.cap
    gen.pfactor = !ismissing(gen.pfactor) ? gen.pfactor : 1
    gen.minload = !ismissing(gen.minload) ? gen.minload : 0
    gen.scap = !ismissing(gen.scap) ? gen.scap : gen.cap
    gen.wcap = !ismissing(gen.wcap) ? gen.wcap : gen.cap
    gen.shut2loadtime = !ismissing(gen.shut2loadtime) ? gen.shut2loadtime : " "
    gen.status = !ismissing(gen.status) ? gen.status : " "
end

function pfacwavg(plant::AbstractArray{Generator})
    on = filter(g -> g.status in ["OP", "SP"], plant)
    if length(on) > 0
        caps = [g.cap for g in on]
        pfacts = [g.pfactor for g in on]
        pfactor = sum(pfacts .* caps)
        return pfactor / sum(caps)
    else
        return 0
    end
end

"""
    place_gens_from_data!(
        grid::Grid,
        plants::Array,
        pfactorfunc::Function = pfacwavg
    )

Create generation buses based on power plant data.

# Arguments:
* grid: Grid instance which will be populated.
* plants: Power plant information already parsed.
* pfactorfunc: Function for computing bus power factor based on some rule.
"""
function place_gens_from_data!(
    grid::Grid,
    plants::Array,
    pfactorfunc::Function = pfacwavg,
)
    # Creating generators
    genbs = sizehint!(Vector{Vector{Generator}}(), length(plants))
    for p in plants
        push!(
            genbs,
            [
                Generator(
                    (g["coords"][1], g["coords"][2]),
                    g["volt"],
                    g["tech"],
                    g["cap"],
                    g["pfactor"],
                    g["minload"],
                    g["scap"],
                    g["wcap"],
                    g["shut2loadtime"],
                    g["status"]
                ) for g in p
            ]
        )
    end
    # Creating generation buses
    for p in genbs
        op = filter(g -> g.status in ["OP", "SP"], p)
        push!(buses(grid), GenBus(length(buses(grid)) + 1, p, op, pfactorfunc))
    end
end

"""
    place_gens_from_data!(
        grid::Grid;
        latlim::Tuple = (-Inf, Inf),
        longlim::Tuple = (-Inf, Inf),
        jsonpath = GENJSONPATH,
        pfactorfunc::Function = pfacwavg,
	)

Create generation buses based on power plant data.

# Arguments:
* grid: Grid instance which will be populated.
* latlim: Tuple of the form (minimum latitude, maximum latitude).
* longlim: Tuple of the form (minimum longitude, maximum longitude).
* jsonpath: Path to the JSON file with the generator data.
* pfactorfunc: Function for computing bus power factor based on some rule.

REFERENCE: https://www.eia.gov/electricity/data/eia860/index.html
"""
function place_gens_from_data!(
    grid::Grid;
    latlim::Tuple = (-Inf, Inf),
    longlim::Tuple = (-Inf, Inf),
    jsonpath = GENJSONPATH,
    pfactorfunc::Function = pfacwavg,
)
    #Loading data from json
    plants = JSON.parsefile(jsonpath)
    #Filtering per latitude and longitude
    plants = [
        p for p in plants if p[1]["coords"][1] > latlim[1] &&
        p[1]["coords"][1] < latlim[2] &&
        p[1]["coords"][2] > longlim[1] &&
        p[1]["coords"][2] < longlim[2]
    ]
    place_gens_from_data!(grid, plants, pfactorfunc)
end

"""
    place_gens_from_data!(
        grid::Grid,
        geolim::Function;
        jsonpath = GENJSONPATH,
        pfactorfunc::Function = pfacwavg,
	)

Create generation buses based on power plant data.

# Arguments:
* grid: Grid instance which will be populated.
* geolim: Function that receives a latitude-longitude pair and returns `true` or `false`.
Specifies the region of interest.
* jsonpath: Path to the JSON file with the generator data.
* pfactorfunc: Function for computing bus power factor based on some rule.

REFERENCE: https://www.eia.gov/electricity/data/eia860/index.html
"""
function place_gens_from_data!(
    grid::Grid,
    geolim::Function;
    jsonpath = GENJSONPATH,
    pfactorfunc::Function = pfacwavg,
)
    #Loading data from json
    plants = JSON.parsefile(jsonpath)
    #Filtering per latitude and longitude
    plants = [p for p in plants if geolim((p[1]["coords"][1], p[1]["coords"][2]))]
    place_gens_from_data!(grid, plants, pfactorfunc)
end

# Here we will adopt transmission lines of type 'Line' with predefined values. We need
# their predefined values because we do not have all properties required to create a
# 'Line'. The use of 'DC Line' is precluded by the fact that, in pandapower, it is only
# capable of unidirectional flow.
const SAFE_LOAD_PERCENT = 100 # Maximum load allowed at lines and transformers.
function add_line(pgrid::PyObject, b1::LoadBus, b2::LoadBus)
    pp.create_line(
        pgrid,
        from_bus = (b1.id - 1), # Adopting zero indexing because this is for Python.
        to_bus = (b2.id - 1), # Adopting zero indexing because this is for Python.
        length_km = haversine(b1, b2),
        std_type = "305-AL1/39-ST1A 110.0",
        max_loading_percent = SAFE_LOAD_PERCENT,
    )
end
function add_line(pgrid::PyObject, b1::LoadBus, b2::GenBus)
    pp.create_transformer(
        pgrid,
        hv_bus = (b2.id - 1), # Adopting zero indexing because this is for Python.
        lv_bus = (b1.id - 1), # Adopting zero indexing because this is for Python.
        std_type = "160 MVA 380/110 kV",
        max_loading_percent = SAFE_LOAD_PERCENT,
    )
end
function add_line(pgrid::PyObject, b1::GenBus, b2::LoadBus)
    pp.create_transformer(
        pgrid,
        hv_bus = (b1.id - 1), # Adopting zero indexing because this is for Python.
        lv_bus = (b2.id - 1), # Adopting zero indexing because this is for Python.
        std_type = "160 MVA 380/110 kV",
        max_loading_percent = SAFE_LOAD_PERCENT,
    )
end
function add_line(pgrid::PyObject, b1::GenBus, b2::GenBus)
    pp.create_line(
        pgrid,
        from_bus = (b1.id - 1), # Adopting zero indexing because this is for Python.
        to_bus = (b2.id - 1), # Adopting zero indexing because this is for Python.
        length_km = haversine(b1, b2),
        std_type = "490-AL1/64-ST1A 380.0",
        max_loading_percent = SAFE_LOAD_PERCENT,
    )
end

"""
    to_pandapower(grid::Grid)

Export a 'grid' to pandapower format.

Important:
Currently, grid voltages and line properties are ignored. This function places all loads at
100kV and all generation at 380kV. Transmission lines use standard values for 110kV and
380kV lines. Transformers are automatically placed between buses that operate at different
voltage values.
"""
function to_pandapower(grid::Grid,)
    LOAD_VOLT = 0.11 # mV, this is a crude approximation
    GEN_VOLT = 0.38 # mV, this is a crude approximation
    voltage(bus::LoadBus) = LOAD_VOLT
    voltage(bus::GenBus) = GEN_VOLT

    # Create empty pandapower network
    pgrid = pycall(pp.create_empty_network, PyObject)

    # Create buses
    for i in 1:length(buses(grid))
        pp.create_bus(
            pgrid,
            voltage(buses(grid)[i]),
            index = (i - 1), # Adopting zero indexing because this is for Python.
            geodata = (buses(grid)[i].coords.lat, buses(grid)[i].coords.lon)
        ) # Despite what the documentation says, 'create_bus' requires voltage.
    end

    # Create loads and generators
    for i in 1:length(buses(grid))
        if isa(buses(grid)[i], LoadBus)
            pp.create_load(
                pgrid,
                bus = (i-1), # Adopting zero indexing because this is for Python.
                p_mw = buses(grid)[i].load,
                controllable = false # for OPF
            )
        else
        # The type of the generator is being stored in 'name' instead
        # of in 'type' because Julia seems to have problems with a parameter called 'type'.
            for gen in buses(grid)[i].gens
                pp.create_gen(
                    pgrid,
                    bus = (i-1), # Adopting zero indexing because this is for Python.
                    p_mw = -gen.cap,  # Negative for generation
                    name = gen.tech,
                    max_p_mw = -gen.minload,  # Min inverts with max because negative
                    min_p_mw = -gen.cap,  # Negative for genereation
                    controllable = true, # for OPF
                )
            end
        end
    end

    for ln in grid.trans_lines
        buses = unique(ln.connecting)
        add_line(pgrid, buses[1], buses[2])
    end
    pp.create_ext_grid(pgrid, 0, vm_pu=1, max_p_kw=0, min_p_kw=0) # Just so pandapower does
                                                    # not throw an error when doing and OPF.
    return pgrid
end

"""
    to_pandapower(grid::Grid, filename::AbstractString)D

Export a 'grid' to pandapower format and saves it as 'filename'. Pandapower requires
'filename' to be a '.p' file.

Important:
Currently, grid voltages and line properties are ignored. This function places all loads at
100kV and all generation at 380kV. Transmission lines use standard values for 110kV and
380kV lines. Transformers are automatically placed between buses that operate at different
voltage values.
"""
function to_pandapower(grid::Grid, filename::AbstractString)
    pgrid = to_pandapower(grid)
    save_pp_grid(pgrid, filename)
    return pgrid
end

"""
    save_pp_grid(pgrid::PyObject, filename::AbstractString)

Save a pandapower grid as `filename`. Pandapower requires `filename` to be a `.p` file.
"""
function save_pp_grid(pgrid::PyObject, path::AbstractString)
    pp.to_pickle(pgrid, path)
end

"""
    load_pp_grid(path)

Load a pandapower grid from `path`.
"""
function load_pp_grid(path)
    return pycall(pp.from_pickle, PyObject,path)
end

# WARNING: exportJLD and loadJLD functions are not working properly yet.
# They should be disregarded.
function save(grid::Grid, outname)
    # Breaking object links in order to avoid StackOverflowError
    dummytuple = (LoadBus(-1, (-1,-1), -1, -1, -1), LoadBus(-1, (-1,-1), -1, -1, -1))
    for ln in grid.trans_lines
        ln.connecting = dummytuple
    end
    for bus in buses(grid)
        bus.connected_to = []
        bus.connections = []
    end
    for sub in substations(grid)
        sub.connected_to = []
        # sub.grouping = []
    end
    #Storing object
    jldopen(outname, "w") do f
        # addrequire(f, SyntheticGrids)
        write(f, "grid", grid)
    end
end

function load_grid(inname)
    # Load grid
    grid = load(inname)
    # Reconnect grid
    connect_buses!(grid)
    count = 1
    for j in 1:(length(buses(grid)) - 1), i in (j + 1):length(buses(grid))
        if grid.bus_conn[i, j] == 1
            grid.trans_lines[count].connecting = (buses(grid)[i], buses(grid)[j])
            push!(buses(grid)[i].connections, grid.trans_lines[count])
            push!(buses(grid)[j].connections, grid.trans_lines[count])
            count += 1
        end
    end
    if length(substations(grid)) > 0
        connect_subs!(grid)
    end
    return grid
end

"""
    save(grid::Grid, outfile::AbstractString)

Save a `grid` in a format that can be recovered later. The format adopted is a JSON file
with a particular structure.
"""
function save(grid::Grid, outfile::AbstractString)
    out_dict = Dict()
    out_dict["seed"] = grid.seed
    bus_vec = []
    for b in buses(grid)
        if isa(b, LoadBus)
            push!(bus_vec, Dict(
                "type" => "LoadBus",
                "coords" => b.coords,
                "load" => b.load,
                "voltage" => b.voltage,
                "population" => b.population
            ))
        elseif isa(b, GenBus)
            push!(bus_vec, Dict(
                "type" => "GenBus",
                "coords" => b.coords,
                "generation" => b.generation,
                "voltage" => b.voltage,
                "tech_type" => b.tech_type,
                "pfactor" => b.pfactor,
                "summgen" => b.summgen,
                "wintgen" => b.wintgen,
                "gens" => [Dict(
                    "coords" => g.coords,
                    "volt" => g.volt,
                    "tech" => g.tech,
                    "cap" => g.cap,
                    "pfactor" => g.pfactor,
                    "minload" => g.minload,
                    "scap" => g.scap,
                    "wcap" => g.wcap,
                    "shut2loadtime" => g.shut2loadtime,
                    "status" => g.status
                ) for g in b.gens]
            ))
        else
            warn("Bus $b is not of a known type. Its type is $(typeof(b))")
        end
    end
    out_dict["buses"] = bus_vec
    out_dict["trans_lines"] = [Dict(
        "impedance" => t.impedance,
        "capacity" => t.capacity,
        "connecting" => [t.connecting[1].id, t.connecting[2].id]
    ) for t in trans_lines(grid)]
    out_dict["substations"] = [Dict(
        "coords" => s.coords,
        "voltages" => s.voltages,
        "load" => s.load,
        "generation" => s.generation,
        "population" => s.population,
        "grouping" => [b.id for b in s.grouping]
    ) for s in substations(grid)]
    out_dict["bus_conn"] = adjacency(grid)

    open(outfile, "w") do f
        JSON.print(f, out_dict)
    end
end

"""
    load_grid(filepath::AbstractString)

Load a synthetic grid from a file previously saved via `save(grid::Grid)`. Due to
conversions during the dumping and the reading of the file, floats may have some noise
added. That should be within machine precision, but may affect operations such as `==`.
"""
function load_grid(filepath::AbstractString)
    in_grid = JSON.parsefile(filepath)
    out_grid = Grid(in_grid["seed"])
    out_grid.bus_conn = hcat(in_grid["bus_conn"]...)
    for id in 1:length(in_grid["buses"])
        if in_grid["buses"][id]["type"] == "LoadBus"
            push!(
                out_grid.buses,
                LoadBus(
                    id,
                    LatLon(
                        in_grid["buses"][id]["coords"]["lat"],
                        in_grid["buses"][id]["coords"]["lon"]
                    ),
                    in_grid["buses"][id]["load"],
                    in_grid["buses"][id]["voltage"],
                    in_grid["buses"][id]["population"],
                )
            )
        elseif in_grid["buses"][id]["type"] == "GenBus"
            push!(
                out_grid.buses,
                GenBus(
                    id,
                    LatLon(
                        in_grid["buses"][id]["coords"]["lat"],
                        in_grid["buses"][id]["coords"]["lon"]
                    ),
                    in_grid["buses"][id]["generation"],
                    in_grid["buses"][id]["voltage"],
                    in_grid["buses"][id]["tech_type"],
                    Set(),
                    Set(),
                    in_grid["buses"][id]["pfactor"],
                    in_grid["buses"][id]["summgen"],
                    in_grid["buses"][id]["wintgen"],
                    [Generator(
                        LatLon(g["coords"]["lat"], g["coords"]["lon"]),
                        g["volt"],
                        g["tech"],
                        g["cap"],
                        g["pfactor"],
                        g["minload"],
                        g["scap"],
                        g["wcap"],
                        g["shut2loadtime"],
                        g["status"]
                    ) for g in in_grid["buses"][id]["gens"]]
                )
            )
        end
    end
    connect_buses!(out_grid)
    for id in 1:length(in_grid["substations"])
        push!(
            out_grid.substations,
            Substation(
                id,
                LatLon(
                    in_grid["substations"][id]["coords"]["lat"],
                    in_grid["substations"][id]["coords"]["lon"]
                ),
                in_grid["substations"][id]["voltages"],
                in_grid["substations"][id]["load"],
                in_grid["substations"][id]["generation"],
                in_grid["substations"][id]["population"],
                Set(),
                [out_grid.buses[i] for i in in_grid["substations"][id]["grouping"]]
            )
        )
    end
    connect_subs!(out_grid)
    for t in in_grid["trans_lines"]
        push!(out_grid.trans_lines, TransLine(
            (out_grid.buses[t["connecting"][1]], out_grid.buses[t["connecting"][2]]),
            t["impedance"],
            t["capacity"]
        ))
        push!(out_grid.buses[t["connecting"][1]].connections, out_grid.trans_lines[end])
        push!(out_grid.buses[t["connecting"][2]].connections, out_grid.trans_lines[end])
    end
    return out_grid
end
