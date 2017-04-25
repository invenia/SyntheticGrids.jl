const CENSUSPATH = joinpath(dirname(@__FILE__), "..", "data", "Census_data.dat")
const GENCOORDPATH = joinpath(dirname(@__FILE__), "..", "data", "Generator_coord.dat")
const GENDATAPATH = joinpath(dirname(@__FILE__), "..", "data", "Generator_data.dat")
const GENJSONPATH = joinpath(dirname(@__FILE__), "..", "data", "GenData.json")

function zipcode_builder(datapath = CENSUSPATH)
    df = CSV.read(datapath, delim='\t')
    zipcodes = Vector(size(df, 1))
    for i in 1:size(df, 1)
        zip = [
            get(df[i, 1]),
            get(df[i, 2]),
            get(df[i, 8]),
            parse(Float64, get(df[i, 9]))
        ]
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

    srand(grid.seed + 23) # Just so we don't keep returning to the same point.
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

    srand(grid.seed + 23) # Just so we don't keep returning to the same point.
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

function get_plant_data(coordpath = GENCOORDPATH)
    df = CSV.read(coordpath)
    pcodes = sizehint!(Int[], size(df, 1))
    pcoords = sizehint!([], size(df, 1))
    pvolts = sizehint!([], size(df, 1))
    for i in 1:size(df, 1)
        push!(pcodes, get(df[i, Symbol("Plant Code")]))
        push!(
            pcoords,
            (get(df[i, Symbol("Latitude")]), get(df[i, Symbol("Longitude")]))
        )
        vs = Real[]
        for s in [
            Symbol("Grid Voltage (kV)"),
            Symbol("Grid Voltage 2 (kV)"),
            Symbol("Grid Voltage 3 (kV)")
        ]
            if get(df[i, s]) != " "
                push!(vs, parse(Float64, get(df[i, s])))
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
    col_types = [
        Int,
        AbstractString,
        AbstractString,
        AbstractString,
        AbstractString,
        AbstractString,
        AbstractString,
        AbstractString,
        AbstractString,
    ]
    df = CSV.read(datapath, types=Dict(zip(keys2keep, col_types)))
    tempdata = sizehint!([], size(df, 1))
    for i in 1:size(df, 1)
        gdata = []
        for k in keys2keep[2:end]
            push!(gdata, get(df[i, Symbol(k)], " "))
        end
        for j in 1:length(pcodes)
            if pcodes[j] == get(df[i, Symbol(keys2keep[1])])
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
    gen["cap"] = parse(Float64, gen["cap"])
    gen["pfactor"] = length(gen["pfactor"])>1 ? parse(Float64, gen["pfactor"]) : 1
    gen["minload"] = length(gen["minload"])>1 ? parse(Float64, gen["minload"]) : 0
    gen["scap"] = length(gen["scap"])>1 ? parse(Float64, gen["scap"]) : gen["cap"]
    gen["wcap"] = length(gen["wcap"])>1 ? parse(Float64, gen["wcap"]) : gen["cap"]
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
    keys = ["coords", "volt", "tech", "cap", "pfactor", "minload", "scap",
            "wcap", "shut2loadtime", "status"]
    for p in plants
        push!(fplants, [Dict(k => v for (k,v) in zip(keys,g)) for g in p])
        for g in fplants[end]
            format_gen!(g)
        end
    end
    open(outfile, "w") do f
        JSON.print(f, fplants)
    end
end

function format_gen!(gen::Generator)
    gen.cap = parse(Float64, gen.cap)
    gen.pfactor = length(gen.pfactor) > 1 ? parse(Float64, gen.pfactor) : 1
    gen.minload = length(gen.minload) > 1 ? parse(Float64, gen.minload) : 0
    gen.scap = length(gen.scap) > 1 ? parse(Float64, gen.scap) : gen.cap
    gen.wcap = length(gen.wcap) > 1 ? parse(Float64, gen.wcap) : gen.cap
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
                    (g["coords"][1],g["coords"][2]),
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

"""
    place_gens_from_data_old!(grid::Grid, latlim=(-Inf,Inf),
                      longlim=(-Inf,Inf), datapath="Generator_data.dat",
                      coordpath = "Generator_coord.dat", pfactorfunc = pfacwavg)

OBSOLETE. Create generation buses based on power plant data.
Obsolete now that we have a standalone function for creating a json with the
data. Reading in the json file via place_gens_from_data!() is considerably faster.

# Arguments:
* grid: Grid instance which will be populated.
* latlim: Tuple of the form (minimum latitude, maximum latitude).
* longlim: Tuple of the form (minimum longitude, maximum longitude).
* datapath: Path to the text file containing detailed generator information.
* coordpath: Path to the text file containing power plant information.
* pfactorfunc: Function for computing bus power factor based on some rule.

REFERENCE: https://www.eia.gov/electricity/data/eia860/index.html
"""
function place_gens_from_data_old!(
    grid::Grid;
    latlim = (-Inf, Inf),
    longlim = (-Inf, Inf),
    datapath = GENDATAPATH,
    coordpath = GENCOORDPATH,
    pfactorfunc = pfacwavg,
)
    plants = get_gen_data(datapath, coordpath)
    plants = [
                p for p in plants if p[1][1][1] > latlim[1] &&
                p[1][1][1] < latlim[2] && p[1][1][2] > longlim[1] &&
                p[1][1][2] < longlim[2]
            ]
    # Creating generators
    genbs = Vector(plants)
    for p in plants
        push!(
            genbs,
            [
                Generator(
                    g[1],
                    g[2],
                    g[3],
                    g[4],
                    g[5],
                    g[6],
                    g[7],
                    g[8],
                    g[9],
                    g[10]
                ) for g in p
            ]
        )
    end
    # Formatting generators
    for p in genbs
        for g in p
            format_gen!(g)
        end
    end
    # Creating generation buses
    for p in genbs
        op = filter(g -> g.status in ["OP", "SP"], p)
        push!(buses(grid), GenBus(length(buses(grid)) + 1, p, op, pfactorfunc))
    end
end

# Here we will adopt transmission lines of type 'Line' with predefined values. We need
# their predefined values because we do not have all properties required to create a
# 'Line'. The use of 'DC Line' is precluded by the fact that, in pandapower, it is only
# capable of unidirectional flow.
const SAFE_LOAD_PERCENT = 100 # Maximum load allowed at lines and transformers.
function add_line(pgrid::PyObject, b1::LoadBus, b2::LoadBus)
    pp[:create_line](
        pgrid,
        from_bus = (b1.id - 1), # Adopting zero indexing because this is for Python.
        to_bus = (b2.id - 1), # Adopting zero indexing because this is for Python.
        length_km = haversine(b1, b2),
        std_type = "305-AL1/39-ST1A 110.0",
        max_loading_percent = SAFE_LOAD_PERCENT,
    )
end
function add_line(pgrid::PyObject, b1::LoadBus, b2::GenBus)
    pp[:create_transformer](
        pgrid,
        hv_bus = (b2.id - 1), # Adopting zero indexing because this is for Python.
        lv_bus = (b1.id - 1), # Adopting zero indexing because this is for Python.
        std_type = "160 MVA 380/110 kV",
        max_loading_percent = SAFE_LOAD_PERCENT,
    )
end
function add_line(pgrid::PyObject, b1::GenBus, b2::LoadBus)
    pp[:create_transformer](
        pgrid,
        hv_bus = (b1.id - 1), # Adopting zero indexing because this is for Python.
        lv_bus = (b2.id - 1), # Adopting zero indexing because this is for Python.
        std_type = "160 MVA 380/110 kV",
        max_loading_percent = SAFE_LOAD_PERCENT,
    )
end
function add_line(pgrid::PyObject, b1::GenBus, b2::GenBus)
    pp[:create_line](
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
    const LOAD_VOLT = 110 # kV # This is a crude approximation
    const GEN_VOLT = 380 # kV # This is a crude approximation
    voltage(bus::LoadBus) = LOAD_VOLT
    voltage(bus::GenBus) = GEN_VOLT
    # Create empty pandapower network
    pgrid = pycall(pp[:create_empty_network], PyObject)
    # Create buses
    for i in 1:length(buses(grid))
        pp[:create_bus](
            pgrid,
            voltage(buses(grid)[i]),
            index = (i - 1) # Adopting zero indexing because this is for Python.
        ) # Despite what the documentation says, 'create_bus' requires voltage.
    end
    # Create loads and generators
    for i in 1:length(buses(grid))
        if isa(buses(grid)[i], LoadBus)
            pp[:create_load](
                pgrid,
                bus = (i-1), # Adopting zero indexing because this is for Python.
                p_kw = 1000*buses(grid)[i].load, # MW to kW
                controllable = false # for OPF
            )
        else
        # We are using static generators because the other type demands more attributes
        # we don't have atm. The type of the generator is being stored in 'name' instead
        # of in 'type' because Julia seems to have problems with a parameter called 'type'.
            for gen in buses(grid)[i].gens
                pp[:create_sgen](
                    pgrid,
                    bus = (i-1), # Adopting zero indexing because this is for Python.
                    p_kw = -1000 * gen.cap, # negative by pandapower convention
                    name = gen.tech,
                    max_p_kw = -1000 * gen.minload, # min inverts with max because negative
                    min_p_kw = -1000 * gen.cap,
                    controllable = true, # for OPF
                )
            end
        end
    end
    # Create lines and transformers
    line_types = ( # All pre-defined overhead line types. Seem to be for lower voltages.
        "149-AL1/24-ST1A 10.0",
        "149-AL1/24-ST1A 110.0",
        "149-AL1/24-ST1A 20.0",
        "15-AL1/3-ST1A 0.4",
        "184-AL1/30-ST1A 110.0",
        "184-AL1/30-ST1A 20.0",
        "24-AL1/4-ST1A 0.4",
        "243-AL1/39-ST1A 110.0",
        "243-AL1/39-ST1A 20.0",
        "305-AL1/39-ST1A 110.0",
        "48-AL1/8-ST1A 0.4",
        "48-AL1/8-ST1A 10.0",
        "48-AL1/8-ST1A 20.0",
        "490-AL1/64-ST1A 220.0",
        "490-AL1/64-ST1A 380.0",
        "94-AL1/15-ST1A 0.4",
        "94-AL1/15-ST1A 10.0",
        "94-AL1/15-ST1A 20.0",
    )
    lines_110 = (
        "149-AL1/24-ST1A 110.0",
        "184-AL1/30-ST1A 110.0",
        "243-AL1/39-ST1A 110.0",
        "305-AL1/39-ST1A 110.0",
    )
    for ln in grid.trans_lines
        buses = unique(ln.connecting)
        add_line(pgrid, buses[1], buses[2])
    end
    pp[:create_ext_grid](pgrid, 0, vm_pu=1, max_p_kw=0, min_p_kw=0) # Just so pandapower does
                                                    # not throw an error when doing and OPF.
    return pgrid
end

"""
        to_pandapower(grid::Grid, filename::AbstractString)

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
    pp[:to_pickle](pgrid, path)
end

"""
    load_pp_grid(path)

Load a pandapower grid from `path`.
"""
function load_pp_grid(path)
    return pycall(pp[:from_pickle], PyObject,path)
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
