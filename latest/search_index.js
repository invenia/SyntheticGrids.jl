var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#SyntheticGrids-v0.5-1",
    "page": "Home",
    "title": "SyntheticGrids v0.5",
    "category": "section",
    "text": "June - 2017Power grid research requires testing in realistic, large-scale,  electric  networks.   However,  in  light  of  security threats,  most  information  on  the  actual  power  grids  is considered  sensitive  and  therefore  not  available  to  the general  public.   So  far,  most  power  transmission  studies have been carried using a few publicly available test grids.  Still,  these test grids are too small to capture the  complexity  of  real  grids.   With  this  in  mind,  there has recently been a strong concentrated effort in developing methodologies for building realistic synthetic grids, based only on publicly available information.  These synthetic grids are supposed to be based on some real example  and  to  present  analogous  properties  —  such  as geographic  load/generation  distribution,  total  load  and generator types — while not actually presenting potentially sensitive information about the real grid.This module provides an open source suite for generating synthetic grids based on real data openly available to the public. Power grids constructed via the SyntheticGrids module can be easily exported to pandapower for running optimum power flow calculations. Currently, information is limited to the USA region, but the framework can be readily applied to any other region, provided there are data sources available. We leverage the works published by Overbye's group and Soltan and Zussman, providing a direct implementation of their methods. For more details on the approaches adopted, please see Model.REFERENCES:Birchfield, Adam B., et al. \"Grid Structural Characteristics as Validation Criteria for Synthetic Networks.\" IEEE Transactions on Power Systems (2016).\nGegner, Kathleen M., et al. \"A methodology for the creation of geographically realistic synthetic power flow models.\" Power and Energy Conference at Illinois (PECI), 2016 IEEE. IEEE, 2016.\nBirchfield, Adam B., et al. \"Statistical considerations in the creation of realistic synthetic power grids for geomagnetic disturbance studies.\" IEEE Transactions on Power Systems 32.2 (2017): 1502-1510.\nSoltan, Saleh, and Gil Zussman. \"Generation of synthetic spatially embedded power grid networks.\" Power and Energy Society General Meeting (PESGM), 2016. IEEE, 2016."
},

{
    "location": "index.html#Detailed-documentation-1",
    "page": "Home",
    "title": "Detailed documentation",
    "category": "section",
    "text": "Types\nFunctions\nModel"
},

{
    "location": "index.html#Current-functionalities-1",
    "page": "Home",
    "title": "Current functionalities",
    "category": "section",
    "text": "Implements basic types.\nBuilds networks from real-world data.\nBuilds connection topology from nodes.\nBuilds transmission lines.\nClusters nodes into substations.\nProvides basic checks for the graph structure.\nInterfaces with pandapower for exporting networks."
},

{
    "location": "index.html#Basic-Usage-1",
    "page": "Home",
    "title": "Basic Usage",
    "category": "section",
    "text": "grid = SynGrid()Create new (empty) gridplace_loads_from_zips!(grid, latitude limits, longitude limits)Build load buses from zipcodesplace_gens_from_data!(grid, latitude limits, longitude limits)Build generation buses from data (buses may also be placed manually)connect_grid!(grid)Generate node connectionscluster_buses!(grid, nloads, nboth, ngens)Build substations by clustering nodes (optional)create_lines!(grid)Build transmission lines from connection topology"
},

{
    "location": "index.html#Exporting-network-to-pandapower-1",
    "page": "Home",
    "title": "Exporting network to pandapower",
    "category": "section",
    "text": "This module makes use of PyCall in order to interface with pandapower. Since the reference charts in pandapower do not contain transmission lines and transformer parameters for several voltages, currently, this module places all loads at 110kV and all generation at 380kV in order to have the proper parameters. Transformers are automatically placed every time two connected buses operate at different voltages. In order to obtain a pandapower object from a SynGrid instance, simply run the command:to_pandapower(grid)This returns a PyObject. A path can be optionally passed as well in order to obtain a file saved in the native pandapower format. These files can later be imported via the loadppgrid command.Note: the exported grid is zero-indexed, as it is being passed to Python routines."
},

{
    "location": "index.html#Current-issues-1",
    "page": "Home",
    "title": "Current issues",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Limitations:-1",
    "page": "Home",
    "title": "Limitations:",
    "category": "section",
    "text": "No type has been implemented for transformers.\nConnections between buses ignore differences in voltages.\nTransmission line properties are still oversimplified.\nThere are missing quantities for AC OPF runs (DC should be fine).\nThese limitations only pertain the SynGrids objects, not the grids exported to pandapower."
},

{
    "location": "index.html#Data:-1",
    "page": "Home",
    "title": "Data:",
    "category": "section",
    "text": "Generator data has inconsistencies.\n2015 survey data has 1337 power plants without any reported generator.\n2014 survey data has 1403 power plants without any reported generator.\nThere are unsited plants and plants with clearly wrong coordinates (the obvious ones were manually corrected)."
},

{
    "location": "Types.html#",
    "page": "Types",
    "title": "Types",
    "category": "page",
    "text": ""
},

{
    "location": "Types.html#Custom-Types-1",
    "page": "Types",
    "title": "Custom Types",
    "category": "section",
    "text": "The Synthetic Grids package implements several custom types."
},

{
    "location": "Types.html#Bus-1",
    "page": "Types",
    "title": "Bus",
    "category": "section",
    "text": ""
},

{
    "location": "Types.html#description:-1",
    "page": "Types",
    "title": "description:",
    "category": "section",
    "text": "Bus is an abstract type that acts as parent to the LoadBus, GenBus and Substation types."
},

{
    "location": "Types.html#LoadBus-1",
    "page": "Types",
    "title": "LoadBus",
    "category": "section",
    "text": ""
},

{
    "location": "Types.html#description:-2",
    "page": "Types",
    "title": "description:",
    "category": "section",
    "text": "Objects of this type represent the electrical loads in a network."
},

{
    "location": "Types.html#attributes:-1",
    "page": "Types",
    "title": "attributes:",
    "category": "section",
    "text": "coords: Geodesy.LatLon. Geocoordinates of the bus.\nload: Real. Load value at bus. Positive MW by convention.\nvoltage: Real. Bus operating voltage. Positive kV by convention.\npopulation: Integer. Bus population. Used for estimating load.\nconnected_to: Set{Bus}. Buses to which it is connected.\nconnections: Set{TransLine}. Transmission lines connected to bus."
},

{
    "location": "Types.html#GenBus-1",
    "page": "Types",
    "title": "GenBus",
    "category": "section",
    "text": ""
},

{
    "location": "Types.html#description:-3",
    "page": "Types",
    "title": "description:",
    "category": "section",
    "text": "GenBus is the type that represents power plants."
},

{
    "location": "Types.html#attributes:-2",
    "page": "Types",
    "title": "attributes:",
    "category": "section",
    "text": "coords: Geodesy.LatLon. Geocoordinates of the bus.\ngeneration: Real. Nameplate capacity value at bus. Positive MW by convention.\nvoltage: Vector{Real}. Bus operating voltages. Multiple values allowed since power plants may have generators that operate at different voltage levels. Positive kV by convention.\ntech_type: Vector{AbtractString}. Technology types of the generators at the power plant.\nconnected_to: Set{Bus}. Buses to which it is connected.\nconnections: Set{TransLine}. Transmission lines connected to bus.\npfactor: Real. Power factor of the power plant.\nsummgen: Real. Summer generation capacity. Positive MW by convention.\nwintgen: Real. Winter generation capacity. Positive MW by convention.\ngens: Vector{Generators}. Generators that operate at the power plant."
},

{
    "location": "Types.html#Substation-1",
    "page": "Types",
    "title": "Substation",
    "category": "section",
    "text": ""
},

{
    "location": "Types.html#description:-4",
    "page": "Types",
    "title": "description:",
    "category": "section",
    "text": "A substation groups several other buses into an aggregate structure."
},

{
    "location": "Types.html#attributes:-3",
    "page": "Types",
    "title": "attributes:",
    "category": "section",
    "text": "coords: Geodesy.LatLon. Geocoordinates of the bus.\ngeneration: Real. Nameplate capacity value at substation. Positive MW by convention.\nvoltage: Vector{Real}. Substation operating voltages. Multiple values allowed since different buses may operate at different voltage levels. Positive kV by convention.\nload: Real. Load value at substation. Positive MW by convention.\npopulation: Integer. Substation population. Used for estimating load.\nconnected_to: Set{Substation}. Substations to which it is connected.\ngrouping: Vector{Bus}. Buses grouped by substation."
},

{
    "location": "Types.html#Generator-1",
    "page": "Types",
    "title": "Generator",
    "category": "section",
    "text": ""
},

{
    "location": "Types.html#description:-5",
    "page": "Types",
    "title": "description:",
    "category": "section",
    "text": "This type models each specific generator that is part of a power plant (GenBus)."
},

{
    "location": "Types.html#attributes:-4",
    "page": "Types",
    "title": "attributes:",
    "category": "section",
    "text": "coords: Geodesy.LatLon. Geocoordinates of the generator.\nvolt:  Vector{Real}. Generator operating voltages. Multiple values allowed since the generator survey provides voltage levels grouped by power plant (see Model). Positive kV by convention.\ntech: AbstractString. Technology type.\ncap: Real. Nameplate capacity.\npfactor: Real. Power Factor.\nminload: Real. Minimum load value. Positive MW by convention.\nscap: Real. Summer generation capacity.\nwcap: Real. Winter generation capacity.\nshut2loadtime: AbstractString. Time for going from cold shutdown to full load. Currently grabbed from the generator survey without parsing (see Model).\nstatus: AbstractString. \"OP\" - operational; \"SB\" - standby; \"OA\" and \"OS\" - out of service."
},

{
    "location": "Types.html#TransLine-1",
    "page": "Types",
    "title": "TransLine",
    "category": "section",
    "text": ""
},

{
    "location": "Types.html#description:-6",
    "page": "Types",
    "title": "description:",
    "category": "section",
    "text": "This type models transmission lines that connect different buses."
},

{
    "location": "Types.html#attributes:-5",
    "page": "Types",
    "title": "attributes:",
    "category": "section",
    "text": "connecting: Tuple{Bus,Bus}. Buses that are being connected by the transmission line.\nimpedance: Real. Impedance of the line in Ohms.\ncapacity: Real. Current carrying capacity of the line in Amperes."
},

{
    "location": "Types.html#Grid-1",
    "page": "Types",
    "title": "Grid",
    "category": "section",
    "text": ""
},

{
    "location": "Types.html#description:-7",
    "page": "Types",
    "title": "description:",
    "category": "section",
    "text": "Represents an entire power grid."
},

{
    "location": "Types.html#attributes:-6",
    "page": "Types",
    "title": "attributes:",
    "category": "section",
    "text": "seed: Integer. Seed number to allow reproducibility of all stochastic procedures.\nbuses: Vector{Buses}. Vector grouping all LoadBus and GenBus objects in a network.\ntrans_lines: Vector{TransLine}. Vector grouping all transmission lines in a network.\nsubstations: Vector{Substations}. Vector grouping all substations in a network.\nbus_conn: AbstractMatrix{Bool}. Boolean matrix indicating which buses are connected to which.\nsub_conn: AbstractMatrix{Int}. Matrix indicating the number of connections between each sbustation."
},

{
    "location": "Model.html#",
    "page": "Model",
    "title": "Model",
    "category": "page",
    "text": ""
},

{
    "location": "Model.html#Model-1",
    "page": "Model",
    "title": "Model",
    "category": "section",
    "text": "The idea behind the SyntheticGrids module is to provide a standalone and easily expandable framework for generating synthetic grids. The goal is to be able to use different levels of approximation when building grids. Interfacing with other packages should be possible through export functions.A synthetic grid is constituted basically by buses connected via transmission lines. There are several relevant points when automatically building a network with this package. We discuss these on a per type basis below."
},

{
    "location": "Model.html#Bus-1",
    "page": "Model",
    "title": "Bus",
    "category": "section",
    "text": "Bus is an abstract type that acts as parent to the LoadBus, GenBus and Substation types. Pragmatically, the one common feature to all of these is their geographical location being specified by latitude and longitude values. This is a better choice than adopting cartesian coordinates, as that would require some kind of projection scheme, such as the Mercator or the UTM one. The issue with adopting a projection scheme is that none of the available ones are adequate for being applied over very large areas, such as the territory span of the USA. Therefore, coordinates are kept in a latitude-longitude format and, whenever we need to compute real distances, we adopt the haversine formula [1] with a radius of 6369.783km, which corresponds to Earth's radius at latitude 38.8 (the center of the region we model).Bus siting and sizing is done according to the approach presented in reference [2]."
},

{
    "location": "Model.html#LoadBus-1",
    "page": "Model",
    "title": "LoadBus",
    "category": "section",
    "text": "The main properties of a LoadBus are geographical location, load and operating voltage values. Since none of this information is public, we have to adopt an approach based on some kind of public data. Overbye et al. [2] have devised the approach which we adopt in the module.Geographical location: load buses are created so as to coincide with each zip code. Official census data [3] provides us with latitude and longitude and population of each zip code in the USA.\nWe adopt the census data from 2010. There are more recent censi available, but these do not offer population data (whose importance will be explained in a moment).\nLoad: we have no public information on load for specific zip codes. Overbye et al. [2] have shown that load is approximately linear with the population size. Therefore, they propose adopting a scaling factor for computing load from population size.\nWe adopt their approximation of 2kW load per capita.\nThe code is written in a way that a function is a passed into the zip code builder for computing the loads. Thus, it is extremely straight-forward to adopt any other form of computing loads.\nOperating voltages: there is no public information on operating voltages of load buses. At this point we can either adopt some constant value or randomly draw from a set of allowed values.\nFor a network with several operating voltages, it is noted that usually loads are connected to the lower voltages."
},

{
    "location": "Model.html#GenBus-1",
    "page": "Model",
    "title": "GenBus",
    "category": "section",
    "text": "GenBus is the type that represents power plants. The reason for creating a bus for each power plant, instead of a bus for each Generator, is the fact that the data we use for placing generation buses is organized in this way. There is an official public yearly survey on all generators that operate in the USA [4] that provides us with most of the important information on power plants.Geographical location: all power plants operating in the USA are listed with their address, latitude and longitude, plant code and other information. A separate file provides information on each of the generators that operate inside a power plant.\nThis information is not completely clean. There are a few plants marked as unsited (which provide no information) and a few with wrong latitude-longitude values (the obvious ones were manually corrected).\nThe most concerning issue is that several power plants have zero operating generators reported. This may be due to incomplete survey data or due to those not having any generator that is not already retired or still just planned. These plants are ignored when creating grids.\nThe data on planned and retired generators is completely ignored when we read in the information.\nSome generators also have a future data as a retirement prediction. This is not taken into account.\nGeneration Capacity: this is obtained by directly summing the nominal capacities of all operational and standby generators in a given plant.\nVoltages: up to three different voltage values are attributed to each power plant. These are likely the voltages of different generators, but are provided on a plant-basis, not on a generator-basis, so we are not able to break them down.\nTechnology type: each plant may have multiple generator types (e.g. nuclear, wind, solar etc.), according to data discriminated in the survey.\nPower factor: the power factor of a power plant is computed as a nominal capacity-weighted average of the power factors of each generator. The code accepts any function for computing the power factor of power plant, so it is trivial to adopt different approximations if necessary.\nSummer and Winter generation capacities: these are taken as the sum of the corresponding quantities of each generator."
},

{
    "location": "Model.html#Generators-1",
    "page": "Model",
    "title": "Generators",
    "category": "section",
    "text": "All data on generators is directly extracted from the power plant survey [4], which provides us with:Geographical location: every generator is directly associated with a power plant. Thus, we adopt the power plant's latitude and longitude as the generator's location.\nOperating voltages: although a given generator probably has a single operating voltage, this information is given to us on a plant-basis, not on a generator-basis. Therefore, each generator will be attributed the whole set of voltages from the parent plant.\nTechnology type: directly obtained from survey data. Generators with different technology may be located at the same plant.\nGeneration capacity: directly obtained from survey data.\nPower factor: directly obtained from survey data. In case it is not given, we assume it equals 1.0.\nMinimum load: directly obtained from survey data. In case it is not given, we assume it equals 0.0.\nSummer and Winter generation capacity: directly obtained from survey data. In case it is not given, we assume it equals the nominal generation capacity.\nTime from cold shutdown to full load: currently we grab it without even parsing, but it is kept due to probably being useful in the future. Several generators do not offer this information.\nStatus: Informs us if the generator was operational, on standby or out of service for most of the survey year."
},

{
    "location": "Model.html#Substations-1",
    "page": "Model",
    "title": "Substations",
    "category": "section",
    "text": "A substation groups several other buses into an aggregate structure. We are directly implementing the clustering algorithm presented by Birchfield et al. in reference [5]. Substations are divided in three kinds: pure load, pure generation and both load and generation. Transmission substations are ignored since we are not using connection buses."
},

{
    "location": "Model.html#TransLines-1",
    "page": "Model",
    "title": "TransLines",
    "category": "section",
    "text": "Transmission lines are placed at every connection between two buses as determined by our connection topology builder (explained ahead). Currently we only attribute two features to the lines:Impedance: the transmission line builder can receive any impedance function as a parameter.\nCurrently we treat this as simple resistance (no AC effects) and adopt a function that is linear on the length of the line.\nThe line length is computed from the geographic distance between buses, ignoring the fact that real lines arc under gravity. Considering the resistivity of the material is already strongly approximated, this should not be too relevant (but can be easily accounted for if necessary).\nWe adopt a single resistivity value of 0.045 Ohms/km. This is a rough estimate based on data sheets for real systems [6],[7].\nWe also provide the non-default option of randomly drawing the resistivity value from a range that is in agreement with the observed in real grids.\nCurrent carrying capacity: the transmission line builder can receive any current capacity function as a parameter.\nBy looking at the charts from [6], we see that capacities tend to increase with increasing line voltage. Based on the values therein, we adopt a capacity of 14 A/kV.\nThe reference voltage is chosen as the one operating voltage that is present in both buses which the line connects. In case there are multiple common voltages, one is randomly drawn. In case there are no common voltages, one is randomly drawn from the set that combines all voltages from those buses."
},

{
    "location": "Model.html#Building-connection-topologies:-1",
    "page": "Model",
    "title": "Building connection topologies:",
    "category": "section",
    "text": "The connections between buses are built by leveraging the method developed and explained in detail by Soltan and Zussman [8]. This method is highly stochastic and aims to reproduce the natural expansion of a power grid. There are two steps in this approach: first, a tunable weight spanning tree is generated, producing a connected network; second, a reinforcement procedure is executed, in order to increase the robustness of the resulting network, in accordance with real cases. Results obtained via this algorithm have been shown to be very realistic [8]."
},

{
    "location": "Model.html#References:-1",
    "page": "Model",
    "title": "References:",
    "category": "section",
    "text": "[1]: Haversine formula: https://en.wikipedia.org/wiki/Haversine_formula[2]: Gegner, Kathleen M., et al. \"A methodology for the creation of geographically realistic synthetic power flow models.\" Power and Energy Conference at Illinois (PECI), 2016 IEEE. IEEE, 2016.[3]: Census data: https://www.census.gov/geo/maps-data/data/gazetteer2010.html[4]: Generator survey: https://www.eia.gov/electricity/data/eia860/index.html[5]: Birchfield, Adam B., et al. \"Grid Structural Characteristics as Validation Criteria for Synthetic Networks.\" IEEE Transactions on Power Systems (2016).[6]: LaForest, J. J. Transmission-line reference book. 345 kV and above. No. EPRI-EL-2500. General Electric Co., Pittsfield, MA (USA). Large Transformer Div.; General Electric Co., Schenectady, NY (USA). Electric Utility Systems Engineering Dept., 1981.[7]: Glover, J. Duncan, Mulukutla S. Sarma, and Thomas Overbye. Power System Analysis & Design, SI Version. Cengage Learning, 2012.[8]: Soltan, Saleh, and Gil Zussman. \"Generation of synthetic spatially embedded power grid networks.\" Power and Energy Society General Meeting (PESGM), 2016. IEEE, 2016."
},

{
    "location": "Functions.html#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "Functions.html#Functions-1",
    "page": "Functions",
    "title": "Functions",
    "category": "section",
    "text": ""
},

{
    "location": "Functions.html#SyntheticGrids.place_loads_from_zips!",
    "page": "Functions",
    "title": "SyntheticGrids.place_loads_from_zips!",
    "category": "Function",
    "text": "place_loads_from_zips!(\n    grid::Grid;\n    latlim = (-Inf,Inf),\n    longlim = (-Inf,Inf),\n    datapath = CENSUSPATH,\n    loadfunc = linearload,\n    voltfunc = randvolt\n)\n\nCreate load buses based on zipcode information.\n\nArguments\n\ngrid: The synthetic grid instance that will be altered.\nlatlim: Tuple of the form (minimum latitude, maximum latitude).\nlonglim: Tuple of the form (minimum longitude, maximum longitude).\ndatapath: Path to the text file containing census information.\nloadfunc: Function for computing bus load based on some rule.\nvoltfunc: Function for computing bus voltage based on some rule.\n\nREFERENCE: https://www.census.gov/geo/maps-data/data/gazetteer2010.html\n\n\n\nplace_loads_from_zips!(\n    grid::Grid,\n    geolim::Function;\n    datapath = CENSUSPATH,\n    loadfunc = linearload,\n    voltfunc = randvolt\n)\n\nCreate load buses based on zipcode information.\n\nArguments\n\ngrid: The synthetic grid instance that will be altered.\ngeolim: Function that receives a latitude-longitude pair and returns true or false.\n\nSpecifies the region of interest.\n\ndatapath: Path to the text file containing census information.\nloadfunc: Function for computing bus load based on some rule.\nvoltfunc: Function for computing bus voltage based on some rule.\n\nREFERENCE: https://www.census.gov/geo/maps-data/data/gazetteer2010.html\n\n\n\n"
},

{
    "location": "Functions.html#SyntheticGrids.place_gens_from_data!",
    "page": "Functions",
    "title": "SyntheticGrids.place_gens_from_data!",
    "category": "Function",
    "text": "place_gens_from_data!(\n    grid::Grid,\n    plants::Array,\n    pfactorfunc::Function = pfacwavg\n)\n\nCreate generation buses based on power plant data.\n\nArguments:\n\ngrid: Grid instance which will be populated.\nplants: Power plant information already parsed.\npfactorfunc: Function for computing bus power factor based on some rule.\n\n\n\nplace_gens_from_data!(\n    grid::Grid;\n    latlim::Tuple = (-Inf, Inf),\n    longlim::Tuple = (-Inf, Inf),\n    jsonpath = GENJSONPATH,\n    pfactorfunc::Function = pfacwavg,\n)\n\nCreate generation buses based on power plant data.\n\nArguments:\n\ngrid: Grid instance which will be populated.\nlatlim: Tuple of the form (minimum latitude, maximum latitude).\nlonglim: Tuple of the form (minimum longitude, maximum longitude).\njsonpath: Path to the JSON file with the generator data.\npfactorfunc: Function for computing bus power factor based on some rule.\n\nREFERENCE: https://www.eia.gov/electricity/data/eia860/index.html\n\n\n\nplace_gens_from_data!(\n    grid::Grid,\n    geolim::Function;\n    jsonpath = GENJSONPATH,\n    pfactorfunc::Function = pfacwavg,\n)\n\nCreate generation buses based on power plant data.\n\nArguments:\n\ngrid: Grid instance which will be populated.\ngeolim: Function that receives a latitude-longitude pair and returns true or false.\n\nSpecifies the region of interest.\n\njsonpath: Path to the JSON file with the generator data.\npfactorfunc: Function for computing bus power factor based on some rule.\n\nREFERENCE: https://www.eia.gov/electricity/data/eia860/index.html\n\n\n\n"
},

{
    "location": "Functions.html#SyntheticGrids.connect!",
    "page": "Functions",
    "title": "SyntheticGrids.connect!",
    "category": "Function",
    "text": "connect!(grid::Grid; k=2.5, m=-1, a=1, b=3.2, g=2.5, t=2, N=10)\n\nConnect buses in an electric grid.\n\nArguments:\n\nk: weight of the spanning tree (see REFERENCE)\nm: total number of connections (computed from the # of nodes if default)\na, b, g, t: paramaters of the model (see REFERENCE)\nN: number of nearest neighbors in the average distance computation\n\nREFERENCE: Soltan, Saleh, and Gil Zussman. \"Generation of synthetic spatially embedded power grid networks.\" Power and Energy Society General Meeting (PESGM), 2016. IEEE, 2016.\n\n\n\n"
},

{
    "location": "Functions.html#SyntheticGrids.cluster!",
    "page": "Functions",
    "title": "SyntheticGrids.cluster!",
    "category": "Function",
    "text": "cluster!(grid::Grid, nloads, nboth, ngens)\n\nCluster all nodes of a grid into 'nloads' load substations, 'ngens' generation substations and 'nboth' substations with both load and generation.\n\nREFERENCE: Birchfield, Adam B., et al. \"Grid Structural Characteristics as Validation Criteria for Synthetic Networks.\" IEEE Transactions on Power Systems (2016).\n\n\n\n"
},

{
    "location": "Functions.html#SyntheticGrids.create_lines!",
    "page": "Functions",
    "title": "SyntheticGrids.create_lines!",
    "category": "Function",
    "text": "create_lines!(grid::Grid; impedfunc = linear_imped, capfunc = volt_cap)\n\nCreate transmission lines for a synthetic grid. Line impedancies are determined by impedfunc and line capacities are calculated via capfunc. This function will only work after the grid has been connected through the connect!() function.\n\n\n\n"
},

{
    "location": "Functions.html#SyntheticGrids.add_load!",
    "page": "Functions",
    "title": "SyntheticGrids.add_load!",
    "category": "Function",
    "text": "add_load!(grid::Grid, args...; reconnect = false)\n\nAdd LoadBus to grid by calling the LoadBus(coords, load, volt, pop, connected_to = Set(), connections = Set()) method with args... as arguments. If reconnect = true, all connections will be remade.\n\n\n\n"
},

{
    "location": "Functions.html#SyntheticGrids.add_gen!",
    "page": "Functions",
    "title": "SyntheticGrids.add_gen!",
    "category": "Function",
    "text": "add_gen!(grid::Grid, args...; reconnect = false)\n\nAdd Genbus to grid by calling the GenBus(coords, gen, volt, tech, connected_to = Set(), connections = Set(), pfactor = -1, summgen = -1, wintgen = -1, gens = []) method with args... as arguments. If reconnect = true, all connections will be remade.\n\n\n\n"
},

{
    "location": "Functions.html#SyntheticGrids.add_substation!",
    "page": "Functions",
    "title": "SyntheticGrids.add_substation!",
    "category": "Function",
    "text": "add_substation!(grid::Grid, args...; reconnect = false)\n\nAdd Substation to grid by calling the Substation(coords, volts, load, gen, pop, con = Set(), group = []) method with args... as arguments.\n\n\n\n"
},

{
    "location": "Functions.html#Main-Functions-1",
    "page": "Functions",
    "title": "Main Functions",
    "category": "section",
    "text": "place_loads_from_zips!\n\nplace_gens_from_data!\n\nconnect!\n\ncluster!\n\ncreate_lines!\n\nadd_load!\n\nadd_gen!\n\nadd_substation!"
},

{
    "location": "Functions.html#Export/Import-Functions-1",
    "page": "Functions",
    "title": "Export/Import Functions",
    "category": "section",
    "text": "to_pandapower\n\nsave_pp_grid\n\nload_pp_grid"
},

{
    "location": "Functions.html#Topology-Building-Functions-1",
    "page": "Functions",
    "title": "Topology Building Functions",
    "category": "section",
    "text": "TWST!\n\nreinforce!\n\nconnect_buses!\n\nconnect_subs!\n\ncluster_loads!\n\ncluster_load_gen!\n\ncluster_gens!"
},

{
    "location": "Functions.html#Checking-Functions-1",
    "page": "Functions",
    "title": "Checking Functions",
    "category": "section",
    "text": "total_links\n\nmean_node_deg\n\ncluster_coeff\n\ncluster_coeff_degw\n\ntest_connectivity\n\nmean_short_path_hop\n\nmean_short_path_real_bus\n\nmean_short_path_real_sub\n\nrobustness_line\n\nrobustness_node"
},

]}
