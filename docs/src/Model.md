# Model

The idea behind the `SyntheticGrids` module is to provide a standalone and easily expandable framework for generating synthetic grids. The goal is to be able to use different levels of approximation when building grids. Interfacing with other packages should be possible through export functions.

A synthetic grid is constituted basically by buses connected via transmission lines. There are several relevant points when automatically building a network with this package. We discuss these on a *per* type basis below.

## Bus
`Bus` is an abstract type that acts as parent to the `LoadBus`, `GenBus` and `Substation` types. Pragmatically, the one common feature to all of these is their geographical location being specified by latitude and longitude values. This is a better choice than adopting cartesian coordinates, as that would require some kind of projection scheme, such as the Mercator or the UTM one. The issue with adopting a projection scheme is that none of the available ones are adequate for being applied over very large areas, such as the territory span of the USA. Therefore, coordinates are kept in a latitude-longitude format and, whenever we need to compute real distances, we adopt the haversine formula [^1] with a radius of 6369.783km, which corresponds to Earth's radius at latitude 38.8 (the center of the region we model).

Bus siting and sizing is done according to the approach presented in reference [^2].

## LoadBus
The main properties of a `LoadBus` are geographical location, load and operating voltage values. Since none of this information is public, we have to adopt an approach based on some kind of public data. Overbye *et al.* [^2] have devised the approach which we adopt in the module.

* **Geographical location:** load buses are created so as to coincide with each zip code. Official census data [^3] provides us with latitude and longitude and population of each zip code in the USA.
  - We adopt the census data from 2010. There are more recent censi available, but these do not offer population data (whose importance will be explained in a moment).
* **Load:** we have no public information on load for specific zip codes. Overbye *et al.* [^2] have shown that load is approximately linear with the population size. Therefore, they propose adopting a scaling factor for computing load from population size.
  - We adopt their approximation of 2kW load per capita.
  - The code is written in a way that a function is a passed into the zip code builder for computing the loads. Thus, it is extremely straight-forward to adopt any other form of computing loads.
* **Operating voltages:** there is no public information on operating voltages of load buses. At this point we can either adopt some constant value or randomly draw from a set of allowed values.
  - For a network with several operating voltages, it is noted that usually loads are connected to the lower voltages.

## GenBus
`GenBus` is the type that represents power plants. The reason for creating a bus for each power plant, instead of a bus for each `Generator`, is the fact that the data we use for placing generation buses is organized in this way. There is an official public yearly survey on all generators that operate in the USA [^4] that provides us with most of the important information on power plants.

* **Geographical location:** all power plants operating in the USA are listed with their address, latitude and longitude, plant code and other information. A separate file provides information on each of the generators that operate inside a power plant.
  - This information is not completely clean. There are a few plants marked as *unsited* (which provide no information) and a few with wrong latitude-longitude values (the obvious ones were manually corrected).
  - The most concerning issue is that several power plants have zero operating generators reported. This may be due to incomplete survey data or due to those not having any generator that is not already retired or still just planned. These plants are ignored when creating grids.
  - The data on planned and retired generators is completely ignored when we read in the information.
  - Some generators also have a future data as a retirement prediction. This is not taken into account.
* **Generation Capacity:** this is obtained by directly summing the nominal capacities of all operational and standby generators in a given plant.
* **Voltages:** up to three different voltage values are attributed to each power plant. These are likely the voltages of different generators, but are provided on a plant-basis, not on a generator-basis, so we are not able to break them down.
* **Technology type:** each plant may have multiple generator types (e.g. nuclear, wind, solar etc.), according to data discriminated in the survey.
* **Power factor:** the power factor of a power plant is computed as a nominal capacity-weighted average of the power factors of each generator. The code accepts any function for computing the power factor of power plant, so it is trivial to adopt different approximations if necessary.
* **Summer and Winter generation capacities:** these are taken as the sum of the corresponding quantities of each generator.

## Generators

All data on generators is directly extracted from the power plant survey [^4], which provides us with:

* **Geographical location:** every generator is directly associated with a power plant. Thus, we adopt the power plant's latitude and longitude as the generator's location.
* **Operating voltages:** although a given generator probably has a single operating voltage, this information is given to us on a plant-basis, not on a generator-basis. Therefore, each generator will be attributed the whole set of voltages from the parent plant.
* **Technology type:** directly obtained from survey data. Generators with different technology may be located at the same plant.
* **Generation capacity:** directly obtained from survey data.
* **Power factor:** directly obtained from survey data. In case it is not given, we assume it equals 1.0.
* **Minimum load:** directly obtained from survey data. In case it is not given, we assume it equals 0.0.
* **Summer and Winter generation capacity:** directly obtained from survey data. In case it is not given, we assume it equals the nominal generation capacity.
* **Time from cold shutdown to full load:** currently we grab it without even parsing, but it is kept due to probably being useful in the future. Several generators do not offer this information.
* **Status:** Informs us if the generator was operational, on standby or out of service for most of the survey year.

## Substations
A substation groups several other buses into an aggregate structure. We are directly implementing the clustering algorithm presented by Birchfield *et al.* in reference [^5]. Substations are divided in three kinds: pure load, pure generation and both load and generation. Transmission substations are ignored since we are not using connection buses.

## TransLines
Transmission lines are placed at every connection between two buses as determined by our connection topology builder (explained ahead). Currently we only attribute two features to the lines:

* **Impedance:** the transmission line builder can receive any impedance function as a parameter.
  - Currently we treat this as simple resistance (no AC effects) and adopt a function that is linear on the length of the line.
  - The line length is computed from the geographic distance between buses, ignoring the fact that real lines arc under gravity. Considering the resistivity of the material is already strongly approximated, this should not be too relevant (but can be easily accounted for if necessary).
  - We adopt a single resistivity value of 0.045 Ohms/km. This is a rough estimate based on data sheets for real systems [^6],[^7].
  - We also provide the non-default option of randomly drawing the resistivity value from a range that is in agreement with the observed in real grids.
* **Current carrying capacity:** the transmission line builder can receive any current capacity function as a parameter.
  - By looking at the charts from [^6], we see that capacities tend to increase with increasing line voltage. Based on the values therein, we adopt a capacity of 14 A/kV.
  - The reference voltage is chosen as the one operating voltage that is present in both buses which the line connects. In case there are multiple common voltages, one is randomly drawn. In case there are no common voltages, one is randomly drawn from the set that combines all voltages from those buses.

## Building connection topologies:
The connections between buses are built by leveraging the method developed and explained in detail by Soltan and Zussman [^8]. This method is highly stochastic and aims to reproduce the natural expansion of a power grid. There are two steps in this approach: first, a tunable weight spanning tree is generated, producing a connected network; second, a reinforcement procedure is executed, in order to increase the robustness of the resulting network, in accordance with real cases. Results obtained via this algorithm have been shown to be very realistic [^8].

## References:
[^1]: Haversine formula: https://en.wikipedia.org/wiki/Haversine_formula

[^2]: Gegner, Kathleen M., et al. "A methodology for the creation of geographically realistic synthetic power flow models." Power and Energy Conference at Illinois (PECI), 2016 IEEE. IEEE, 2016.

[^3]: Census data: https://www.census.gov/geo/maps-data/data/gazetteer2010.html

[^4]: Generator survey: https://www.eia.gov/electricity/data/eia860/index.html

[^5]: Birchfield, Adam B., et al. "Grid Structural Characteristics as Validation Criteria for Synthetic Networks." IEEE Transactions on Power Systems (2016).

[^6]: LaForest, J. J. Transmission-line reference book. 345 kV and above. No. EPRI-EL-2500. General Electric Co., Pittsfield, MA (USA). Large Transformer Div.; General Electric Co., Schenectady, NY (USA). Electric Utility Systems Engineering Dept., 1981.

[^7]: Glover, J. Duncan, Mulukutla S. Sarma, and Thomas Overbye. Power System Analysis & Design, SI Version. Cengage Learning, 2012.

[^8]: Soltan, Saleh, and Gil Zussman. "Generation of synthetic spatially embedded power grid networks." arXiv:1508.04447 [cs.SY], Aug. 2015.
