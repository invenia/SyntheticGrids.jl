# Custom Types

The Synthetic Grids package implements several custom types.

## Bus

### description:
`Bus` is an abstract type that acts as parent to the `LoadBus`, `GenBus` and `Substation` types.


## LoadBus

### description:
Objects of this type represent the electrical loads in a network.

### attributes:

* **coords**: `Geodesy.LatLon`. Geocoordinates of the bus.
* **load**: `Real`. Load value at bus. Positive MW by convention.
* **voltage**: `Real`. Bus operating voltage. Positive kV by convention.
* **population**: `Integer`. Bus population. Used for estimating load.
* **connected_to**: `Set{Bus}`. Buses to which it is connected.
* **connections**: `Set{TransLine}`. Transmission lines connected to bus.


## GenBus

### description:
`GenBus` is the type that represents power plants.

### attributes:

* **coords**: `Geodesy.LatLon`. Geocoordinates of the bus.
* **generation**: `Real`. Nameplate capacity value at bus. Positive MW by convention.
* **voltage**: `Vector{Real}`. Bus operating voltages. Multiple values allowed since power plants may have generators that operate at different voltage levels. Positive kV by convention.
* **tech_type**: `Vector{AbtractString}`. Technology types of the generators at the power plant.
* **connected_to**: `Set{Bus}`. Buses to which it is connected.
* **connections**: `Set{TransLine}`. Transmission lines connected to bus.
* **pfactor**: `Real`. Power factor of the power plant.
* **summgen**: `Real`. Summer generation capacity. Positive MW by convention.
* **wintgen**: `Real`. Winter generation capacity. Positive MW by convention.
* **gens**: `Vector{Generators}`. Generators that operate at the power plant.

## Substation

### description:
A substation groups several other buses into an aggregate structure.

### attributes:
* **coords**: `Geodesy.LatLon`. Geocoordinates of the bus.
* **generation**: `Real`. Nameplate capacity value at substation. Positive MW by convention.
* **voltage**: `Vector{Real}`. Substation operating voltages. Multiple values allowed since different buses may operate at different voltage levels. Positive kV by convention.
* **load**: `Real`. Load value at substation. Positive MW by convention.
* **population**: `Integer`. Substation population. Used for estimating load.
* **connected_to**: `Set{Substation}`. Substations to which it is connected.
* **grouping**: `Vector{Bus}`. Buses grouped by substation.


## Generator

### description:
This type models each specific generator that is part of a power plant (GenBus).

### attributes:
* **coords**: `Geodesy.LatLon`. Geocoordinates of the generator.
* **volt**:  `Vector{Real}`. Generator operating voltages. Multiple values allowed since the generator survey provides voltage levels grouped by power plant (see [Model](Model.md)). Positive kV by convention.
* **tech**: `AbstractString`. Technology type.
* **cap**: `Real`. Nameplate capacity.
* **pfactor**: `Real`. Power Factor.
* **minload**: `Real`. Minimum load value. Positive MW by convention.
* **scap**: `Real`. Summer generation capacity.
* **wcap**: `Real`. Winter generation capacity.
* **shut2loadtime**: `AbstractString`. Time for going from cold shutdown to full load. Currently grabbed from the generator survey without parsing (see [Model](Model.md)).
* **status**: `AbstractString`. "OP" - operational; "SB" - standby; "OA" and "OS" - out of service.


## TransLine

### description:
This type models transmission lines that connect different buses.

### attributes:
* **connecting**: `Tuple{Bus,Bus}`. Buses that are being connected by the transmission line.
* **impedance**: `Real`. Impedance of the line in Ohms.
* **capacity**: `Real`. Current carrying capacity of the line in Amperes.


## Grid

### description:
Represents an entire power grid.

### attributes:
* **seed**: `Integer`. Seed number to allow reproducibility of all stochastic procedures.
* **buses**: `Vector{Buses}`. Vector grouping all LoadBus and GenBus objects in a network.
* **trans_lines**: `Vector{TransLine}`. Vector grouping all transmission lines in a network.
* **substations**: `Vector{Substations}`. Vector grouping all substations in a network.
* **bus_conn**: `AbstractMatrix{Bool}`. Boolean matrix indicating which buses are connected to which.
* **sub_conn**: `AbstractMatrix{Int}`. Matrix indicating the number of connections between each sbustation.
