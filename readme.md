# PyGeodesy
## What is this?
Geodesy is graphical tool meant to mesure distances and paths along the Earth via several methods

## Graphical interface
The interface allow to set either a pair of position on Earth surface or a starting point and an intended route.

Coordiantes and Azimuths are intended in degrees. Distance unit can be cycled through rightmost button. Just insert plain numbers!

Three different methods are implemented in order to evaluate the geodesic scenario: Vincenty (the most used), Sodano (a non-iterative method, alternative to Vincenty) and shperical Earth.

A database file, originated from [this generator](https://github.com/RossWorks/NAV-DB/), starting from Arinc 424 files may be loaded to allow to search for coded points via their ICAO identificator. The FAA releases a file covering the USA [here](https://www.faa.gov/air_traffic/flight_info/aeronav/digital_products/cifp/download/), available to the general public.


## Implementation details
The main computation code is included in the Geodesy.py module. All functions accepts arguments given in mks system: distances in meters, coordinates in radians. In the complete tool, conversion between units is handled by the GUI.