# SOT.jl

The code in this package processes seismic and acoustic data to infer temperature change in the ocean. The principal functions to be performed are:

1. Obtain an earthquake catalog from ISC.
2. Enhance the catalog through template matching at reference stations.
3. Find *P*-wave pairs at the reference stations. [Done.]
4. Cross-correlate the waveforms of the identified *P*-wave pairs at the *T*-wave station. [Done.]
5. Invert the measured travel time changes for travel time anomalies and apply a cycle skipping correction. [Done.]

To run, enter the package manager with `]` and type `activate .`. Then exit the package manager with a backspace and execute a script, e.g., `include("examples/nias_isc_DGAR.jl")`.
