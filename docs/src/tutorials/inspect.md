# Kernels Inspection

This tutorial will walk you through the basic features and interfaces that allow you to inspect the properties of binary ephemeris kernels.

## Available Times 

### Timescale
SPK and PCK segments support two types of ephemeris timescales, namely, TDB and TCB. The timescale of the loaded kernels can be retrieved as follows: 

```julia
using Ephemerides 

# Load the kernel
eph = EphemerisProvider("kernel.bsp")

# Retrieve the ID of the kernel timescales
id = ephem_timescale_id(eph)
```
The retrieved ID is 1 for TDB and 2 for TCB. A value of -1 is returned if the kernels are empty. 

!!! note 
    Only one timescale is admissed within a single `EphemerisProvider` object, 

### Timespan

To retrieve the first and last available time in the ephemeris files associated to a provider 
object, two functions are available to distinguish between SPK and PCK data: 

```
ephem_spk_timespan(eph)
ephem_pck_timespan(eph)
```

where `eph` is an `EphemerisProvider` instance. Both functions return the minimum and maximum 
available time in TDB seconds since J2000, as well as a continuity parameter defined as 
follows: 

- **0** no SPK or PCK data is available.
- **1** the quantities of all bodies are available for any time between the first and last time.
- **2** the quantities of some bodies are available on discontinuous time intervals between the 
    first and last time.
- **3** the quantities of each body are available on a continuous time interval between the first 
    and the last time, but not available for any time between the first and last time.

## Available Points and Axes

To retrieve the list of NAIF IDs with the points or axes that have available ephemeris data, 
these function should be called: 
```
ephem_get_points(eph)
ephem_get_axes(eph)
```

## Segment Records

Position and orientation metadata relative to the records loaded in the ephemeris kernels can be retrieved with the following two functions, respectively:

```
ephem_spk_records(eph)
ephem_pck_records(eph)
```

Both functions return a vector of `EphemRecordSPK` or `EphemRecordPCK` ordered by priority, i.e., they use the highest priority records when there are multiple records that could satisfy the same target, center pair for a given epoch.

In particular SPK records contain the following information: 

- target: NAIF ID of the target object
- center: NAIF ID of the center object
- axes: NAIF ID of the reference axes
- t_start: start times of each sub-window, in TDB seconds since J2000
- t_end: final times of each sub-window, in TDB seconds since J2000

Similarly, PCK records contain these information:

- target: NAIF ID of the target axes
- center: NAIF ID of the reference axes
- t_start: start times of each sub-window, in TDB seconds since J2000
- t_end: final times of each sub-window, in TDB seconds since J2000

For PCK records, the reference axes ID is set into the `center` field. Please notice that 
whenever a gap between the data of a given pair of (center, target) objects is present, 
`t_start` and `t_end` will store the start and end times of each window with available data, respectively.