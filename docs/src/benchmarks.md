# Performance Benchmarks

The performance of this package have been tested against both CALCEPH and SPICE, two of the most-popular open-source ephemeris readers used in the space industry.

The results show that `Ephemerides.jl` largely outperforms SPICE as well as CALCEPH for most SPK segment types. For example, for state vector computations (i.e., position and velocity) the mean execution times are the following:

```@raw html
<p align="center">
<img src="https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl/assets/85893254/e54ac790-7421-47ff-9b68-d35bdea74de5" width="512"/>
</p>
```

Additionally, it is better optimised to compute higher order derivatives (i.e., acceleration and jerk) with respect to CALCEPH. 

```@raw html
<p align="center">
<img src="https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl/assets/85893254/ec2df247-9dde-44e9-82f4-4e8372317b8e" width="512"/>
</p>
```

!!! note 
    These time benchmarks have been obtained on an Intel Core i7-6700 CPU @ 3.40 GHz with 16 GB of RAM