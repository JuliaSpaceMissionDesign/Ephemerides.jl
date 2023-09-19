# [Low-level API](@id ephemerides_api)

This functions are not meant to be used outside of the package. They are documented 
only to aid future developments of the package.

## Ephemeris Provider 

```@docs
Ephemerides.get_daf 
Ephemerides.spk_links 
Ephemerides.pck_links

```

## DAF Routines 

### DAF Header
```@docs 
Ephemerides.DAFHeader 

Ephemerides.initial_record 
Ephemerides.final_record 
Ephemerides.free_address
Ephemerides.endian 
Ephemerides.filename 
Ephemerides.summary_size
```

### DAF Descriptor 
```@docs
Ephemerides.DAFSegmentDescriptor 

Ephemerides.segment_type 
Ephemerides.initial_time
Ephemerides.final_time
Ephemerides.center
Ephemerides.target
Ephemerides.axes
Ephemerides.initial_address
Ephemerides.final_address

Ephemerides.parse_spk_segment_descriptor
Ephemerides.parse_pck_segment_descriptor
```

```@docs
Ephemerides.DAF

Ephemerides.comment
Ephemerides.header
Ephemerides.array
Ephemerides.descriptors
Ephemerides.segment_list
Ephemerides.filepath
Ephemerides.is_spk
Ephemerides.is_pck

Ephemerides.parse_daf_comment
Ephemerides.parse_daf_summaries

Ephemerides.initialise_segments!
Ephemerides.create_spk_segment
```

## SPK Links 
```@docs 
Ephemerides.SPKLink 
Ephemerides.SPKLinkTable

Ephemerides.descriptor 
Ephemerides.file_id
Ephemerides.list_id
Ephemerides.element_id
Ephemerides.factor
Ephemerides.reverse_link
Ephemerides.create_linktables
Ephemerides.add_spklinks!

```

## SPK Segment Types 

### SPK Type 1 and 21
```@docs 
Ephemerides.SPKSegmentHeader1
Ephemerides.SPKSegmentCache1 
Ephemerides.SPKSegmentType1
```

### SPK Type 2 and 3
```@docs 
Ephemerides.SPKSegmentHeader2
Ephemerides.SPKSegmentCache2
Ephemerides.SPKSegmentType2
```

### SPK Type 8 and 12
```@docs 
Ephemerides.SPKSegmentHeader8
Ephemerides.SPKSegmentCache8
Ephemerides.SPKSegmentType8
```

### SPK Type 9 and 13
```@docs 
Ephemerides.SPKSegmentHeader9
Ephemerides.SPKSegmentCache9
Ephemerides.SPKSegmentType9
```

### SPK Type 18 and 19
```@docs 
Ephemerides.SPKSegmentHeader18
Ephemerides.SPKSegmentCache18
Ephemerides.SPKSegmentHeader19
Ephemerides.SPKSegmentCache19
Ephemerides.SPKSegmentType19
```

### SPK Type 20
```@docs 
Ephemerides.SPKSegmentHeader20
Ephemerides.SPKSegmentCache20
Ephemerides.SPKSegmentType20
```


## Interpolating Functions

### Chebyshev Polynomials 

```@docs
Ephemerides.chebyshev
Ephemerides.∂chebyshev
Ephemerides.∂²chebyshev
Ephemerides.∂³chebyshev
Ephemerides.∫chebyshev
```

### Lagrange Polynomials 

```@docs 
Ephemerides.lagrange
Ephemerides.∂lagrange
Ephemerides.∂²lagrange
```

### Hermite Polynomials 

```@docs 
Ephemerides.hermite
Ephemerides.∂hermite
Ephemerides.∂²hermite
Ephemerides.∂³hermite
```




