# Creation of travel time tables for LocSAT

The script make-locsat-tables-iasp91.py creates travel time tables for the IASP91 model in LocSAT format to be used in SeisComP. It is provided here for documentation and to make it easier to add new tables.

The script uses SeisComP and therefore can generate travel-time tables for all phases available in the libtau version used in SeisComP. The default phase set includes approximately the phases you might want to use in event location using LocSAT, including some phases like PP, which are normally not used in event location but which we want to associate and thus compute proper residuals for.

## Invocation

Invocation is as simple as:

```
~/seiscomp/bin/seiscomp exec seiscomp-python make-locsat-tables-iasp91.py
```
or for a specific phase, e.g. `Pg`:
```
~/seiscomp/bin/seiscomp exec seiscomp-python make-locsat-tables-iasp91.py Pg
```



## Remarks:

### Depth phases

We produce a depth phase table for depth of 0 by simply adopting the P (S)
table for a depth of 0 as for that depth no depth phases are generated. We do need depth phase tables for a depth of 0 in order to cover the entire depth range.


### Crustal phases

The IASP91 velocity model has a Conrad discontinuity at a depth of 20 km. This implies that for the IASP91 no Pg exists for focal depths larger than 20 km. While formally OK, in practice this may create problems if Pg is used for event location, which is usually the case in local to regional setups. The difficulty is that while the IASP91 model has a sharp Conrad discontinuity at 20 km depth, the crust in the vicinity of the network may not have and thus Pg may be observed and picked at distances/depths where there is formally no Pg in IASP91.

In order to ensure continuous Pg/Sg travel time tables spanning the distance range from 0 to 10 degrees and event depth range 0 to 35 km, we compute travel times using a constant IASP91 upper-crustal velocity of 5.8 km/s (vp) and 3.36 km/s (vs) slightly increasing with depth to match the libtau times where available. For many networks, 5.8 (3.36) km/s might not be the best choice [1] and by changing the velocities the Pg/Sg tables can easily be customized. But make sure you don't call them "iasp91" to avoid confusion.

Please note that the Pg/Sg tables are defaults for IASP91 and obviously don't give the most accurate locations for regions in which the crustal velocities deviate significantly from IASP91.


[1] https://doi.org/10.1111/j.1365-246X.2009.04109.x
