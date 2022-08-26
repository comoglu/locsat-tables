#!/usr/bin/env python

# Can be either "iasp91" or "ak135".
# Other models are NOT supported by the libtau version used by SC.
model = "iasp91"
prefix = model # for filename generation

# These depths are the same for iasp91 and ak135
ConradDepth = 20.
MohoDepth = 35.

# If True then phases like Pg, Pn, P, Pdiff, PKPdf are combined into
# a single table
combinePhases = True

# If True then a LocSAT incompatible variant of the format will be
# generated that has a few more features. But this is highly
# experimental and just an idea. might also be useful for debugging.
extendedFormat = False

# Add extra comments incl. depth *and* distance before *each* travel
# time line, for debugging etc.
extraComments = False


#################################################################
# No config below
#################################################################

import sys
import os
from numpy import arange, array, concatenate

# Use SC to compute the travel times
from seiscomp.seismology import TravelTimeTable

ttt = TravelTimeTable()
ttt.setModel(model)


class Arrival:
    def __str__(self):
        return "%-12s %8.3f %8.5f %8.5f" % (arr.phase, arr.time, arr.dtdd, arr.toff)

def computeTravelTimesSC(delta, depth):
    delta, depth = float(delta), float(depth)
    arrivals = ttt.compute(0, 0, depth, 0, delta, 0, 0)
    return arrivals

regional_distances = [
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,
    1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5 ]

distances = {
    "P":      regional_distances + list(arange(0, 181,  1.)),
    "Pg":     regional_distances + list(arange(0, 10.2, 0.2)),
    "Pb":     regional_distances + list(arange(0,  9.5, 0.5)),
    "Pn":     regional_distances + list(arange(0, 20.5, 0.5)),
    "S":      regional_distances + list(arange(0, 116,  1.)),
    "Sb":     regional_distances + list(arange(0,  9.5, 0.5)),
    "Sg":     regional_distances + list(arange(0, 10.2, 0.2)),
    "Sn":     regional_distances + list(arange(0, 20.5, 0.5)),

    "pP":     arange( 20, 105, 1),
    "sP":     arange( 20, 105, 1),
    "sS":     arange( 20, 105, 1),

    "PP":     arange( 30, 181, 1),
    "SS":     arange( 30, 181, 1),

    "PcP":    arange( 25,  61, 1),
    "ScS":    arange( 25,  61, 1),
    "ScP":    arange( 25,  61, 1),

    "PKP":    arange( 90, 181, 1),
    "PKPdf":  arange( 90, 181, 1),
    "PKPab":  arange(140, 181, 1),
    "PKPbc":  arange(140, 161, 1),
    "SKPdf":  arange( 90, 181, 1),
    "pPKPdf": arange( 90, 181, 1),
    "pPKPab": arange(140, 181, 1),
    "pPKPbc": arange(140, 161, 1),
    "sPKPdf": arange( 90, 181, 1),
    "sPKPab": arange(140, 181, 1),
    "sPKPbc": arange(140, 161, 1),
}


crustal_depths = [ 0, 5, 10, 15, 20, 25, 30, 35 ]

teleseismic_depths = crustal_depths + [
     40,  50,  75, 100, 150, 200, 250, 300, 350, 400,
    450, 500, 550, 600, 650, 700, 750, 800 ]

depths = {
    "P":      teleseismic_depths,
    "pP":     teleseismic_depths,
    "sP":     teleseismic_depths,
    "Pg":     crustal_depths,
    "Pb":     crustal_depths,
    "Pn":     crustal_depths + [ 40, 50, 75, 100, 150, 200, 250, 300, 400 ],
    "S":      teleseismic_depths,
    "sS":     teleseismic_depths,
    "Sb":     crustal_depths,
    "Sg":     crustal_depths,
    "Sn":     crustal_depths + [ 40, 50, 75, 100, 150, 200, 250, 300, 400 ],

    "PP":     teleseismic_depths,
    "SS":     teleseismic_depths,

    "PcP":    teleseismic_depths,
    "ScS":    teleseismic_depths,
    "ScP":    teleseismic_depths,

    "PKP":    teleseismic_depths,
    "PKPdf":  teleseismic_depths,
    "PKPab":  teleseismic_depths,
    "PKPbc":  teleseismic_depths,
    "SKPdf":  teleseismic_depths,
    "pPKPdf": teleseismic_depths,
    "pPKPab": teleseismic_depths,
    "pPKPbc": teleseismic_depths,
    "sPKPdf": teleseismic_depths,
    "sPKPab": teleseismic_depths,
    "sPKPbc": teleseismic_depths,
}


distanceLimits = {
    "Pdiff": (   0, 116 ),
    "PKPdf": ( 118, 180 ),
    "Sdiff": (   0, 116 ),
    "SKSdf": ( 118, 180 ),
}


phaseCombinations = {
    # Make P a generic phase of the first-arriving P wave.
    # Note that there are gaps in the slope at the 
    # transitions Pg - Pb - Pn (not P - Pdiff)
    "P"    :  [ "Pg", "Pb", "Pn", "P", "Pdiff", "PKPdf" ],
    "PKP"  :  [ "PKPab", "PKPbc", "PKPdf" ],

    "pP"   :  [ "pP", "pPn", "pPdiff" ],
    "sP"   :  [ "sP", "sPn", "sPdiff" ],
    "PP"   :  [ "PP", "PnPn" ],

    # Generic S phase. See P for comments.
    "S"    :  [ "S","Sn","Sg","Sb", "Sdiff" ],
    "sS"   :  [ "sS", "sSg", "sSb", "sSn", "sSdiff" ],
    "SS"   :  [ "SS", "SnSn" ],

    # Crustal phases
    # NOTE: Pg and Sg are computed independently

    "Pg"   :  [ "Pg", "PgPg"],
    "Sg"   :  [ "Sg", "SgSg"],
}


def print_as_10_columns(xx):
    output = ""
    line = ""
    for x in xx:
        try:
            line += "%8.2f" % x
        except:
            print(x, xx)
        if len(line) >= 80:
            output += "%s\n" % line
            line = ""
    output += "%s\n" % line
    return output


def create_table(phase):

    if extendedFormat:
        output = "%s\n" % phase
    else:
        output = "# %s travel-time tables\n" % phase
    output += "%d    # number of depth samples\n" % len(depths[phase])
    output += print_as_10_columns(list(depths[phase]))
    output += "%d    # number of distance samples\n" % len(distances[phase])
    output += print_as_10_columns(list(distances[phase]))

    for z in depths[phase]:

        ph = phase
        # special case of depth phase at zero depth
        if z == 0 and ph[0] in ["p", "s"]:
            ph = ph[1:]

        if extendedFormat:
            output += "# z = %g km\n" % z
        else:
            output += "# Travel time for z = %g km\n" % z

        for d in distances[phase]:
            arrivals = computeTravelTimesSC(d, z)
            found = False
            for arr in arrivals:
                if arr.phase in distanceLimits:
                    dmin, dmax = distanceLimits[arr.phase]
                    if not dmin <= d <= dmax:
                        continue

                if ph in phaseCombinations:
                    if arr.phase in phaseCombinations[ph]:
                        tt = arr.time
                        found = True
                        break
                else:
                    if arr.phase == ph:
                        tt = arr.time
                        found = True
                        break

            if not found:
                arr = None

            if phase == "Pg":
                # We approximate Pg time using a nearly constant velocity
                # of 5.8 km/s, which is the IASP91 upper-crust velocity,
                # and letting it increase very slightly with source
                # depth The advantage of this approximation is the full
                # coverage of the whole crustal depth range with Pg
                # arrivals and smooth travel times.

                # The CTBTO uses a constant vp of 6.1 km/s for Pg at
                # distances above 1.5 degrees, but this is not in agreement
                # with IASP91. It probably works well in many regions,
                # though.
                #
                # See also https://doi.org/10.1111/j.1365-246X.2009.04109.x

                arr = Arrival()
                v = 5.8 + 0.01*z + 0.01*d
                v = 5.8 + 0.0006*z
                arr.phase = "Pg"
                arr.time = ((111.195*d)**2+z**2)**0.5/v
                arr.dtdd = -1
                arr.dtdh = -1

            if phase == "Sg":
                arr = Arrival()
                # We approximate Sg time using a nearly constant velocity
                # of 3.36 km/s, which is the IASP91 upper-crust velocity.
                # Only very small depth gradient gives very good agreement
                # with IASP91 even without using the libtau routines.
                v = 3.36 + 0.00037*z
                arr.phase = "Sg"
                arr.time = ((111.195*d)**2+z**2)**0.5/v
                arr.dtdd = -1
                arr.dtdh = -1

            if extraComments:
                output += "# depth: %8.3f    distance: %.3f\n" % (z,d)

            if arr:
                if z==0. and d==0.:
                    arr.time = 0.
                if extendedFormat:
                    # t, p=dt/dx, q=dt/dz
                    # i.e.: time, horizontal slowness, vertical slowness
                    output += "%12.3f %12.6f %12.6f %s\n" % (arr.time, arr.dtdd, arr.dtdh, arr.phase)
                else:
                    output += "%12.3f\n" % arr.time
            else:
                if extendedFormat:
                    output += "%12.3f %12.6f %12.6f\n" % (-1, -1, -1)
                else:
                    output += "%12.3f\n" % -1

    return output


defaultPhaseSet = """
    P Pg Pb Pn S Sg Sb Sn
    PP SS
    pP sP sS
    PKP PKPab PKPbc PKPdf pPKPab pPKPbc pPKPdf sPKPab sPKPbc sPKPdf SKPdf
    PcP ScS ScP
"""

if sys.argv[1:]:
    phaseSet = sys.argv[1:]
else:
    phaseSet = defaultPhaseSet.strip().split()

try:
    os.makedirs("tables")
except FileExistsError:
    pass

for phase in phaseSet:
    table = create_table(phase)
    filename = os.path.join("tables", prefix + "." + phase)
    with open(filename, "w") as f:
        f.write(table)
