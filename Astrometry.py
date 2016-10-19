#################################################################
# Name:     Astrometry.py                                       #
# Author:   Yuan Qi Ni                                          #
# Date:     June 22, 2016                                       #
# Function: Program contains functions that perform essential   #
#           astrometric tasks.                                  #
#################################################################

#essential modules
import numpy as np
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_moon

#function: converts isot time to day of year float
def isot_day(time):
    #create astropy time object
    time = Time(time, format='isot', scale='utc')
    #set reference time
    t_ref = "2015-01-01T00:00:00.000"
    t_ref = Time(t_ref, format='isot', scale='utc')
    #return day of year
    return float((time - t_ref).value)

#function: converts day of year float to isot time
def day_isot(day, year):
    #create astropy time object
    t_ref = str(year)+"-01-01T00:00:00.000"
    t_ref = Time(t_ref, format='isot', scale='utc')
    #create astropy time difference object
    t_diff = TimeDelta(day, format='jd')
    #return isot time
    return (t_ref+t_diff).value

#function: return RA and DEC of the moon at utc isot time, location
def moonEQC(time, loc):
    t = Time(time, format='isot', scale='utc')
    l = EarthLocation.from_geodetic(*loc)
    eq = get_moon(t,l,'de430')
    aa = eq.transform_to(AltAz(obstime=t, location=l))
    RA, DEC = eq.ra.degree, eq.dec.degree
    return RA, DEC

#function: return Alt and Az of the moon at utc isot time, location
def moonLC(time, loc):
    t = Time(time, format='isot', scale='utc')
    l = EarthLocation.from_geodetic(*loc)
    eq = get_moon(t,l,'de430')
    lc = eq.transform_to(AltAz(obstime=t, location=l))
    ALT, AZ = lc.alt.degree, lc.az.degree
    return ALT, AZ

#function: return separation angle between sky coordinates
def sepAngle(coord1, coord2):
    coord1 = SkyCoord(*coord1, unit='deg')
    coord2 = SkyCoord(*coord2, unit='deg')
    return coord1.separation(coord2).degree
