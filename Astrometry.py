#################################################################
# Name:     Astrometry.py                                       #
# Author:   Yuan Qi Ni                                          #
# Version:  August 25, 2016                                     #
# Function: Program contains functions that perform essential   #
#           astrometric tasks.                                  #
#################################################################

#essential modules
import numpy as np

#function converts degrees to sexagesimal (hms, dms)
def deg_sex(ra, dec):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    return c.to_string('hmsdms')

#function converts ra, dec to galactocentric
def Hel_toGal(ra, dec):
    import astropy.coordinates as coord
    import astropy.units as u
    
    c1 = coord.SkyCoord(ra=ra*u.degree, dec=dec*u.degree,frame='icrs')
    return c1.galactic.l.value, c1.galactic.b.value #(l, b)

#function: converts KSPTime to isot time
def ksp_isot(time):
    '''
    ################################################################
    # Desc: Converts KSPTime to isot time.                         #
    # ------------------------------------------------------------ #
    # Imports:                                                     #
    # ------------------------------------------------------------ #
    # Input                                                        #
    # ------------------------------------------------------------ #
    # time: str time format in KSP time, YYMMDD_HHMM               #
    # ------------------------------------------------------------ #
    # Output                                                       #
    # ------------------------------------------------------------ #
    # time: str time format in ISOT time, YYYY-MM-DDTHH:MM:SS.SSS  #
    ################################################################
    '''
    return "20"+time[:2]+"-"+time[2:4]+"-"+time[4:6]+"T"+time[7:9]+":"+time[9:11]+":00.000"

#function: converts isot time to day of year float
def isot_day(time, year):
    '''
    ################################################################
    # Desc: Converts isot time to day of year float.               #
    # ------------------------------------------------------------ #
    # Imports: astropy.time.Time                                   #
    # ------------------------------------------------------------ #
    # Input                                                        #
    # ------------------------------------------------------------ #
    # time: str time format in ISOT UTC, YYYY-MM-DDTHH:MM:SS.SSS   #
    # ------------------------------------------------------------ #
    # Output                                                       #
    # ------------------------------------------------------------ #
    # day: float time in days since start of year YYYY             #
    ################################################################
    '''
    
    from astropy.time import Time
    
    #create astropy time object
    time = Time(time, format='isot', scale='utc')
    #set reference time
    t_ref = str(year)+"-01-01T00:00:00.000"
    t_ref = Time(t_ref, format='isot', scale='utc')
    #return day of year
    return float((time - t_ref).value)

#function: converts day of year float to isot time
def day_isot(day, year):
    '''
    ################################################################
    # Desc: Converts day of year float to isot time.               #
    # ------------------------------------------------------------ #
    # Imports: astropy.time.(Time, TimeDelta)                      #
    # ------------------------------------------------------------ #
    # Input                                                        #
    # ------------------------------------------------------------ #
    #  day: float time in days since start of year YYYY            #
    # year: int reference year YYYY                                #
    # ------------------------------------------------------------ #
    # Output                                                       #
    # ------------------------------------------------------------ #
    # time: str time format in ISOT UTC, YYYY-MM-DDTHH:MM:SS.SSS   #
    ################################################################
    '''
    
    from astropy.time import Time, TimeDelta
    
    #create astropy time object
    t_ref = str(year)+"-01-01T00:00:00.000"
    t_ref = Time(t_ref, format='isot', scale='utc')
    #create astropy time difference object
    t_diff = TimeDelta(day, format='jd')
    #return isot time
    return (t_ref+t_diff).value

def day_mjd(day, year):
    '''
    ################################################################
    # Desc: Converts day of year float to mjd time.               #
    # ------------------------------------------------------------ #
    # Imports: astropy.time.(Time, TimeDelta)                      #
    # ------------------------------------------------------------ #
    # Input                                                        #
    # ------------------------------------------------------------ #
    #  day: float time in days since start of year YYYY            #
    # year: int reference year YYYY                                #
    # ------------------------------------------------------------ #
    # Output                                                       #
    # ------------------------------------------------------------ #
    # time: float time format in mjd                                 #
    ################################################################
    '''
    
    from astropy.time import Time, TimeDelta
    
    #create astropy time object
    t_ref = str(year)+"-01-01T00:00:00.000"
    t_ref = Time(t_ref, format='isot', scale='utc')
    #create astropy time difference object
    t_diff = TimeDelta(day, format='jd')
    #return isot time
    return (t_ref+t_diff).mjd

#function: return RA and DEC of the moon at utc isot time, location
def moonEQC(time, loc):
    '''
    #################################################################
    # Desc: Return RA, DEC of the moon at utc isot time, location.  #
    # ------------------------------------------------------------- #
    # Imports: astropy.time.Time                                    #
    #          astropy.coordinates.(Earthlocation, AltAz, get_moon) #
    # ------------------------------------------------------------- #
    # Input                                                         #
    # ------------------------------------------------------------- #
    # time: str time format in ISOT, UTC                            #
    #  loc: iterable floats position coordinate in [long, lat, alt] #
    # ------------------------------------------------------------- #
    # Output                                                        #
    # ------------------------------------------------------------- #
    # RA, DEC: float equatorial coordinate position in degree       #
    #################################################################
    '''

    from astropy.time import Time
    from astropy.coordinates import EarthLocation, AltAz, get_moon
    
    t = Time(time, format='isot', scale='utc')
    l = EarthLocation.from_geodetic(*loc)
    eq = get_moon(t,l,'de440')
    aa = eq.transform_to(AltAz(obstime=t, location=l))
    RA, DEC = eq.ra.degree, eq.dec.degree
    return RA, DEC

#function: return Alt and Az of the moon at utc isot time, location
def moonLC(time, loc):
    '''
    #################################################################
    # Desc: Return ALT, AZ of the moon at utc isot time, location.  #
    # ------------------------------------------------------------- #
    # Imports: astropy.time.Time                                    #
    #          astropy.coordinates.(Earthlocation, AltAz, get_moon) #
    # ------------------------------------------------------------- #
    # Input                                                         #
    # ------------------------------------------------------------- #
    # time: str time format in ISOT, UTC                            #
    #  loc: iterable floats position coordinate in [long, lat, alt] #
    # ------------------------------------------------------------- #
    # Output                                                        #
    # ------------------------------------------------------------- #
    # ALT, AZ: float alt-az coordinate position in degree           #
    #################################################################
    '''

    from astropy.time import Time
    from astropy.coordinates import EarthLocation, AltAz, get_moon
    
    t = Time(time, format='isot', scale='utc')
    l = EarthLocation.from_geodetic(*loc)
    eq = get_moon(t,l,'de440')
    lc = eq.transform_to(AltAz(obstime=t, location=l))
    ALT, AZ = lc.alt.degree, lc.az.degree
    return ALT, AZ

#function: return separation angle between sky coordinates
def sepAngle(coord1, coord2):
    '''
    #################################################################
    # Desc: Return separation angle between spherical coordinates.  #
    # ------------------------------------------------------------- #
    # Imports: astropy.coordinates.SkyCoord                         #
    # ------------------------------------------------------------- #
    # Input                                                         #
    # ------------------------------------------------------------- #
    # coord1: (X, Y) angular coordinate 1                           #
    # coord2: (X, Y) angular coordinate 2                           #
    # ------------------------------------------------------------- #
    # Output                                                        #
    # ------------------------------------------------------------- #
    # sep: float degree angular distance between coord1 and coord2  #
    #################################################################
    '''

    from astropy.coordinates import SkyCoord
    
    coord1 = SkyCoord(*coord1, unit='deg')
    coord2 = SkyCoord(*coord2, unit='deg')
    return coord1.separation(coord2).degree

#function: returns separation angle on single plate (<1deg)
def smallAngle(coord1, coord2):
    '''
    #################################################################
    # Desc: Returns separation angle on single plate (<1deg).       #
    # ------------------------------------------------------------- #
    # Imports:                                                      #
    # ------------------------------------------------------------- #
    # Input                                                         #
    # ------------------------------------------------------------- #
    # coord1: (X, Y) angular coordinate 1                           #
    # coord2: (X, Y) angular coordinate 2                           #
    # ------------------------------------------------------------- #
    # Output                                                        #
    # ------------------------------------------------------------- #
    # sep: float degree angular distance between coord1 and coord2  #
    #################################################################
    '''
    
    da = coord1[0] - coord2[0]
    dd = coord1[1] - coord2[1]
    d = 0.5*(coord1[1] + coord2[1])
    return np.sqrt(np.square(da*np.cos(d*np.pi/180.0))+np.square(dd))
