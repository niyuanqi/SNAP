#################################################################
# Name:     MakeReg.py                                          #
# Author:   Yuan Qi Ni                                          #
# Date:     July 7, 2016                                        #
# Function: Program makes a reg file from coordinates.          #
#################################################################

def padzero(string, num):
    if len(string)<num:
        return padzero("0"+string, num)
    else:
        return string

def makeReg(name, ra, dec):
    f = open(name, 'w')
    rahms = deg_toHMS(ra)
    deghms = deg_toDMS(dec)
    circle1 = "fk5; circle "+rahms+"  "+deghms+" 5\""
    circle2 = "fk5; circle "+rahms+"  "+deghms+" 30\""
    comment = "#  "+str(ra)+"  "+str(dec)
    f.write(circle1+'\n')
    f.write(circle2+'\n')
    f.write(comment)

#RA conversion
def deg_toHMS(deg):
    hour = 24.*abs(deg)/360.
    minute = (hour%1)*60
    hour = int(hour - hour%1)
    second = (minute%1)*60
    minute = int(minute - minute%1)
    second = round(second*1000)/1000.0
    frac = int(second%1*1000)
    second = int(second - second%1)
    
    hour = padzero(str(hour),2)
    minute = padzero(str(minute),2)
    second = padzero(str(second),2)
    frac = str(frac)

    HMS = hour+":"+minute+":"+second+"."+frac
    if deg<0:
        HMS = "-"+HMS
    return HMS

#DEC conversion
def deg_toDMS(deg):
    minute = (abs(deg)%1)*60
    degree = int(abs(deg) - abs(deg)%1)
    second = (minute%1)*60
    minute = int(minute - minute%1)
    second = round(second*10)/10.0
    frac = int(second%1*10)
    second = int(second - second%1)
    
    degree = str(degree)
    minute = padzero(str(minute),2)
    second = padzero(str(second),2)
    frac = str(frac)

    DMS = degree+":"+minute+":"+second+"."+frac
    if deg<0:
        DMS = "-"+DMS
    return DMS

if __name__=='__main__':

    import argparse
    
    parser = argparse.ArgumentParser(description="Create ds9 Reg File")
    parser.add_argument("name", type=str, help="name of reg file")
    parser.add_argument("position", type=str, help="RA:DEC as deg:deg")
    args = parser.parse_args()

    #extract RA, DEC from position argument
    RA, DEC = [float(coord) for coord in args.position.split(':')]

    #write reg file
    makeReg(args.name, RA, DEC) 
