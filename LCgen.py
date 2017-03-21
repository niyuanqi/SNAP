#################################################################
# Name:     LCgen.py                                            #
# Author:   Yuan Qi Ni                                          #
# Date:     March, 21, 2017                                     #
# Function: Program uses MagCalc routine and Source Extractor   #
#           generated light curves to generate light curve file #
#           with both magnitudes for cross comparison.          #
#################################################################

#essential modules
from glob import glob
import numpy as np
import math
import os

from .Processing.MagCalc import *

print "hw"
