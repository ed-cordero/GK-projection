# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from math import pi
from math import cos
from math import sin
from math import sqrt
from math import radians

import pandas as pd

#Module to calcuale distances over earth surface using the Geodetic model

# Parameters model as described by IGAC (Agustin Codazzi Geographic 
# Institute of Colombia)

M_axis = 6378137                # Semi-major axis of the ellipsoid (radius at the equator) in meters
f_1 = 298.257222101             # Ellipsoid reciprocal flattening (IGAC -2004)
m_axis = (1-1/f_1)*M_axis       # Semi-minor axis length in meters
e2_1 = 0.0066943800229          # First eccentricity
e2_2 = 6.73949677547816E-03     # Second eccentricity
    
def radius(lat_1, lat_2):
    # This fucntion calculates the equivalent radius
    # for a projection in a given region.

    # Medium latitud in radians    
    lat = radians((lat_1 + lat_2)/2)
    
    # Computing the radius at that latitude.
    radius = (M_axis**2*cos(lat))**2 +(m_axis**2*sin(lat))**2
    radius = radius / ((M_axis*cos(lat))**2+(m_axis*sin(lat))**2)
    radius = sqrt(radius)
    
    return radius






    