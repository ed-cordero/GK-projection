# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 20:32:17 2023

@author: corde
"""

import csv
from MercatorProjection import Coordinates 

print("Loading program\n")
print("Welcome to Geodesic Project - Use this program to convert from geodetic data to\n\
      Gauss-Krüger plane coordinates for Colombia") 
while True:
    f_opt= input("-Is your data on a CSV file? [y/n]:\n")
    f_opt = f_opt.lower()
    if f_opt in ['y', 'yes', 'n', 'no']:
        break
    else:
        print("{} is not a valid option".format(f_opt))
#with open()