# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 20:32:17 2023

@author: corde
"""

import csv
from MercatorProjection import Coordinates 

def convert():
    t = Coordinates()


datums = {'1' : "Bogotá - MAGNA", '2': "Este Cetral - MAGNA", '3': "Este Este - MAGNA", '4': "Oeste - MAGNA", '5':"Oeste Oeste - MAGNA"}
print("Loading program\n")
print("""
Welcome to Geodesic Project - Use this program to convert from geodetic data to 
Gauss-Krüger plane coordinates for Colombia.
Geodesic supports conversion from and to Geodetic coordinates (WGS84).
In order to convert a datum must be selected. All datums referred to the MAGNA-SIRGAS system.
The datum by default is Bogotá - MAGNA SIRGAS.
""") 

while True:
    f_opt= input("-Is your data on a CSV file? [y/n]:\n")
    f_opt = f_opt.lower()

    if f_opt not in ['y', 'yes', 'n', 'no']:
        print("{} is not a valid option".format(f_opt))

    elif f_opt in ['n', 'no']:
        datum = ''

        while datum not in datums.keys():
            datum = input(
"""
Select a datum to perform the conversion:
1. Bogotá
2. East - Central
3. East - East
4. West
5. West - West
"""
            )
            if datum not in datums.keys():
                print(f"{datum} is not a valid option")
        
    else:
        break
        
#with open()