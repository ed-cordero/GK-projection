# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 20:32:17 2023

@author: corde
"""

import csv
from MercatorProjection import Coordinates 

# Global variables
datums = {'1' : "Bogotá - MAGNA", '2': "Este Cetral - MAGNA", '3': "Este Este - MAGNA", '4': "Oeste - MAGNA", '5':"Oeste Oeste - MAGNA"}

def convert(conv_type, coord, datum_selected):
    if conv_type == 1:
        t = Coordinates(coor_ellip=coord, datum=datum_selected)
        t.to_Gauss_Kruger()
        return t.North, t.East
    else:
        t = Coordinates(coor_gauss_kruger=coord, datum=datum_selected)
        t.to_geodetic()
        return t.Lat, t.Long    

def get_datum():
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
             return datum

def get_type_convesion():
     while True:
        conv_type = input("""
Select type of conversion:
1. From Geodetic data
2. From Plane Coordinates
""")
        if conv_type in ['1', '2']:
            return conv_type
        else:
            print("{} is not a valid option".format(conv_type))

def get_coord ():
    coord = input("-Please enter the coordinates separated by a comma (N,E) -- Just numbers:")
    return list(map(int, coord.strip().split(',')))

def main(): 
    print("Loading program\n")
    print("""
    Welcome to Geodesic Project - Use this program to convert from geodetic data to 
    Gauss-Krüger plane coordinates for Colombia.
    Geodesic supports conversion from and to Geodetic coordinates (WGS84).
    In order to convert a datum must be selected. All datums referred to the MAGNA-SIRGAS system.
    The datum by default is Bogotá - MAGNA SIRGAS.
    """)
    
    datum = get_datum()
    conv_type = get_type_convesion()

    while True:
        f_opt= input("-Is your data on a CSV file? [y/n]:\n")
        f_opt = f_opt.lower()

        if f_opt not in ['y', 'yes', 'n', 'no']:
            print("{} is not a valid option".format(f_opt))

        elif f_opt in ['n', 'no']:
            coord = get_coord()
            converted = convert(conv_type, coord, datums[datum])
            
        else:
            file = input("-please enter the name of the file: ")
            with open(file, 'r') as coor_file:
                coord_read = csv.DictReader(coor_file, fieldnames=['latitude', 'longitude'])

main()