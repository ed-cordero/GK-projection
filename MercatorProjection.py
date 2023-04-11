# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 19:02:17 2023

@author: corde
"""

import numpy as np

"""
Parameters model as described by IGAC (Agustin Codazzi Geographic 
Institute of Colombia)
"""
M_axis = 6378137                # Semi-major axis of the ellipsoid (radius at the equator) in meters
f_1 = 298.257222101             # Ellipsoid reciprocal flattening (IGAC -2004)
m_axis = (1-1/f_1)*M_axis       # Semi-minor axis length in meters
e1_2 = 0.0066943800229          # First eccentricity (squared)
e2_2 = 6.73949677547816E-03     # Second eccentricity (squared)
n = 1/(2*f_1-1)                 # Third flattening
Q_0 = M_axis/(1 + n)            # Multiplier in the meridian arc calculation [Deankin & Hunter, 2010]

del f_1

"""
# # Coefficients for meridian arc calculations
# coef_00 = np.array([1,   0 , 1/4  ,   0   ,   1/64 ,     0    ,    1/256 ,      0     ,    25/16384  ])
# coef_04 = np.array([0,   0 , 15/16,   0   ,  -15/64,     0    ,  -75/2048,      0     ,   -105/8192  ])
# coef_08 = np.array([0,   0 ,  0   ,   0   , 315/512,     0    , -441/2048,      0     ,   1323/32768 ])
# coef_12 = np.array([0,   0 ,  0   ,   0   ,    0   ,     0    , 1001/2048,      0     ,   -1573/8192 ])
# coef_16 = np.array([0,   0 ,  0   ,   0   ,    0   ,     0    ,     0    ,      0     , 109395/262144])

# coef_02 = np.array([0, -3/2,  0   ,  3/16 ,    0   ,    3/128 ,     0    ,    15/2048 ,       0      ])
# coef_06 = np.array([0,   0 ,  0   , -35/48,    0   ,  175/768 ,     0    ,   245/6144 ,       0      ])
# coef_10 = np.array([0,   0 ,  0   ,   0   ,    0   , -693/1280,     0    ,  2079/10240,       0      ])
# coef_14 = np.array([0,   0 ,  0   ,   0   ,    0   ,     0    ,     0    , -6435/14336,       0      ])

# coef_matrix = np.array([coef_00, coef_02, coef_04, coef_06, coef_08, coef_10, coef_12, coef_14, coef_16])
"""
n_powers_8 = np.array([n**i for i in range(9)])

# # Array of coefficients in the series for calculating the meridian arc.
# # Obtain by multiplying coef_matrix and n_powers in that order []
# coef_meridian = coef_matrix @ n_powers_8

# Coefficients using the rectified latitude for arc meridian length. IGAC uses that method for Gauss-Krüger projection
# Rectified radius (A in Deakin & Hunter, alpha in IGAC)
Q_1 = Q_0 * (16 + 3*n**2)/(16 - n**2)
# Coefficients to compute mu, meridian Arc (M = A*mu)
coef_00 = np.array([1,   0 ,  0   ,   0    ,    0   ,     0    ,     0    ,      0     ,       0      ])
coef_02 = np.array([0, -3/2,  0   ,  9/16  ,    0   ,    3/32  ,     0    ,      0     ,       0      ])
coef_04 = np.array([0,   0 , 15/16,   0    ,  -15/32,     0    ,     0    ,      0     ,       0      ])
coef_06 = np.array([0,   0 ,  0   , -35/48 ,    0   ,  105/256 ,     0    ,      0     ,       0      ])
coef_08 = np.array([0,   0 ,  0   ,    0   , 315/512,     0    ,     0    ,      0     ,       0      ])

coef_matrix = np.array([coef_00, coef_02, coef_04, coef_06, coef_08])
coef_m_rectified = coef_matrix @ n_powers_8

# Coefficients to compute phi_f
coef_02 = np.array([0, 3/2,  0   , -27/32 ,    0    ,  269/512,     0    ,      0     ,       0      ])
coef_04 = np.array([0,  0 , 21/16,   0    ,  -55/32 ,    0    ,     0    ,      0     ,       0      ])
coef_06 = np.array([0,  0 ,  0   , 151/96 ,    0    , -417/128,     0    ,      0     ,       0      ])
coef_08 = np.array([0,  0 ,  0   ,    0   , 1097/512,     0   ,     0    ,      0     ,       0      ])

coef_matrix = np.array([coef_00, coef_02, coef_04, coef_06, coef_08])
coef_phi_f = coef_matrix @ n_powers_8

del coef_00,coef_02, coef_04, coef_06, coef_08, coef_matrix, n_powers_8
# del coef_10, coef_12, coef_14, coef_16

# Coefficients of the polynomial in t and eta in the IGAC paper for calculating the North

p0 = np.array([     1    ,      0     ,    0   ,     0    ,     0   ,  0  ,    0    , 0, 0, 0])
p1 = np.array([   5/12   ,    -1/12   ,   9/12 ,     0    ,     0   , 4/12,    0    , 0, 0, 0])
p2 = np.array([  61/360  ,   -58/360  , 270/360,   1/360  , -330/360,  0  ,    0    , 0, 0, 0])
p3 = np.array([1385/20160, -3111/20160,    0   , 543/20160,     0   ,  0  , -1/20160, 0, 0, 0])

p_matrix = np.array([p0, p1, p2, p3])

# Coefficients of the polynomial in t and eta in the IGAC paper for calculating the East

p1 = np.array([   1/6    ,   -1/6   ,   1/6  ,     0    ,     0   ,  0  ,    0   , 0, 0, 0])
p2 = np.array([   5/120  ,  -18/120 ,  14/120,   1/120  , -58/120 ,  0  ,    0   , 0, 0, 0])
p3 = np.array([  61/5040 , -479/5040,    0   , 179/5040 ,     0   ,  0  , -1/5040, 0, 0, 0])

q_matrix = np.array([p0, p1, p2, p3])

# Coefficients for latitude (IGAC)

p0 = np.array([   -1   ,     0    ,    -1  ,     0    ,  0 ,   0 ,    0   ,  0 ,   0 , 0])
p1 = np.array([   5/6  ,    1/2   ,     1  ,     0    , -1 , -1/2,    0   ,  0 , -3/2, 0])
p2 = np.array([ -61/90 ,     -1   , -170/90,   -1/2   , 9/5,   0 ,    0   , 1/2,   0 , 0])
p3 = np.array([ 277/502, 3633/2510,     0  , 2048/1255,  0 ,   0 , 315/502,  0 ,   0 , 0])

g_matrix = np.array([p0, p1, p2, p3])

# Coefficients for longitude (IGAC)

p0 = np.array([    1    ,     0    ,   0  ,    0   ,  0  , 0,   0  , 0, 0, 0])
p1 = np.array([  -1/6   ,   -1/3   , -1/6 ,    0   ,  0  , 0,   0  , 0, 0, 0])
p2 = np.array([   1/24  ,    7/30  ,  1/20,   1/5  , 1/15, 0,   0  , 0, 0, 0])
p3 = np.array([ -61/5040, -331/2520,   0  , -61/252,  0  , 0, -1/7 , 0, 0, 0])

h_matrix = np.array([p0, p1, p2, p3])

del p0, p1, p2, p3

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def Geo_f_to_Dec(coord):
    ''' 
     Geo_f_to_Dec: Converts a geographical position (latitude or longitude) in sexagesimal format into decimal
     coord: the coordinates (either latitude/longitude) in sexagesimal format -- string
     Conventions: string format "dd°mm'ss.sss"H"
                  H is hemisphere either: N(north), E (east), W (west), S (south)
                  mm can be a float number, in which case ss.sss must be left out.
    return: The corresponding value in decimal. Negative for west/south hemispheres
    '''
    point = coord.strip()
    
    i = point.find("°")
    j = point.find("\'")
    k = point.find("\"")
    
    dec = float(point[:i])+float(point[i+1:j])/60
    
    if k != -1:
        dec += float(point[j+1:k])/3600
      
    if point[-1] =="W" or point[-1] =="S":
        dec = -dec
    return dec

class Datum:
    # Datum defines the parameters for the different Datums in Colombia
    
    def __init__(self, datum_select = "Bogotá - MAGNA"):
        
        # False Northing and Easting according to IGAC for the Bogotá-MAGNA (in meters)
        self.N_0 = 1000000
        self.E_0 = 1000000
        self.k_0 = 0.9996
        self.N_1 = 491767.5344
        
        # Ellipsoidal coordinates according to Datum selected
        self.Lat_0 = Geo_f_to_Dec("4°35\'46.3215\"N")
        self.Long_0 = Geo_f_to_Dec("74°04\'39.0285\"W")
        if datum_select == "Bogotá - MAGNA":
            return
        elif datum_select == "Este Cetral - MAGNA":
            self.Long_0 = self.lambda_0 + 3
        elif datum_select == "Este Este - MAGNA":
            self.Long_0 = self.lambda_0 + 6
        elif datum_select == "Oeste - MAGNA":
            self.Long_0 = self.lambda_0 - 3
        elif datum_select == "Oeste Oeste - MAGNA":
            self.Long_0 = self.lambda_0 - 6
        else:
            raise Exception("Datum not defined")

class Coordinates:
    
    def __init__(self, coor_ellip = [], coor_gauss_kruger = [], datum = None):
        """ 
        coor_ellip and coor_gaus_kruger must be 2d-lists
        coor_ellip in the format (latitiude, longitude)
        coor_gauss_kruger in the format (north, east)
        datum a string according to class Datum
        """
        self.Lat = None
        self.Long = None
        self.North = None
        self.East = None
        
        if coor_ellip != []:
            self.Lat = coor_ellip[0]
            self.Long = coor_ellip[1]
        if coor_gauss_kruger != []:
            self.North = coor_gauss_kruger[0]
            self.East = coor_gauss_kruger[1]
        if datum == None:
            self.datum = Datum()
        else:
            self.datum = Datum(datum)
        self.meridian_datum = self.meridian_arc(datum = True)
    
    def set_coor_ellip(self, coor_ellip):
        self.Lat = coor_ellip[0]
        self.Long = coor_ellip[1]
    
    def meridian_arc(self, datum = False):
        """
        This function calculates the length of the meridian arc using Helmert's method 
        according to (Deakin & Hunter, 2010)
        """
        if self.datum != None and datum:
            phi = np.radians(self.datum.Lat_0)
            
        elif self.Lat != None:
            phi = np.radians(self.Lat)
        
        # Calculates de meridian arc using Helmert's formula
        phi_array = np.array([phi if i==0 else np.sin(2*i*phi) for i in range(5)])
        meridian_length = Q_1*(coef_m_rectified @ phi_array)
        
        return meridian_length
    
    def poly_xy(self, x=0, y=0, grad=0):
        """
        
        Parameters
        ----------
        x : TYPE, optional
            DESCRIPTION. The default is 0.
        y : TYPE, optional
            DESCRIPTION. The default is 0.
        grad : TYPE, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        poly : TYPE
            DESCRIPTION.

        """
        poly = []

        for k in range(0,grad+1):
            for i in range(k+1):
                poly.append((y)**i * (x**2)**(k-i))
                # if k-i == 3:
                #     break
        return poly
    
    def to_Gauss_Kruger(self):
        
        if self.Lat == None:
            raise Exception("No geodetic coordinates defined")
            return
        # Parameters for finding the North coordinate:
        l = np.radians(self.Long-self.datum.Long_0)
           
        phi = np.radians(self.Lat)
        t = np.tan(phi)
        eta_2 = e2_2*np.cos(phi)**2                     # eta squared
        
        # Curvature radius
        M = M_axis/np.sqrt(1-e1_2*np.sin(phi)**2)
    
        A = (t*M*(l*np.cos(phi))**2)/2
        
        # Organise expansion terms in an array for easy computation
        cos_powers = np.array([(l*np.cos(phi))**i for i in range(0, 7,2)])
        poly_t_eta = np.array(self.poly_xy(t, eta_2, 3))
        
        # Coefficients of terms in formula
        term_coef = np.matmul(p_matrix, poly_t_eta)
        
        #Caluculate North
        self.North = (self.meridian_arc()-self.meridian_datum)+\
            A*np.dot(term_coef, cos_powers)+self.datum.N_0
        
        # Similar procedure for calculationg East
        B = np.sqrt(2*M*A/t)
        term_coef = np.matmul(q_matrix, poly_t_eta)
        
        self.East = B*np.dot(term_coef, cos_powers)+self.datum.E_0
        
    def to_geodetic(self):
        """
        This function calculates the latitude and longitude from the North and
        east coordinate in the Gauss-Krüger projection according to IGAC docs.
        The default datum is the MAGNA-SIRGAS Bogotá
        """
        if self.North == None:
            raise Exception("No plane coordinates defined")
            return
        
        # Parameters for finding the latitude coordinate:
        N = self.North - self.datum.N_1
        E = self.East - self.datum.E_0
        
        # Computing the reference point latitude (IGAC)
        theta = N / Q_1
        theta_hat = np.array([theta if i == 0 else np.sin(2*i*theta) for i in range (5)])
        phi_f = coef_phi_f @ theta_hat
        
        # Parameters which depend on phi_f
        t_f = np.tan(phi_f)
        eta_f_2 = e2_2*np.cos(phi_f)**2                 # eta_f squared
        
        # Curvature radius
        M = M_axis/np.sqrt(1-e1_2*np.sin(phi_f)**2)                 # N_f in IGAC paper
        
        # Computing latitude
        x = 0.5*(E/M)**2                                            # Our variable for the series
        x_powers = np.array([x**(i+1) for i in range(4)])           # Powers of x
        poly_t_eta = np.array(self.poly_xy(t_f, eta_f_2, 3))        # Polynomial of powers of eta_f and t_f
        terms = g_matrix @ poly_t_eta                               # g polynomials of eta_f and t_f
        self.Lat = np.degrees(phi_f + t_f*np.dot(terms, x_powers))

        # Computing longitude
        y = E/M
        y_powers = np.array([y**(i*2) for i in range(4)])           # Powers of y
        terms = h_matrix @ poly_t_eta                               # h polynomials of eta_f and t_f
        self.Long = self.datum.Long_0 + \
                    np.degrees(y/np.cos(phi_f)*np.dot(terms, y_powers))
        
# Test code:
    
t = Coordinates()
t.Lat = 7.325036667
t.Long = -72.48744667
t.to_Gauss_Kruger()
coor_plane = [t.North, t.East]

print("Gauss-Krüger coordinates:\n{} N, {} E".format(t.North, t.East))

u = Coordinates(coor_gauss_kruger = coor_plane)
u.to_geodetic()
print("\nGeodetic coordinates:\nLat {}, Long {}".format(u.Lat, u.Long))


            
        
            
            
        
    
        
    
        