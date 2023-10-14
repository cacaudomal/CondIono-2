#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 22:54:35 2020

@author: Ciaran Beggan (British Geological Survey)

A set of functions to parse the three input options of the IGRF code
         '1 - values at one or more locations & dates'
         '2 - values at yearly intervals at one location'
         '3 - values on a latitude/longitude grid at one date'

"""

import igrf_utils as iut
import numpy as np

degree_sign= u'\N{DEGREE SIGN}'

def option1(itype,LTD,LND,alt,date):
    '''
    Option 1 is the simplest: a single point and time
    
    Parameters
    ----------
    LTD : Latitude in decimal degrees
    LND : Longitude in decimal degrees
    alt : altitude in kilometers
    date : decimal date in years 1900-2025
    
    '''
 
    while 1:
        itype = iut.check_int(itype)
        if itype < 1 or itype > 2: continue
        else:
            break
        
    #print( '\nio_options_clara func_options1: Enter latitude & longitude in decimal degrees',LTD,LND)
    # LTD,LND = input('-> ').rstrip().split(' ') 
    latd = iut.check_float(LTD)
    lond = iut.check_float(LND)
        
    lat, lon = iut.check_lat_lon_bounds(latd,0,lond,0)
    colat = 90-lat
        
    while 1:
        if itype==1:
            #alt = input( 'Enter altitude in km: ').rstrip()
            alt = iut.check_float(alt)
            alt, colat, sd, cd = iut.gg_to_geo(alt, colat)
        else:
            #alt = input( 'Enter radial distance in km (>3485 km): ').rstrip()
            alt = iut.check_float(alt)
            sd = 0; cd = 0
        
        if (itype == 2) and (alt < 3485):
            print('Alt must be greater then CMB radius (3485 km)').rstrip()
            continue
        else:
            break
        
    while 1: 
        #date = input('Enter decimal date in years 1900-2025: ').rstrip()
        date = iut.check_float(date)
        if date < 1900 or date > 2030: continue
        else:
            break   

    return date, alt, lat, colat, lon, itype, sd, cd


def write1(name, date, alt, lat, colat, lon, X, Y, Z, dX, dY, dZ, \
                  dec, hoz, inc, eff, decs, hozs, incs, effs, itype):
    '''
    Write out a single lat/long/alt the main field and SV values to screen or a file
    '''
    if itype == 1:
         alt, lat = iut.geo_to_gg(alt, colat)
         lat = 90-lat
         
    print("entrou",name)
    
    if not name: # Print to screen
        print('\nGeomagnetic field values at: ', str(np.round(lat, decimals=4)) 
            + degree_sign  + ' / ' + str(lon) 
            + degree_sign + ', at altitude ' 
            + str(np.round(alt, decimals=3)) + ' for ' + str(date))
        print('Declination (D):', '{: .3f}'.format(dec), degree_sign)
        print('Inclination (I):', '{: .3f}'.format(inc), degree_sign)
        print('Horizontal intensity (H):', '{: .1f}'.format(hoz), 'nT')
        print('Total intensity (F)     :', '{: .1f}'.format(eff), 'nT')
        print('North component (X)     :', '{: .1f}'.format(X), 'nT')
        print('East component (Y)      :', '{: .1f}'.format(Y), 'nT')
        print('Vertical component (Z)  :', '{: .1f}'.format(Z), 'nT')
        print('Declination SV (D):', '{: .2f}'.format(decs), 'arcmin/yr')
        print('Inclination SV (I):', '{: .2f}'.format(incs), 'arcmin/yr')
        print('Horizontal SV (H):', '{: .1f}'.format(hozs), 'nT/yr')
        print('Total SV (F)     :', '{: .1f}'.format(effs), 'nT/yr')
        print('North SV (X)     :', '{: .1f}'.format(dX), 'nT/yr')
        print('East SV (Y)      :', '{: .1f}'.format(dY), 'nT/yr')
        print('Vertical SV (Z)  :', '{: .1f}'.format(dZ), 'nT/yr')
    
    else: # Print to filename 
        print("entrou file")
        with open(name, 'a') as file: 
            file.writelines(['Geomagnetic field values at: \nLat: ',  str(np.round(lat, decimals=4)) 
               + degree_sign +' \nLon: ' + str(lon) 
               + degree_sign + '\nAltitude: ' 
               + str(np.round(alt, decimals=3)) +'\nYear: ' + str(date) +'\n'])
            file.writelines(['Declination (D):', '{: 5.2f}'.format(dec), degree_sign + '\n'])
            file.writelines(['Inclination (I):', '{: 5.2f}'.format(inc), degree_sign + '\n'])
            file.writelines(['Horizontal intensity (H):', '{: 8.1f}'.format(hoz), 'nT\n'])
            file.writelines(['Total intensity (F)     :', '{: 8.1f}'.format(eff), 'nT\n'])
            file.writelines(['North component (X)     :', '{: 8.1f}'.format(X), 'nT\n'])
            file.writelines(['East component (Y)      :', '{: 8.1f}'.format(Y), 'nT\n'])
            file.writelines(['Vertical component (Z)  :', '{: 8.1f}'.format(Z), 'nT\n'])
            file.writelines(['Declination SV (D) :', '{: 5.2f}'.format(decs), 'arcmin/yr\n'])
            file.writelines(['Inclination SV (I) :', '{: 5.2f}'.format(incs),  'arcmin/yr\n'])
            file.writelines(['Horizontal SV (H)  :', '{: 7.1f}'.format(hozs), 'nT/yr\n'])
            file.writelines(['Total SV (F)       :', '{: 7.1f}'.format(effs), 'nT/yr\n'])
            file.writelines(['North SV (X)       :', '{: 7.1f}'.format(dX), 'nT/yr\n'])
            file.writelines(['East SV (Y)        :', '{: 7.1f}'.format(dY), 'nT/yr\n'])
            file.writelines(['Vertical SV (Z)    :', '{: 7.1f}'.format(dZ), 'nT/yr\n\n'])
            
def write(name, date, alt, lat, colat, lon, X, Y, Z, dX, dY, dZ, \
                  dec, hoz, inc, eff, decs, hozs, incs, effs, itype):
    '''
    Write out a single lat/long/alt the main field and SV values to screen or a file
    '''
    if itype == 1:
         alt, lat = iut.geo_to_gg(alt, colat)
         lat = 90-lat
         
    #("entrou",name)
    
    if not name: # Print to screen
        print('\nGeomagnetic field values at: ', str(np.round(lat, decimals=4)) 
            + degree_sign  + ' / ' + str(lon) 
            + degree_sign + ', at altitude ' 
            + str(np.round(alt, decimals=3)) + ' for ' + str(date))
        print('Declination (D):', '{: .3f}'.format(dec), degree_sign)
        print('Inclination (I):', '{: .3f}'.format(inc), degree_sign)
        print('Horizontal intensity (H):', '{: .1f}'.format(hoz), 'nT')
        print('Total intensity (F)     :', '{: .1f}'.format(eff), 'nT')
        print('North component (X)     :', '{: .1f}'.format(X), 'nT')
        print('East component (Y)      :', '{: .1f}'.format(Y), 'nT')
        print('Vertical component (Z)  :', '{: .1f}'.format(Z), 'nT')
        print('Declination SV (D):', '{: .2f}'.format(decs), 'arcmin/yr')
        print('Inclination SV (I):', '{: .2f}'.format(incs), 'arcmin/yr')
        print('Horizontal SV (H):', '{: .1f}'.format(hozs), 'nT/yr')
        print('Total SV (F)     :', '{: .1f}'.format(effs), 'nT/yr')
        print('North SV (X)     :', '{: .1f}'.format(dX), 'nT/yr')
        print('East SV (Y)      :', '{: .1f}'.format(dY), 'nT/yr')
        print('Vertical SV (Z)  :', '{: .1f}'.format(dZ), 'nT/yr')
    
    else: # Print to filename 
        #print("entrou file")
        with open(name, 'a') as file: 
           
            file.writelines([str(np.round(lat, decimals=4)),",", str(lon),",",str(np.round(alt, decimals=3)),",",str(date),",",'{: 5.2f}'.format(dec),\
                             ",",'{: 5.2f}'.format(inc),",",'{: 8.1f}'.format(hoz),",", '{: 8.1f}'.format(eff),",", '{: 8.1f}'.format(X), \
                             ",",'{: 8.1f}'.format(Y),",",'{: 8.1f}'.format(Z), ",",'{: 5.2f}'.format(decs),",",'{: 5.2f}'.format(incs), \
                             ",",'{: 7.1f}'.format(hozs),",",'{: 7.1f}'.format(effs),",",'{: 7.1f}'.format(dX),",",'{: 7.1f}'.format(dY),\
                             ",",'{: 7.1f}'.format(dZ), '\n'])
            