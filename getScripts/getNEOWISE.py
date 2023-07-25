#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 29 17:46:57 2022

@author: james
"""


import numpy as np
import matplotlib.pyplot as plt
from astroquery.ipac.irsa import Irsa
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation
from miscAstro import * # analysis:ignore

#WISE 5σ photometric sensitivity is estimated to be 0.068, 0.098, 0.86 and 5.4 mJy
#(16.6, 15.6, 11.3, 8.0 Vega mag) at 3.4, 4.6, 12 and 22 μm in unconfused regions on the
#ecliptic plane


# PLEASE READ THIS if you wonder why I change the MJD https://wise2.ipac.caltech.edu/docs/release/neowise/expsup/sec1_2.html#exptime_offset
# new data release every spring https://irsa.ipac.caltech.edu/Missions/wise.html
class NEOWISE(object):
    def queryWise(RAdeg,Decdeg,radius):
        Irsa.ROW_LIMIT = 5000   # 5000 is the new value for row limit here.
        Irsa.TIMEOUT = 200
        t = Irsa.query_region(SkyCoord(RAdeg, Decdeg, unit=(u.deg,u.deg)), spatial="Cone",
                                                            radius=radius * u.arcsec, catalog="neowiser_p1bs_psd",
                                                            selcols="ra,dec,mjd,w1mpro,w1sigmpro,w2mpro,w2sigmpro") # these selcols are ignored only here?
        if (len(t) > 0):
            print("Found neoWISE data: %d measurement(s).\n" %len(t))
            t.write('NEOWISE.csv', overwrite=True)
        else:
            print("No neoWISE data found.\n")


    def saveFilters(RAdeg,Decdeg):
        try:
            table=Table.read('NEOWISE.csv')

            #IMPORTANT NOTE: WISE WAS A SPACECRAFT AND THE STORED TIME WAS MJD.
            # I DO NOT KNOW THE EXACT LOCATION FOR BJD

            # Write W1 data to separate file
            mask1=(table['w1mpro'] > 0)
            mjd=table['mjd'][mask1] + 0.57 # for WISE and neoWISE W1 and W2
            mjd[mjd < 57196.48882]+=1./86400. # correct leap second error
            bjd=miscAstro.jd_corr(mjd,RAdeg,Decdeg,loc=EarthLocation.of_site("lapalma")).value
            w1mag=table['w1mpro'][mask1]
            w1mage=table['w1sigmpro'][mask1]
            np.savetxt("NEOWISE_W1.csv", np.array([bjd, w1mag, w1mage]).T)

            # Write W2 data to separate file
            mask2=(table['w2mpro'] > 0)
            mjd=table['mjd'][mask2] + 0.57
            mjd[mjd < 57196.48882]+=1./86400.
            bjd=miscAstro.jd_corr(mjd,RAdeg,Decdeg,loc=EarthLocation.of_site("lapalma")).value
            w2mag=table['w2mpro'][mask2]
            w2mage=table['w2sigmpro'][mask2]
            np.savetxt("NEOWISE_W2.csv", np.array([bjd, w1mag, w1mage]).T)

        except FileNotFoundError:
            return

    def Plot():
        try:
            mjd,mag,mage=np.loadtxt("NEOWISE1_W1.csv",unpack=True)
            plt.scatter(mjd,mag, c='b', label="WISE1")
            plt.errorbar(mjd,mag,yerr=mage,ls=" ", c='b')
        except: None

        try:
            mjd,mag,mage=np.loadtxt("NEOWISE2_W2.csv",unpack=True)
            plt.scatter(mjd,mag, c='g', label="WISE2")
            plt.errorbar(mjd,mag,yerr=mage,ls=" ", c='g')
        except: None
