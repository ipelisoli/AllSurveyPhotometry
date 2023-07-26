#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 15 23:40:25 2022

@author: james
"""

# general libraries
import os, sys, sqlite3, matplotlib
sys.path.insert(0, os.getcwd()+"/getScripts")
import matplotlib.pyplot as plt
#from termcolor import colored
import numpy as np
from pathlib import Path
from astropy.table import Table
from multiprocessing import Process
from datetime import datetime

# classes used (mine)
from getK2 import K2
from getTESSnew import TESS
from getZTF import ZTF
from getATLASforcedPHOT import getATLASforcedPHOT
from miscAstro import * # analysis:ignore
from get_SDSS import SDSSclass
from get_photometric_SED_CDS import photometricSED
from OverplotGaiaLocation import OverplotGaia
from getPTF import PTF
from getPanstarrs import getPanstarrs
from plotEverything import plotEverything
from getASASSN import ASASSN
from getWISE import WISE
from getNEOWISE import NEOWISE
from getCatalinaData import getCatalinaData
from getCDS import CDS
from checkLocalStars import checkLocalStars
from getpwds import getpwds
from getASASSNweb import getASASSNweb
from distutils.dir_util import copy_tree
from getGaiaDatalink import GaiaDatalink
import yaml  # Read from yml configuration file.

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rc('font', size=14)
#plt.xticks(fontsize=12)
#plt.yticks(fontsize=12)

import warnings
warnings.filterwarnings("ignore")

#address warnings:
#import warnings
#warnings.filterwarnings('error')



def Final_SDSS(RAdeg, Decdeg, rad=1.):
    # get_SDSS
    rad=rad/60. #arcmin
    scale=0.1
    try:
        query=SDSSclass.get_SDSSquery(RAdeg,Decdeg,rad)
        photSDSS = SDSSclass.search_SDSS_phot(query)

        if photSDSS:
            spec=SDSSclass.search_SDSS_spectrum(query,RAdeg,Decdeg,rad)
            if spec is not None:
                SDSSclass.plot_SDSS_spec("spec",spec)
                fspec=True
            else:
                fspec=False

            uv,g,r,i,z=SDSSclass.get_SDSSmagsUGRIZ("ugriz",RAdeg,Decdeg)

            url=SDSSclass.get_findingchart_url(RAdeg,Decdeg,scale)
            SDSSclass.createShortcut(url, str(cwd)+"/")
            SDSSclass.scrapeImageFromShortcut(url)
            return uv,g,r,i,z,fspec
        else:
            print("No SDSS data found.\n")
            return None
    except:
        print("Warning! SDSS query failed.\n")


def Final_phot_SED_CDS(RADec):
    # get_photometric_SED_CDS
    search_radius_SED = 3 # arcsec
    try:
        url=photometricSED.get_url_CDS(RADec, rad=search_radius_SED)
        photometricSED.plot_SED(url)
    except: None


def Final_K2(RADec, radius):
    # getK2
    #split=RADec.split(":")
    #split2=split[2].split(" ")
    # kepler coverage: RA=19h22m40s and Dec=+44˚30'00"(J2000)

    #if int(split[0]) >=18 and int(split[0]) <=20 and int(split2[1]) >= 40 and int(split2[1])<=50: # this is where K2 imaged
    # some issue here, got no data for something I know is in K2
    for time in exptimes:
        K2.get_K2(RADec, exptime=time, radius=radius, ignore_any_dodgyness=False) # radius in arcsec
    #else:
    #    print("Not in the K2 footprint.\n")



def Final_TESS(RADec, radius):
    ########### getTESSnew  NOTE: made a limiting mag for tess of 18
    for time in exptimes:
        print("- %s cadence data:" %time)
        TESS.get_tess(RADec, time=time, radius=radius, ignore_any_dodgyness=False) # ignore_any_dodgyness= True # this command will not process any lightcurve that includes nans at any point
    print("\n")

def Final_ZTF(RAdeg, Decdeg, RA, Dec, radius):
    #getPTF
    radZTF = radius/3600.   #  arcseconds in degrees
    if Decdeg > -34:
        if wantPTF: # this is globally set at the start
            print("#Checking PTF...")
            PTF.queryPTF(RA,Dec,str(radius))
            nPTF = PTF.splitRandG_PTF()
            if nPTF:
                print("Found %d measurement(s).\n"%nPTF)
            else:
                print("No data found.\n")
        # getZTF
        if wantZTF: # this is globally set at the start
            print("#Checking ZTF...")
            urlZTF=ZTF.create_url(None, [RAdeg,Decdeg,radZTF], BAD_CATFLAGS_MASK=True) # radius in deg
            try:
                nZTF = ZTF.save_data(RADec,     ZTF.get_data(urlZTF, (getpwds.ZTF()[0], getpwds.ZTF()[1])))
                if nZTF > 0:
                    print("Found %d measurement(s).\n"%nZTF)
                else:
                    print("No data found.\n")
            except:
                print("No data found.\n")
    else:
        print("#Checking PTF/ZTF...")
        print("Outside of PTF/ZTF footprint!\n")



def Final_ATLAS_forced(RAdeg, Decdeg, RA, Dec, reference_epoch, pmra, pmdec): #come back to this and also don't recompute stuff I already have
    # get ATLAS forced photometry
    # you can change how you want the photometry to be computed - difference imaging or placing an aperture. check the online docs and the script that handles ATLAS things
    # https://fallingstar-data.com/forcedphot/apiguide/
    try:
        if "data.txt" in os.listdir(os.getcwd()) and not "dataATLASalreadyprocessed.txt" in os.listdir(os.getcwd()): #if I manually had to get the file from ATLAS because of a break
            a = dfresult = pd.read_csv("data.txt", delim_whitespace=True)
            if minimumMJD==50000:
                getATLASforcedPHOT.plot_and_save_data(a,"FirstTime")
            else:
                getATLASforcedPHOT.plot_and_save_data(a,"AddToOriginal")
            np.savetxt("dataATLASalreadyprocessed.txt", np.array(["nope"]))

        else:
            if "ATLAS_c.dat" in os.listdir(os.getcwd()) or "ATLAS_o.dat" in os.listdir(os.getcwd()):
                try:
                    MJD_c = np.loadtxt("ATLAS_c.dat", unpack=True, usecols=(0))
                    maxMJD_c=np.amax(MJD_c)+0.1 # +0.1 just to ignore the first measurement
                except: None

                try:
                    MJD_o = np.loadtxt("ATLAS_o", unpack=True, usecols=(0))
                    maxMJD_o=np.amax(MJD_o)+0.1 # +0.1 just to ignore the first measurement
                except: None

                # first get the largest MJD value between both files if they exist
                try: minimumMJD=np.amax(np.array([maxMJD_c,maxMJD_o]))
                except:
                    try: #otherwise get the largest of the _c file
                        minimumMJD=maxMJD_c*-1
                    except: #otherwise get the largest of the _o file
                        minimumMJD=maxMJD_o*-1
            else: minimumMJD=50000 # else there is no prior entry and we ask for all the data

            dt = datetime.today().strftime('%Y-%m-%d')
            today=(list(sqlite3.connect(":memory:").execute("select julianday('" + dt + "')"))[0][0] -2400000.5)

            if abs(today-minimumMJD)>182.5:
                a=getATLASforcedPHOT.ATLAS(getpwds.ATLAS()[0], getpwds.ATLAS()[1], RAdeg, Decdeg,
                                           reference_epoch, pmra, pmdec, minMJD=minimumMJD)#minimumMJD)

                if minimumMJD==50000:
                    nATLAS = getATLASforcedPHOT.plot_and_save_data(a,"FirstTime")
                    print("Found %d measurements.\n" %nATLAS)
                else:
                    getATLASforcedPHOT.plot_and_save_data(a,"AddToOriginal")
                    print("Found %d measurements.\n" %nATLAS)

    except:
        with open("../../list_of/bad_atlas.txt", "a") as atlasfile:
            atlasfile.write(str(RA) + " " + str(Dec)+ "\n")



def Final_Catalina(RAdeg,Decdeg,ref_epoch,pmra,pmdec,radius):
    # get Catalina
    radius_catalina = radius/60  # passed as arcsec, so this is 3 arc seconds
    if not "Catalina.csv" in os.listdir(os.getcwd()): # this is globally set at the start
        try:
            getCatalinaData.getData(RAdeg,Decdeg,ref_epoch,pmra,pmdec,radius_catalina)
            #getCatalinaData.plot()
            nCatalina = getCatalinaData.handleData(RAdeg,Decdeg)
            print("Found %d measurement(s).\n"%nCatalina)
        except:
            print("No Catalina data found.\n")
    else:
        print("Catalina data already present in directory.\n")

def FinalGAIA(RADec, BP_RP, Abs_g, gmag):
    OverplotGaia.plotGaia(RADec, BP_RP, Abs_g, gmag)

def FinalPanstarrs(RAdeg,Decdeg,radius):
    radius_Panstarrs = radius/3600.
    getPanstarrs.getAllData(RAdeg,Decdeg, rad=radius_Panstarrs)
    nPS = getPanstarrs.getPanstarrsLCs(RAdeg,Decdeg)
    if nPS is not None:
        print("Found %d measurement(s).\n" %nPS)
    #getPanstarrs.getPanstarrsMeanMags()

def FinalASASSN(RAdeg, Decdeg, eDR3name="a"):
    radius_ASASSN = 5/3600
    if wantASASSN==True and not ("ASASSNv_lc.dat" in os.listdir(os.getcwd()) or "ASASSNg_lc.dat" in os.listdir(os.getcwd())):

        ASASSN_name = ASASSN.isItInASASSN(eDR3name, RAdeg, Decdeg, radius_ASASSN)
        filenameG, filenameV = ASASSN.getASASSN(ASASSN_name)

        try: ASASSN.plotVband(RAdeg,Decdeg,filenameV)
        except:
            try:
                try: json_link, lightcurve_link = getASASSNweb.get_urls(RAdeg, Decdeg)
                except: np.savetxt("ASSASNweb_empty.txt", np.array([1]))
                getASASSNweb.get_json_info(json_link, lightcurve_link)
                ASASSNweb.get_web_csv(lightcurve_link)
            except: None

        try: ASASSN.plotGband(RAdeg,Decdeg,filenameG)
        except: None # NOTE for the future : I don't know if the catalogue has any G only entries

        apphoturl=getASASSNweb.ASASSN_apphot_url(RAdeg,Decdeg)
        getASASSNweb.createShortcut(apphoturl, str(cwd)+"/")


def FinalWISE(RAdeg, Decdeg, gmag):
    if not "WISE.csv" in os.listdir(os.getcwd()):
        try:
            WISE.queryWise(RAdeg,Decdeg,ref_epoch=2020,pmra=0,pmdec=0)
            nWISE = WISE.saveFilters(RAdeg,Decdeg)
            if nWISE is not None:
                print("Found %d measurement(s).\n" %nWISE)
        except:
            print("Failed!\n")
    else:
        print("WISE data already present in directory.\n")

def FinalNEOWISE(RAdeg, Decdeg, radius):
    try:
        NEOWISE.queryWise(RAdeg, Decdeg, radius)
        NEOWISE.saveFilters(RAdeg, Decdeg)
    except Exception as e: print("error neowise"); print(e)

def FinalCDS(RAdeg, Decdeg):
    url=CDS.CDSurl(RAdeg,Decdeg)
    CDS.createShortcut(url, str(cwd)+"/")






if __name__ == '__main__':
    # don't touch these two
    joinTESS=True
    joinK2=True

    input_file = sys.argv[1]

    # What data do you want?
    with open('flags_photometry.yml') as f:
        phot_flags = yaml.safe_load(f)

    wantSDSS=phot_flags['SDSS']['download']
    radSDSS=float(phot_flags['SDSS']['radius'])

    wantK2=phot_flags['K2']['download']
    radK2=float(phot_flags['K2']['radius'])

    wantTESS=phot_flags['TESS']['download']
    radTESS=float(phot_flags['TESS']['radius'])
    # TESS exposure times to search for? not all are always available
    exptimes=phot_flags['TESS']['exptimes']
    for exptime in exptimes:
        if exptime not in ['fast', 'short', 'long']:
            print("TESS exptime not recognised; should be one or more of: fast, short, long.")
            sys.exit()

    wantZTF=phot_flags['ZTF']['download']
    radZTF=phot_flags['ZTF']['radius']

    # ATLAS might take a couple of minutes as a request is queued to their server
    wantATLASforced=phot_flags['ATLAS']['download']

    wantCatalina=phot_flags['Catalina']['download']
    radCatalina=phot_flags['Catalina']['radius']

    wantPTF=phot_flags['PTF']['download']

    wantPanstarrs=phot_flags['PanSTARRS']['download']
    radPanstarrs=phot_flags['PanSTARRS']['radius']

    wantWISE=phot_flags['WISE']['download']

    # neoWISE can be long to query
    wantNEOWISE=phot_flags['neoWISE']['download']
    radNEOWISE=phot_flags['neoWISE']['radius']

    wantASASSN=phot_flags['ASASSN']['download']

    # What extras are required?
    with open('flags_extras.yml') as f:
        extra_flags = yaml.safe_load(f)

    wantSED=extra_flags['SED']['plot']
    wantGaiaHR=extra_flags['GaiaHR']['plot']
    wantCDS=extra_flags['CDS']['create']
    wantGaiaDatalink=extra_flags['GaiaDatalink']['create']

    remove_old_dir=False  # are you sure? put twice to make sure you are positive. this will junk the contents of the created folders
    remove_old_dir=False # are you sure? put twice to make sure you are positive


    t = Table.read(input_file)


    NameOfNewDir = "Objects"
    try: os.mkdir(NameOfNewDir)
    except: None

    # this can be named anything, but this place has to be a separate folder to dump files into
    os.chdir(NameOfNewDir)
    homedirectory=os.getcwd()


    for count, (RAdeg, Decdeg) in enumerate(zip(t['ra'].value,t['dec'].value)):
        # Many of these things rely off of Gaia metrics. These are a) to handle proper motion in search queries b) plot on the Gaia HR diagram
        # You can set these to different values if you just care about getting the data... but the options are included as they are typical things to plot
        # Generic values could be:
        # Gaia_Gmag = XX  # you need this number! I apply saturation limits and minimum magnitudes for e.g. TESS, WISE
        # ref_epoch = 2020
        # propermRA = 0
        # propermDec = 0
        # probWD = 0
        # BPRP = 0
        # GaiaSourceID = YY # you need this if you want to get Gaia datalink data
        # GaiaABS_G = 0
        # Teff = 0

        Gaia_Gmag = t['phot_g_mean_mag'].value[count]   #   the Gaia magnitude in the Gband
        ref_epoch = t['ref_epoch'].value[count]   #   the Gaia reference epoch (e.g. 2016)
        propermRA = t['pmra'].value[count]     # the Gaia proper motion in RA
        propermDec = t['pmdec'].value[count]   # the Gaia proper motion in Dec
        #probWD = t['Pwd'].value[count]    # This is specific for white dwarfs and is just for the title of a graph. set this equal to 0 if you do not care about WDs. I use it as the probability of a white dwarf
        BPRP = t['bp_rp'].value[count]    # Gaia BP-RP
        GaiaSourceID = t['source_id'].value[count]    # Gaia source ID. Working as of DR3
        parallax = t['parallax'].value[count]
        GaiaABS_G = 5.+5.*np.log10(parallax/1000.)+Gaia_Gmag  # Absolute magnitude in Gaia G
        #Teff = t['teff_H'].value[count]     # Teff of the object



        RA,Dec = miscAstro.ra_dec_deg_to_hr(RAdeg,Decdeg)

        if Decdeg<0:
            RADec=str(RA) + " " + str(Dec)
            RAtemp = str(RA).replace(":", "")
            RAtemp = RAtemp.ljust(9, "0")[:9]
            DEtemp = str(Dec).replace(":", "")
            DEtemp = DEtemp.ljust(9, "0")[:9]
            jname = "J" + RAtemp + DEtemp
        else:
            RADec=str(RA) + " +" + str(Dec)
            RAtemp = str(RA).replace(":", "")
            RAtemp = RAtemp.ljust(9, "0")[:9]
            DEtemp = str(Dec).replace(":", "")
            DEtemp = DEtemp.ljust(8, "0")[:8]
            jname = "J" + RAtemp + "+" + DEtemp

        folder = os.getcwd()+"/"+str(jname)

        print("### WORKING ON OBJECT %s ###\n"%jname)

        print("Gaia ID: %s" %GaiaSourceID)
        print("Gaia G: %5.2f\n" %Gaia_Gmag)

        try:
            if remove_old_dir==True:     miscAstro.remDir(folder)
        except:None

        # make new dir and go to dir. if this breaks, dir already exists and we enter current dir
        try: os.mkdir(str(jname))
        except: None
        os.chdir(folder)

        #try: os.mkdir("files")
        #except: None

        cwd=os.getcwd()

        # bring back all files I neatened
        if remove_old_dir==False:
            for filename2 in os.listdir(cwd+"/"):
                Path(cwd+"/"+filename2).rename(cwd+"/"+filename2)


        #np.savetxt('TargetRADecDegrees.dat', np.array([RAdeg,Decdeg]).T)


        #print("We have ", str(len(t['ra'].value) - count), "out of ", str(len(t['ra'].value)), "remaining (progress = ", str(np.round(100*count/len(t['ra'].value),2)), "% )")
        #print(colored((RA,Dec),'cyan')); print(RAdeg, Decdeg)

        if wantSDSS:
            print("#Checking SDSS...")
            photSDSS=Final_SDSS(RAdeg,Decdeg,radSDSS)
            if photSDSS is not None:
                uv,g,r,i,z,spec = photSDSS
                print("u=%4.1f, g=%4.1f, r=%4.1f, i=%4.1f, z=%4.1f" %(uv,g,r,i,z))
                if spec:
                    print("Found at least one spectrum.\n")
                else:
                    print("No spectra found.\n")

        if wantNEOWISE:
            print("#Checking neoWISE...")
            p0 = Process(target = FinalNEOWISE(RAdeg, Decdeg, radNEOWISE))
            p0.start()

        if wantSED:
            p1 = Process(target = Final_phot_SED_CDS(RADec))
            p1.start()

        try:
            if wantTESS == True or wantK2 == True or wantATLASforced==True:
                obj, star_mag = checkLocalStars.find_star_in_gaia_edr3(RAdeg,Decdeg) #added in case there is a formatting mismatch, otherwise you could use Nicola's catalogue
        except: None

        if wantTESS:
            print("#Checking TESS...")
            returnClause = checkLocalStars.localTESS(obj,star_mag)
            if returnClause:
                p2 = Process(target = Final_TESS(RADec,radTESS))
                p2.start()
                joinTESS=True

        if wantK2:
            print("#Checking K2...")
            returnClause = checkLocalStars.localK2(obj,star_mag)
            if returnClause:
                p3 = Process(target = Final_K2(RADec,radK2))
                p3.start()
                joinK2=True

        if wantZTF or wantPTF:
            if Gaia_Gmag >12.5: # saturation limit
            # note that I put my own quality control cut on the ztf data by airmass and zeropoint rms
                p4 = Process(target = Final_ZTF(RAdeg,Decdeg, RA, Dec, radZTF))
                p4.start()
            else:
                print("Target too bright for ZTF/PTF.\n")

        if wantATLASforced:
            print("#Checking ATLAS...")
            try:
                if Gaia_Gmag >12.5:
                    if Decdeg>=-45: # saturation limit
                        returnClause = checkLocalStars.localATLAS(obj,star_mag)
                        if returnClause:
                            print("Querying ATLAS database, this might take a few minutes.")
                            p5 = Process(target = Final_ATLAS_forced(RAdeg,Decdeg,RA,Dec,reference_epoch=ref_epoch, pmra=propermRA, pmdec=propermDec))
                            #p5 = Process(target = Final_ATLAS_forced(RAdeg,Decdeg,RA,Dec,reference_epoch=2016, pmra=-146.303, pmdec=-155.864))
                            p5.start()
                    else:
                        print("Outside of ATLAS footprint!\n")
                else:
                    print("Target too bright for ATLAS.\n")
            except:
                print("Failed!\n")

        if wantCatalina:
            print("#Checking Catalina...")
            if Gaia_Gmag >= 13: # saturation limit
                p6 = Process(target = Final_Catalina(RAdeg,Decdeg, ref_epoch=ref_epoch, pmra=propermRA, pmdec=propermDec,radius=radCatalina)) # bit larger since not as good astrometric solution
                p6.start()
            else:
                print("Target too bright for Catalina.\n")

        if wantPanstarrs:
            print("Checking PanSTARRS...")
            if Gaia_Gmag >= 12: # saturation limit
                p7 = Process(target = FinalPanstarrs(RAdeg, Decdeg,radPanstarrs))
                p7.start()
            else:
                print("Target too bright for PanSTARRS.\n")

        #if Gaia_Gmag >= 11: # saturation limit
        #    p8 = Process(target = FinalASASSN(RAdeg, Decdeg, eDR3name="EDR3 "+str(GaiaSourceID)))
        #    p8.start()

        if wantWISE:
            print("Checking WISE...")
            if Gaia_Gmag<16.6:
                p9 = Process(target = FinalWISE(RAdeg, Decdeg, gmag=Gaia_Gmag))
                p9.start()
            else:
                print("Target too faint for WISE.\n")

        # plot gaia hr
        if wantGaiaHR:
            print("Making Gaia HR diagram...")
            try:
                p10 = Process(target = FinalGAIA(RADec, BPRP, GaiaABS_G, gmag=Gaia_Gmag))
                p10.start()
                print("Done!\n")
            except:
                print("Failed!\n")

        # get CDS shortcut
        if wantCDS:
            print("Creating CDS shortcut...")
            p11 = Process(target = FinalCDS(RAdeg,Decdeg))
            p11.start()
            print("Done!\n")

        if wantNEOWISE:
            p0.join()
        if wantSED:
            p1.join() # these make it so that the bunch terminates when the final process pX does
        if wantTESS and joinTESS:
            try: p2.join()
            except: None
        if wantK2 and joinK2:
            try: p3.join()
            except: None
        if (wantZTF or wantPTF) and (Gaia_Gmag >12.5):
            try: p4.join()
            except: None
        if wantATLASforced:
            try: p5.join()
            except: None
        if wantCatalina and Gaia_Gmag >= 13: # saturation limit
            try: p6.join()
            except: None
        if wantPanstarrs and Gaia_Gmag >= 12: # saturation limit
            try: p7.join()
            except: None
        #if Gaia_Gmag >= 11: # saturation limit
        #    try: p8.join()
        #    except: None
        if wantWISE and Gaia_Gmag < 16.6:
            p9.join()



        #pX = Process(target=plotEverything.plot());       pX.start()

        if wantGaiaDatalink:
            print("Checking for additional Gaia products...")
        #    try:
            source_id_input=str(GaiaSourceID)
            GaiaDatalink.getData(source_id_input)
        #    except:
#                print("None found!\n")




        # here is an example period search following the bare basics. you need to import the file you want.
        try:
            from astropy.timeseries import LombScargle
            from astropy import units as u
            MJD, mag, mage = np.loadtxt("ZTF_zg.csv", unpack=True)
            # remember that this is MJD!!! you should always convert your system to the barycentric reference frame.
            # see tdb Barycentric Dynamical Time (TDB)   under https://docs.astropy.org/en/stable/time/index.html
            # you can do it all with astropy routines
            frequency, power = LombScargle(MJD*u.d, mag, dy=mage).autopower()
            plt.plot(frequency, power)
            plt.xlabel("Frequency " +str(frequency.unit))
            plt.ylabel("Power")
            plt.savefig("ZTFperiodogram.png")
            plt.clf()


            # fold data at the highest peak... this may not be the true frequency of the system
            best_frequency = frequency[power==np.amax(power)]
            period=1/best_frequency.value[0]
            plt.errorbar((MJD%period)/period, mag, yerr=mage, fmt='.k')
            plt.savefig("ZTFphasefold.png")
            plt.clf()

            # and if you want to get into periodograms and period searching more seriously, PLEASE read these
            # https://ui.adsabs.harvard.edu/abs/2015ApJ...812...18V/abstract
            # https://iopscience.iop.org/article/10.3847/1538-4365/aab766     <- particularly this one. it's fantastic

            # lastly I recommend investigating the BLS search, typical for eclipsing systems/exoplanets

        except: None


        # neaten file list
        exceptions=["files", "GaiaLoc.png", "PhotSED.pdf", "SDSS.png", "ZTFphot.pdf",
                    "AllPhot.png", "CDS_objectofinterest.url"]

        for filename in os.listdir(cwd):
            if filename not in exceptions:
                Path(cwd+"/"+filename).rename(cwd+"/"+filename)



        # delete all memory of variables so I never get confused with next loop. variables defined locally in a function are never saved globally
        try: del respAll
        except: None
        try: del respZTF
        except: None

        del RA; del Dec
        del RAdeg; del Decdeg
        del folder; del RADec

        try: del uv; del g; del r; del i; del z
        except: None

        try: del p1
        except: None
        try: del p2
        except: None
        try: del p3
        except: None
        try: del p4
        except: None
        try: del p5
        except: None
        try: del p6
        except: None
        try: del p7
        except: None
        try: del p8
        except: None
        try: del p9
        except: None
        try: del p10
        except: None
        try: del p11
        except: None


        # back to the start
        os.chdir(homedirectory)



# observing biases:
# asassn - 30s
# goto - 60s, maybe changed in data out of comissioning
# gaia - 4s
# ztf -
# crts -
# asassn -

# todo:
# query cds references, e.g.: http://simbad.u-strasbg.fr/simbad/sim-id-refs?submit=sort+references&Ident=SDSS%20J053332.05%2B020911.5
