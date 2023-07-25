#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 23:21:30 2022

@author: james
"""

import matplotlib.pyplot as plt
import numpy as np
import lightkurve as lk
#from astropy.timeseries import LombScargle

plt.ioff()

# https://docs.lightkurve.org/tutorials/2-creating-light-curves/2-1-combining-multiple-quarters.html
# https://heasarc.gsfc.nasa.gov/docs/tess/Aperture-Photometry-Tutorial.html


#transit_mask = corrected_lc.create_transit_mask(period=[3.2277, 7.3745],
                                                #duration=[0.25, 0.25],
                                                #transit_time=[2385.6635, 2389.9635])
#https://docs.lightkurve.org/tutorials/2-creating-light-curves/2-3-k2-pldcorrector.html?highlight=to_corrector


class TESS(object):
    def get_tess(RADec, time, radius=5, ignore_any_dodgyness=False):
        print("Searching for light curve...")
        # https://docs.lightkurve.org/reference/api/lightkurve.search_lightcurve.html
        search_result_lc=lk.search_lightcurve(RADec, radius=radius,author="TESS", exptime=time)

        if search_result_lc:
            lc_obj=search_result_lc.stitch().remove_outliers().remove_nans()   # the flux here is meant to automatically be PDCSAP
            lcfound=True

        else:
            print("No light curve found, searching for target pixel file...")
            search_result = lk.search_targetpixelfile(RADec, mission="TESS", exptime=time)
            if search_result:
                # inspect tpf: https://docs.lightkurve.org/tutorials/1-getting-started/interactively-inspecting-data.html
                print("Target pixel file found, performing photometry...")
                all_tpf = search_result.download_all("hard")
                all_lcs=[]

                for tpf in all_tpf:
                    bkg = tpf.get_bkg_lightcurve()
                    for z in range(3):
                        medianbkg=np.median(bkg.flux).value
                        sigmabkg=np.std(bkg.flux).value
                        mask=((bkg.flux.value > (medianbkg-3.5*sigmabkg)) & (bkg.flux.value < (medianbkg+2*sigmabkg)))
                        bkg=bkg[mask]
                        tpf=tpf[mask]

                    aper = tpf.create_threshold_mask()
                    raw_lc = tpf.to_lightcurve(aperture_mask=aper)
                    raw_lc = raw_lc.remove_nans()
                    try:
                        lc=raw_lc.to_corrector(method="sff").correct() #cbv #pld
                    except:
                        try:
                            lc=raw_lc.to_corrector(method="cbv").correct()
                        except:
                            None
                    try:
                        lc=lc.remove_outliers(sigma=5)
                        lc=lc.normalize()

                        all_lcs.append(lc)
                    except: None

                try:
                    print("Stitching light curve...")
                    lc_obj=lk.LightCurveCollection(all_lcs).stitch()
                    #ax=lc_obj.scatter(title=time)
                    #ax.figure.savefig('AllLC'+str(time)+'.png')
                    #plt.close()
                    lcfound=True
                    print("Success!")
                except:
                    print("Failed!")
                    lcfound=False
            else:
                print("No target pixel file found.")
                lcfound=False


            if lcfound:
                np.savetxt("TESS_"+str(time)+".csv", np.array([lc_obj.time + 2457000- 2400000.5, lc_obj.flux, lc_obj.flux_err]).T, fmt="%s")
                # time is output as BJD







#old code here to do periodograms with lightkurve
                #try: # this will only break if lc_obj does not exist
                    ##ts=lc_obj.time
                    ##fluxes=lc_obj.flux
                    ##fluxes_e=lc_obj.flux_err

                    ##times_to_sample=np.linspace(1/86400,1/250,1000)*u.Hz
                    ##ls = LombScargle(ts,fluxes,fluxes_e)
                    ##power=ls.power(times_to_sample)


                    ##plt.plot((1/times_to_sample)/86400, power)
                    ##plt.show()




                    #periodogram=lc_obj.to_periodogram("bls", minimum_period=250/86400, maximum_period=2)
                    #plt.clf()
                    #print("I GOT HERE")
                    #periodogram.plot().figure.savefig('Periodogram_'+str(time)+'.png')
                    #period=periodogram.period_at_max_power
                    #folded_lc=lc_obj.fold(period)
                    #folded_lc.scatter(title=str(period) + "   " + str(time)).savefig('folded_lc_'+str(time)+'.png')
                    #plt.close()
                    ##folded_lc.plot_river(method='sigma')
                #except Exception as exc:
                #    print(exc)





            # for kepler data, https://docs.lightkurve.org/tutorials/2-creating-light-curves/2-3-k2-sffcorrector.html


            # limitation: might struggle with crowded fields... you can use difference imaging with lightkurve https://docs.lightkurve.org/tutorials/3-science-examples/periodograms-verifying-the-location-of-a-signal.html




#TESS.get_tess("23:53:00.9 -38:51:47.67", "short")
