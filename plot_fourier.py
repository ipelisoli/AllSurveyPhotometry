import numpy as np

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['xtick.top']=True
plt.rcParams['ytick.right']=True

from astropy.timeseries import LombScargle
from astropy.stats import sigma_clip

from scipy import optimize

import sys

import yaml

def sort_arrays(ar1, ar2):
    '''Sort arrray1 according to array2'''
    ar1 = [ar1 for ar2, ar1 in sorted(zip(ar2, ar1))]
    return np.array(ar1)

def chi_sq(guess, x, y, err, factor):
    a, b = guess
    model = np.median(y)+a*np.sin(factor*2.0*np.pi*x + b)
    var = np.array(err)*np.array(err)
    chisq = np.sum((model - y) * (model - y)/var)
    return chisq

obj = sys.argv[1]

with open('flags_fourier.yml') as f:
    fourier_flags = yaml.safe_load(f)

fmin = fourier_flags['general']['fmin']
fmax = fourier_flags['general']['fmax']
period = fourier_flags['general']['period']
interac = fourier_flags['general']['interactive']
sigma = fourier_flags['general']['sigmaclip']

TESS = fourier_flags['TESS']['include']
ZTF = fourier_flags['ZTF']['include']
ATLAS = fourier_flags['ATLAS']['include']
Catalina = fourier_flags['Catalina']['include']
PTF = fourier_flags['PTF']['include']
PanSTARRS = fourier_flags['PanSTARRS']['include']
WISE = fourier_flags['WISE']['include']
neoWISE = fourier_flags['neoWISE']['include']
ASASSN = fourier_flags['ASASSN']['include']
Gaia = fourier_flags['Gaia']['include']

files = []

if neoWISE:
    files.append("NEOWISE_W1")
    files.append("NEOWISE_W2")

if TESS:
    files.append("TESS_long")

if PTF:
    files.append("PTF_r")
    files.append("PTF_g")

if ZTF:
    files.append("ZTF_zg")
    files.append("ZTF_zr")
    files.append("ZTF_zi")

if ATLAS:
    files.append("ATLAS_c")
    files.append("ATLAS_o")

if Catalina:
    files.append("Catalina")

if PanSTARRS:
    files.append("PANSTARRS_g")
    files.append("PANSTARRS_r")
    files.append("PANSTARRS_i")
    files.append("PANSTARRS_y")
    files.append("PANSTARRS_z")

if WISE:
    files.append("WISE_W1")
    files.append("WISE_W2")
    files.append("WISE_W3")
    files.append("WISE_W4")

#if ASASSN:
#    files.append()

if Gaia:
    files.append("Gaia_BP")
    files.append("Gaia_RP")
    files.append("Gaia_G")


for file in files:
    print("\nWorking on %s..." %file)
    try:
        mjd, m, me = np.loadtxt("Objects/%s/%s.dat"%(obj,file), unpack=True)
        try:
            nor = len(mjd)
        except TypeError:
            print("Only one detection in %s for %s, not enough for Fourier analysis.\n"%(file, obj))
            continue
        # Remove upper limits
        mjd = mjd[me > 0]
        m = m[me > 0]
        me = me[me > 0]
        n = len(mjd)
        if n < nor:
            print("%d measurements with null uncertainties (presumably upper limits) removed."%(nor-n))
        if len(mjd) < 30:
            print("Less than 30 detections in %s for %s, not enough for Fourier analysis.\n"%(file, obj))
            continue
    except OSError:
        print("File %s.dat not found for object %s." %(file, obj))
        continue

    # Sigma clipping
    if sigma>0:
        clip = sigma_clip(m, sigma=sigma, maxiters=3, masked=True)
        mask = ~(clip.mask)
        mjd = mjd[mask]
        m = m[mask]
        me = me[mask]
        print("%d measurement(s) removed with sigma-clipping."%(n-len(mjd)))

    ls = LombScargle(mjd, m, me)

    freq, power = ls.autopower(minimum_frequency=fmin, maximum_frequency=fmax,
                               samples_per_peak=5)

    if (period == 0):
        best_f = freq[np.argmax(power)]
        period = 1.0/best_f #period from the LS periodogram
        fap_p = ls.false_alarm_probability(power.max())

    fap_01 = ls.false_alarm_level(0.01)

    phase = np.mod((mjd-mjd[0])/period, 1)
    mphase = sort_arrays(m, phase)
    mephase = sort_arrays(me, phase)
    sphase = np.sort(phase)

    initial = [0.1, 0.]
    factor = 1.0
    solution = optimize.minimize(chi_sq, initial,
                                 args=(sphase, mphase, mephase, factor))

    phase_fit = np.arange(0,2,0.05)
    flux_fit = np.median(mphase) + solution.x[0] * np.sin(factor*2.*np.pi*phase_fit + solution.x[1])

    fig = plt.figure(figsize=(20,12))

    plt.subplot(221)
    plt.title("%s" %obj)
    plt.xlabel("Period [days]")
    plt.ylabel("Power")
    plt.xlim(1./fmax, 1./fmin)
    try:
        plt.ylim(0, np.max([1.1*np.max(power), 1.1*fap_01]))
    except:
        print("Could not set axis range.")
    plt.plot(1./freq, power, 'k')
    plt.axhline(fap_01, c='r', ls='--')

    plt.subplot(222)
    plt.title("%s" %file)
    plt.xlabel("MJD")
    plt.ylabel("Magnitude")
    plt.gca().invert_yaxis()
    plt.errorbar(mjd, m, me, fmt='k.')

    plt.subplot(212)
    plt.gca().invert_yaxis()
    plt.xlim(0, 2)
    plt.ylabel("Magnitude")
    plt.xlabel("Phase")
    plt.errorbar(phase, m, me, fmt='.k', alpha=0.5)
    plt.errorbar(1+phase, m, me, fmt='.k', alpha=0.5)
    plt.plot(phase_fit, flux_fit, 'r--', lw=5, alpha=0.5)
    plt.title("Period = %5.2f hours" %(24.*period))

    if interac:
        plt.tight_layout()
        plt.show()
    else:
        fig.savefig("%s_%s_ft.png"%(obj,file), bbox_inches='tight')

    plt.close()

    period = 0
