### Note about this fork
This was my attempt at making AllSurveyPhotometry a bit more user-friendly in terms of input/output. Whether I have been successful or not is for you to decide. Credits for the heavy-lifting part of the script go entirely to James Munday.

# All Survey Time-Series Photometry
Are you sick of people saying after talks "have you checked XXX for data"? Look no further! This package will obtain all time-series photometry associated with the object that is publicly available. All you need is the RA and Dec!

# What do you need to do?
- Download this repository and copy getScripts/getpwds_example.py to getScripts/getpwds.py. Include there your [ATLAS](https://fallingstar-data.com/forcedphot/) login details and your [IRSA](https://irsa.ipac.caltech.edu/Missions/ztf.html) account details.
- Open the flags_photometry.yml to define which photometry surveys you require by setting download to either True or False. When a radius can be set, that should be given in arcsec. For TESS, three cadences are available: fast, short, long. Only the latter seems to be working as expected.
- Open the flags_extras.yml to define which extras you required. SDSS will check for both photometry and spectroscopy; SED will download and plot SED from VizieR; GaiaHR will plot the HR diagram (white dwarf focused, requires parallax and bp_rp as input); CDS will create a link for the CDS object view; GaiaDataLink will check for additional data in Gaia (epoch photometry, BPRP spectra, RV spectra).
- Create a list with your target input parameters. Currently required are source_id, ra, dec, ref_epoch, pmra, pmdec (proper motion correction is applied for querying some surveys), phot_g_mean_mag. If wantGaiaHR = True, parallax and bp_rp are also needed. This list can be in any format readable by astropy.table (e.g. fits, csv). There is an example file that can be used for testing.
- Run the script! If you have all the required packages, it should work. The directory "Objects" should be created and you will see folders for each object. You can enter each folder to see what the outputs look like. Light curves are created as ".dat" files with a name corresponding to the survey+filter (or cadence in the case of TESS/K2). Other extensions are data taken from the survey "as is" (useful if you need to double check any parameters). There will also be a log file reporting e.g. how many measurements were found for each photometric survey.

# Current list of implemented surveys:  

| Survey      | Function  | Comments     |  Output Time Format (UTC, unless special)  |
| :---        |    :----   |    :----       | :----       |
| ATLAS        |    Forced photometry from the ATLAS survey     |    Only query if there are no close contaminants       | MJD |
| ASASSN        |    Find light curves in the full variable star catalogue     |    Requires download of the full catalogue (>40 GB), off by default  | BJD  |
| ASASSN Web        |    Autogenerate web search criteria for forced photometry   |           |
| Catalina/CRTS        |    All epoch photometry   | |    BJD       |
| CDS        |    Obtain photometric SED for any search radius. Clickable link to CDS for the object   |           |
| Gaia        |    All epoch photometry/spectra/RVS   |   | Gaia units... BJD but slightly different       |
| Kepler        |    Query the K2 field, extract lightcurves from tpf files    |       Only query if there are no close contaminants   (in progress, I recommend using your own scripts) | Units from Kepler |
| NEOWISE        |    Query all photometry with RA/Dec entries around the object   |          The processor to WISE |  BJD |
| Panstarrs        |    All epoch photometry from DR1   | |     BJD      |
| PTF        |    All epoch photometry   |          (iPTF to be included) |  MJD   |
| SDSS        |    Spectra, finding charts, mean star magnitudes   |           |
| TESS        |    Query TESS, extract lightcurves from tpf files   |        Only query if there are no close contaminants  (in progress, I recommend using your own scripts) | BJD |
| WISE        |    Query all photometry with RA/Dec entries around the object  | |    BJD       |
| ZTF        |    All epoch photometry  | |    MJD       |

### Extras:  
Gaia HR plotter  
checkLocalStars.py - performs a search for close stars using Gaia to not bother obtaining TESS/ATLAS/Kepler photometry. The search radius conditions reflect the pixel scale of the detector (read it)  

# A note on science usage
We recommend using this package to simply to check what data is available on various surveys perform preliminary analysis. We take no responsibility for supplying data for published work!

# Where might the queries break?
- Is your proper motion huge? Not all surveys have a proper motion option that can be accounted for. You should increase your search radius... but be prepared for some junk.
- Remember that ground based telescopes will have their declination limits - your source might not be available everywhere

## Required (non standard) packages:

| Package      | Tested with version  |
| :---        |    :----   |
| Python | 3.7/3.9 |
|[Multiprocess](https://pypi.org/project/multiprocess/)|  
|[Astropy](https://docs.astropy.org/en/stable/install.html)| 5.1 |   
|[Beautiful Soup](https://pypi.org/project/beautifulsoup4/)| 4.11.1 |  
|[PIL](https://pypi.org/project/Pillow/)| 9.2.0 |
|[Lightkurve](https://docs.lightkurve.org/about/install.html)| 2.2.0 |  
|[Astroquery](https://astroquery.readthedocs.io/en/latest/)|  0.4.6 |  
|[jdcal](https://pypi.org/project/jdcal/)| 1.4.1 |  

# Missing a source that you want included?
Push a script to access the data and/or supply basic details of how the data is obtained.

# Bugs/improvements?
If you find any, let us know! Email/make an issue.

# To-do list
Include ASAS  
Include APASS
Improve TESS/Kepler  
Neaten many bits of code
