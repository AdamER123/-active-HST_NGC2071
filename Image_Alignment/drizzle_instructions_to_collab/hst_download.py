#always start in the right env
#can run full script or work in command line, ipython, jupyter, etc
source activate astroconda 
ipython



'''
Make sure to update header as we go

Can check alignment with drizzle pac ipynbs [add some update note needs from research logbook]

For alignment examples, can see 
or https://github.com/spacetelescope/notebooks/tree/master/notebooks/DrizzlePac
or https://github.com/spacetelescope/drizzlepac

one idea is can try out https://github.com/spacetelescope/WFC3Library/blob/master/notebooks/calwf3_recalibration/calwf3_recal_tvb.ipynb see what happens
after this, run astrodrizzle (for further refinement, see https://github.com/spacetelescope/notebooks/blob/master/notebooks/DrizzlePac/align_sparse_fields/align_sparse_fields.ipynb)
will also need to apply a tweakreg with updatehdr = True then...
one last option if this does not improve things is to use https://github.com/spacetelescope/wfc3tools/blob/master/docs/notebooks/wf3persist.ipynb 
or https://github.com/spacetelescope/WFC3Library/blob/master/notebooks/persistence/wfc3_ir_persistence.ipynb (whatever is more recent)
https://github.com/spacetelescope/gaia_alignment/blob/master/Gaia_alignment.ipynb
'''


# All packages for entire work here
import glob
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

# Astropy packages we'll need
from astropy import units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.units import Quantity
from astropy.visualization import wcsaxes
from astropy.wcs import WCS

# Astroquery packages used for queries
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from astroquery.skyview import SkyView

# Drizzle related packages we'll need
from drizzlepac import astrodrizzle
from drizzlepac import tweakreg
from stsci.tools import teal
from stwcs import updatewcs

# Other handy parts
from ginga.util import zscale
from multiprocessing import Pool


# Some extra configuration 
SkyView.TIMEOUT = 15
mpl.rcParams['xtick.labelsize'] = 10
plt.rcParams.update({'axes.titlesize' : '18',
                     'axes.labelsize' : '14',
                     'xtick.labelsize' : '14',
                     'ytick.labelsize' : '14'})


# Option - Determine coordinates from data (not sure where this came from)
from matplotlib.patches import Polygon
import matplotlib.cm as cm



#querying and downloading observations
from astroquery.mast import Observations

# Query the ACS and WFC3 obervations HST program 11548 (HOPS 335, 361, F16*) and 16493
# I believe this will change what ht folders look like, but you can test that
# You can also specify filters with the parameter filters='...'
# see https://mast.stsci.edu/api/v0/_c_a_o_mfields.html for potential search criteria
obs_table_megeath = Observations.query_criteria(target_name='HOPS-335', obs_collection='HST', proposal_pi="Megeath*")
# Show the results of the observations query
print(obs_table_megeath)
# Repeat for other sets of observations
obs_table_karnath = Observations.query_criteria(target_name='HOPS361*', obs_collection='HST', proposal_pi="Karnath*")
print(obs_table_karnath)

# joining the two obs_Tables
obs_table_list = [obs_table_megeath, obs_table_karnath]

# This performs the data download
download_reduced_set = False #if you would like to filter the data...

# Create download directory
import os
if not os.path.exists('data_download/'):
    os.mkdir('data_download/')

for table in obs_table_list:
    # Get the information about the individual products (files) for each of the observations returned in the query
    products = Observations.get_product_list(table)
    # Show the results of the products query
    products

    # Filter the products to keep only the *flc.fits files
    filtered_products = Observations.filter_products(products,mrp_only=False,productSubGroupDescription='FLT')
    # Show the filtered products to ensure we kept the right products
    filtered_products

      
    if download_reduced_set:
        inds = [0, 1, 8, 9, 16, 17, 24, 25]
        filtered_products = filtered_products[inds]

    Observations.download_products(filtered_products,download_dir='data_download/',mrp_only=False)

path = './proper_motions_default/' #change this to copy from the download path to some new folder

#moving files to working directory, added option to make new one in code with the lines below
if not os.path.exists(path):
    os.mkdir(path)

# if not os.path.exists('./simult_tweakreg/ngc2071_megeath_11548/'):
#     os.mkdir('./simult_tweakreg/ngc2071_megeath_11548/')

# if not os.path.exists('./simult_tweakreg/ngc2071_karnath_16493/'):
#     os.mkdir('./simult_tweakreg/ngc2071_karnath_16493/')

import shutil
for f in glob.glob('data_download/mastDownload/HST/*11548*/*flt.fits'):
    filename = os.path.split(f)[-1]
    if not os.path.exists(filename):
        shutil.copy(f, path) #ngc2071_karnath_16493/')

import shutil
for f in glob.glob('data_download/mastDownload/HST/*16493*/*flt.fits'):
    filename = os.path.split(f)[-1]
    if not os.path.exists(filename):
        shutil.copy(f, path) #ngc2071_megeath_11548/')

#cd to working directory and then work on files further...

