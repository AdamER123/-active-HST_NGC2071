#always start in the right env
source activate astroconda 
ipython



'''
Make sure to update header as we go

Can check alignment with drizzle pac ipynbs [add some update note needs from research logbook]

For alignment, try https://github.com/spacetelescope/gaia_alignment/blob/master/Gaia_alignment.ipynb
or https://github.com/spacetelescope/notebooks/tree/master/notebooks/DrizzlePac
or https://github.com/spacetelescope/drizzlepac

one idea is can try out https://github.com/spacetelescope/WFC3Library/blob/master/notebooks/calwf3_recalibration/calwf3_recal_tvb.ipynb see what happens
after this, run astrodrizzle (for further refinement, see https://github.com/spacetelescope/notebooks/blob/master/notebooks/DrizzlePac/align_sparse_fields/align_sparse_fields.ipynb)
will also need to apply a tweakreg with updatehdr = True then...
one last option if this does not improve things is to use https://github.com/spacetelescope/wfc3tools/blob/master/docs/notebooks/wf3persist.ipynb 
or https://github.com/spacetelescope/WFC3Library/blob/master/notebooks/persistence/wfc3_ir_persistence.ipynb (whatever is more recent)

Use https://docs.astropy.org/en/stable/coordinates/velocities.html
https://allendowney.github.io/AstronomicalData/03_motion.html
Or cross correlation for proper motions from https://iopscience.iop.org/article/10.1086/339837/pdf
'''

# All packages for entire notebook here, but also imported in first relevant cell


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
from astropy.table import vstack
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


# Option 3- Determine coordinates from data
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import wcsaxes

from matplotlib.patches import Polygon
import matplotlib.cm as cm



#querying and downloading observations
from astroquery.mast import Observations

# Query the ACS and WFC3 obervations HST program 11548 (HOPS 335, 361, F16*)
# see https://mast.stsci.edu/api/v0/_c_a_o_mfields.html for potential search criteria
obs_table_megeath = Observations.query_criteria(target_name='HOPS-335', obs_collection='HST', proposal_pi="Megeath*")
# Show the results of the observations query
obs_table_megeath
# Repeat for other sets of observations
obs_table_karnath = Observations.query_criteria(target_name='HOPS361*', obs_collection='HST', proposal_pi="Karnath*")
obs_table_karnath

# joining the two obs_Tables
obs_table_list = [obs_table_megeath, obs_table_karnath ]

# This cell performs the data download
download_reduced_set = False

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


#moving files to convenient directory, create new one is likely
if not os.path.exists('./simult_tweakreg/'):
    os.mkdir('./simult_tweakreg/')

# if not os.path.exists('./simult_tweakreg/ngc2071_megeath_11548/'):
#     os.mkdir('./simult_tweakreg/ngc2071_megeath_11548/')

# if not os.path.exists('./simult_tweakreg/ngc2071_karnath_16493/'):
#     os.mkdir('./simult_tweakreg/ngc2071_karnath_16493/')

import shutil
for f in glob.glob('data_download/mastDownload/HST/*11548*/*flt.fits'):
    filename = os.path.split(f)[-1]
    if not os.path.exists(filename):
        shutil.copy(f, './simult_tweakreg/') #ngc2071_karnath_16493/')

import shutil
for f in glob.glob('data_download/mastDownload/HST/*16493*/*flt.fits'):
    filename = os.path.split(f)[-1]
    if not os.path.exists(filename):
        shutil.copy(f, './simult_tweakreg/') #ngc2071_megeath_11548/')

#cd to alignment tests and then work on files further...




# ----------------------------------------------------------------------------------------------------------

# 
def get_footprints(im_name):
    """Calculates positions of the corners of the science extensions of some image 'im_name' in sky space"""
    footprints = []
    hdu = fits.open(im_name)
    
    flt_flag = 'flt.fits' in im_name or 'flc.fits' in im_name
    
    # Loop ensures that each science extension in a file is accounted for.  This is important for 
    # multichip imagers like WFC3/UVIS and ACS/WFC
    for ext in hdu:
        if 'SCI' in ext.name:
            hdr = ext.header
            wcs = WCS(hdr, hdu)
            footprint = wcs.calc_footprint(hdr, undistort=flt_flag)
            footprints.append(footprint)
    
    hdu.close()
    return footprints

# ----------------------------------------------------------------------------------------------------------
def bounds(footprint_list):
    """Calculate RA/Dec bounding box properties from multiple RA/Dec points"""
    
    # flatten list of extensions into numpy array of all corner positions
    merged = [ext for image in footprint_list for ext in image]
    merged = np.vstack(merged)
    ras, decs = merged.T
    
    # Compute width/height
    delta_ra = (max(ras)-min(ras))
    delta_dec = max(decs)-min(decs)

    # Compute midpoints
    ra_midpt = (max(ras)+min(ras))/2.
    dec_midpt = (max(decs)+min(decs))/2.
    

    return ra_midpt, dec_midpt, delta_ra, delta_dec
# ----------------------------------------------------------------------------------------------------------

def plot_footprints(footprint_list, axes_obj=None, fill=True):
    """Plots the footprints of the images on sky space on axes specified by axes_obj """
    
    if axes_obj != None: 
        ax = axes_obj
    
    else: # If no axes passed in, initialize them now
        merged = [ext for image in footprint_list for ext in image] # flatten list of RA/Dec
        merged = np.vstack(merged)
        ras, decs = merged.T
        
        # Calculate aspect ratio
        delta_ra = (max(ras)-min(ras))*np.cos(math.radians(min(np.abs(decs))))
        delta_dec = max(decs)-min(decs)
        aspect_ratio = delta_dec/delta_ra
    
        # Initialize axes
        fig = plt.figure(figsize=[8,8*aspect_ratio])
        ax = fig.add_subplot(111)
        ax.set_xlim([max(ras),min(ras)])
        ax.set_ylim([min(decs),max(decs)])
       
        # Labels
        ax.set_xlabel('RA [deg]')
        ax.set_ylabel('Dec [deg]')
        ax.set_title('Footprint sky projection ({} images)'.format(len(footprint_list)))
        
        ax.grid(ls = ':')
    
        
    colors = cm.rainbow(np.linspace(0, 1, len(footprint_list)))
    alpha = 1./float(len(footprint_list)+1.)+.2
    
    if not fill:
        alpha =.8

    for i, image in enumerate(footprint_list): # Loop over images
        for ext in image: # Loop over extensions in images
            if isinstance(ax, wcsaxes.WCSAxes): # Check axes type
                rect = Polygon(ext, alpha=alpha, closed=True, fill=fill, 
                               color=colors[i], transform=ax.get_transform('icrs'))
            else:
                rect = Polygon(ext, alpha=alpha, closed=True, fill=fill, color=colors[i])

            ax.add_patch(rect)
    
    return ax

# ----------------------------------------------------------------------------------------------------------

# Cut the sources

def get_error_mask(catalog, max_error):
    """Returns a mask for rows in catalog where RA and Dec error are less than max_error"""
    ra_mask = catalog['ra_error']< max_error
    dec_mask = catalog['dec_error'] < max_error
    mask = ra_mask & dec_mask
#     print('Cutting sources with error higher than {}'.format(max_error))
#     print('Number of sources befor filtering: {}\nAfter filtering: {}\n'.format(len(mask),sum(mask)))
    return mask
mask = get_error_mask(r, 10.)
# Plot RA/Dec Positions after clipping 

fig = plt.figure(figsize=[10,20])
ax1 = fig.add_subplot(211)
plt.scatter(ras[mask],decs[mask],c=mags[mask],alpha=.5,s=10,vmin=14,vmax=20)
ax1.set_xlim(max(ras),min(ras))
ax1.set_ylim(min(decs),max(decs))
ax1.grid(ls = ':')
ax1.set_xlabel('RA [deg]')
ax1.set_ylabel('Dec [deg]')
ax1.set_title('Source location (error < 10. mas)')
cb = plt.colorbar()
cb.set_label('G Magnitude')


ax2 = fig.add_subplot(212)
for err_threshold in [40., 10., 5., 2.]:
    mask = get_error_mask(r, err_threshold)
    hist, bins, patches = ax2.hist(mags[mask],bins=20,rwidth=.925,
                                   range=(10,20),label='max error: {} mas'.format(err_threshold))
ax2.grid(ls = ':')
ax2.set_xlabel('G Magnitude')
ax2.set_ylabel('N')
ax2.set_title('Photometry Histogram (No Log Scale)')
# ax2.set_yscale("log")
legend = ax2.legend(loc='best')
plt.savefig('source_locations.png')

# ----------------------------------------------------------------------------------------------------------

from astropy.units import Quantity
from astroquery.gaia import Gaia

width = Quantity(delta_ra, u.deg)
height = Quantity(delta_dec, u.deg)

images_megeath = glob.glob('*megeath*/*flt.fits')
images_karnath_164 = glob.glob('*karnath*/f164n/*flt.fits')
images_karnath_167 = glob.glob('*karnath*/f167n/*flt.fits')
image_list = [images_megeath , images_karnath_164,images_karnath_167]
image_labels = ['megeath', 'karnath_164', 'karnath_167']

for i in range(len(image_list)):
    images = image_list[i]
    footprint_list = list(map(get_footprints, images))

    # # If that's slow, here's a version that runs it in parallel:
    # from multiprocessing import Pool
    # p = Pool(8)
    # footprint_list = list(p.map(get_footprints, images))
    # p.close()
    # p.join()

    ra_midpt, dec_midpt, delta_ra, delta_dec = bounds(footprint_list)

    coord = SkyCoord(ra=ra_midpt, dec=dec_midpt, unit=u.deg)
    print(coord)

    plot_footprints(footprint_list)
    plt.savefig('footprints_'+image_labels[i]+'.png')

    # Perform the query!
    #pip install astroquery --upgrade
    r = Gaia.query_object_async(coordinate=coord, width=width, height=height)

    # Print the table
    print(r)

    ras = r['ra']
    decs = r['dec']
    mags = r['phot_g_mean_mag']
    ra_error = r['ra_error']
    dec_error = r['dec_error']

    fig = plt.figure(figsize=[15,15])

    # Plot RA and Dec positions, color points by G magnitude
    ax1 = fig.add_subplot(221)
    plt.scatter(ras,decs,c=mags,alpha=.5,s=6,vmin=14,vmax=20)
    ax1.set_xlim(max(ras),min(ras))
    ax1.set_ylim(min(decs),max(decs))
    ax1.grid(ls = ':')
    ax1.set_xlabel('RA [deg]')
    ax1.set_ylabel('Dec [deg]')
    ax1.set_title('Source location')
    cb = plt.colorbar()
    cb.set_label('G Magnitude')

    # Plot photometric histogram
    ax2 = fig.add_subplot(222)
    hist, bins, patches = ax2.hist(mags,bins=15,rwidth=.925)
    ax2.grid(ls = ':')
    ax2.set_xlabel('G Magnitude')
    ax2.set_ylabel('N')
    ax2.set_title('Photometry Histogram')
    ax2.set_yscale("log")


    ax3a = fig.add_subplot(425)
    hist, bins, patches = ax3a.hist(ra_error,bins=40,rwidth=.9)
    ax3a.grid(ls = ':')
    ax3a.set_title('RA Error Histogram')
    ax3a.set_xlabel('RA Error [mas]')
    ax3a.set_ylabel('N')
    ax3a.set_yscale("log")

    ax3b = fig.add_subplot(427)
    hist, bins, patches = ax3b.hist(dec_error,bins=40,rwidth=.9)
    ax3b.grid(ls = ':')
    ax3b.set_title('Dec Error Histogram')
    ax3b.set_xlabel('Dec Error [mas]')
    ax3b.set_ylabel('N')
    ax3b.set_yscale("log")


    ax4 = fig.add_subplot(224)
    plt.scatter(ra_error,dec_error,alpha=.2,c=mags,s=1)
    ax4.grid(ls = ':')
    ax4.set_xlabel('RA error [mas]')
    ax4.set_ylabel('Dec error [mas]')
    ax4.set_title('Gaia Error comparison')
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    cb = plt.colorbar()
    cb.set_label('G Magnitude')
    plt.tight_layout()
    plt.savefig('radec_hist_'+image_labels[i]+'.png')


    # cutting sources

    from astropy.table import Table

    tbl = Table([ras, decs]) # Make a temporary table of just the positions
    cat = 'gaia_'+image_labels[i]+'.cat'
    tbl.write(cat, format='ascii.fast_commented_header') # Save the table to a file.  The format argument ensures
                                                            # the first line will be commented out.

    thresh = 10.
    mask = get_error_mask(r, thresh)

    tbl_filtered = Table([ras[mask], decs[mask]]) 
    tbl.write('gaia_filtered_{}_mas_'+image_labels[i]+'.cat'.format(thresh), format='ascii.fast_commented_header')


#ls *cat #can use as a check
#if updating wcs needed (?)
#from stwcs import updatewcs
#derp = list(map(updatewcs.updatewcs, input_images))

# # Parallelized option
# p = Pool(8)
# derp = p.map(updatewcs.updatewcs, input_images)
# p.close()
# p.join()

wcsname ='GAIA'
teal.unlearn('tweakreg')
teal.unlearn('imagefindpars')

cw = 3.5 # psf width measurement (2*FWHM).  Use 3.5 for WFC3/UVIS and ACS/WFC and 2.5 for WFC3/IR
thresh = 100

# image_folders = ['*megeath*/', '*karnath*/f164n/', '*karnath*/f167n/']
# image_labels = ['megeath', 'karnath_164', 'karnath_167']

# for i in range(len(image_folders)):
#     input_images = sorted(glob.glob(image_folders[i]+'*flt.fits')) 

    # tweakreg.TweakReg(input_images, # Pass input images
    #               updatehdr=False, # update header with new WCS solution
    #               imagefindcfg={'threshold':250.,'conv_width':cw},# Detection parameters, threshold varies for different data
    #               separation=0.0, # Allow for very small shifts
    #               refcat=cat, # Use user supplied catalog (Gaia)
    #               clean=True, # Get rid of intermediate files
    #               interactive=False,
    #               see2dplot=False,
    #               shiftfile=True, # Save out shift file (so we can look at shifts later)
    #               outshifts='shift_'+image_labels[i]+'flt_cw6.txt'
    #               wcsname=wcsname, # Give our WCS a new name
    #               reusename=True,
    #               fitgeometry='general', # Use the 6 parameter fit
    #               minobj=10)

    # for i in range(len(image_labels)):
    #     shift_table = Table.read('defaultshift_'+image_labels[i]+'_thresh'+str(thresh)+'_cw'+str(cw)+'_flt.txt', format='ascii.no_header',
    #                 names=['file', 'dx', 'dy', 'rot', 'scale', 'xrms', 'yrms'])
    #     formats = ['.2f', '.2f', '.3f', '.5f', '.2f', '.2f']
    #     for i, col in enumerate(shift_table.colnames[1: ]):
    #         shift_table[col].format = formats[i]
    #     print(shift_table)


#or
input_images = glob.glob('*flt.fits')


 
# a simpler option
tweakreg.TweakReg(input_images, 
     imagefindcfg={'threshold':thresh, 'conv_width': cw}, 
     shiftfile=True, 
     outshifts='shifts_allepochs_thresh'+str(thresh)+'_cw'+str(cw)+'_flt.txt', 
     updatehdr=True, #ALWAYS CHECK TO CHANGE THIS!!!
     interactive=False, 
     wcsname='simult_epochs')
         
# Give the 'fit residual plots' a unique name for comparison with subsequent tests.
residual_pngs = glob.glob("residual*png")
for png in residual_pngs: 
    path = os.path.abspath(os.path.join(os.curdir, png))
    new_path = os.path.abspath(os.path.join(os.curdir, 'test1_{}'.format(png)))
    os.rename(path, new_path)

#Inspect the shift file for Test1
shift_table = Table.read('shifts_allepochs_thresh'+str(thresh)+'_cw'+str(cw)+'_flt.txt', format='ascii.no_header', names=['file', 'dx', 'dy', 'rot', 'scale', 'xrms', 'yrms'])
formats = ['.2f', '.2f', '.3f', '.5f', '.2f', '.2f']
for i, col in enumerate(shift_table.colnames[1: ]):
    shift_table[col].format = formats[i]
print(shift_table)

astrodrizzle.AstroDrizzle(input_images,
    output='thresh'+str(thresh)+'_cw'+str(cw)+'_combined',
    preserve=False,
    driz_sep_bits='64,16',
    driz_cr_corr=True,
    final_bits='64,16',
    clean=False,
    configobj=None,
    build=True)

plt.figure(figsize = (10, 10))
drc_dat = fits.open('thresh'+str(thresh)+'_cw'+str(cw)+'_combined_drz.fits')['SCI', 1].data #final drizzled image in SCI,1 extension
z1, z2 = zscale.zscale(drc_dat)
plt.imshow(drc_dat, origin='lower', vmin=z1, vmax=z2, cmap='Greys')
plt.title('thresh'+str(thresh)+'_cw'+str(cw)+' drizzled image', fontsize=30)
plt.savefig('drizz.png')




#yet another option...
#simbad astroquery if you want to try this as an option
from astroquery.simbad import Simbad
import pandas as pd

def all_stars(table, keys):
    star_list = []

    for i in table:
        name, ra, dec = i['MAIN_ID'], i['RA'], i['DEC']
        for j in keys:
             if name.find(j) != -1:
                star_list.append([name, SkyCoord(ra, dec, unit=(u.hourangle, u.deg))])   
    return star_list

def qual_all_stars(table, keys):
    star_list = []

    for i in table:
        name, ra, dec, qual = i['MAIN_ID'], i['RA'], i['DEC'], i['COO_QUAL']
        for j in keys:
             if qual.find(j) != -1:
                star_list.append([name, SkyCoord(ra, dec, unit=(u.hourangle, u.deg))])   
    return star_list

# result_table = Simbad.query_region("[WRA93] IRS 7", radius=1.2/60. * u.deg)                    
coord_list = ['5h46m57.7803s +0d20m41.872s', '5h47m04.6822s +0d22m24.126s', '5h47m10.0882s +0d23m44.782s']
radius_list = ['1.26100m', '1.16885m', '1.24114m']
result_table = vstack([Simbad.query_region(SkyCoord(coord, frame='icrs'), radius=rad) for coord,rad in zip(coord_list, radius_list)])
star_keys = ['GAIA', '2MASS', 'HD'] #alternative option is to filter by coo_err_maja and coo_err_mina of approx less than 250 mas? that's what this basically grabs
star_coords = all_stars(result_table, star_keys)
qual_star_coords = qual_all_stars(result_table, ['A'])

#now searching diff images for stars, composing cat
from photutils.centroids import centroid_sources
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs.utils import pixel_to_skycoord
path_list = glob.glob('../astrodrizzle_tests/*drz*')[:2]
hdu_list = [fits.open(i) for i in path_list]
primaryhdu_list = [hdu[0].header for hdu in hdu_list]
wcs_header_list = [hdu[1].header for hdu in hdu_list]
data_list = [hdu[1].data for hdu in hdu_list]

#grabbing coords and saving into cat file
for i in range(len(hdu_list)):
    refstar_pix_list = []
    stars_in_field = []

    for refstar in qual_star_coords:
        refstar_sky = refstar[1]
        refstar_pix = np.array(skycoord_to_pixel(refstar_sky, WCS(wcs_header_list[i])))
        refstar_pix_list.append(refstar_pix)
        refstar_centroids_pix = centroid_sources(data_list[i], refstar_pix[0], refstar_pix[1])
        refstar_centroids_sky = pixel_to_skycoord(refstar_centroids_pix[0], refstar_centroids_pix[1], WCS(wcs_header_list[i]))
        stars_in_field.append(np.array([refstar_centroids_sky.ra.value[0], refstar_centroids_sky.dec.value[0]]))
    df = pd.DataFrame(stars_in_field, columns = ['X','Y'])
    df.to_csv(primaryhdu_list[i]['FILENAME'].split('.')[0] + '.cat', encoding='ascii',index=False, header=False, sep=' ')

#splitting images up to be tweakreg'd together per epoch based on a reference catalog
id_list = ['*11548*', '*16493*']
cat_name_list = glob.glob('*.cat')
cw = 3.5 # psf width measurement (2*FWHM).  Use 3.5 for WFC3/UVIS and ACS/WFC and 2.5 for WFC3/IR
thresh = 100
for i in range(len(id_list)):
    input_images = glob.glob(id_list[i])

    tweakreg.TweakReg(input_images, # Pass input images
              updatehdr=True, # update header with new WCS solution
              imagefindcfg={'threshold':thresh,'conv_width':cw},# Detection parameters, threshold varies for different data
              separation=0.0, # Allow for very small shifts, this is default though
              refcat=cat_name_list[i], # Use user supplied catalog
              clean=True, # Get rid of intermediate files
              interactive=False,
              see2dplot=False,
              shiftfile=True, # Save out shift file (so we can look at shifts later)
              outshifts='shift_'+cat_name_list[i].split('.')[0]+'searchrad5.txt',
              wcsname=cat_name_list[i].split('.')[0]+'searchrad5_wcs', # Give our WCS a new name,
              fitgeometry='shift',
              minobj=3,
              searchrad=5, #in arcsec?
              reusename=True)    

for line in open('shift_'+cat_name_list[i].split('.')[0]+'searchrad5.txt').readlines():
    print(line)

lam_list = ['126n', '128n', '130n', '160', '164n', '167n']
for i in range(len(lam_list)):
    input_images = glob.glob('*' + lam_list[i] + '*.fits')

    astrodrizzle.AstroDrizzle(input_images,
        output=lam_list[i] + '_combined',
        preserve=False,
        driz_sep_bits='64,16',
        driz_cr_corr=True,
        final_bits='64,16',
        clean=False,
        configobj=None,
        build=True,
        driz_sep_pixfrac=0.5,
        final_pixfrac=0.5)

    plt.figure(figsize = (10, 10))
    drc_dat = fits.open(lam_list[i]+'_combined_drz.fits')['SCI', 1].data #final drizzled image in SCI,1 extension
    z1, z2 = zscale.zscale(drc_dat)
    plt.imshow(drc_dat, origin='lower', vmin=z1, vmax=z2, cmap='Greys')
    plt.title(lam_list[i]+' drizzled image', fontsize=30)
    plt.savefig('drizz.png')