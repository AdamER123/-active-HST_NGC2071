
# ----------------------------------------------------------------------------------------------------------
# using simbad astroquery if you want to gather many catalogs
from astroquery.simbad import Simbad
import pandas as pd

# a function to grab all the reference stars from simbad
def all_stars(table, keys):
    star_list = []

    for i in table:
        name, ra, dec = i['MAIN_ID'], i['RA'], i['DEC'] #can also use other metrics, like cooqual, etc at http://simbad.u-strasbg.fr/simbad/sim-fsam
        for j in keys:
             if name.find(j) != -1: #this is the filter I'm using, which is by checking if the IDs match any star_keys I'm interested in (mainly 2MASS)
                star_list.append([name, SkyCoord(ra, dec, unit=(u.hourangle, u.deg))])   
    return star_list

result_table = Simbad.query_region("[WRA93] IRS 7", radius=1.2/60. * u.deg)                    
star_keys = ['GAIA', '2MASS', 'V*', 'EM*']
star_coords = all_stars(result_table, star_keys)

#here is an alternative example
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
qual_star_coords = qual_all_stars(result_table, ['A'])


#now searching diff images for stars, composing cat
from photutils.centroids import centroid_sources
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs.utils import pixel_to_skycoord

#you'll need a sample WCS for this to work. I just use one from a default astrodrizzle image that I like, depends if you want absolute or relative alignment
path_list = glob.glob('../astrodrizzle_tests/*drz*')[:2] 
hdu_list = [fits.open(i) for i in path_list]
primaryhdu_list = [hdu[0].header for hdu in hdu_list]
wcs_header_list = [hdu[1].header for hdu in hdu_list]
data_list = [hdu[1].data for hdu in hdu_list]

#grabbing coords and saving into cat file...this gets messy as you'll see...
for i in range(len(hdu_list)):
    refstar_pix_list = []
    stars_in_field = []

    for refstar in star_coords:
        refstar_sky = refstar[1] #reference star's skycoords from simbad
        refstar_pix = np.array(skycoord_to_pixel(refstar_sky, WCS(wcs_header_list[i]))) #converting them to simbad with our template header
        refstar_pix_list.append(refstar_pix)
        refstar_centroids_pix = centroid_sources(data_list[i], refstar_pix[0], refstar_pix[1]) #we needed that because astropy's centroid function only takes pixel coordinates (?)
        refstar_centroids_sky = pixel_to_skycoord(refstar_centroids_pix[0], refstar_centroids_pix[1], WCS(wcs_header_list[i])) #then we convert BACK to skycoords so we can check everything is consistent by hand and with drizzlepac
        stars_in_field.append(np.array([refstar_centroids_sky.ra.value[0], refstar_centroids_sky.dec.value[0]]))
    df = pd.DataFrame(stars_in_field, columns = ['X','Y'])
    df.to_csv(primaryhdu_list[i]['FILENAME'].split('.')[0] + '.cat', encoding='ascii',index=False, header=False, sep=' ') #and so we made a cat

#splitting images up to be tweakreg'd together per epoch based on a reference cat
id_list = ['*11548*', '*16493*']
cat_name_list = glob.glob('*.cat')
cw = 3.5 # psf width measurement (2*FWHM).  Use 3.5 for WFC3/UVIS and ACS/WFC and 2.5 for WFC3/IR
thresh = 100
for i in range(len(id_list)):
    input_images = glob.glob(id_list[i])


    # THE SAME APPLIES AS IN ALIGN_DEFAULT.PY. Please check your shift files, THEN set updatehdr to True
    tweakreg.TweakReg(input_images, # Pass input images
              updatehdr=False, # update header with new WCS solution
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

#outputting shift files
for line in open('shift_'+cat_name_list[i].split('.')[0]+'searchrad5.txt').readlines():
    print(line)

#need to split everything up by wavelength, then we astrodrizzle. Here I was experimenting with drop size too
lam_list = ['126n', '128n', '130n', '160', '164n', '167n']
for i in range(len(lam_list)):
    input_images = glob.glob('*' + lam_list[i] + '*.fits')

    astrodrizzle.AstroDrizzle(input_images,
        output='thresh'+str(thresh)+'_cw'+str(cw)+'_f'+lam_list[i] + '_combined',
        preserve=False,
        driz_sep_pixfrac=0.75,
        final_wcs=True,
        final_outnx=2229,
        final_outny=1982,
        final_refimage=newrefimg,                      
        driz_sep_bits='512', #look into raw image for DQ... try removing 16 (which is hot pixels)
        driz_cr_corr=True,
        final_bits='512',
        clean=False,
        configobj=None,
        build=True)

    plt.figure(figsize = (10, 10))
    drc_dat = fits.open(lam_list[i]+'_combined_drz.fits')['SCI', 1].data #final drizzled image in SCI,1 extension
    z1, z2 = zscale.zscale(drc_dat)
    plt.imshow(drc_dat, origin='lower', vmin=z1, vmax=z2, cmap='Greys')
    plt.title(lam_list[i]+' drizzled image', fontsize=30)
    plt.savefig('drizz.png')


# ----------------------------------------------------------------------------------------------------------
# Implementing an example that uses GAIA

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


from astropy.units import Quantity
from astroquery.gaia import Gaia

width = Quantity(delta_ra, u.deg)
height = Quantity(delta_dec, u.deg)

images_megeath = glob.glob('*megeath*/*flt.fits')
images_karnath_164 = glob.glob('*karnath*/f164n/*flt.fits')
images_karnath_167 = glob.glob('*karnath*/f167n/*flt.fits')
image_list = [images_megeath , images_karnath_164,images_karnath_167]
image_labels = ['megeath', 'karnath_164', 'karnath_167']

from astropy.table import Table
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

image_folders = ['*megeath*/', '*karnath*/f164n/', '*karnath*/f167n/']
image_labels = ['megeath', 'karnath_164', 'karnath_167']

for i in range(len(image_folders)):
    input_images = sorted(glob.glob(image_folders[i]+'*flt.fits')) 

    tweakreg.TweakReg(input_images, # Pass input images
                  updatehdr=False, # update header with new WCS solution
                  imagefindcfg={'threshold':250.,'conv_width':cw},# Detection parameters, threshold varies for different data
                  separation=0.0, # Allow for very small shifts
                  refcat=cat, # Use user supplied catalog (Gaia)
                  clean=True, # Get rid of intermediate files
                  interactive=False,
                  see2dplot=False,
                  shiftfile=True, # Save out shift file (so we can look at shifts later)
                  outshifts='shift_'+image_labels[i]+'flt_cw6.txt'
                  wcsname=wcsname, # Give our WCS a new name
                  reusename=True,
                  fitgeometry='general', # Use the 6 parameter fit
                  minobj=10)

    for i in range(len(image_labels)):
        shift_table = Table.read('defaultshift_'+image_labels[i]+'_thresh'+str(thresh)+'_cw'+str(cw)+'_flt.txt', format='ascii.no_header',
                    names=['file', 'dx', 'dy', 'rot', 'scale', 'xrms', 'yrms'])
        formats = ['.2f', '.2f', '.3f', '.5f', '.2f', '.2f']
        for i, col in enumerate(shift_table.colnames[1: ]):
            shift_table[col].format = formats[i]
        print(shift_table)


