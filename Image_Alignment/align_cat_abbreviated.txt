
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