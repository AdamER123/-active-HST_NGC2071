


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

#ls *cat #can use as a check
#if updating wcs needed (?)
#from stwcs import updatewcs
#derp = list(map(updatewcs.updatewcs, input_images)) #this is their code, not mine, I think

# # Parallelized option
# p = Pool(8)
# derp = p.map(updatewcs.updatewcs, input_images)
# p.close()
# p.join()