#this is a pseudo-script combining some notebooks by Adam Rubinstein (advised by Tom Megeath, Sam Federman, Dan Watson, Joel Green)
#purpose is to align HST/WFC3 images from NGC 2071

#generally I run this in ipython or as a shell script...
#always start in the right env
# source activate astroconda 
# ipython

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
import numpy as np
# import matplotlib.pyplot as plt
import os
import drizzlepac.tweakreg as tweakreg  # from Drizzle packages we'll need
import drizzlepac.astrodrizzle as astrodrizzle
from astropy.io import ascii    # from Astropy packages we'll need
from astropy.io import fits
import sys

#### here we try a basic example to make sure tweakreg and astrodrizzle work!

#references to original data file and processing we're trying to reproduce
refimg = 'hops335.fits' #the first epoch image created by marina kounkel
hdu=fits.open('hops335.fits')
header1=hdu[1].header
# fits.writeto('hops335_simple.fits', hdu[1].data, header=header1, overwrite=True) # or data1=fits.getdata('hops335.fits')
newrefimg='hops335_simple.fits' #needed for first step of astrodrizzle (?)

#setting up for a loop over images with different wavelengths
# lambda_sets = ['130n'] #if just testing one example
# image_sets = glob.glob('*'+lambda_sets[0]+'*_flt.fits')
#note that f160w is ID 11548, while the rest are ID 16493 
path = '../original_center_frames/'
lambda_sets = ['126n', '128n', '130n', '164n', '167n', '160w'] #this set only uses the epoch 2 frames
image_sets = [glob.glob(path + '*'+i+'*') for i in lambda_sets]

#other parameters
cw = 2.5 # psf width measurement (2*FWHM).  Use 3.5 for WFC3/UVIS and ACS/WFC and 2.5 for WFC3/IR
thresh_list = [100, 100, 100, 100, 100, 100] #feel free to change this
min_obj_list = [3, 3, 3, 3, 3, 3]
combine_list = ['imedian', 'imedian', 'imedian', 'imedian', 'imedian', 'imedian'] #recommendation is minmed for <6 frames, or median and mean for >4 frames (but you need to change another paramter if you use it)
# blot_param = 'linear' #controls blotting and interpolating onto the WCS


#looping through the different images
for i in range(len(image_sets)):
    # # tweakreg time!
    # '''
    # Goal: Rerun just tweakreg until you find the shiftfile to be satisfactory
    # when this is the case, turn updatehdr to True, and it will begin to tweak registration for the images (why I copy them to a new folder)
    # you can check the shifts with the shift_table and/or residual_pngs below (should be mostly automated)
    # '''
    # outshifts_path = 'shifts_' + lambda_sets[i] + '_thresh'+str(thresh_list[i])+'_minobj'+str(min_obj_list[i])+'_flt.txt' #change this as needed...

    # tweakreg.TweakReg(image_sets[i], 
    #     imagefindcfg={'threshold':thresh_list[i], 'conv_width': cw}, 
    #     refimagefindcfg={'threshold':thresh_list[i], 'conv_width': cw},
    #     #enforce_user_order=False,
    #     #catfile=cat_file,
    #     #expand_refcat=True,
    #     refimage=refimg,
    #     nclip=3,
    #     sigma=3,
    #     minobj=min_obj_list[i],
    #     searchrad=3,
    #     shiftfile=True, 
    #     outshifts=outshifts_path, 
    #     updatehdr=True, #ALWAYS CHECK TO CHANGE THIS!!! When you want to continue to astrodrizzle
    #     interactive=False, 
    #     wcsname='epoch2_mod_' + lambda_sets[i])

    # #Inspect the shift file for Test1
    # shift_table = ascii.read(outshifts_path, format='no_header', names=['file', 'dx', 'dy', 'rot', 'scale', 'xrms', 'yrms'])
    # formats = ['.2f', '.2f', '.3f', '.5f', '.2f', '.2f']
    # for j, col in enumerate(shift_table.colnames[1: ]):
    #     shift_table[col].format = formats[j]
    # print(shift_table)

    '''
    astrodrizzle, next part!
    For default, I don't think many changes should be needed here other than making sure you're using the same input_images
    '''

    #need to find out parameters for putting onto the grid
    astrodrizzle.AstroDrizzle(image_sets[i],
        output=lambda_sets[i] + '_thresh'+str(thresh_list[i])+'_minobj'+str(min_obj_list[i])+'_'+combine_list[i]+'_realigned'+'_dropsize1p0_finalpixfrac1_nocr',
        # preserve=False,
        driz_sep_pixfrac=1.0,
        # blot_interp = blot_param,
        final_pixfrac=1.0,
        final_wcs=True,
        final_refimage=newrefimg,                      
        final_outnx=2229,
        final_outny=1982,
        final_rot = 0.0, #needed to point everything to north as up on the image...
        final_crpix1 = 1114.788726376634,
        final_crpix2 = 991.4713893880687,
        final_ra = 86.7698350000,
        final_dec = 0.369595643698,
        combine_type = combine_list[i],
        # driz_sep_bits='512', #look into raw image for DQ... try removing 16 (which is hot pixels)
        # final_bits='16 64 512',
        # driz_cr=True,
        # driz_cr_corr=True,
        clean=False,
        configobj=None,
        build=True)