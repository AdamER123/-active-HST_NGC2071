#here we try a basic example to make sure tweakreg and astrodrizzle work!

#ls *cat #can use as a check
#if updating wcs needed (?)
#from stwcs import updatewcs
#derp = list(map(updatewcs.updatewcs, input_images))

# # Parallelized option
# p = Pool(8)
# derp = p.map(updatewcs.updatewcs, input_images)
# p.close()
# p.join()

#I think this is to get the drizzlepac commands working? forgot...
teal.unlearn('tweakreg')
teal.unlearn('imagefindpars')
cw = 3.5 # psf width measurement (2*FWHM).  Use 3.5 for WFC3/UVIS and ACS/WFC and 2.5 for WFC3/IR
thresh = 100 #feel free to change this

# get image files in a list
# please change this to make sure epoch1 and epoch2 are handled separately
# this is because tweakreg must be applied to each visit, though you can do multiple wavelengths at once
# astrodrizzle, though, must be done for each wavelength separately
input_images_epoch1 = glob.glob('*flt.fits')
input_images_epoch2 = glob.glob('*flt.fits')

# tweakreg time
'''
Rerun just tweakreg until you find the shiftfile to be satisfactory
when this is the case, turn updatehdr to True, and it will begin to tweak registration for the images (why I copy them to a new folder)
you can check the shifts with the shift_table and residual_pngs below (should be mostly automated)
'''
outshifts_path = 'shifts_allepochs_thresh'+str(thresh)+'_cw'+str(cw)+'_flt.txt' #change this as needed...

tweakreg.TweakReg(input_images_epoch1, 
     imagefindcfg={'threshold':thresh, 'conv_width': cw}, 
     shiftfile=True, 
     outshifts=outshifts_path, 
     updatehdr=False, #ALWAYS CHECK TO CHANGE THIS!!! When you want to continue to astrodrizzle
     interactive=False, 
     wcsname='simult_epochs')
         
# Give the 'fit residual plots' a unique name for comparison with subsequent tests.
residual_pngs = glob.glob("residual*png")
for png in residual_pngs: 
    path = os.path.abspath(os.path.join(os.curdir, png))
    new_path = os.path.abspath(os.path.join(os.curdir, 'test1_{}'.format(png)))
    os.rename(path, new_path)

#Inspect the shift file for Test1
shift_table = Table.read(outshifts_path, format='ascii.no_header', names=['file', 'dx', 'dy', 'rot', 'scale', 'xrms', 'yrms'])
formats = ['.2f', '.2f', '.3f', '.5f', '.2f', '.2f']
for i, col in enumerate(shift_table.colnames[1: ]):
    shift_table[col].format = formats[i]
print(shift_table)


'''
astrodrizzle, next part!
For default, I don't think many changes should be needed here other than making sure you're using the same input_images
'''
astrodrizzle.AstroDrizzle(input_images_epoch1,
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
