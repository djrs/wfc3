import numpy as np
import pyfits
import os

##################################################################
#                                                                #
#  Routines to identify cosmic rays from Dolphot output and to   #
#  subsequently mask them out in the original F606W FLT image.   #
#                                                                #
#  Variables Used:                                               #
#     peak_threshold - Minimum value of cosmic ray to mask out   #
#     background_threshold - All values above this limit         #
#          that are adjacent to the cosmic ray will be masked    #
#     recursion_depth - Maximum recursion depth for searching    #
#          adjacent pixels                                       #
#     wcs_path - Path to the astrometry.net routines             #
#                                                                #
#  - D. Radburn-Smith, April 2013                                #
#                                                                #
##################################################################

peak_threshold=300
background_threshold=90
recursion_depth=30
wcs_path='/usr/local/astrometry/bin/'


def MaskOut(x,y,data,mask,shape,depth):
    ''' Recursively search for bright pixels adjacent to peak detection '''
    if (depth<recursion_depth):
        for xc in [x,x-1,x+1]:
            for yc in [y,y-1,y+1]:
                if (xc>0) & (xc<2050) & (yc>0) & (yc<4096):
                    if (mask[xc,yc]==0) & (data[xc,yc]>background_threshold):
                        mask[xc,yc]=1
                        depth+=1
                        mask=MaskOut(xc,yc,data,mask,shape,depth)
    return mask


def CosmicSearch(x,y,data):
    ''' Search a 7x7 pixel block around the probable location for the peak '''
    block=data[x-3:x+4,y-3:y+4]
    ind=np.indices(np.shape(block))
    xval=x-3+ind[0][(block==np.max(block))]
    yval=y-3+ind[1][(block==np.max(block))]
    return xval[0],yval[0]


def ApplyMask(path,fltfiles):
    ''' Apply the mask images to the original F606W image '''
    print 'Applying masks to originial F606W FLT image'
    mask1=pyfits.getdata('mask1.fits')
    mask2=pyfits.getdata('mask2.fits')
    image_name=fltfiles[0].replace('chip1.fits','fits') # Name of original, unprocessed FLT file
    os.system('cp %s/../hst/%s %s/'%(path,image_name,path)) # Copy original to dph directory
    sourcefile='%s/%s'%(path,image_name)
    hdu=pyfits.open(sourcefile,mode='update')
    tmpdata=hdu[3].data
    sel=(mask1==1)
    hdu[3].data[sel]=4096 # Just tag as cosmic ray regardless of previous tags
    tmpdata=hdu[6].data
    sel=(mask2==1)
    hdu[6].data[sel]=4096 # Just tag as cosmic ray regardless of previous tags
    hdu.flush()
    hdu.close()
    os.system('rm mask1.fits mask2.fits')
    print 'Running Dolphot preprocessing steps'
    os.system('wfc3mask %s/%s'%(path,image_name)) # Run Dolphot tasks on masked FLT image
    os.system('splitgroups %s/%s'%(path,image_name))
    os.system('calcsky %s/%s 15 35 4 2.25 2.00'%(path,image_name.replace('.fits', '.chip1')))
    os.system('calcsky %s/%s 15 35 4 2.25 2.00'%(path,image_name.replace('.fits', '.chip2')))


def ParseInfoFile(path,prefix):
    ''' search the [prefix].phot.info file for the F606W chip images '''
    with open('%s/%s.phot.info'%(path,prefix),'r') as infile:
        filenames=[]
        f606w_images=[]
        startrecording=False
        for line in infile:
            elements=line.split()
            if (len(elements)>0):
                if (len(elements)>1):
                    if (elements[1]=='sets'):
                        startrecording=True
                    if (len(elements)>3):
                        if (elements[1]=='image') & (elements[3]=='F606W'):
                            f606w_images=np.append(f606w_images,int(elements[2][:-1]))
                if (startrecording) & (line.find('chip')>0):
                    filenames.append(line.strip())
                if (elements[0]=='EXTENSION'):startrecording=False
    f606w_names=[]
    for element in f606w_images: f606w_names.append('%s.fits'%filenames[int(element)-1])
    return f606w_names


def CreateMask(path,prefix):
    ''' Main Routine: path is the path to the dph photometry directory, prefix is the galaxy-field combination '''
    print 'Selecting cosmic rays from Dolphot output'
    dat=np.loadtxt('%s/%s.phot.gz'%(path,prefix)) # Path to the raw photometry
    sel=((dat[:,28]>25.75) & (dat[:,15]<26.0) & (dat[:,15]-dat[:,28]<-1.0)) | (dat[:,28]>50) & (dat[:,15]<26.0) # Select Cosmics based on magnitudes and color
    #sel = (sel) & (dat[:,10]<=2)# & (dat[:,18]<5.5) & (dat[:,20]>-0.09) & (dat[:,20]<0.35) & (dat[:,21]>-1.0) & (dat[:,21]<1.0) & (dat[:,22]<0.1) & (dat[:,23]<1)
    # Originaly selected only good stars, but dolphot repeatedly used cosmics that were just less well defined
    xpos=dat[sel,2]
    ypos=dat[sel,3]
    flt_files=np.sort(ParseInfoFile(path,prefix))
    
    print 'Converting coordinates'
    data_array1=pyfits.getdata('%s/%s'%(path,flt_files[0])) # F606W['SCI'] chip 1 image
    data_array2=pyfits.getdata('%s/%s'%(path,flt_files[1])) # F606W['SCI'] chip 2 image
    mask1=np.zeros(np.shape(data_array1))
    mask2=np.zeros(np.shape(data_array2))
    col1=pyfits.Column(name='x',format='E',array=xpos)
    col2=pyfits.Column(name='y',format='E',array=ypos)
    if os.path.exists('xy.fits'):os.system('rm xy.fits')
    tbhdu=pyfits.new_table(pyfits.ColDefs([col1,col2])).writeto('xy.fits') # Create Fits table with x, y positions of detections
    if os.path.exists('rd.fits'):os.system('rm rd.fits')
    os.system('%s/wcs-xy2rd -w %s/%s-f814w_drz.chip1.fits -i xy.fits -o rd.fits'%(wcs_path,path,prefix)) # Convert x, y to RA, Dec
    os.system('rm xy.fits')
    if os.path.exists('xy1.fits'):os.system('rm xy1.fits')
    if os.path.exists('xy2.fits'):os.system('rm xy2.fits')
    os.system('%s/wcs-rd2xy -w %s/%s -i rd.fits -o xy1.fits'%(wcs_path,path,flt_files[0])) # Convert RA, Dec to x, y in chip images
    os.system('%s/wcs-rd2xy -w %s/%s -i rd.fits -o xy2.fits'%(wcs_path,path,flt_files[1]))
    os.system('rm rd.fits')
    tbdata=pyfits.getdata('xy1.fits')
    x=tbdata['X']
    y=tbdata['Y']
    os.system('rm xy1.fits')
    sel=(x>4) & (x<4095) & (y>4) & (y<2046) # Don't get too close to the chip edges
    x=x[sel];y=y[sel]
    print 'Computing mask for chip 1'
    for i,xpos in enumerate(x):
        ypos=y[i]
        xp,yp=CosmicSearch(ypos,xpos,data_array1) # Search around the predicted position for the peak brightness as the x, y / RA, dec conversions may be off
        if (data_array1[xp,yp]>peak_threshold):
            mask1=MaskOut(xp,yp,data_array1,mask1,np.shape(mask1),0) # Mask out bright pixels adjacent to the peak

    print 'Computing mask for chip 2'
    tbdata=pyfits.getdata('xy2.fits')
    x=tbdata['X']
    y=tbdata['Y']
    os.system('rm xy2.fits')
    sel=(x>4) & (x<4095) & (y>4) & (y<2046)
    x=x[sel];y=y[sel]
    for i,xpos in enumerate(x):
        ypos=y[i]
        xp,yp=CosmicSearch(ypos,xpos,data_array2)
        if (data_array2[xp,yp]>peak_threshold):
            mask2=MaskOut(xp,yp,data_array2,mask2,np.shape(mask2),0)
    
    if os.path.exists('mask1.fits'):os.system('rm mask1.fits mask2.fits')
    hdu=pyfits.PrimaryHDU(mask1) # Write out masks 0 = good, 1 = bad
    hdu.writeto('mask1.fits')
    hdu=pyfits.PrimaryHDU(mask2)
    hdu.writeto('mask2.fits')

    ApplyMask(path,flt_files) # Apply mask to the ['DQ'] HDU of the original file


if __name__=='__main__':
    CreateMask('/Users/djrs/ghosts/data/ngc3031/field10/dph','ngc3031-field10')
    
