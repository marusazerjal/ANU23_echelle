import os,sys,string,glob
from astropy.io import fits as pyfits
from numpy import *

def average_header_fits(folder,input_header,output_file):
    print('ARGS, folder,input_header,output_file', folder,input_header,output_file)
    fitslist = sort(glob.glob(folder+"*.fits"))
    headermatch = []
    #print('fitslist', fitslist)
    for fits in fitslist:
        if pyfits.getheader(fits)["OBJECT"] == input_header:
        #~ if pyfits.getheader(fits)["object"] == input_header:
            print('object fits', fits)
            headermatch.append(pyfits.getdata(fits))

    headermatch = array(headermatch)
    print('headermatch', headermatch)
    headermatch = average(headermatch,axis=0)

    hdu = pyfits.PrimaryHDU(headermatch)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(output_file, overwrite=True)

    return headermatch


if __name__ == "__main__":
    folder = "/data/yjzhou0/ANU23/echelle/20180702/"
    input_header = "bias"
    output_file = folder+"temp/masterbias.fits"

    average_header_fits(folder,input_header,output_file)
    
    
