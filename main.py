#####!/usr/bin/python

"""
marusa@mash:~>source venv_echelle/bin/activate
"""


#activate_this = '/data/mash/marusa/reductions_echelle/venv_echelle/bin/activate'
#execfile(activate_this, dict(__file__=activate_this))

#print('venv_echelle activated.')

import os,sys,string,glob,pickle
from astropy.io import fits as pyfits
from numpy import *
import config_file
import calibrate
import mask_orders
#import thar_straighten
import extract_order
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import average_adjacent_obs
import wavecal_all
import copy

import imp
config_filename = sys.argv[1]
print('CONFIG FILENAME', config_filename)
config_file = imp.load_source(config_filename.replace('.py', ''), config_filename)
config = config_file.set_config()

try:
    do_main_loop = sys.argv[2]
    if do_main_loop=='True':
        do_main_loop=True
    elif do_main_loop=='False':
        do_main_loop=False
except:
    do_main_loop=True

print(do_main_loop)


binning = config["binning"]
ccdsec_min = 53
ccdsec_max = 2095

def return_obslist(folder):
    fitslist = sort(glob.glob(folder+"*.fits"))

    nonobjectlist = [
        config["bias"],
        config["flat"],
        config["dark"],
        config["arc"],
        "test",
        "sky",
        "Sky"]

        
    obslist = [[],[]]
    for fits in fitslist:
        if not pyfits.getheader(fits)["OBJECT"] in nonobjectlist:
            obslist[0].append(fits)
            obslist[1].append(pyfits.getheader(fits)[config["jdheader"]])

    return obslist
    
def return_tharlist(folder):
    fitslist = sort(glob.glob(folder+"*.fits"))
    tharlist = [[],[]]
    for fits in fitslist:
        if pyfits.getheader(fits)["OBJECT"] == config["arc"]:
            tharlist[0].append(fits)
            tharlist[1].append(pyfits.getheader(fits)[config["jdheader"]])
        #~ else:
            #~ print(pyfits.getheader(fits)["OBJECT"])
    return tharlist

def loadthar(jd,tharlist,masterbias):
    tharjd = array(tharlist[1])
    indexes = argsort(abs(tharjd-jd))

    tdiff = [abs(tharjd-jd)[indexes[0]]]
    thar = []
    ### get the nearest one
    thar.append(pyfits.getdata(tharlist[0][indexes[0]])[:,ccdsec_min:ccdsec_max]-masterbias)
    # if abs(tharjd-jd)[indexes[1]] < 1/24.: ### if the second thar is less than 1 hr diff, get it
    #     thar.append(pyfits.getdata(tharlist[0][indexes[1]])[:,ccdsec_min:ccdsec_max]-masterbias)
    #     tdiff.append(abs(tharjd-jd)[indexes[1]])


    thar.append(pyfits.getdata(tharlist[0][indexes[1]])[:,ccdsec_min:ccdsec_max]-masterbias)
    tdiff.append(abs(tharjd-jd)[indexes[1]])

        

    tdiff = array(tdiff)
    weights = 1/tdiff
    weights /= sum(weights)

    print("time diff between obs and thar",tdiff*24,"hr")

    thar_averaged = thar[0]*weights[0] + thar[1]*weights[1]

    save(config["folder"]+"/temp/thar.npy",thar_averaged)
    
    return thar_averaged
    
        
def main():
    # Check if RV_standards includes only one line per a unique object.
    # E.g. if someone adds a line with the same object name, the code
    # later crashes.
    tmp = loadtxt('RV_standard.dat', dtype=str)
    l=list(tmp[:,0])
    ex=False
    for x in tmp[:,0]:
        if l.count(x)>1:
            print('%s is listed more than once in the RV_standard.dat. Please keep only one!'%x)
            ex=True
    if ex:
        exit()
    
    
    if do_main_loop:
        print("making folders")
        print(config)
        try:
            os.stat(os.path.join(config["folder"], "temp/"))
        except:
            print(os.path.join(config["folder"], "temp/"))
            os.mkdir(os.path.join(config["folder"], "temp/"))  
        try:
            os.stat(os.path.join(config["folder"], "reduced/"))
        except:
            os.mkdir(os.path.join(config["folder"], "reduced/"))

        if config["delete_prev"]:
            os.system("rm -rf "+config["folder"]+"temp/*")
            os.system("rm -rf "+config["folder"]+"reduced/*")
        
        if not os.path.exists(os.path.join(config["folder"], "temp/masterbias.fits")):
            print("create master bias")
            print('ccdsec_min:ccdsec_max', ccdsec_min, ccdsec_max)
            #test=calibrate.average_header_fits(config["folder"], config["bias"], os.path.join(config["folder"], "temp/masterbias.fits"))
            #print (test)
            print('LEN', (calibrate.average_header_fits(config["folder"], config["bias"], os.path.join(config["folder"], "temp/masterbias.fits"))))
            masterbias = calibrate.average_header_fits(config["folder"], config["bias"], os.path.join(config["folder"], "temp/masterbias.fits"))[:,ccdsec_min:ccdsec_max]

        else:
            masterbias = pyfits.getdata(os.path.join(config["folder"], "temp/masterbias.fits"))[:,ccdsec_min:ccdsec_max]

        # MASTERFLAT
        if not os.path.exists(os.path.join(config["folder"], "temp/masterflat.fits")):
            print("create master flat")
            print(ccdsec_min, ccdsec_max, type(ccdsec_min), type(ccdsec_max))
            masterflat = calibrate.average_header_fits(config["folder"], config["flat"], os.path.join(config["folder"], "temp/masterflat.fits"))[:,ccdsec_min:ccdsec_max]
            print('flat before bias', masterflat)
            masterflat -= masterbias
            print('flat AFTER bias', masterflat)

        else:
            masterflat = pyfits.getdata(os.path.join(config["folder"], "temp/masterflat.fits"))[:,ccdsec_min:ccdsec_max]
            print('flat before bias', masterflat)
            masterflat -= masterbias
            print('flat AFTER bias', masterflat)

        savetxt('masterflat.dat', masterflat)

        print("determine order masks")
        if os.path.exists(os.path.join(config["folder"], "temp/order_masks_LETS_DISABLE_THIS.pkl")):
            order_masks = pickle.load(open(os.path.join(config["folder"], "temp/order_masks_LALALA.pkl"), "rb"))
        else:
            print('MASTERFLAT', masterflat)
            plt.imshow(masterflat, aspect="auto")
            order_masks = mask_orders.return_masks(masterflat, config=config)
            print 'ORDER MASKS DETERMINED', order_masks
        plt.show()

        masterflat_extracted = extract_order.mask_order(masterflat,order_masks, binning=binning)

    ### determine a list of observations that need extracting
    tharlist = return_tharlist(config["folder"])
    obslist = return_obslist(config["folder"])
    
    #~ #"""
    if do_main_loop:
        for i in range(len(obslist[0])):
            fitsname = os.path.basename(obslist[0][i])
            print(fitsname)
            #print('&&&&&&& fitsname', fitsname, os.path.exists(os.path.join(config["folder"], "temp/", fitsname+".spec.pkl")))
            if not os.path.exists(os.path.join(config["folder"], "temp/", fitsname+".spec.pkl")): ### check whether obs was reduced
                print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
                print("Reducing",obslist[0][i])
                thar = loadthar(obslist[1][i],tharlist,masterbias)

                #print("Calculating shear"
                #shear_list = thar_straighten.find_straighten_function(thar,order_masks)
                #print("Applying shear"
                #thar_shear = thar_straighten.shear_obs(extract_order.mask_order(thar,order_masks),shear_list)
                thar_shear = extract_order.mask_order(thar,order_masks, binning=binning)
                fits = pyfits.getdata(obslist[0][i])[:,ccdsec_min:ccdsec_max]
                fits = fits.astype(float64) # Marusa
                #print fits
                #print(type(fits[0][0]), type(masterbias[0][0]))

                fits -= masterbias

                # Mask order
                fits = extract_order.mask_order(fits,order_masks, binning=binning)

                # import copy
                fits_noflat = copy.deepcopy(fits)
                print("Dividing by flat")
                for order in range(len(fits)):
                    flat_order = masterflat_extracted[order][0]
                    flat_order /= nanmax(flat_order.flatten())

                    fits[order][0] /= flat_order
                    
                    
                #fits_shear = thar_straighten.shear_obs(fits,shear_list)
                #fits_noflat_shear = thar_straighten.shear_obs(fits_noflat,shear_list)
                fits_shear = fits
                fits_noflat_shear = fits_noflat

                #pickle.dump(thar_shear,open(config["folder"]+"/temp/"+fitsname+".thar.pkl","wb"))
                #pickle.dump(fits_shear,open(config["folder"]+"/temp/"+fitsname+".shear.pkl","wb"))

                # Tracing orders
                trace_array = extract_order.find_trace(fits_noflat_shear, binning=binning)
                plt.savefig(os.path.join(config["folder"], "temp/", fitsname+"trace.pdf"))
                plt.clf()
                pickle.dump(trace_array, open(os.path.join(config["folder"], "temp/", fitsname+".trace.pkl"),"wb"))

                print("extracting spectra")
                print(fits_shear)
                print(len(fits_shear))
                print(len(trace_array))
                spectrum,background = extract_order.extract_trace(fits_shear,trace_array)
                spectrum_noflat,background_noflat = extract_order.extract_trace(fits_noflat_shear,trace_array)
                tharspec,bk = extract_order.extract_trace(thar_shear,trace_array)
                for i in range(len(tharspec)):
                    tharspec[i] += bk[i]

                pickle_filename = os.path.join(config["folder"], "temp/", fitsname+".spec.pkl")
                print 'DUMP pickle', pickle_filename
                pickle.dump([spectrum,background,tharspec,spectrum_noflat,background_noflat], open(pickle_filename, "wb"))
                
                #sys.exit()
        #"""
    print("Creating fits files") # Produces ANU..fits files in the temp/ folder.
    print('obslist', obslist)
    print('tharlist', tharlist)
    average_adjacent_obs.average_adjacent_obs(obslist,tharlist,config["folder"])
    
    print('\n\ndo wavecal_all')
    wavecal_all.main(config["folder"], config=config)

    if config["run_analysis"]:

        print("Running spectral analysis")
        
        import run_analysis
        # Photometry to estimate temperature to find the template. Photometry should be a dictionary with objectnames and BP-RP
        #~ run_analysis.main(config["folder"], photometry)
        run_analysis.main(config["folder"])
    
if __name__ == "__main__":

    main()
    
    print('DONE.')
