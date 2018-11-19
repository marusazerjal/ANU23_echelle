import os,sys,string,pyfits,glob,pickle
from numpy import *
import config_file
import calibrate
import mask_orders
#import thar_straighten
import extract_order
config = config_file.set_config()
import matplotlib.pyplot as plt
import average_adjacent_obs
import wavecal_all

binning = config["binning"]
ccdsec_min = 53
ccdsec_max = 2095

def return_obslist(folder):
    fitslist = sort(glob.glob(folder+"*.fits"))

    nonobjectlist = [
        config["bias"],
        config["flat"],
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

    print "time diff between obs and thar",tdiff*24,"hr"

    thar_averaged = thar[0]*weights[0] + thar[1]*weights[1]

    save(config["folder"]+"/temp/thar.npy",thar_averaged)
    
    return thar_averaged
    
        
def main():
    
    print "making folders"
    try:
        os.stat(config["folder"]+"/temp/")
    except:
        os.mkdir(config["folder"]+"/temp/")    
    try:
        os.stat(config["folder"]+"/reduced/")
    except:
        os.mkdir(config["folder"]+"/reduced/")

    if config["delete_prev"]:
        os.system("rm -rf "+config["folder"]+"/temp/*")
        os.system("rm -rf "+config["folder"]+"/reduced/*")
    
    if not os.path.exists(config["folder"]+"/temp/masterbias.fits"):
        print "create master bias"
        masterbias = calibrate.average_header_fits(config["folder"],config["bias"],config["folder"]+"/temp/masterbias.fits")[:,ccdsec_min:ccdsec_max]

    else:
        masterbias = pyfits.getdata(config["folder"]+"/temp/masterbias.fits")[:,ccdsec_min:ccdsec_max]

    if not os.path.exists(config["folder"]+"/temp/masterflat.fits"):
        print "create master flat"
        masterflat = calibrate.average_header_fits(config["folder"],config["flat"],config["folder"]+"/temp/masterflat.fits")[:,ccdsec_min:ccdsec_max]
        masterflat -= masterbias

    else:
        masterflat = pyfits.getdata(config["folder"]+"/temp/masterflat.fits")[:,ccdsec_min:ccdsec_max]
        masterflat -= masterbias

    print "determine order masks"
    if os.path.exists(config["folder"]+"/temp/order_masks.pkl"):
        order_masks = pickle.load(open(config["folder"]+"/temp/order_masks.pkl","rb"))
    else:
        order_masks = mask_orders.return_masks(masterflat)

    masterflat_extracted = extract_order.mask_order(masterflat,order_masks)


    ### determine a list of observations that need extracting
    tharlist = return_tharlist(config["folder"])
    obslist = return_obslist(config["folder"])
    
    for i in range(len(obslist[0])):
        fitsname = os.path.basename(obslist[0][i])
        if not os.path.exists(config["folder"]+"/temp/"+fitsname+".spec.pkl"): ### check whether obs was reduced
            print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
            print "Reducing",obslist[0][i]
            thar = loadthar(obslist[1][i],tharlist,masterbias)

            #print "Calculating shear"
            #shear_list = thar_straighten.find_straighten_function(thar,order_masks)
            #print "Applying shear"
            #thar_shear = thar_straighten.shear_obs(extract_order.mask_order(thar,order_masks),shear_list)
            thar_shear = extract_order.mask_order(thar,order_masks)
            fits = pyfits.getdata(obslist[0][i])[:,ccdsec_min:ccdsec_max]
            fits -= masterbias
            fits = extract_order.mask_order(fits,order_masks)

            import copy
            fits_noflat = copy.deepcopy(fits)
            print "Dividing by flat"
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

            trace_array = extract_order.find_trace(fits_noflat_shear)
            plt.savefig(config["folder"]+"/temp/"+fitsname+"trace.pdf")
            plt.clf()
            pickle.dump(trace_array,open(config["folder"]+"/temp/"+fitsname+".trace.pkl","wb"))

            print "extracting spectra"
            spectrum,background = extract_order.extract_trace(fits_shear,trace_array)
            spectrum_noflat,background_noflat = extract_order.extract_trace(fits_noflat_shear,trace_array)
            tharspec,bk = extract_order.extract_trace(thar_shear,trace_array)
            for i in range(len(tharspec)):
                tharspec[i] += bk[i]

            pickle.dump([spectrum,background,tharspec,spectrum_noflat,background_noflat],open(config["folder"]+"/temp/"+fitsname+".spec.pkl","wb"))
            
            #sys.exit()

    print "Creating fits files"
    average_adjacent_obs.average_adjacent_obs(obslist,tharlist,config["folder"])
    wavecal_all.main(config["folder"])

    if config["run_analysis"]:

        print "Running spectral analysis"
        
        import run_analysis
        run_analysis.main(config["folder"])
    
if __name__ == "__main__":

    main()
    
