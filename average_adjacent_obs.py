import os,sys,string,pickle,glob
from numpy import *
import matplotlib.pyplot as plt
import astropy.time
from PyAstronomy import pyasl

from astropy.io import fits as pyfits
import matplotlib
matplotlib.use('Agg')




# Convert HH:MM:SS.SSS into Degrees :
def convHMS(ra):
   try :
      sep1 = ra.find(':')
      hh=int(ra[0:sep1])
      sep2 = ra[sep1+1:].find(':')
      mm=int(ra[sep1+1:sep1+sep2+1])
      ss=float(ra[sep1+sep2+2:])
   except:
      raise
   else:
      pass
   
   return(hh*15.+mm/4.+ss/240.)

# Convert +DD:MM:SS.SSS into Degrees :
def convDMS(dec):

   Csign=dec[0]
   if Csign=='-':
      sign=-1.
      off = 1
   elif Csign=='+':
      sign= 1.
      off = 1
   else:
      sign= 1.
      off = 0

   try :
      sep1 = dec.find(':')
      deg=int(dec[off:sep1])
      sep2 = dec[sep1+1:].find(':')
      arcmin=int(dec[sep1+1:sep1+sep2+1])
      arcsec=float(dec[sep1+sep2+2:])
   except:
      raise
   else:
      pass

   return(sign*(deg+(arcmin*5./3.+arcsec*5./180.)/100.))




def combine_order(order_specs):
    medvals = nanmedian(order_specs,axis=1)
    normed_spec = order_specs.copy()/transpose([medvals])
    init_median = median(normed_spec,axis=0)

    
    errors = sqrt(order_specs)/order_specs
    for i in range(len(normed_spec[0])):

        ### if val deviates 2x from photon error, reject and replace with median
        mask = abs(normed_spec[:,i]-init_median[i]) > 2*errors[:,i]        
        order_specs[:,i][mask] = init_median[i]*medvals[mask]
        
    return sum(order_specs,axis=0)

def return_header_numbers(header_list,header_name):
    header_number = []
    for header in header_list:
        header_number.append(header[header_name])

    return header_number

def combine_header(header_list):
    new_header = header_list[0]

    new_header["MJD-OBS"] = mean(return_header_numbers(header_list,"MJD-OBS"))
    new_header["AIRMASS"] = mean(return_header_numbers(header_list,"AIRMASS"))
    new_header["EXPTIME"] = sum(return_header_numbers(header_list,"EXPTIME"))

    jd = new_header["MJD-OBS"]+2400000.5
    t = astropy.time.Time(jd, format='jd', scale='utc')
    new_header["DATE-OBS"] = str(t.fits)

    # Coordinates of European Southern Observatory
    # (Coordinates of UT1)
    longitude = new_header["LONG-OBS"]
    latitude = new_header["LAT-OBS"]
    altitude = new_header["ALT-OBS"]

    # Coordinates of HD 12345 (J2000)
    ra2000 = convHMS(new_header["RA"])
    dec2000 = convDMS(new_header["DEC"])

    # Calculate barycentric correction (debug=True show
    # various intermediate results)
    corr, hjd = pyasl.helcorr(longitude, latitude, altitude, \
                              ra2000, dec2000, jd, debug=False)
    
    print jd,corr,hjd
    
    new_header.set('HJD',hjd)
    new_header.set('JD',jd)
    new_header.set('BCORR',corr)
    new_header.set('AP0',"Spectrum FFdiv+BKsub")
    new_header.set('AP1',"Background FFdiv")
    new_header.set('AP2',"ThAr")
    new_header.set('AP3',"Spectrum BKsub")
    new_header.set('AP4',"Background")
    new_header.set('AP5',"Wavelength")

    return new_header,string.replace(string.replace(str(t.fits),"(UTC)",""),":","-")
    

def average_adjacent_obs(obslist,tharlist,folder):

    obslist_jd = array(obslist[1])
    
    for i in range(len(tharlist[0])-1):
        obsindx = []
        objname = []
        for j in range(len(obslist[1])):
            fitsname = os.path.basename(obslist[0][j])

            ### return spectra taken between these two thars
            
            if obslist[1][j] > tharlist[1][i] and obslist[1][j] < tharlist[1][i+1] and os.path.exists(folder+"/temp/"+fitsname+".spec.pkl"):
                objname.append(pyfits.getheader(obslist[0][j])["OBJECT"])
                obsindx.append(j)


        if len(obsindx) > 0:

            objectlist = unique(objname)
            print objectlist

            for obj in objectlist:
                
                spectrum_list = []
                background_list = []
                thar_list = []
                spectrum_noflat_list = []
                background_noflat_list = []
                header_list = []
                
                for i in range(len(objname)): ### read in all the spectra of that object
                    if objname[i] == obj:
                        fitsname = os.path.basename(obslist[0][obsindx[i]])
                        spec = pickle.load(open(folder+"/temp/"+fitsname+".spec.pkl","rb"))
                        header_list.append(pyfits.getheader(obslist[0][obsindx[i]]))


                        spectrum_list.append(spec[0])
                        background_list.append(spec[1])
                        thar_list.append(spec[2])
                        spectrum_noflat_list.append(spec[3])
                        background_noflat_list.append(spec[4])

                header,fitsTIME = combine_header(header_list)
                print obj,len(spectrum_list)
                        
                if len(spectrum_list) == 1:
                    spectrum_master,background_master,thar_master,spectrum_noflat_master,background_noflat_master = spectrum_list[0],background_list[0],thar_list[0],spectrum_noflat_list[0],background_noflat_list[0]

                else:
                    spectrum_master,background_master,thar_master,spectrum_noflat_master,background_noflat_master = [],[],[],[],[]
                    for order in range(len(spectrum_list[0])):

                        order_specs = []
                        for s in background_list:
                            order_specs.append(s[order])
                            #plt.plot(s[order])
                            #plt.show()
                        bk_summed = sum(array(order_specs),axis=0)
                        background_master.append(bk_summed)


                        order_specs = []
                        for s in spectrum_list:
                            order_specs.append(s[order])

                        spec_summed = combine_order(array(order_specs))
                        spectrum_master.append(spec_summed)#-bk_summed)
            
        
                        order_specs = []
                        for s in thar_list:
                            order_specs.append(s[order])
                        spec_summed = mean(array(order_specs),axis=0)
                        thar_master.append(spec_summed)


                        order_specs = []
                        for s in background_noflat_list:
                            order_specs.append(s[order])
                        bk_summed = sum(array(order_specs),axis=0)
                        background_noflat_master.append(bk_summed)
                        
                        order_specs = []
                        for s in spectrum_noflat_list:
                            order_specs.append(s[order])

                        spec_summed = combine_order(array(order_specs))
                        spectrum_noflat_master.append(spec_summed)#-bk_summed)
            
        
                ### perform checks on spectrum for nan and infs
                def nanchecks(spec):
                   for i in range(len(spec)):
                      mask = spec[i] != spec[i]
                      mask += abs(spec[i]) == inf
                      
                      try:
                          print i, len(spec[i]), spec[i], mask, ~any(mask)
                          #~ print
                          indx = arange(len(spec[i]))
                          
                          if sum(mask) > 0 and sum(mask)<len(spec[i]):
                             print 'sum(mask)', sum(mask)
                             for j in indx[mask]:
                                adjacent_indx = abs(indx-j)
                                adjacent_indx_mask = adjacent_indx > 0
                                #~ print adjacent_indx
                                #~ print len(adjacent_indx)
                                #~ print any(invert(mask))
                                #~ print adjacent_indx[invert(mask)]
                                adjacent_indx_mask *= adjacent_indx < min(adjacent_indx[invert(mask)]) + 10 # this line crashes
                                fixval = nanmean(spec[i][adjacent_indx_mask])
                                if fixval != fixval or abs(fixval) == inf:
                                   fixval = nanmedian(spec[i])
                                spec[i][j] = fixval
                          elif sum(mask)==len(spec[i]):
                                print 'sum(mask)', sum(mask), 'setting to 0'
                                spec[i] = zeros(len(spec[i]))
                      except:
                          print 'len0', 'setting to 0'
                          spec[i] = zeros(len(spec[i]))
                      print


                   return spec


                  
                spectrum_master = nanchecks(spectrum_master)
                background_master = nanchecks(background_master)
                thar_master = nanchecks(thar_master)
                spectrum_noflat_master = nanchecks(spectrum_noflat_master)
                background_noflat_master = nanchecks(background_noflat_master)

                

                hduspec = pyfits.PrimaryHDU(array(spectrum_master),header=header)
                hdubk = pyfits.ImageHDU(array(background_master))
                hduthar = pyfits.ImageHDU(array(thar_master))
                hduspec_nf = pyfits.ImageHDU(array(spectrum_noflat_master))
                hdubk_nf = pyfits.ImageHDU(array(background_noflat_master))
                hdulist = pyfits.HDUList([hduspec,hdubk,hduthar,hduspec_nf,hdubk_nf])

                output_name = folder+"/temp/"+"ANU23e_"+obj+"_"+fitsTIME+".fits"
                os.system("rm "+output_name)
                hdulist.writeto(output_name,overwrite=True)

                
if __name__ == "__main__":
    import main

    folder = "/media/Onion/Data/ANU23echelle/20190121/"
    obslist = main.return_obslist(folder)
    print obslist
    tharlist = main.return_tharlist(folder)
    print tharlist
    average_adjacent_obs(obslist,tharlist,folder)
