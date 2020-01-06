import os,sys,glob,string#,pyfits
from astropy.io import fits as pyfits
from numpy import *
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord

# source /data/mash/marusa/reductions_echelle/venv_echelle/bin/activate

import cc_rv
import lsd
import plot_spectra_anl
import spectype
#from astroquery.vizier import Vizier
from astroquery.gaia import Gaia

import prepare_synthetic_templates

import warnings
warnings.filterwarnings("ignore")

#~ speclib = "/media/Onion/Data/spectral_library/lib/"
speclib = "templates/"


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





def roundnearest(val,r):
    val = val / r
    val = round(val)
    val = val*r
    return val

# def getteff(c):
#     gaia = Vizier(columns=["*"], catalog="I/345")
#     radius = 5
#     found = False
#     while not found:
#        results = gaia.query_region(c, radius=str(radius)+"s")
#        if len(results) > 0:
#           if len(results[0]) > 0:
#               print "Nstars found",len(results[0])
#               gmag = []
#               for i in range(len(results[0])):
#                   gmag.append(results[0][i]["Gmag"])
#               print gmag
#               print max(gmag)

#               if min(gmag) < 15:

#                  gmag = array(gmag)
#                  indx = argmin(gmag)

#                  results = results[0][indx]

#                  teff = results["Teff"]
#                  rstar = results["Rad"]
#                  gmag = results["Gmag"]
#                  found = True
#               else:
#                  radius += 5
#           else:
#               radius += 5
                 

#        else:
#           radius += 5
              
#     print teff,rstar,gmag 
#     return teff,rstar,gmag


def estimate_teff_from_bp_rp(bp_rp):
    # Estimate temperature from the colour
    z = [33.14962128, -476.14587496, 2546.45692104, -6360.12338185, 9466.73218724]
    p = poly1d(z)
    teff = p(bp_rp)
    return teff

def getteff(c,cone=10):

    G = nan
    while G != G and cone < 180:
       print "cone search radius",cone
       #if True:
       try:
           width = u.Quantity(cone, u.arcsec)
           height = u.Quantity(cone, u.arcsec)
           result = Gaia.cone_search_async(c,width)
           result =  result.get_results()
           #print result.columns

           useindx = 0
           print "len result",len(result)
           if len(result) > 1:
               Gmag = []
               for i in range(len(result)):
                   Gmag.append(result[i]["phot_g_mean_mag"])
               Gmag = array(Gmag)
               print Gmag
               useindx = argmin(Gmag)

           G =  result[useindx]["phot_g_mean_mag"]
           teff =  result[useindx]["teff_val"]
           rstar = result[useindx]["radius_val"]

       except IndexError:
           G = nan
           teff = nan
           rstar = nan

       cone += 20
        
    return teff,rstar,G






def writerv(rvlist,delim=","):
    o = ""
    for i in rvlist:
        o += str(i)+delim

    return o

def measure_snr(fitsname):
   fits = pyfits.open(fitsname)
   data = fits[3].data
   snrlist = []
   for order in data:
      snr_order = nanmedian(order[len(order)/2-100:len(order)/2+100])
      snr_order = snr_order/sqrt(snr_order)
      snr_order *= 2.5
      snrlist.append(snr_order)

   return nanmax(array(snrlist))

 
#~ def main(folder, photometry):
def main(folder):
    """
    Photometry should be a dictionary with objectnames and BP-RP
    """
    fitslist = sort(glob.glob(folder+"/reduced/*.fits"))
    print('folder', folder+"/reduced/*.fits")
    print('FITSLIST', fitslist)
    for fits in fitslist:
    #for fits in [fitslist[1]]:

        if not os.path.exists(fits+".rv"):
           print fits

           fitsheader = pyfits.getheader(fits)
           objectname = fitsheader["OBJECT"]
           ra,dec = fitsheader["RA"],fitsheader["DEC"]
           hjd = fitsheader["HJD"]
           bcorr = fitsheader["BCORR"]
           exptime = fitsheader["EXPTIME"]
           print 'objectname', objectname

           c = SkyCoord(ra+" "+dec, frame='icrs', unit=(u.hourangle, u.deg))

           # Find BP-RP
           #~ try:
               #~ bp_rp = photometry[objectname]
           #~ except:
               #~ bp_rp = None
               # Now what?

           try:
               teff,rstar,gmag = getteff(c)
               #~ teff, _, _ = estimate_teff_from_bp_rp(bp_rp)
               #~ teff = estimate_teff_from_bp_rp(bp_rp)
               print("Gaia teff,rstar,gmag",teff,rstar,gmag)
               #~ print("Estimated photometric teff",teff#,rstar,gmag
               teff = int(roundnearest(teff,250))
               rstar = int(roundnearest(rstar,0.5))
               gmag = float(gmag)

           except (ValueError,IndexError):
               print "Error on teff"
               teff = 6000
               rstar = 1.0
               gmag = -99.


           #~ if teff < 4000:
               #~ teff = 4000
           #~ if teff > 10000:
               #~ teff = 10000

            
           # Create a template with this temperature
           #~ template = '../template.dat'
           #~ w, flux = prepare_synthetic_templates.get_spectrum(teff)
           #~ fle = open(template, 'wb')
           #~ print w
           #~ print flux
           #~ print len(w)
           #~ print len(flux)
           #~ for ww, ff in zip(w, flux):
               #~ fle.write('%f %f\n'%(ww, ff))
           #~ fle.close()

           template = speclib+"template_"+str(teff)+"_4.0_0.0.dat"

           snr = measure_snr(fits)
           ccf,vsini,lsdshift = lsd.run_spectrum(fits,template)

           #~ if rstar < 3:
              #~ logg = 3.5
           #~ else:
              #~ logg = 2.0 ### initial logg minimum. If rstar < 3, allow all logg, else only allow dwarfish logg
              
           logg = 5.0 # MZ hack

           teff,logg,feh,vsini,lsdshift = spectype.get_best_template(fits, os.path.join(folder, "reduced/lsd_"+os.path.basename(fits)+".pkl"),teff,logg)
           template = speclib+"template_"+str(int(teff))+"_"+str(logg)+"_"+str(feh)+".dat"
           
           #~ template = os.path.join(speclib, 'template_%d.dat'%int(teff)) # MZ. Only temperature for now
           
           #~ print 'TEMPLATE', template
           
           if not os.path.exists(template):
              print("Template doesn't exist, trying a logg=4.5,feh=0.0 template")
              feh,logg = 0.0,4.5
              template = speclib+"template_"+str(int(teff))+"_"+str(logg)+"_"+str(feh)+".dat"
              print("using template",template)

           rvlist,ccfrv = cc_rv.main(fits,template,vsini=vsini)
           rvlist[:,1] += bcorr
           ccfrv += bcorr
           telrv = cc_rv.measure_telluric_rv(fits,template,vsini,lsdshift)

           print("CCF RV",ccfrv)

           print("diagnostic",fits,template,int(teff),logg,feh,vsini,lsdshift)

           plot_spectra_anl.main(fits,template,int(teff),logg,feh,vsini,lsdshift)
           rvlist = [objectname,fits,hjd,convHMS(ra),convDMS(dec),gmag,exptime,teff,logg,feh,snr,bcorr,telrv,vsini,lsdshift+bcorr,ccfrv]+list(rvlist[:,1])
           rvlist = writerv(rvlist)+"\n"
           out = open(fits+".rv","w")
           out.write(rvlist)
           out.close()


           plt.clf()

        
if __name__ == "__main__":

    #~ p=loadtxt('../observed_id_color.dat', dtype=str)

    #~ photometry=dict()
    #~ for x in p:
        #~ photometry[x[2]]=float(x[0])


    import imp
    config_filename = sys.argv[1]
    print('CONFIG FILENAME', config_filename)
    config_file = imp.load_source(config_filename.replace('.py', ''), config_filename)

    #~ import config_file
    config = config_file.set_config()

    folder = config["folder"]#
    #folder = "/media/Onion/Data/ANU23echelle/20181115/bin2/"
    #~ main(folder, photometry)
    main(folder)
