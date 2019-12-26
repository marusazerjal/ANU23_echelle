import os,sys,glob,string
from astropy.io import fits as pyfits
from numpy import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import config_file
#config = config_file.set_config()

import identify_order_solarspec
#import wavesolve
import wavesolve2 as wavesolve
import fitarc_chebyshev

import pandas
rv_standards = pandas.read_csv("RV_standard.dat",delim_whitespace=True)

def create_filelist(folder, config=None):
    #print('path', os.path.join(folder, "temp/ANU23e*.fits"))
    raw_fitslist = sort(glob.glob(os.path.join(folder, "temp/ANU23e*.fits")))
    print('raw_fitslist', raw_fitslist)
    print('init_wave_obj', config["init_wave_obj"])
    for fits in raw_fitslist:
        print('THE FILE', pyfits.getheader(fits)["OBJECT"], 'INIT_WAVE_OBJ', config["init_wave_obj"])
        if pyfits.getheader(fits)["OBJECT"].replace(' ', '') == config["init_wave_obj"]:
            init_wave_fits = fits

    try:
        print('HHHHHH', raw_fitslist,init_wave_fits)
    except:
        print('HHHHHH cant print (wavecall_all)')

    try:
        return raw_fitslist,init_wave_fits
    except:
        print("Error, no observation suitable for inital wave calibration, check your init_wave_obj setting in config file")
        sys.exit()

def main(folder, config=None):
    raw_fitslist,init_wave_fits = create_filelist(folder, config=config)

    print("estimating wavelength solution based on stellar spectrum")

    if not os.path.exists(folder+"/temp/order_initial_solutions"):
    #if True:
    
        #~ print rv_standards
    
        init_wave_spec = pyfits.open(init_wave_fits)[0].data
        objectname = pyfits.getheader(init_wave_fits)["OBJECT"].replace(' ', '')
        print 'do wavecall objectname', objectname
        bcorr = pyfits.getheader(init_wave_fits)["BCORR"]
        mask = rv_standards["Star"] == objectname
        if len(rv_standards[mask])>1:
            print('WARNING: %s appears in the RVstandards list more than once!'%objectname)
        print('V_r', rv_standards[mask]["V_r"], 'END')
        truerv = float(rv_standards[mask]["V_r"])
        print("rv standard",objectname,truerv,bcorr)
        
        initial_solutions = identify_order_solarspec.iterate_whole_spectrum(init_wave_spec)
        savetxt(folder+"/temp/init_solutions",array(initial_solutions))
        initial_solutions = loadtxt(folder+"/temp/init_solutions")
        initial_solutions = identify_order_solarspec.fit_echelle_solution_recc(initial_solutions,config["norder"],init_wave_spec,wshift=truerv-bcorr)
        savetxt(folder+"/temp/order_initial_solutions",initial_solutions,fmt="%.3f")

        #os.system("mv initial_wave_solution.pdf "+folder+"/temp/")
        #plt.savefig(folder+"/temp/echelle_equation_fit.pdf")
        #plt.clf()




        
    initial_solutions = loadtxt(folder+"/temp/order_initial_solutions")
    
    for fitsname in raw_fitslist:
        print(fitsname)
        fits = pyfits.open(fitsname)
        arc = fits[2].data
        peaklist = wavesolve.run_spectrum(arc,initial_solutions) ### identifies a list of ThAr peaks
        savetxt(fitsname+".tharpeaks",peaklist,fmt="%.5f")

        #peaklist = loadtxt("thar_peaks")
        
        # MZ:  mask = peaklist[:,0] != 0
        # IndexError: too many indices for array (this happens sometimes)
        # So I added a try-except statement
        try:
            mask = peaklist[:,0] != 0
            mask *= peaklist[:,1] != 0
            peaklist = peaklist[mask]
        except:
            pass
        
        ### ndeg-1 is the degree of the polynomial (counting from 0)
        ndeg_wave = 5 ### degree of cheb for the wave axis
        ndeg_order = 4 ### degree of cheb for the order axis 

        x0 = ones((ndeg_wave,ndeg_order))*0.5
        x0[0,0] = 5200.
        x0[0,1] = -1200.

        x0 = x0.flatten()    
        x0 = fitarc_chebyshev.fit_chebyshev_lstsq(x0,peaklist,(ndeg_wave,ndeg_order))

        plt.clf()
        try:
            fitarc_chebyshev.nice_plot(x0,peaklist,(ndeg_wave,ndeg_order))
        except:
            pass
        plt.savefig(fitsname+".wavecal.pdf")

        wave_hdu = fitarc_chebyshev.apply_solution(x0,fits[0].data,(ndeg_wave,ndeg_order))

        #~ output_name = folder+"/reduced/"+os.path.basename(fitsname)
        output_name = os.path.join(folder, 'reduced', os.path.basename(fitsname))
        os.system("rm "+output_name)
        fits.append(wave_hdu)
        fits.writeto(output_name)


        #sys.exit()
        
if __name__ == "__main__":

    #folder = "/media/Onion/Data/ANU23echelle/20180705/"
    folder = config["folder"]
    main(folder)
