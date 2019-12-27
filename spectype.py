import os,sys,string,pickle#,pyfits
from astropy.io import fits
import emcee
from numpy import *
from scipy import interpolate,optimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# Commented: Marusa
# import config_file
#config = config_file.set_config()

#~ library_path = os.path.join(config["spectral_library"]+"lib/")
library = pickle.load(open('synthetic_library_for_echelle_rv.pkl', "rb"))


c = 3.*10**5

def make_rot_prof(v,vsini=100,epsilon=0.6451,scale=1,vshift=0):
    
    vrot = 2*(1-epsilon)*(1-((v-vshift)/vsini)**2)**(0.5) + 0.5*pi*epsilon*(1-((v-vshift)/vsini)**2)
    #vrot /= pi*vsini*(1-epsilon/3)


    mask = vrot != vrot
    vrot[mask] = 0

    vrot /= max(vrot)
    vrot *= scale
    
    return vrot

def vgauss(v,vrot,macro=0,scale=1):
    gauss = exp(-v**2/(2*macro**2))
    gauss /= sum(gauss)
    vrot = convolve(vrot,gauss,mode="same")
    #vrot /= max(vrot)
    #vrot *= scale

    return vrot

def make_vrot(x0,vel):
    vrot = make_rot_prof(vel,vsini=x0[0],scale=x0[1],vshift=x0[2])
    vrot = vgauss(vel,vrot,macro=13.04/2.355,scale=1.)
        
    return vrot

def fitlsd(lsd,vsini_init=40):
    lsd = interpolate.splrep(lsd[:,0],lsd[:,1])
    vel = arange(-500,500,0.1)
    lsd = interpolate.splev(vel,lsd)
    lsd = transpose(array([vel,lsd]))
    

    x0 = [vsini_init,max(lsd[:,1]),lsd[:,0][argmax(lsd[:,1])]]
    def minfunc(x0):
        f = make_vrot(x0,lsd[:,0])
        lstsq = sum((f-lsd[:,1])**2) 
        #print x0,lstsq
        return lstsq
    x0 = optimize.fmin(minfunc,x0)
    return x0

def calc_noise(vel,ccf,vsini):
    v0 = vel[argmax(ccf)]
    mask = abs(vel-v0) > vsini
    return std(ccf[mask])

def fitlsd_mcmc(lsd,vsini_init=40,nrun=500):
    print('LSD0', lsd)
    lsd = interpolate.splrep(lsd[:,0],lsd[:,1])
    vel = arange(-500,500,0.1)
    lsd = interpolate.splev(vel,lsd,ext=1)
    lsd = transpose(array([vel,lsd]))

    def vrot_minfunc(x0,vel,ccf,noise):
        if x0[0] > 500 or x0[1] < 0 or abs(x0[2]) > 500:
            return -1*inf
        else:
            
            f = make_vrot(x0,vel)
            #plt.plot(vel,ccf)
            #plt.plot(vel,f)
            #plt.show()
            chisq = sum(((f-ccf)/noise)**2)

            if chisq != chisq:
                chisq = -1*inf
            #print chisq,x0
            return -0.5*chisq

    noise = calc_noise(lsd[:,0],lsd[:,1],vsini_init)
    print('LSD', lsd)
    maxheight = max(lsd[:,1])
    nwalkers = 20

    vmean = lsd[:,0][argmax(lsd[:,1])]
    
    p0 = []
    for i in range(nwalkers):
        print('fitlsd_mcmc', random.normal(vsini_init,0.05*vsini_init))
        print('second', maxheight, random.normal(maxheight,0.1*maxheight))
        print('third', random.normal(vmean,1.))
        x0_i = [random.normal(vsini_init,0.05*vsini_init),random.normal(maxheight,0.1*maxheight),random.normal(vmean,1.)]
        p0.append(x0_i)
    ndim = len(p0[0])

    sampler = emcee.EnsembleSampler(nwalkers,ndim,vrot_minfunc,args=(lsd[:,0],lsd[:,1],noise),threads=1)

    print('Run emcee...')
    pos,prob,state = sampler.run_mcmc(p0,nrun*2)
    sampler.reset()

    # print "******************************"
    chain = []
    sampler.run_mcmc(pos,nrun,storechain=True)


    chain = sampler.flatchain
    # print "Length of chain",len(chain)

    x0 = []

    
    for i in range(len(chain[0])):
        #hist,bins = histogram(chain[:,i],bins=50)
            
        chain_i = sort(chain[:,i])
        chain_mode = median(chain_i)
        low = chain_i[int(0.18*float(len(chain_i)))]
        high = chain_i[int(0.82*float(len(chain_i)))]

        #plt.axvline(x=chain_mode)
        #plt.axvline(x=low)
        #plt.axvline(x=high)
        low = chain_mode-low
        high = high-chain_mode

        print(chain_mode,low,high)
        
        #plt.hist(chain[:,i],bins=200,histtype="step",color="k")
        #plt.show()
        
        x0.append(chain_mode)


    return x0

def match_template(template_array,spectrum_array,lsd):

    lstsq = 0
    for i in range(len(spectrum_array[0])):

        mask = spectrum_array[1][i] > 5100
        mask *= spectrum_array[1][i] < 5600
        #mask = spectrum_array[1][i] > 5150
        #mask *= spectrum_array[1][i] < 5300
        spectrum_i = transpose(array([spectrum_array[1][i][mask],spectrum_array[0][i][mask]]))
        spectrum_i = spectrum_i[5:-5]

        if len(spectrum_i) > 100:
        
            template_mask = template_array[:,0] > min(spectrum_i[:,0])-10
            template_mask *= template_array[:,0] < max(spectrum_i[:,0])+10
            template_i = template_array[template_mask]

            vel = c*(template_i[:,0]-median(template_i[:,0]))/template_i[:,0]
            lsd_interp = interpolate.splrep(lsd[:,0],lsd[:,1])
            lsd_interp = interpolate.splev(vel,lsd_interp,ext=1)
            lsd_interp /= sum(lsd_interp)

            template_i[:,1] = convolve(template_i[:,1],lsd_interp,mode="same")

            template_i_interp = interpolate.splrep(template_i[:,0],template_i[:,1])
            template_i_interp = interpolate.splev(spectrum_i[:,0],template_i_interp,ext=1)
            mask = template_i_interp != 0

            offset = spectrum_i[:,1]/template_i_interp
            try:
                fit = polyfit(spectrum_i[:,0][mask],offset[mask],10)
                fit = polyval(fit,spectrum_i[:,0])
                template_i_interp *= fit
            except:
                print("bad normalisation")

            # plt.subplot(211)
            # plt.plot(spectrum_i[:,0],spectrum_i[:,1])
            # plt.plot(spectrum_i[:,0],template_i_interp)
            # plt.subplot(212)
            # plt.plot(spectrum_i[:,0],spectrum_i[:,1]-template_i_interp)
            # plt.show()

            diff = spectrum_i[:,1]-template_i_interp
            lstsq += sum(diff**2)

    return lstsq

def plot_template(template_array,spectrum_array,lsd):

    lstsq = 0
    for i in range(len(spectrum_array[0])):

        mask = spectrum_array[1][i] > 5100
        mask *= spectrum_array[1][i] < 5800
        #mask = spectrum_array[1][i] > 5150
        #mask *= spectrum_array[1][i] < 5300
        spectrum_i = transpose(array([spectrum_array[1][i][mask],spectrum_array[0][i][mask]]))
        spectrum_i = spectrum_i[5:-5]
        if len(spectrum_i) > 100:
        
            template_mask = template_array[:,0] > min(spectrum_i[:,0])-10
            template_mask *= template_array[:,0] < max(spectrum_i[:,0])+10
            template_i = template_array[template_mask]

            vel = c*(template_i[:,0]-median(template_i[:,0]))/template_i[:,0]
            lsd_interp = interpolate.splrep(lsd[:,0],lsd[:,1])
            lsd_interp = interpolate.splev(vel,lsd_interp,ext=1)
            lsd_interp /= sum(lsd_interp)


            template_i[:,1] = convolve(template_i[:,1],lsd_interp,mode="same")
            plt.plot(template_i[:,0],template_i[:,1])
            plt.show()

            template_i_interp = interpolate.splrep(template_i[:,0],template_i[:,1])
            template_i_interp = interpolate.splev(spectrum_i[:,0],template_i_interp,ext=1)
            mask = template_i_interp != 0


            plt.plot(spectrum_i[:,0],template_i_interp)
            plt.show()

            offset = spectrum_i[:,1]/template_i_interp

            fit = polyfit(spectrum_i[:,0][mask],offset[mask],10)
            fit = polyval(fit,spectrum_i[:,0])
            template_i_interp *= fit

            plt.subplot(211)
            plt.plot(spectrum_i[:,0],spectrum_i[:,1])
            plt.plot(spectrum_i[:,0],template_i_interp)
            plt.subplot(212)
            plt.plot(spectrum_i[:,0],spectrum_i[:,1]-template_i_interp)
            plt.show()

            diff = spectrum_i[:,1]-template_i_interp
            lstsq += sum(diff**2)

    return lstsq

def get_best_template(spectrum, lsd, teffinit, logginit):
    lsd_list, lsd_master, vsini, shift = pickle.load(open(lsd,"rb"))
    spectrum = fits.open(spectrum)

    spectrum_array = [spectrum[0].data, spectrum[5].data]

    #print "vsini,shift", vsini,shift
    
    print('LSD MASTER', lsd_master)

    x0 = fitlsd_mcmc(lsd_master, vsini_init=vsini)
    print("lsdprof fit", x0)
    lsd_interp = make_rot_prof(lsd_master[:,0], vsini=x0[0], scale=x0[1], vshift=x0[2])
    lsd_interp = vgauss(lsd_master[:,0], lsd_interp, macro=13.04/2.355, scale=1.)
    lsd_interp = transpose(array([lsd_master[:,0], lsd_interp]))

    library_mask = library[0][:,0] >= teffinit - 2000
    library_mask *= library[0][:,0] <= teffinit + 2000
    library_mask *= library[0][:,1] >= logginit
    #~ #library_mask = library[0][:,0] == teffinit 

    bestfit = []
    for i in range(len(library[0][library_mask])):
        template_array = library[1][library_mask][i]
        template_array = transpose(array([library[2], template_array]))
            
        lstsq = match_template(template_array,spectrum_array, lsd_interp)
        #lstsq = match_template(template_array,spectrum_array,lsd_master)
        #print library[0][library_mask][i],lstsq
        bestfit.append([lstsq]+list(library[0][library_mask][i]))


    spectrum.close() # Added by MZ

    bestfit = array(bestfit)
    bestfit = bestfit[argmin(bestfit[:,0])]
    print("best fit spectype", bestfit)
    return bestfit[1],bestfit[2], bestfit[3], x0[0], x0[2]
    #~ return 0,0,0,x0[0],x0[2] # MARUSA

    
if __name__ == "__main__":
    #~ change='marusa'
    spectrum = '/data/mash/marusa/2m3data/echelle/20190514/reduced/ANU23e_13172883+2024199_2019-05-14T12-01-33.542.fits'
    lsd = '/data/mash/marusa/2m3data/echelle/20190514/reduced/lsd_ANU23e_13172883+2024199_2019-05-14T12-01-33.542.fits.pkl'

    teff,logg,feh,vsini,vshift =  get_best_template(spectrum,lsd,7000,3.5)
    print('teff,logg,feh,vsini,vshift', teff,logg,feh,vsini,vshift)
    
    lsd_list,lsd_master,vsini,shift = pickle.load(open(lsd,"rb"))
    spectrum = fits.open(spectrum)

    spectrum_array = [spectrum[0].data,spectrum[5].data]

    #print "vsini,shift", vsini,shift

    x0 = fitlsd_mcmc(lsd_master,vsini_init=vsini)
    print("lsdprof fit", x0)
    lsd_interp = make_rot_prof(lsd_master[:,0],vsini=x0[0],scale=x0[1],vshift=x0[2])
    lsd_interp = vgauss(lsd_master[:,0],lsd_interp,macro=13.04/2.355,scale=1.)
    lsd_interp = transpose(array([lsd_master[:,0],lsd_interp]))


    #teff,logg,feh = 5750,4.5,0.0
    #vsini,shift = 6.740170462845786,34.74406961975362

    
    library_mask = library[0][:,0] == teff
    library_mask *= library[0][:,1]  == logg
    library_mask *= library[0][:,2]  == feh

    
    template_array = library[1][library_mask][0]
    template_array = transpose(array([library[2],template_array]))
            

    plot_template(template_array,spectrum_array,lsd_interp)
    print('DONE')
