from numpy import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import optimize
import os,sys,string
import emcee
import pyfits

def chebyshev(x0,peaklist,shape):
    x0 = reshape(x0,shape)
    
    peaklist[:,0] += 34.145
    
    ### Follow buchhave's thesis

    xmax = 1000
    xmin = -1000 ### pixels
    xnorm = (2*peaklist[:,1] - (xmax+xmin)) / (xmax-xmin) ### normalise the x-axis

    order_min = 0+34
    order_max = 20+34
    order_norm = (2*peaklist[:,0] - (order_max+order_min)) / (order_max-order_min) ### normalise the order axis

 
    Pm_list = [1,xnorm]
    Pn_list = [1,order_norm]
    
    for m in range(2,len(x0)):
        Pm_list.append(2*xnorm*Pm_list[m-1]-Pm_list[m-2]) ### Chebyshev
        #Pm_list.append(((2*(m-1)+1)*xnorm*Pm_list[m-1] - m*Pm_list[m-2])/(m)) ### Legendre
        
    for n in range(2,len(x0[0])):
        Pn_list.append(2*order_norm*Pn_list[n-1]-Pn_list[n-2]) ### Chebyshev
        #Pn_list.append(((2*(n-1)+1)*order_norm*Pn_list[n-1] - n*Pn_list[n-2])/(n)) ### Legendre
                       
        
    ### compute polynomial
    f = zeros(len(peaklist))
    for m in range(len(x0)):
        for n in range(len(x0[0])):
            f += x0[m,n] * Pm_list[m] * Pn_list[n]
            

    ### finally, compute to wavelength
    wave = f / peaklist[:,0]


    return wave



def fit_chebyshev_lstsq(x0_init,peaklist,shape):

    # x0 = optimize.least_squares(minfunc,x0_init).x


    niter = 5
    clip = 0.9
    i = 0
    mask = peaklist[:,0]==peaklist[:,0]
    while i < niter:
        def minfunc(x0):
            wave = chebyshev(x0,peaklist.copy()[mask],shape)
            lstsq = (wave-peaklist[:,2][mask])**2
            return lstsq
        
        x0 = optimize.least_squares(minfunc,x0_init).x

        wave = chebyshev(x0,peaklist.copy()[mask],shape)
        lstsq = (wave-peaklist[:,2][mask])**2
        
        lstsq_percentile = sort(lstsq)[int(clip*len(lstsq))]

        wave = chebyshev(x0,peaklist.copy(),shape)
        lstsq = (wave-peaklist[:,2])**2
        mask = lstsq < lstsq_percentile ### remove outliers
        print len(lstsq[mask])
        i += 1

    
    savetxt("chebyshev_solution",reshape(x0,shape))
    
    return x0
    

def nice_plot(x0,peaklist,shape):

    plt.figure(figsize=(8,10))
    
    gs = gridspec.GridSpec(4, 1)
    ax = plt.subplot(gs[:3, 0])
    
    plt.scatter(peaklist[:,1],peaklist[:,2])

    orders = arange(min(peaklist[:,0]),max(peaklist[:,0])+1)

    for o in orders:
        xpos = arange(-1000,1000)
        input_list = transpose(array([o*ones(len(xpos)),xpos]))
        wave = chebyshev(x0,input_list,shape)
        plt.plot(xpos,wave,"r-")


    plt.xlabel("Pixel",fontsize=15,weight="black")
    plt.ylabel("Wavelength (A)",fontsize=15,weight="black")

    plt.xlim(-1000,1000)
    plt.ylim(min(peaklist[:,2])-500,max(peaklist[:,2])+500)

    ax = plt.subplot(gs[3, 0])

    model = chebyshev(x0,peaklist,shape)
    plt.scatter(peaklist[:,2],peaklist[:,2]-model,color="k",s=5)

    plt.text(0.05,0.8,"$\sigma="+str(round(std(model-peaklist[:,2]),5))+"$",fontsize=18,transform=ax.transAxes,ha="left")

    plt.xlabel("Wavelength",fontsize=15,weight="black")
    plt.ylabel("Wavelength Resdual (A)",fontsize=15,weight="black")

    plt.ylim(-0.2,0.2)
    plt.xlim(min(peaklist[:,2]),max(peaklist[:,2]))

def apply_solution(x0,input_spectrum,shape):

    wave_solution = []
    for o in range(len(input_spectrum)):
        spectrum = input_spectrum[0]
        xpos = arange(len(spectrum))-len(spectrum)/2
        input_list = transpose(array([o*ones(len(xpos)),xpos]))
        wave = chebyshev(x0,input_list,shape)
        wave_solution.append(wave)

    hdu_wave = pyfits.ImageHDU(array(wave_solution))

    return hdu_wave
        
    
if __name__ == "__main__":
    peaklist = loadtxt("/media/Onion/Data/ANU23echelle/20180705/temp/RAWSPEC_HD100623_2018-07-05T08-51-31.844.fits.tharpeaks")
    #peaklist = loadtxt("/media/Onion/Data/ANU23echelle/20180703/temp/RAWSPEC_WASP167_2018-07-03T13-36-30.516.fits.tharpeaks")
    #peaklist = loadtxt("thar_peaks")
    mask = peaklist[:,0] != 0
    mask *= peaklist[:,1] != 0
    #mask *= peaklist[:,2] > 4000
    peaklist = peaklist[mask]

    ### ndeg-1 is the degree of the polynomial (counting from 0)
    ndeg_wave = 5 ### degree of cheb for the wave axis
    ndeg_order = 5 ### degree of cheb for the order axis 

    x0 = ones((ndeg_wave,ndeg_order))*0.5
    x0[0,0] = 5200.
    x0[0,1] = -1200.

    x0 = x0.flatten()    
    x0 = fit_chebyshev_lstsq(x0,peaklist,(ndeg_wave,ndeg_order))
    
    nice_plot(x0,peaklist,(ndeg_wave,ndeg_order))
    plt.show()
    ### Apply solution to spectrum
    #apply_solution(x0,"HD189625_spec.fits",(ndeg_wave,ndeg_order))
