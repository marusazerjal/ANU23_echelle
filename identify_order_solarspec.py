import os,sys,string
from numpy import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from scipy import interpolate,optimize,signal,stats
from PyAstronomy import pyasl

c = 3*10**5.

solar_template = loadtxt("sky_template.dat")
template_interp = interpolate.splrep(solar_template[:,0],solar_template[:,1],k=1) ### interpolate 


def echelle_equation(x0,orders):
    ### x0[0] is the 1st order blaze lambda
    ### x0[1] is the offset of our extracted order vs real order 
    return x0[0]/(orders+x0[1])

def echelle_equation_fit_sigclip(x,y,order=2,clip=1,niter=3): ### clip in sigma
    
    def minfunc(x0,x,y):
        f = echelle_equation(x0,x)
        return sum((f-y)**2)

    ### do a first fit with a median filter
    fit = optimize.fmin(minfunc,[100000,10.],args=(x,y))
    
    #print "length of original array",len(x)
    mask = abs(echelle_equation(fit,x)-y) < clip * std(echelle_equation(fit,x)-y)

    for i in range(niter):
        fit = optimize.fmin(minfunc,fit,args=(x[mask],y[mask]))
        mask *= abs(echelle_equation(fit,x)-y) < clip * std(echelle_equation(fit,x[mask])-y[mask])

        #print "length of poly fit clipped array",len(x[mask])
    fit = optimize.fmin(minfunc,fit,args=(x[mask],y[mask]))

    stdev = std(echelle_equation(fit,x[mask])-y[mask])
    return fit,stdev


    

    

def normalise(spec,niter=2,sigma_low = 0.05,deg=5):

    ### normalise the spectrum to 1
    
    #~ print 'spec', spec
    #~ x = arange(len(spec)) # MZ: commented out
    mask = spec == spec
    spec_iter = spec[mask]
    print 'mask', mask
    print 'spec_iter', spec_iter
    x = arange(len(spec_iter)) # MZ: added this line
    x_iter = x.copy()

    i = 0
    while i < niter:
        if len(x_iter)>0: # MZ
            print 'identify_order_solarspec: x_iter', x_iter, any(x_iter), len(x_iter)
            print 'identify_order_solarspec: spec_iter', spec_iter, any(spec_iter), len(spec_iter)
            fit = polyfit(x_iter,spec_iter,deg)
            fit = polyval(fit,x_iter)
            
            #~ for x, y in zip(x_iter, spec_iter):
                #~ print x, y

            print 'fit', fit, any(fit)
            mask = spec_iter - fit > sigma_low * std(spec_iter-fit)
            #~ print 'mask', any(mask)
            spec_iter = spec_iter[mask]
            x_iter = x_iter[mask]
        else:
            pass
        i += 1

    fit = polyfit(x_iter,spec_iter,deg)
    fit = polyval(fit,x)
    max_prefit = max(fit)
    spec -= fit
    spec += max_prefit

    mask = spec == spec
    spec_iter = spec[mask]
    x_iter = x[mask]

    fit = polyfit(x_iter,spec_iter,1)
    fit = polyval(fit,x_iter)

    mask = spec_iter - fit > sigma_low * std(spec_iter-fit)
    spec_iter = spec_iter[mask]
    x_iter = x_iter[mask]

    fit = polyfit(x_iter,spec_iter,1)
    fit = polyval(fit,x)
    spec /= fit

    mask = spec > 1.1
    spec[mask] = 1.

    #plt.plot(spec)
    #plt.show()
    
    return spec



def cross_correlation_to_find_sol(spec,testwave_centre=6500,testwave_width=300,testdelta_centre=0.09,testdelta_width=0.005,toplot=False):

    ### Cross correlate order to solar spectrum

    #~ print 'spec1', spec
   
    
    xpos = arange(len(spec))-len(spec)/2
    mask = xpos > min(xpos)+150
    mask *= xpos < max(xpos)-150

    print 'MZ: If this fails then the spectrum might not look like a spectrum but something weird. Does "trace order" work?'

    spec -= nanmin(spec[spec!=-inf])
    #~ print 'spec2', spec
    spec = normalise(spec[mask],deg=10)
    xpos = xpos[mask]

    
    def errfunc(x0,return_ccf=False):
        wave = polyval(x0,xpos)
        spec_i = -1*(spec-1)

        wave_template = arange(min(wave)-150,max(wave)+150,abs(x0[0]))

        template_i = interpolate.splev(wave_template,template_interp,ext=1) ### interpolate the template spectrum
        template_i = normalise(template_i) ### comment out this step if using normalised template spectrum
        template_i = -1*(template_i-1)

        #spec_i /= max(spec_i)
        #template_i /= max(template_i)

        mask = wave_template == wave_template


        drv,cc = pyasl.crosscorrRV(wave, spec_i, wave_template[mask], template_i[mask], -5000, +5000, 1, mode='doppler', skipedge=0, edgeTapering=None)

        peak_loc = drv[argmax(cc)]
        mask = abs(drv-peak_loc) > 1000
        fit = polyfit(drv[mask],cc[mask],2)
        fit = polyval(fit,drv)
        cc -= fit

        snr = (max(cc)-median(cc[mask]))/std(cc[mask])

        if return_ccf:
            return drv,cc
        else:
            print(x0,snr)
            #return -1*max(cc)
            return -1*snr


    ### loop around for a grid search
    wave0 = arange(testwave_centre-testwave_width,testwave_centre+testwave_width+1,100.)
    deltaA = arange(testdelta_centre-testdelta_width,testdelta_centre+testdelta_width+0.001,0.001) ### the tested A/pixel scale

    ww,aa = meshgrid(wave0,deltaA)
    errarray = ones(shape(ww))

    for i in range(len(ww)):
        for j in range(len(ww[0])):
            errarray[i,j] = errfunc([aa[i,j],ww[i,j]])


    mask = errarray == min(errarray.flatten())
    deltaA_init,wave0_init = aa[mask],ww[mask]


    print("minimum found for grid search",deltaA_init,wave0_init,min(errarray.flatten()))

    ### use the cc peak to find the correct wave0_init
    x0 = [deltaA_init,wave0_init]
    
    drv,cc = errfunc(x0,return_ccf=True)
    rvshift = drv[argmax(cc)]
    waveshift = (rvshift/c)*x0[1][0]
    x0[1] -= waveshift
    drv,cc = errfunc(x0,return_ccf=True)
    rvshift = drv[argmax(cc)]
    waveshift = (rvshift/c)*x0[1][0]
    x0[1] -= waveshift

    print("with the correct wavelength shift",x0)


    if toplot:

        plt.figure(figsize=(10,10))
        plt.subplots_adjust(hspace=0.25)
        ax = plt.subplot(211)
        
        drv,cc = errfunc(x0,return_ccf=True)
        plt.plot(drv,cc,"k-")
        plt.xlabel("RV (km/s)",fontsize=15)
        plt.xlim(-5000,5000)
        [i.set_linewidth(3.) for i in ax.spines.itervalues()]

        ax = plt.subplot(212)
        
        wave = polyval(x0,xpos)
        wave_template = arange(min(wave)-50,max(wave)+50,0.075)

        template_i = interpolate.splev(wave_template,template_interp,ext=1)
        template_i = normalise(template_i) ### comment out this step if using normalised template spectrum

        plt.plot(wave_template,template_i,"r-")
        plt.plot(wave,spec,"k-")
        plt.xlim(min(wave),max(wave))
        plt.xlabel("Wavelength (A)",fontsize=15)

        [i.set_linewidth(3.) for i in ax.spines.itervalues()]

        #plt.savefig("echelle_equation_fit.pdf")
        plt.show()
        plt.clf()

    return x0

def iterate_whole_spectrum(spectrum):

    """
    iterate through the entire spectrum to get rough linear solutions to each order
    """

    wave_between_orders = 200 #distance in wavelength between orders

    initial_solutions = zeros((10,3))
    
    ### do the first order
    print("finding initial solution for order 0")
    x0 = cross_correlation_to_find_sol(spectrum[0],testwave_centre=6700,testwave_width=300,testdelta_centre=-0.095,testdelta_width=0.005,toplot=False)
    initial_solutions[0,0] = 0
    initial_solutions[0,1] = x0[0]
    initial_solutions[0,2] = x0[1]

    norders = len(spectrum)
    if norders > len(initial_solutions):
        norders = len(initial_solutions)
    
    
    for order in arange(1,norders):
        print("order",order )
        x0 = cross_correlation_to_find_sol(spectrum[order],testwave_centre=initial_solutions[order-1,2]-wave_between_orders,testwave_width=200,testdelta_centre=initial_solutions[order-1,1],toplot=False)
        initial_solutions[order,0] = order
        initial_solutions[order,1] = x0[0]
        initial_solutions[order,2] = x0[1]


    return initial_solutions


def fit_echelle_solution(initial_solutions,norders):

    fit_w0,sigma = echelle_equation_fit_sigclip(initial_solutions[:,0],initial_solutions[:,2])
    print(fit_w0)

    fit_deltaA = median(initial_solutions[:,1]/echelle_equation(fit_w0,initial_solutions[:,0]))
    print(fit_deltaA)
    
    savetxt("order_initial_solutions",transpose(array([arange(0,norders),echelle_equation(fit_w0,arange(0,norders)),fit_deltaA*echelle_equation(fit_w0,arange(0,norders))])),fmt="%.3f")
    
    plt.subplot(221)
    plt.scatter(initial_solutions[:,0],initial_solutions[:,2])
    plt.plot(arange(0,norders),echelle_equation(fit_w0,arange(0,norders)))

    plt.subplot(222)
    plt.scatter(initial_solutions[:,0],initial_solutions[:,2]-echelle_equation(fit_w0,initial_solutions[:,0]))

    plt.subplot(223)
    plt.scatter(initial_solutions[:,0],initial_solutions[:,1])
    plt.plot(arange(0,norders),fit_deltaA*echelle_equation(fit_w0,arange(0,norders)))

    plt.subplot(224)
    plt.scatter(initial_solutions[:,0],initial_solutions[:,1]-fit_deltaA*initial_solutions[:,2])

    #plt.savefig("echelle_equation_fit.pdf")
    #plt.show()

    return transpose(array([arange(0,norders),echelle_equation(fit_w0,arange(0,norders)),fit_deltaA*echelle_equation(fit_w0,arange(0,norders))]))


def fit_echelle_solution_recc(initial_solutions,norders,spectrum,wshift=0):

    fit_w0,sigma = echelle_equation_fit_sigclip(initial_solutions[:,0],initial_solutions[:,2])
    print(fit_w0)

    fit_deltaA = median(initial_solutions[:,1]/echelle_equation(fit_w0,initial_solutions[:,0]))
    print(fit_deltaA)
    
    initial_solutions = transpose(array([arange(0,norders),echelle_equation(fit_w0,arange(0,norders)),fit_deltaA*echelle_equation(fit_w0,arange(0,norders))]))


    for order in range(len(initial_solutions)):
        x0 = initial_solutions[order]
        
        spec = spectrum[order]
        xpos = arange(len(spec))-len(spec)/2
        mask = xpos > min(xpos)+10
        mask *= xpos < max(xpos)-10
        
        spec -= min(spec)
        spec = normalise(spec[mask],deg=10)
        xpos = xpos[mask]

        wave = polyval([x0[2],x0[1]],xpos)
        spec_i = -1*(spec-1)

        wave_template = arange(min(wave)-150,max(wave)+150,abs(x0[2]))

        template_i = interpolate.splev(wave_template,template_interp,ext=1) ### interpolate the template spectrum
        template_i = normalise(template_i) ### comment out this step if using normalised template spectrum
        template_i = -1*(template_i-1)

        #spec_i /= max(spec_i)
        #template_i /= max(template_i)

        mask = wave_template == wave_template

        try:

            drv,cc = pyasl.crosscorrRV(wave, spec_i, wave_template[mask], template_i[mask], -5000, +5000, 1, mode='doppler', skipedge=0, edgeTapering=None)

            epos = pyasl.quadExtreme(drv, cc, mode='max', dp=(10, 10), exInd=None, fullOutput=False, fullPoint=False)[0]
            epos -= wshift
            wcenter = x0[1]-epos*x0[1]/c
        except:
            print("Could not cross correlate order",order)
        print(wcenter,x0[1])
        initial_solutions[order,1] = wcenter


    fit_w0,sigma = echelle_equation_fit_sigclip(initial_solutions[:,0],initial_solutions[:,1])
    print(fit_w0)

    fit_deltaA = median(initial_solutions[:,2]/echelle_equation(fit_w0,initial_solutions[:,0]))
    print(fit_deltaA)
    
    initial_solutions = transpose(array([arange(0,norders),echelle_equation(fit_w0,arange(0,norders)),fit_deltaA*echelle_equation(fit_w0,arange(0,norders))]))

    
    # plt.subplot(221)
    # plt.scatter(initial_solutions[:,0],initial_solutions[:,2])
    # plt.plot(arange(0,norders),echelle_equation(fit_w0,arange(0,norders)))

    # plt.subplot(222)
    # plt.scatter(initial_solutions[:,0],initial_solutions[:,2]-echelle_equation(fit_w0,initial_solutions[:,0]))

    # plt.subplot(223)
    # plt.scatter(initial_solutions[:,0],initial_solutions[:,1])
    # plt.plot(arange(0,norders),fit_deltaA*echelle_equation(fit_w0,arange(0,norders)))

    # plt.subplot(224)
    # plt.scatter(initial_solutions[:,0],initial_solutions[:,1]-fit_deltaA*initial_solutions[:,2])

    # #plt.savefig("echelle_equation_fit.pdf")
    # #plt.show()

    return initial_solutions



if __name__ == "__main__":

    #initial_solutions = iterate_whole_spectrum(pyfits.getdata("sky_spec.fits"))
    #savetxt("init_solutions",array(initial_solutions))


    initial_solutions = loadtxt("init_solutions")
    fit_echelle_solution(initial_solutions,17)
