import os,sys,string
from numpy import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from scipy import interpolate,optimize,signal,stats


atlas_coarse = loadtxt("atlas_coarse")
#atlas_fine = loadtxt("atlas_fine")

def gaussian(x0,xpos):
    return x0[0]*exp(-(xpos-x0[1])**2/(2*x0[2]**2))+x0[3]

def fitgaussian(xpos,order,method="highest_peak",toplot=False):
    ### method: highest_peak,only_peak,nearest_peak,middle

    

    def minfunc(x0):
        if x0[0] < 0 or x0[2] > 5 or x0[1] < min(xpos)+5 or x0[1] > max(xpos)-5:
            return nan
        else:
            f = gaussian(x0,xpos)
            return sum((f-order)**2)


    if method == "highest_peak":
        x0 = [max(order),xpos[argmax(order)],1.,min(order)]

    if method == "only_peak":
        peaks = signal.find_peaks(order,prominence = 1000)[0]

        if len(peaks) > 1:
            peak_pos = nan
        else:
            
            peak_pos = xpos[argmax(order)]
            
        x0 = [max(order),peak_pos,1.,min(order)]

    if method == "nearest_peak":

        try:

            npeaks = 0
            prominence = 500
            while npeaks == 0 and prominence > 10:
                
                peaks = signal.find_peaks(order,prominence = prominence)[0]
                npeaks = len(peaks)
                prominence /= 2
            #print xpos[peaks]
            peak_pos = xpos[peaks][argmin(abs(xpos[peaks]-median(xpos)))]

            #plt.plot(xpos,order)
            #plt.scatter(xpos[peaks],zeros(len(peaks)))
            #plt.show()
            
            mask = abs(xpos-peak_pos) < 4
            xpos = xpos[mask]
            order = order[mask]            
        except ValueError:
            print("didn't find peaks")
            peak_pos = nan #median(xpos)
        x0 = [max(order),peak_pos,1.,min(order)]

    if method == "middle":
        try:

            npeaks = 0
            prominence = 500
            while npeaks == 0 and prominence > 10:
                
                peaks = signal.find_peaks(order,prominence = prominence)[0]
                npeaks = len(peaks)
                prominence /= 2

            if npeaks >= 1:

                peak_pos = median(xpos)
            else:
                peak_pos = nan

        except ValueError:
            print("didn't find peaks")
            peak_pos = nan #median(xpos)
        x0 = [max(order),peak_pos,1.,min(order)]
        
        



    if x0[1] == x0[1]:
        x0 = optimize.fmin(minfunc,x0,disp=0)

    

    if toplot:
        ### PLOT it
        f = gaussian(x0,xpos)
        plt.plot(xpos,order)
        plt.plot(xpos,f)
        plt.show()

    if x0[0] < 100:
        return nan
    else:
        return x0[1]


    
def polyfit_sigclip(x,y,order=2,clip=2,niter=3): ### clip in sigma
    ### do a first fit with a median filter
    fit = polyfit(signal.medfilt(x,kernel_size=3),y,2)
    #print "length of original array",len(x)
    mask = abs(polyval(fit,x)-y) < clip * std(polyval(fit,x)-y)

    for i in range(niter):
        fit = polyfit(x[mask],y[mask],order)
        mask *= abs(polyval(fit,x)-y) < clip * std(polyval(fit,x[mask])-y[mask])

        #print "length of poly fit clipped array",len(x[mask])

    fit = polyfit(x[mask],y[mask],order)
    stdev = std(polyval(fit,x[mask])-y[mask])
    return fit,stdev,mask
    
        
    

def doorder(order,x0_init,toplot=True):
    # MZ: set stdev so it exists in any case
    stdev=99
    
    xpos = arange(len(order))-len(order)/2
    wave_init = polyval(x0_init,xpos)


    atlas_use = loadtxt("atlas_coarse")
    atlasmask = atlas_use < max(wave_init)
    atlasmask *= atlas_use > min(wave_init)

    # plt.plot(wave_init,order)
    # plt.scatter(atlas_coarse[atlasmask],ones(len(atlas_coarse[atlasmask])),color="r")
    # plt.show()

    extent = 1000
    step = 50
    #midpoint = len(order)/2
    midpoint = len(order)
    mode = "only_peak"
    radius_search = 1
    atlas_use = loadtxt("atlas_coarse")

    while extent < len(order)-10:#midpoint-10:

        atlasmask = atlas_use < max(wave_init[midpoint-extent:])
        atlasmask *= atlas_use > min(wave_init[midpoint-extent:])

        peaklist = zeros([len(atlas_use[atlasmask]),2])

        for i in range(len(atlas_use[atlasmask])):
            line = atlas_use[atlasmask][i]
            linemask = abs(wave_init-line) < radius_search

            if len(xpos[linemask]) > 5:
                line_centre = fitgaussian(xpos[linemask],order[linemask],method=mode,toplot=False)

                peaklist[i,0] = atlas_use[atlasmask][i]
                peaklist[i,1] = line_centre

                
        mask = peaklist[:,1] == peaklist[:,1]
        mask *= peaklist[:,1] != 0
        mask *= peaklist[:,0] != 0
        peaklist = peaklist[mask]
        
        if extent > 1200:
            deg = 2
            mode = "only_peak"
            radius_search = 1.
            atlas_use = loadtxt("atlas_coarse")
            nclip = 5
            niter = 1

        if extent > 1500:
            deg = 3
            mode = "only_peak"
            radius_search = 0.5
            atlas_use = loadtxt("atlas_coarse")
            nclip = 2
            niter = 1

        else:
            deg = 1
            mode = "only_peak"
            radius_search = 2
            atlas_use = loadtxt("atlas_coarse")
            nclip = inf
            niter = 1


        if len(peaklist) > 5:
        
            fit,stdev,mask = polyfit_sigclip(peaklist[:,1],peaklist[:,0],order=deg,clip=nclip,niter=niter)
            # if stdev > 0.13:
            #     deg = 3
            #     nclip = 1
            #     fit,stdev,mask = polyfit_sigclip(peaklist[:,1],peaklist[:,0],order=deg,clip=nclip,niter=niter)
            wave_init = polyval(fit,xpos)

            #print extent
            #plt.plot(wave_init,order)
            #plt.scatter(atlas_coarse[atlasmask],ones(len(atlas_coarse[atlasmask])),color="r")
            #plt.show()

        extent += step




    # ### Do it one last time, but use the predicted wavelength position and a finer atlas
    # atlas_fine = loadtxt("atlas_fine")

    # wave_init = polyval(fit,xpos)

    # atlasmask = atlas_fine > min(wave_init)+1
    # atlasmask *= atlas_fine < max(wave_init)-1

    # peaklist = zeros([len(atlas_fine[atlasmask]),2])
    
    # for i in range(len(atlas_fine[atlasmask])):
    #     line = atlas_fine[atlasmask][i]
    #     min_dist = sort(abs(line-atlas_fine))[1]
    #     linemask = abs(wave_init-line) < min_dist/2.
        
    #     if len(xpos[linemask]) > 5:
    #         line_centre = fitgaussian(xpos[linemask],order[linemask],method="middle",toplot=False)

    #         peaklist[i,0] = atlas_fine[atlasmask][i]
    #         peaklist[i,1] = line_centre

    # mask = peaklist[:,1] == peaklist[:,1]
    # mask *= peaklist[:,1] != 0
    # peaklist = peaklist[mask]

    # fit,stdev,mask = polyfit_sigclip(peaklist[:,1],peaklist[:,0],order=5,clip=5,niter=1)

    try:
        peaklist = peaklist[mask]
    except IndexError:
        print("no peaks found")
        stdev = 99

        


    if toplot:
        
        plt.figure(figsize=(10,15))
        plt.subplots_adjust(hspace=0.25)
        ax = plt.subplot(311)
        wave = polyval(fit,xpos)
        plt.plot(wave,order,color="k")
        plt.scatter(polyval(fit,peaklist[:,1]),zeros(len(peaklist)),color="r")
        plt.xlim(min(wave),max(wave))
        plt.xlabel("Wavelength",fontsize=15)
        plt.ylabel("ThAr",fontsize=15)
                    
        [i.set_linewidth(3.) for i in ax.spines.itervalues()]

        ax = plt.subplot(312)
        plt.scatter(peaklist[:,1],peaklist[:,0],color="k")
        f = polyval(fit,arange(min(xpos),max(xpos)))
        plt.plot(arange(min(xpos),max(xpos)),f,color="r")
        plt.xlim(min(xpos),max(xpos))
        plt.xlabel("Pixel",fontsize=15)
        plt.ylabel("Wavelength (A)",fontsize=15)
        [i.set_linewidth(3.) for i in ax.spines.itervalues()]

        ax = plt.subplot(313)
        plt.scatter(peaklist[:,1],polyval(fit,peaklist[:,1])-peaklist[:,0])
        plt.text(0.95,0.95,"$\sigma = "+str(round(stdev,4))+"\,\AA$",transform=ax.transAxes,ha="right",va="top",fontsize=18)
        if stdev > 0.1:
            plt.text(0.95,0.75,"Wave Cal has failed",transform=ax.transAxes,ha="right",va="top",fontsize=18,color="r",weight="black")
            
        plt.axhline(y=0,color="r")
        plt.xlim(min(xpos),max(xpos))
        plt.xlabel("Pixel",fontsize=15)
        plt.ylabel("Wavelength residual (A)",fontsize=15)
        [i.set_linewidth(3.) for i in ax.spines.itervalues()]

        plt.show()

    return peaklist,stdev

def run_spectrum(arc,initial_solutions):

    all_peaks = []
    norders = min([len(arc),len(initial_solutions)])
    for i in range(0,norders):
        peaklist,stdev = doorder(arc[i],[initial_solutions[i,2],initial_solutions[i,1]],toplot=False)
        print("order",i,stdev,len(peaklist))
        if stdev < 0.13:
            peaklist = list(transpose(array([i*ones(len(peaklist)),peaklist[:,1],peaklist[:,0]])))
            all_peaks += peaklist

    all_peaks = array(all_peaks)
    #savetxt("thar_peaks",all_peaks,fmt="%.5f")
    return all_peaks

if __name__ == "__main__":

    #arc = pyfits.getdata("arc_spec.fits")
    arc = pyfits.open("/media/Onion/Data/ANU23echelle/20181115/bin2/temp/ANU23e_2MASSJ02224418-6022476_2018-11-15T15-08-55.040.fits")
    arc = arc[2].data

    initial_solutions = loadtxt("/media/Onion/Data/ANU23echelle/20181115/bin2/temp/order_initial_solutions")
    run_spectrum(arc,initial_solutions)
