import os,sys,string,pickle,glob
from astropy.io import fits as pyfits
from numpy import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import config_file
from scipy import optimize,signal,interpolate
#config = config_file.set_config()

#binning = config["binning"]

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

def mask_order(fits,order_masks, binning=None):
    extraction = []
    x,y = arange(len(fits[0])),arange(len(fits))
    xx,yy = meshgrid(x,y)
    
    for order in order_masks:
        fits_order = fits.copy()
        fits_order[invert(order[0])] = nan
        
        ymax = nanmax(yy[order[0]].flatten())+10/binning
        ymin = nanmin(yy[order[0]].flatten())-10/binning

        if ymin < 0:
            ymin = 0
        if ymax > max(yy.flatten()):
            ymax = max(yy.flatten())

        # # print "order",ymax,ymin

        #plt.imshow(fits_order[ymin:ymax],aspect="auto")
        #plt.show()
        

        extraction.append([fits_order[ymin:ymax],xx[ymin:ymax],yy[ymin:ymax]])
    return extraction


def gaussian(x0,x):
    return x0[0]*exp(-(x-x0[1])**2/(2*x0[2]**2))+x0[3]

def fitgaussian(x,y):
    print(x, y)
    mask = y == y
    x,y = x[mask],y[mask]
    x = x[2:-2]
    y = y[2:-2]

    def minfunc(x0):
        f = gaussian(x0,x)
        return sum((f-y)**2)

    x0 = [max(y)-min(y),x[argmax(y)],2.5,min(y)]
    x0 = optimize.fmin(minfunc,x0,disp=0)

    # plt.plot(x,y)
    # plt.plot(x,gaussian(x0,x))
    # plt.show()
    return x0

def shiftgaussian(x,y,sigma): ### shift the same gaussian
    mask = y == y
    x,y = x[mask],y[mask]
    x = x[2:-2]
    y = y[2:-2]

    def minfunc(x0):
        if x0[1] > min(x) and x0[1] < max(x):
            f = gaussian([x0[0],x0[1],sigma,x0[2]],x)
            return sum((f-y)**2)
        else:
            return nan
    
    x0 = [max(y)-min(y),x[argmax(y)],min(y)]
    x0 = optimize.fmin(minfunc,x0,disp=0)
    
    #plt.plot(x,y)
    #plt.plot(x,gaussian([x0[0],x0[1],sigma,x0[2]],x))
    #plt.show()

    return x0


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
    
def find_trace(fits_extracted, binning=None):
    print("tracing orders")
    plt.figure(figsize=(20,20))

    trace_array = []
    for order in range(len(fits_extracted)):
    #for order in [19,20,21]:

        
        fits_extracted_order,xx,yy = fits_extracted[order]
        x,y = xx[0],yy[:,0]
        mask = fits_extracted_order == 0
        fits_extracted_order[mask] = nan

        ### fit the middle to get a width
        center = mean(fits_extracted_order[:,len(fits_extracted_order)/2-10:len(fits_extracted_order)/2+10],axis=1)
        width = fitgaussian(y,center)[2]
        if abs(width) > 5.0/binning:
            width = 5.0/binning
        if width < 1.5:
            width = 1.5
        print("width of extraction",width)

        trace = []
        for xpos in x[100:-100]:
            try:
                f = nanmedian(fits_extracted_order[:,xpos-50:xpos+50],axis=1)
                x0 = shiftgaussian(y,f,width)
                trace.append([xpos,x0[1]])
            except ValueError:
                print("Value error while gaussian fitting")
                print(y,f,width)

        trace = array(trace)
        trace[:,1] = signal.medfilt(trace[:,1],5)
        fit,dymmy,dummy = polyfit_sigclip(trace[:,0],trace[:,1],order=2,clip=2)
        fit = polyval(fit,x)
        #plt.scatter(trace[:,0],trace[:,1],s=2,color="r")
        plt.plot(x,fit,color="k",lw=2)
        plt.imshow(fits_extracted_order,extent=(min(x),max(x),min(y),max(y)),aspect="auto",origin="lowerleft")
        #plt.show()
        #plt.clf()
        trace_array.append(transpose(array([x,fit,width*ones(len(x))])))

    plt.xlim(0,2048)
    plt.ylim(0,2048/binning)
    #plt.savefig("trace.pdf")
    #plt.show()
    #plt.clf()
    #plt.close()
            
    return array(trace_array)


def extract_trace(fits_extracted,trace_array):
    spectrum = []
    background = []
    for order in range(len(fits_extracted)):

        spectrum_order = []
        background_order = []
        
        fits_extracted_order,xx,yy = fits_extracted[order]
        #plt.imshow(fits_extracted_order,aspect="auto",interpolation="nearest")
        #plt.show()
        
        x,y = xx[0],yy[:,0]
        mask = fits_extracted_order == 0
        fits_extracted_order[mask] = nan
    
        for i in range(len(trace_array[order])):
            x0 = [1,trace_array[order,i,1],trace_array[order,i,2],0]
            weights = gaussian(x0,y)
            f = fits_extracted_order[:,x[i]]
            

            mask = f == f
            f,weights = f[mask],weights[mask]
            spectrum_order.append(sum(f*weights))
            
            bk = abs(y[mask]-trace_array[order,i,1]) > 3*trace_array[order,i,2]
            #print weights
            try:

                #plt.plot(y[mask],f,color="b",alpha=0.1)
                
                bkfit,dummy,dummy = polyfit_sigclip(y[mask][bk],f[bk],order=1,clip=2,niter=3)
                bk = polyval(bkfit,y[mask])

                #plt.plot(y[mask],bk,color="r",alpha=0.1)
                #plt.show()
            
            except:
                print("Polyfit error"                )
                bk = nanmedian(f[bk])*ones(len(f))

            background_order.append(sum(bk*weights))
            

        #plt.show()

        background_order = array(background_order)
        spectrum_order = array(spectrum_order)
        
        if median(background_order) != median(background_order):
            background_order = 0
        spectrum.append(spectrum_order-background_order)
        background.append(background_order)

        #plt.clf()
        #plt.plot(spectrum_order-background_order)
        #plt.plot(background_order)
        #plt.show()

    return spectrum,background



if __name__ == "__main__":

    ccdsec_min = 53
    ccdsec_max = 2095

    folder = "/media/Onion/Data/ANU23echelle/20181129/"
    fitsname = "T2M3Ec-20181129.110305-0204.fits"
    bias = pyfits.getdata(folder+"/temp/masterbias.fits")
    order_masks = pickle.load(open(folder+"/temp/order_masks.pkl","rb"))
    fits = pyfits.getdata(folder+fitsname)-bias
    fits = fits[:,ccdsec_min:ccdsec_max]
    fits_extracted = mask_order(fits,order_masks)

    #fits_extracted = pickle.load(open(folder+fitsname+".shear.pkl","rb"))
    trace_array = find_trace(fits_extracted)
    pickle.dump(trace_array,open(folder+"/temp/"+fitsname+".trace.pkl","wb"))
    trace_array = pickle.load(open(folder+"/temp/"+fitsname+".trace.pkl","rb"))
    extract_trace(fits_extracted,trace_array)
