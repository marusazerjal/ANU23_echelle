#!/usr/bin/env python
# This line is needed for matplotlib to work
import os,sys,glob,string,pickle
from astropy.io import fits as pyfits
from numpy import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import config_file
from scipy import signal,interpolate,optimize




def polyfit_sigclip(x,y,order=2,clip=2,niter=1): ### clip in sigma
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
    return fit#,stdev,mask


def return_masks(masterflat,toplot=False, config=None):
    """
    MZ (guessing): Find orders. Return masks that include orders.
    """
    #config = config_file.set_config() # Marusa

    masterflat.astype(float)

    ybin = config["binning"]
    print("assuming a ybin of",ybin)

    print("creating a set of steps for the centres of the orders")

    steps = mean(masterflat[:,int(len(masterflat[0])/2-10):int(len(masterflat[0])/2+10)],axis=1)
    #steps = steps[:-100/ybin]

    ### normalise the flats
    norm = []
    j = int(200/ybin)
    i = 0
    while i < len(steps):
        norm.append([i+j/2, min(steps[i:i+j]), max(steps[i:i+j]) - min(steps[i:i+j])])
        i += j

    norm = array(norm)
    x = arange(len(steps))
    norm_min = interpolate.splrep(norm[:,0],norm[:,1],k=1)
    norm_min = interpolate.splev(x,norm_min)
    norm_max = interpolate.splrep(norm[:,0],norm[:,2],k=1)
    norm_max = interpolate.splev(x,norm_max)

    # plt.plot(x,norm_min)
    # plt.plot(x,norm_max+norm_min)
    # plt.plot(steps)
    # plt.show()

    steps -= norm_min
    steps /= norm_max

    from scipy.signal import find_peaks_cwt
    indexes = find_peaks_cwt(steps,arange(20,70,10)/ybin,min_snr=0.2)
    indexes = sort(indexes) ### just to make sure
    indexes = indexes[:19]
    x = arange(len(indexes))

    # plt.plot(steps)
    # plt.scatter(indexes,zeros(len(indexes)))
    # plt.show()

    def minfunc(x0):
        f = x0[0]*log(x-x0[1])+x0[2]
        return (f-indexes)**2
    x0 = optimize.least_squares(minfunc,[-1.,-1.,1.,1.]).x
    x = arange(config["norder"]+1)
    f = x0[0]*log((x-x0[1])*x0[3])+x0[2]

    #f = polyfit(x,indexes,3)
    #x = arange(config["norder"]+1)
    #f = polyval(f,x)

    # plt.scatter(arange(len(indexes)),indexes)
    # plt.plot(x,f)
    # plt.show()

    indexes = f.astype(int)

    if len(indexes)-1 < config["norder"]:
        print("!!!!!!!!!!! Ah shit, did not identify enough orders, check your flats !!!!!!!!!!!!!!!!!")
    
    # plt.scatter(indexes,zeros(len(indexes)))
    # plt.plot(steps,color="r")
    # plt.show()
    print("finding edges for each order")

    order_masks = []
    stepsize = 50
    cutoff = 0.5
    
    for order in range(len(indexes)-1):
    #for order in [22]:
        if order < config["norder"]:
            if order == 0:
                trace_top = []
                
                ### initiate trace
                xpos = arange(10/ybin,indexes[order])
                edge = mean(masterflat[int(10/ybin):indexes[order], int(len(masterflat[0])/2) - int(10/ybin) : int(len(masterflat[0])/2+10)],axis=1)
                mask = edge-min(edge) > cutoff*(max(edge)-min(edge))
                edge = min(xpos[mask])
                #print('edge0', len(edge), len(xpos))
                edge0 = edge
                width = indexes[order]-edge
                trace_top.append([len(masterflat[0])/2,edge])

                ### look left
                x = int(len(masterflat[0])/2-10)
                while x > 100:
                    xpos = arange(10/ybin,edge+width)
                    edge = mean(masterflat[int(10/ybin):int(edge+width), int(x-stepsize):int(x+stepsize)],axis=1)
                    mask = edge-min(edge) > cutoff*(max(edge)-min(edge))
                    edge = min(xpos[mask])
                    #print('edge1', len(edge), len(xpos))
                    trace_top.append([x,edge])

                    x -= stepsize

                zpt = indexes[order]-(indexes[order]-indexes[order-1])/2
                edge = edge0
                
                ### look right
                x = len(masterflat[0])/2+10
                while x < len(masterflat[0])-100:
                    xpos = arange(10/ybin,edge+width)
                    edge = mean(masterflat[int(10/ybin):int(edge+width),int(x-stepsize):int(x+stepsize)],axis=1)
                    mask = edge-min(edge) > cutoff*(max(edge)-min(edge))
                    edge = min(xpos[mask])
                    #print('edge2', len(edge), len(xpos))
                    trace_top.append([x,edge])

                    x += stepsize

                    


            else:
                trace_top = []
                zpt = indexes[order]-(indexes[order]-indexes[order-1])/2
                ### initiate trace
                xpos = arange(zpt,indexes[order])
                edge = mean(masterflat[int(zpt):indexes[order],int(len(masterflat[0])/2)-10:int(len(masterflat[0])/2)+10],axis=1)
                mask = edge-min(edge) > cutoff*(max(edge)-min(edge))
                edge = min(xpos[mask])
                width = indexes[order]-edge
                trace_top.append([len(masterflat[0])/2,edge])
                edge0 = edge
                zpt = edge-5

                ### look left
                x = len(masterflat[0])/2-10
                while x > 100:
                    xpos = arange(zpt+(edge-edge0),edge+width)
                    edge = mean(masterflat[int(zpt+(edge-edge0)):int(edge+width),int(x-stepsize):int(x+stepsize)],axis=1)
                    #plt.plot(xpos,edge)

                    emin,eminpos = min(edge),xpos[argmin(edge)]
                    if min(xpos) < eminpos:
                        xpos = xpos[argmin(edge):]
                        edge = edge[argmin(edge):]


                    mask = (edge-min(edge)) > cutoff*(max(edge)-min(edge))
                    edge = min(xpos[mask])
                    #plt.axvline(x=edge)
                    #plt.show()
                    
                    trace_top.append([x,edge])

                    x -= stepsize

                #zpt = indexes[order]-(indexes[order]-indexes[order-1])/2
                edge = edge0
                    
                ### look right
                x = len(masterflat[0])/2+10
                while x < len(masterflat[0])-100:
                    xpos = arange(zpt+(edge-edge0),edge+width)
                    edge = mean(masterflat[int(zpt+(edge-edge0)):int(edge+width),int(x-stepsize):int(x+stepsize)],axis=1)
                    #plt.plot(xpos,edge)

                    emin,eminpos = min(edge),xpos[argmin(edge)]
                    if min(xpos) < eminpos:
                        xpos = xpos[argmin(edge):]
                        edge = edge[argmin(edge):]


                    
                    mask = edge-min(edge) > cutoff*(max(edge)-min(edge))
                    edge = min(xpos[mask])
                    #plt.axvline(x=edge)
                    #plt.show()
                    
                    trace_top.append([x,edge])

                    x += stepsize


            trace_top = array(trace_top)-2.
            trace_top_fit = polyfit_sigclip(trace_top[:,0],trace_top[:,1])
            #trace_top_fit = polyval(trace_top_fit,arange(len(masterflat)))
                               



            ################# 


            trace_bottom = []
            zpt = indexes[order]+(indexes[order+1]-indexes[order])/2
            #~ if order == 0:
                 #~ zpt+= 1 # Marusa: +1
            ### initiate trace
            xpos = arange(indexes[order],zpt)
            if order ==0:
                #~ edge = mean(masterflat[int(indexes[order]):int(zpt)+1,int(len(masterflat[0])/2)-10:int(len(masterflat[0])/2)+10],axis=1) # Marusa: added +1 in zpt+1
                edge = mean(masterflat[int(indexes[order]):int(zpt),int(len(masterflat[0])/2)-10:int(len(masterflat[0])/2)+10],axis=1)
            else:
                edge = mean(masterflat[int(indexes[order]):int(zpt),int(len(masterflat[0])/2)-10:int(len(masterflat[0])/2)+10],axis=1)
            print(len(xpos), len(edge), order)
            print(indexes[order])
            print(zpt)


            mask = edge-min(edge) > cutoff*(max(edge)-min(edge))
            edge = max(xpos[mask])
            width = edge-indexes[order]
            trace_bottom.append([len(masterflat[0])/2,edge])
            edge0 = edge
            zpt = edge+5

            ### look left
            x = len(masterflat[0])/2-10
            while x > 100:
                xpos = arange(edge-width,zpt+(edge-edge0))
                edge = mean(masterflat[int(edge-width):int(zpt+(edge-edge0)),int(x-stepsize):int(x+stepsize)],axis=1)

                emin,eminpos = min(edge),xpos[argmin(edge)]
                if max(xpos) > eminpos:
                    xpos = xpos[:argmin(edge)]
                    edge = edge[:argmin(edge)]

                #plt.plot(xpos,edge)
                
                mask = edge-min(edge) > cutoff*(max(edge)-min(edge))
                edge = max(xpos[mask])

                #plt.axvline(x=edge)
                #plt.show()
                
                trace_bottom.append([x,edge])

                x -= stepsize

            #zpt = indexes[order]+(indexes[order+1]-indexes[order])/2
            edge = edge0

            ### look right
            x = len(masterflat[0])/2+10
            while x < len(masterflat[0])-100:
                xpos = arange(edge-width,zpt+(edge-edge0))
                edge = mean(masterflat[int(edge-width):int(zpt+(edge-edge0)),int(x-stepsize):int(x+stepsize)],axis=1)
                emin,eminpos = min(edge),xpos[argmin(edge)]
                if max(xpos) > eminpos:
                    xpos = xpos[:argmin(edge)]
                    edge = edge[:argmin(edge)]

                #plt.plot(xpos,edge)
                
                mask = edge-min(edge) > cutoff*(max(edge)-min(edge))
                print('MZ: mask edge, ', cutoff, edge, min(edge), max(edge))
                edge = max(xpos[mask])

                #plt.axvline(x=edge)
                #plt.show()

                
                trace_bottom.append([x,edge])

                x += stepsize


            trace_bottom = array(trace_bottom)+2.
            trace_bottom_fit = polyfit_sigclip(trace_bottom[:,0],trace_bottom[:,1])

            if toplot:
                                
                plt.scatter(trace_bottom[:,0],trace_bottom[:,1])
                plt.plot(arange(len(masterflat[0])),polyval(trace_bottom_fit,arange(len(masterflat[0]))))
 
                plt.scatter(trace_top[:,0],trace_top[:,1])
                plt.plot(arange(len(masterflat[0])),polyval(trace_top_fit,arange(len(masterflat[0]))))
                plt.show()
            
            x = arange(len(masterflat[0]))
            y = arange(len(masterflat))
            xx,yy = meshgrid(x,y)
            mask = masterflat == masterflat
            mask *= yy >  polyval(trace_top_fit,xx) 
            mask *= yy < polyval(trace_bottom_fit,xx)

            order_masks.append([mask,trace_top_fit,trace_bottom_fit])

    pickle.dump(order_masks, open(os.path.join(config["folder"], "temp/", "order_masks.pkl"), "wb"))
        
            
    return order_masks
            

if __name__ == "__main__":

    ccdsec_min = 53
    ccdsec_max = 2095

    masterflat = pyfits.getdata("/media/Onion/Data/ANU23echelle/20181115/bin2/temp/masterflat.fits")
    masterflat = masterflat[:,ccdsec_min:ccdsec_max]
    masks = return_masks(masterflat,toplot=True)
