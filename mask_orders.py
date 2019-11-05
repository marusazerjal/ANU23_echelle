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
    MZ (guessing): Find orders from the masterflat. Return masks that include orders.
    """
    
    print('MASTERFLAT', masterflat)
    #config = config_file.set_config() # Marusa
    #~ print('config', config)
    
    # MZ: A HACK assuming that all flats are the same and stable
    #~ filename = os.path.join(config['folder'], 'temp', 'order_masks.pkl')
    #~ if os.path.isfile(filename):
        #~ f=open(filename, 'rb')
        #~ order_masks = pickle.load(f)
        #~ f.close()
        #~ return order_masks
    #~ else:
        #~ pass

    masterflat.astype(float)

    ybin = config["binning"]
    print("assuming a ybin of",ybin)

    print("creating a set of steps for the centres of the orders")

    # Steps: Mean values of signal in the middle of the x axis for each row
    steps = mean(masterflat[:,int(len(masterflat[0])/2-10):int(len(masterflat[0])/2+10)], axis=1)
    #steps = steps[:-100/ybin]

    #~ print('steps0')
    #~ steps0=1.0*steps
    #~ for x in steps:
        #~ print(x)

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

    #~ fig=plt.figure()
    #~ ax=fig.add_subplot(111)
    #~ ax.plot(x,norm_min)
    #~ ax.plot(x,norm_max+norm_min)
    #fig.plot(steps)
    #~ plt.show()

    steps -= norm_min
    steps /= norm_max

    #~ print('steps1')
    #~ for x0, x1 in zip(steps0, steps):
        #~ print(x0, x1)

    from scipy.signal import find_peaks_cwt
    indexes = find_peaks_cwt(steps,arange(20,70,10)/ybin,min_snr=0.2)
    indexes = sort(indexes) ### just to make sure
    indexes = indexes[:19]
    x = arange(len(indexes))

    #~ print('indexes', len(indexes), indexes)

    # plt.plot(steps)
    # plt.scatter(indexes,zeros(len(indexes)))
    # plt.show()

    # WHAT HAPPENS HERE?
    # MZ: I think this is to fit how indices increase along the y-axis (orders)
    # Probably because the last orders have worse signal, and they come very close to each other.
    def minfunc(x0):
        f = x0[0]*log(x-x0[1])+x0[2]
        return (f-indexes)**2
    x0 = optimize.least_squares(minfunc,[-1.,-1.,1.,1.]).x
    x = arange(config["norder"]+1)
    f = x0[0]*log((x-x0[1])*x0[3])+x0[2] # x0[3] is not fitted!!

    #f = polyfit(x,indexes,3)
    #x = arange(config["norder"]+1)
    #f = polyval(f,x)

    # plt.scatter(arange(len(indexes)),indexes)
    # plt.plot(x,f)
    # plt.show()

    indexes = f.astype(int)
    #~ print('INDEXES')
    #~ print(len(indexes), indexes)

    if len(indexes)-1 < config["norder"]:
        print("!!!!!!!!!!! Ah shit, did not identify enough orders, check your flats !!!!!!!!!!!!!!!!!")
    
    # plt.scatter(indexes,zeros(len(indexes)))
    # plt.plot(steps,color="r")
    # plt.show()
    print("finding edges for each order")
    #~ print(indexes)

    order_masks = []
    stepsize = 50
    cutoff = 0.5
    
    for order in range(len(indexes)-1):
    #for order in [22]:
        if order < config["norder"]:
            if order == 0:
                trace_top = []
                
                ### initiate trace (in the middle of the frame). MZ: It seems that the edges found above were just initial guesses. Now do this for real here.
                xpos = arange(10/ybin,indexes[order])
                #~ print('xpos', xpos)
                # edge: average y-values for each row from the bottom of the image (10/ybin) to the y-value where order has been detected in the middle
                edge = mean(masterflat[int(10/ybin):indexes[order], int(len(masterflat[0])/2) - int(10/ybin) : int(len(masterflat[0])/2+10)], axis=1)
                #~ print('edge', edge)
                mask = edge-min(edge) > cutoff*(max(edge)-min(edge))
                edge = min(xpos[mask]) # Lower y-value where the trace begins in the middle of the image
                #~ print('edge0', edge, len(xpos))
                edge0 = edge
                width = indexes[order]-edge
                #~ print('width', width)
                trace_top.append([len(masterflat[0])/2,edge])

                ### look left
                x = int(len(masterflat[0])/2-10)
                while x > 100:
                    xpos = arange(10/ybin,edge+width)
                    # Find the edge in the same manner as before - the only difference is that in each iteration you shift leftward for 50 (stepsize) pixels.
                    edge = mean(masterflat[int(10/ybin):int(edge+width), int(x-stepsize):int(x+stepsize)], axis=1)
                    mask = edge-min(edge) > cutoff*(max(edge)-min(edge))
                    edge = min(xpos[mask])
                    #print('edge1', len(edge), len(xpos))
                    trace_top.append([x,edge])

                    x -= stepsize

                # ??
                #~ print('order', order, indexes[order-1])
                #~ zpt = indexes[order]-(indexes[order]-indexes[order-1])/2 # indexes[order-1] here means indexes[-1]!!
                #~ print('zpt', zpt, trace_top)
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
                    edge_marusa=edge
                    edge = mean(masterflat[int(zpt+(edge-edge0)):int(edge+width),int(x-stepsize):int(x+stepsize)],axis=1)
                    #plt.plot(xpos,edge)
                    print('edgggge', int(zpt+(edge_marusa-edge0)), int(edge_marusa+width), int(x-stepsize), int(x+stepsize), masterflat.shape, edge, masterflat)

                    # ADDED by MARUSA:
                    if len(xpos)>len(edge):
                        #~ print('##########CUTTING')
                        xpos=xpos[:len(edge)]

                    emin,eminpos = min(edge),xpos[argmin(edge)]
                    if min(xpos) < eminpos:
                        print('OLD', len(xpos), len(edge))
                        xpos = xpos[argmin(edge):]
                        edge = edge[argmin(edge):]
                        print('NEW', len(xpos), len(edge))


                    print('LEN', len(xpos), len(edge), int(zpt+(edge_marusa-edge0)), edge_marusa+width, masterflat.shape, x)

                    

                    if len(xpos)>1: # MZ
                        mask = (edge-min(edge)) > cutoff*(max(edge)-min(edge))
                        edge = min(xpos[mask]) # THE CODE CRASHES HERE!!
                        #plt.axvline(x=edge)
                        #plt.show()
                    else:
                        edge=edge[0]     # MZ
                    trace_top.append([x,edge])
                    print('edge', edge)
                    print('\n')
                    x -= stepsize

                #zpt = indexes[order]-(indexes[order]-indexes[order-1])/2
                edge = edge0
                    
                ### look right
                x = len(masterflat[0])/2+10
                while x < len(masterflat[0])-100:
                    # xpos: rows in a y-direction
                    xpos = arange(zpt+(edge-edge0),edge+width) # George
                    #~ xpos = arange(int(zpt+(edge-edge0)),int(edge+width)) # MZ
                    edge_marusa = edge
                    edge = mean(masterflat[int(zpt+(edge-edge0)):int(edge+width),int(x-stepsize):int(x+stepsize)],axis=1)
                    #~ print('XXX', order, x, len(masterflat[0])-100, len(xpos), len(edge), masterflat.shape, 'zpt', int(zpt+(edge_marusa-edge0)), int(edge_marusa+width))
                    #~ if int(edge_marusa+width)>masterflat.shape[0]:
                        #~ print('*************', int(edge_marusa+width), masterflat.shape[0])
                    #plt.plot(xpos,edge)

                    # ADDED by MARUSA:
                    if len(xpos)>len(edge):
                        #~ print('##########CUTTING')
                        xpos=xpos[:len(edge)]

                    emin,eminpos = min(edge),xpos[argmin(edge)]
                    #~ print('emin, eminpos', emin, eminpos, min(xpos))
                    if min(xpos) < eminpos: # George
                    #~ if min(xpos) <= eminpos: # Marusa
                        #~ print('OLD', len(xpos), len(edge))
                        xpos = xpos[argmin(edge):]
                        edge = edge[argmin(edge):]
                        #~ print('NEW', len(xpos), len(edge))


                    
                    mask = edge-min(edge) > cutoff*(max(edge)-min(edge))
                    #~ print('ANY(MASK)', any(mask), len(mask), len(xpos))
                    edge = min(xpos[mask]) # THIS IS WHERE THE CODE CRASHES!!
                    #plt.axvline(x=edge)
                    #plt.show()
                    
                    trace_top.append([x,edge])
                    #~ print([x, edge], edge0)
                    #~ print('\n')

                    x += stepsize


            trace_top = array(trace_top)-2.
            trace_top_fit = polyfit_sigclip(trace_top[:,0],trace_top[:,1]) # (x,y) pixels
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
            #~ print(len(xpos), len(edge), order)
            #~ print(indexes[order])
            #~ print(zpt)


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
                #~ print('MZ: mask edge, ', cutoff, edge, min(edge), max(edge))
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
            #~ print('mask', mask.shape)
            order_masks.append([mask,trace_top_fit,trace_bottom_fit])
    #~ print(order_masks)
    pickle.dump(order_masks, open(os.path.join(config["folder"], "temp/", "order_masks.pkl"), "wb"))
        
            
    return order_masks
            

if __name__ == "__main__":

    #~ ccdsec_min = 53
    #~ ccdsec_max = 2095

    #~ masterflat = pyfits.getdata("/media/Onion/Data/ANU23echelle/20181115/bin2/temp/masterflat.fits")
    #~ masterflat = masterflat[:,ccdsec_min:ccdsec_max]
    #~ masks = return_masks(masterflat,toplot=True)

    import imp
    config_filename = sys.argv[1]
    print('CONFIG FILENAME', config_filename)
    config_file = imp.load_source(config_filename.replace('.py', ''), config_filename)
    config = config_file.set_config()

    masterflat = pyfits.getdata("../masterflat0719.fits")
    masks = return_masks(masterflat, toplot=True, config=config)
