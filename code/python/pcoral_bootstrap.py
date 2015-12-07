# This script computes the seasonal range of a specific file
import cdms2, cdutil, cdtime
import bootstrap, bandpass
import numpy as np
import numpy.ma as ma

#  Routine to compute and isolate the seasonal cycle
def seasonal_cycle(Xb):
    nb,nt = Xb.shape
    ny    = int(nt/12)
    clim  = np.empty((nb,12))
    anom  = Xb*0
    
    for i in range(12):
        clim[:,i]=Xb[:,i::12].mean(axis=1)
    
    anom = Xb - np.tile(clim,(1,ny))
    return clim, anom


def computer(name, start_lon, end_lon, start_lat, end_lat, Nb=200, Lb=24, windows = [20,30,40,50,75,100]):
    '''
    Where the hell is my docstring?
    '''
    # filtering parameters    
    fs = 1; f_hi = 1/(12*2.0); f_lo = fs/(12*7.0)
    
    # open file
    f = cdms2.open(name, 'r')
    start_time = f.getAxis('time').asRelativeTime()[0]
    end_time = f.getAxis('time').asRelativeTime()[-1]

    # extract variable of interest in east pacific area
    coral = f('pseudocoral',latitude=(start_lat,end_lat),longitude=(start_lon,end_lon))
    # print 'coral'
    # print coral
    f.close()

    # compute spatial mean
    cdutil.setTimeBoundsMonthly(coral,stored=0)
    spatial_mean = cdutil.averager(coral,axis='xy')
    
    # generate boostrap samples
    Xb = bootstrap.block_bootstrap_ET(spatial_mean, Lb, Nb)
    #print 'spatial_mean_bootstrap'
    #print spatial_mean_bootstrap
    nw = len(windows) # number of windows
   
    
    seasonal_amp = np.empty((nw,Nb))
    variance     = np.empty((nw,Nb))
    
    index = 0  # loop over windows
    for i in windows:
        Xw =  Xb[:,:i*12]  # sample over window
        clim, anom = seasonal_cycle(Xw)  # isolate seasonal cycle
        # compute seasonal amplitude
        smax = np.nanmax(clim, axis=1)
        smin = np.nanmin(clim, axis=1)
        seasonal_amp[index,:] = smax - smin
        
        # compute ENSO variance
        anom2_7 = np.empty(anom.shape)
        for b in range(Nb):
            # apply bandpass filter        
            anom2_7[b,:] = bandpass.butter_bandpass_filter(anom[b,:],f_lo,f_hi,fs)
        # compute variance per se     
        variance[index,:]  = np.var(anom2_7,axis=1)
        index +=1  # update index


    return (variance, seasonal_amp)


# DEBUG: plot the distributions and scatterplots
#s = seasonal_amp[4,:]
#v = variance[4,:]
#plt.scatter(s,v,alpha=0.3)
