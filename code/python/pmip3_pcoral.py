# this script compute the pseudocoral based on lon, lat, sos, and tos
# total means no seasonal cycle is removed from the computation

# Import modules
from __future__ import print_function
import numpy as np
import numpy.ma as ma
from os import listdir, chdir
from os.path import isfile, join
import string
import sys

# Open data file
fns = 'KCM_midHolocene_SSS_500_999.nc'
fns = 'KCM_midHolocene_SST_500_999.nc'
fns = 'KCM_midHolocene_pcoral_500_999.nc'
sname = 'sosaline'
tname = 'sosstsst'

coral_sensor_apply(fnt,fns,fnc,tname='tos',sname='sos'):


def coral_sensor_field(latArray,lonArray,sst,sss):
'''  This function implements the bivariate model of [1] to SST and SSS fields

INPUTS:  
    - latArray,lonArray, numpy 1D arrays
    - SST (in K or degC), masked array
    - SSS (in psu*), masked array
OUTPUTS
    - coral, the pseudocoral at the same locations as SST, SSS
    = tosContri, the thermal contribution
    - sosContri, the hydrological contribution
    
    
NB: assumes that latArray,lonArray are regular grids (no curvilinear shenanigans, please).
  * assumes SSS in psu, so need to convert if this is not the case 

[1] Thompson, D. M. , T. R. Ault , M. N. Evans , J. E. Cole , and J. Emile-Geay (2011), Comparison of observed and simulated tropical climate trends using a forward model of coral δ18O, Geophys. Res. Lett., 38, L14706, doi:10.1029/2011GL048224.  
'''
    import numpy as np
    import numpy.ma as ma
    # center the fields
    nt,ny,nx = sss.shape 
    sss_m = ma.mean(sss,axis=0)
    sss_c = sss - np.tile(sss_m,(nt,1,1)) 
    sst_m = ma.mean(sst,axis=0)
    sst_c = sst - np.tile(sst_m,(nt,1,1)) 
    
    #assign different b values based on location
    a = -0.22
    b1 = 0.3007062
    b2 = 0.1552032
    b3 = 0.2619054
    b4 = 0.436509
    
    # possibly the least effficient way of doing this, but it works (we think)
    b = np.empty((len(latArray),len(lonArray)))
    for lat in range(len(latArray)):
        for lon in range(len(lonArray)):
            #Red sea
            if lonArray[lon]>=32.83 and lonArray[lon]<=43.5 and latArray[lat]>=12.38 and latArray[lat]<=28.5:
                b[lat][lon]=b1
            #Indian ocean
            elif lonArray[lon]<=120:
                b[lat][lon]=b2
            #Tropical Pacific
            elif latArray[lat]>= -5 and latArray[lat]<=13:
                b[lat][lon]=b3
            #South Pacific
            elif latArray[lat]< -5:
                b[lat][lon]=b4
            #Default: Tropical Pacific
            else:
                b[lat][lon]=b3

    # store coordinates of four b values seperately
    b1_index = np.where(b == b1)
    b2_index = np.where(b == b2)
    b3_index = np.where(b == b3)
    b4_index = np.where(b == b4)

    # create a new array with the same shape as IPsos and compute coral
    coral = np.empty_like(sss)
    tosContri = np.empty_like(sst)
    sosContri = np.empty_like(sss)
    # hydrological contribution
    for b_index, b in ((b1_index,b1), (b2_index,b2), (b3_index,b3), (b4_index,b4)):
        sosContri[:,b_index[0],b_index[1]] = b * sss_c[:,b_index[0],b_index[1]]
    # thermal comtribution    
    tosContri = a * sst_c
    # total contribution
    coral = sosContri + tosContri
    # export all three
    return (coral, tosContri, sosContri)


def coral_sensor_apply(fnt,fns,fnc,tname='tos',sname='sos'):
'''
This function converts CMIP5/PMIP3 output to pseudocoral, according to the bivariate model of [1] 
and writes the result to a separate file

INPUTS:  
    - fnt: filename for SST field, with variable name tname [default = 'tos']
    - fns: filename object field, with variable name sname [default = 'sos']
    - fnc: filename for pseudocoral file
    
NB: assumes that SST and SSS are on identical, regular grids (no curvilinear shenanigans, please)

[1] Thompson, D. M. , T. R. Ault , M. N. Evans , J. E. Cole , and J. Emile-Geay (2011), Comparison of observed and simulated tropical climate trends using a forward model of coral δ18O, Geophys. Res. Lett., 38, L14706, doi:10.1029/2011GL048224.
'''
    import cdms2, cdutil
    import numpy as np
    import numpy.ma as ma

    print ("Working on sss " + fns[:-3] + " and sst" + fnt[:-3])
    # open files    
    fs = cdms2.open(fns)
    ft = cdms2.open(fnt)

    # get the start and end time steps
    start_time_sos = fs.getAxis('time').asRelativeTime()[0]
    end_time_sos   = fs.getAxis('time').asRelativeTime()[-1]
    start_time_tos = ft.getAxis('time').asRelativeTime()[0]
    end_time_tos   = ft.getAxis('time').asRelativeTime()[-1]

    # extract Indo-Pacific region data
    IPsosVar = fs(sname, latitude=(-30.,30.),longitude=(30.,300.))
    IPtosVar = ft(tname, latitude=(-30.,30.),longitude=(30.,300.))
    # define missing values
    ma.set_fill_value(IPsosVar, 1e20)
    ma.set_fill_value(IPtosVar, 1e20)
    
    # load into arrays
    IPsos = IPsosVar.getValue()
    IPtos = IPtosVar.getValue()
    
    # get the values for computations
    sos_ma = ma.masked_equal(IPsos, 1e20)
    sos_ma = ma.array(sos_ma, mask=np.isnan(sos_ma))
    tos_ma = ma.masked_equal(IPtos, 1e20)
    tos_ma = ma.array(tos_ma, mask=np.isnan(tos_ma))
    #print ('sos_ma')
    #print (sos_ma)
    #print ('tos_ma')
    #print (tos_ma)

    # get the means map
    sos_mean = ma.mean(sos_ma,axis=0)
    tos_mean = ma.mean(tos_ma,axis=0)

    # get total mean
    sos_mean_total = ma.mean(sos_ma)
    if sos_mean_total <= 1:   #  NEED MUCH MORE SOPHISTICATED EXCEPTION HANDLING HERE
        print ('times sos by 1000')
        sos_ma = sos_ma * 1000
       
    # apply coral sensor model 
    
    coral,tosContri,sosContri = coral_sensor_field(latArray,lonArray,tos_ma,sos_ma)

    # print some diagnostics
    print ('coral before mask')
    print (coral)
    coral = ma.masked_equal(coral, 1e20)
    print ('coral after mask')
    print (coral)
    tosContri = ma.masked_equal(tosContri, 1e20)
    sosContri = ma.masked_equal(sosContri, 1e20)

    print ('the ratio of thermal to hydrological variance is ')
    print (ma.mean(ma.var(tosContri,axis=0))/ma.mean(ma.var(sosContri,axis=0)))


       # prepare a netCDF file to store the values
    ##############test###########
    if "Indo-Pacific_regridded" in fc:
        fc = cdms2.createDataset(fnc[9:fnc.find('_regridded')]+'_coral_total.nc')
    else:
        fc = cdms2.createDataset(fnc[9:-3]+'_Indo-Pacific_coral_total.nc')

    
    ########test ends############
    tobj = IPsosVar.getTime()
    timeArray = tobj.getValue()
    tobj2 = fc.copyAxis(tobj)
    tobj2.attributes = tobj.attributes
    latobj = IPsosVar.getLatitude()
    lonobj = IPsosVar.getLongitude()

    # detect whether the variable is in curvilinear grid by detecting the shape of its latitude array. If so, get only the demension needed for computation
    # yes it is in curvilinear
    if len(IPsosVar.getLatitude().getValue().shape) == 2:
        latArray = latobj.getValue()[:,0]
        lonArray = lonobj.getValue()[0,:]
        lonArray.sort()
    # nope it is in rectangular
    else:
        latArray = latobj.getValue()
        lonArray = lonobj.getValue()
        lonArray.sort()

    # latobj2 = fc.createAxis('latitude',latArray)
    latobj2 = fc.copyAxis(latobj)
    latobj2.attributes = latobj.attributes
    # lonobj2 = fc.createAxis('longitude',lonArray)
    lonobj2 = fc.copyAxis(lonobj)
    lonobj2.attributes = lonobj.attributes

    # save coral to new netCDF's variable
    coralVar = fc.createVariable("pseudocoral",np.float,(tobj2,latobj2,lonobj2))
    coralVar[:] = coral
    coralVar.units = "permil"
    coralVar.long_name = "pseudocoral d18O"

    # save temperature contribution to new netCDF's variable
    tosContriVar = fc.createVariable("tosContri",np.float,(tobj2,latobj2,lonobj2))
    tosContriVar[:] = tosContri
    tosContriVar.units = "K"
    tosContriVar.long_name = "pseudocoral d18O temperature contribution"

    # save coral to new netCDF's variable
    sosContriVar = fc.createVariable("sosContri",np.float,(tobj2,latobj2,lonobj2))
    sosContriVar[:] = sosContri
    sosContriVar.units = "permil"
    sosContriVar.long_name = "pseudocral d18O salinity contribution"

    ft.close()
    fs.close()
    fc.close()

print ("Done!")
