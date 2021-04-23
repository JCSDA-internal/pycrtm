#!/usr/bin/env python3
import os, h5py, sys 
import numpy as np
from matplotlib import pyplot as plt
from pyCRTM import pyCRTM, profilesCreate
 
def main(sensor_id):
    thisDir = os.path.dirname(os.path.abspath(__file__))
    cases = os.listdir( os.path.join(thisDir,'data') )
    cases.sort()
    # create 4 profiles for each of the 4 cases
    profiles = profilesCreate( 4, 92, nAerosols=0, nClouds = 0 )
    storedTb = []
    storedEmis = []
    # populate the cases, and previously calculated Tb from crtm test program.    
    for i,c in enumerate(cases):
        h5 = h5py.File(os.path.join(thisDir,'data',c) , 'r')
        profiles.Angles[i,0] = h5['zenithAngle'][()]
        profiles.Angles[i,1] = 999.9 
        profiles.Angles[i,2] = 100.0  # 100 degrees zenith below horizon.
        profiles.Angles[i,3] = 0.0 # zero solar azimuth 
        profiles.Angles[i,4] = h5['scanAngle'][()]
        profiles.DateTimes[i,0] = 2001
        profiles.DateTimes[i,1] = 1
        profiles.DateTimes[i,2] = 1
        profiles.Pi[i,:] = np.asarray(h5['pressureLevels'] )
        profiles.P[i,:] = np.asarray(h5['pressureLayers'][()])
        profiles.T[i,:] = np.asarray(h5['temperatureLayers'])
        profiles.Q[i,:] = np.asarray(h5['humidityLayers'])
        profiles.O3[i,:] = np.asarray(h5['ozoneConcLayers'])
        profiles.climatology[i] = h5['climatology'][()]
        profiles.surfaceFractions[i,:] = h5['surfaceFractions']
        profiles.surfaceTemperatures[i,:] = h5['surfaceTemperatures']
        profiles.Salinity[i] = 33.0 
        profiles.windSpeed10m[i] = 5.0
        profiles.LAI[i] = h5['LAI'][()]
        profiles.windDirection10m[i] = h5['windDirection10m'][()]
        # land, soil, veg, water, snow, ice
        profiles.surfaceTypes[i,0] = h5['landType'][()]
        profiles.surfaceTypes[i,1] = h5['soilType'][()]
        profiles.surfaceTypes[i,2] = h5['vegType'][()]
        profiles.surfaceTypes[i,3] = h5['waterType'][()]
        profiles.surfaceTypes[i,4] = h5['snowType'][()]
        profiles.surfaceTypes[i,5] = h5['iceType'][()]
        storedTb.append(np.asarray(h5['Tb_atms_no_clouds'][0:22]))
        storedEmis.append(np.asarray(h5['emissivity_atms_no_clouds'][0:22]))
        h5.close()

    crtmOb = pyCRTM()
    crtmOb.profiles = profiles
    crtmOb.sensor_id = sensor_id
    crtmOb.nThreads = 4

    crtmOb.loadInst()

    crtmOb.runDirect()
    forwardTb = crtmOb.Bt
    forwardEmissivity = crtmOb.surfEmisRefl[0,:]
    crtmOb.surfEmisRefl = []

    crtmOb.runK()
    kTb = crtmOb.Bt
    kEmissivity = crtmOb.surfEmisRefl[0,:]

    if ( all( np.abs( forwardTb.flatten() - np.asarray(storedTb).flatten() ) <= 1e-5)  and all( np.abs( kTb.flatten() - np.asarray(storedTb).flatten() ) <= 1e-5) ):
        print("Yay! all values are close enough to what CRTM test program produced!")
    else: 
        print("Boo! something failed. Look at atms no cloud plots")
        wavenumbers = np.zeros([4,22])
        wavenumbers[0:4,:] = np.linspace(1,23,22)
        
        plt.figure()
        plt.plot(wavenumbers.T,forwardTb.T-np.asarray(storedTb).T ) 
        plt.legend(['1','2','3','4'])
        plt.savefig(os.path.join(thisDir,'atms_no_clouds'+'_spectrum_forward.png'))
        plt.figure()
        plt.plot(wavenumbers.T,forwardEmissivity.T-np.asarray(storedEmis).T)
        plt.savefig(os.path.join(thisDir,'atms_no_clouds'+'_emissivity_forward.png')) 
    
        plt.figure()
        plt.plot(wavenumbers.T,kTb.T-np.asarray(storedTb).T)
        plt.savefig(os.path.join(thisDir,'atms'+'_spectrum_k.png'))
        plt.figure()
        plt.plot(wavenumbers.T,kEmissivity.T-np.asarray(storedEmis).T)
        plt.savefig(os.path.join(thisDir,'atms'+'_emissivity_k.png')) 
        sys.exit("Boo! didn't pass tolerance with CRTM test program.")


if __name__ == "__main__":
    sensor_id = 'atms_npp'
    main(sensor_id)
 
