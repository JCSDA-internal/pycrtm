# pyCRTM - python interface to CRTM.

## Bryan M. Karpowicz, Ph.D. - USRA/GESTAR/NASA 610.1 Global Modeling and Assimilation Office, with Contributions from Patrick Stegmann, Dr.-Ing. - JCSDA

This is a basic python interface to CRTM v2.4.0. 

The user interface is designed to be very similar to the python RTTOV interface. So, the user sets profiles, passes them to an object for a desired sensor, runs either the forward model/K-matrix, and pulls the brightness temperature/Jacobian/transmission/emissivity out of the object.  


This `README` has 4 parts:

1. Installation -- installing this library
2. Test/Examples -- describing test in testCases subdirectory
3. Importing -- how to use this library in a project.
4. Using the interface -- HOWTO/run through on how to use this interface

- Bryan Karpowicz -- March 20, 2021
---------------------------------------------------------------------------------------- 

## 1. Installation:
- For the novice that doesn't care about where or how this installs, look at the crtm-bundle and run kickstart_pyCRTM.sh Otherwise...
- Dependencies CRTM, h5py, numpy and scikit-build (install those first, if you don't have them). Note crtm must be built with the static option (`ecbuild --static`) 
- Configuration
First modify `setup.cfg` to point to the crtm install location (path underneath should contain `lib/libcrtm.a`). 
```
[Setup]
# Specify the location of the crtm install (ecbuild install ONLY)
crtm_install = /discover/nobackup/bkarpowi/github/JCSDA_crtm/crtm-bundle/crtm/build
# Download Coefficients
# Controls whether coefficients are downloaded
download = True
#This will move the coefficients with the package install.
coef_with_install = True
[Coefficients]
# Use to specify alternative coefficient file location where Little Endian Coefficient files are stored.
# If user desires coefficients to be stored with the installed package, leave this alone.
# If user selects coef_with_install = False, this must be specified.
# set argument below (path) to the full path of the coefficients.
path = /discover/nobackup/projects/gmao/obsdev/bkarpowi/tstCoef/
```
In the example above the coefficients will be included with the pycrtm install. To change this, set `coef_with_install` and set `path` to the location where you would like crtm coefficients stored. If you already have a directory with coefficients, you can set `download` and `coef_with_install` to False, and set `path` to that location. The pycrtm configuration will then point to the location in `path`.  

- Installation 
There are two recommended ways to install. The first, if the user has full write access to their python distribution, it may be installed globally using:
```
python3 setup.py install 
```
This will take some time as it will download coefficients, move them around, compile the pycrtm module, and link against th crtm library.

The second, if the user doesn't have full write access to their python distribution is to first build a wheel, and install using pip:
```
python3 setup.py bdist_wheel
```
This will take some time as it will download coefficients, move them around, compile the pycrtm module, and link against the crtm library. Once the wheel has been built, it may be installed locally using pip:
```
pip install dist/pyCRTM_JCSDA*.whl --target /discover/nobackup/projects/gmao/obsdev/bkarpowi/pythonModules/
```
paired with appending `/discover/nobackup/projects/gmao/obsdev/bkarpowi/pythonModules/` to the `PYTHONPATH` environment variable in your .bashrc or .cshrc.

For Bash this is:
```
export PYTHONPATH="${PYTHONPATH}:/discover/nobackup/projects/gmao/obsdev/bkarpowi/pythonModules/"
```
For Tcsh/csh:
```
setenv PYTHONPATH ${PYTHONPATH}:/discover/nobackup/projects/gmao/obsdev/bkarpowi/pythonModules
```

Compiler options are handled autmoatically through cmake. On HPC systems this means loading the right set of modules. For example, if you would like pycrtm compiled with intel, you would load the same intel modules you used to build crtm. 

---------------------------------------------------------------------------------------- 

## 2. Tests/Examples:

A few test cases have been developed using input/output grabbed from the CRTM test program.
Two basic scripts which will perform cases 1-4 (stored in $PWD/testCases/data/case[n].h5) from the CRTM test program on 4 OpenMP threads: 
* `$PWD/testCases/test_atms.py`
* `$PWD/testCases/test_cris.py`
These *should* just say Yay, and not produce any plots if successful. 

The following scripts will do the same thing, only this time load up the same 4 profiles multiple times to further test threading with 10 threads (turn your laptop into a space heater more or less).
* `$PWD/testCases/test_atms_threads.py`
* `$PWD/testCases/test_cris_threads.py`
These *should* just say Yay, and not produce any plots if successful. 


The following scripts will run CRTM without aerosols or clouds:
* `$PWD/testCases/test_atms_no_clouds.py`
* `$PWD/testCases/test_cris_no_clouds.py`

For those Jupyter notebook fans, there is even Jupyter notebook example simulating ATMS:
* `$PWD/testCases/test_atms.ipynb`

## 3. Importing 

```Python
from pyCRTM import pyCRTM, profilesCreate
```
---------------------------------------------------------------------------------------- 

## 4. Using the interface (designed to be pretty much like the RTTOV equivalent python library):

Create a profiles data structure using `profilesCreate(nprofiles, nlayers)` which will generate an object with user specified number of profiles, and number of layers in the profiles provided.
```Python
profiles = profilesCreate(4, 92) # will generate an empty object with 4 profiles each with 92 layers. 
```
Once initialized, the user will need to provide values for the desired profiles (see example scripts). Next, the user initializes a crtm instance, set desired parameters, and passes profiles to the CRTM:

```Python
crtmOb = pyCRTM()
crtmOb.coefficientPath = pathInfo['CRTM']['coeffs_dir']
crtmOb.sensor_id = sensor_id
crtmOb.nThreads = 4
crtmOb.profiles = profiles
```
Note that while the python interface allows for threads to be set, if your environment has the environment variable `OMP_NUM_THREADS` set, the number of threads will be set by the environment variable, and the argument passed through pyCRTM will be ignored.

Next, the instrument is loaded/number of channels in the output structure are initialized:

```Python
crtmOb.loadInst()
```

Next, the user can run either the forward model (runDirect), or the K-matrix (runK) 
```Python
crtmOb.runDirect()
crtmOb.runK()
```
Finally, the user can pull out desired parameters such as brightness temperatures, Jacobians, or Transmission along path (TauLevels - the derivative will be the weighting function):
```Python
# brightness temperature (nprofiles, nchan):
brightnessTemperature = crtmOb.Bt 

#Transmission (to compute weighting functions) ( nprofiles, nchan, nlayers)
Tau = crtmOb.TauLevels 

#Temperature, Water Vapo[u]r, and Ozone Jacobians ( npforfiles, nchan, nlayers)
O3_Jacobian = crtmOb.O3K
Water_Vapor_Jacobian = crtmOb.QK
Temperature_Jacobian = crtm.TK

#Emissivity (nprofiles, nchan)
Emissivity = crtmOb.surfEmisRefl
```
Futher detail on setting profiles can be seen in the testCases directory, however, for quick reference and explanation of profile object:
```Python
        profiles.Angles[i,0] =  # instrument zenith Angle
        profiles.Angles[i,1] =  # instrument azimuth Angle (optional)
        profiles.Angles[i,2] =  # Solar zenith Angle 
        profiles.Angles[i,3] =  # Solar Azimuth Angle 
        profiles.Angles[i,4] =  # Instrument scan angle (see CRTM documentation. e.g., https://ftp.emc.ncep.noaa.gov/jcsda/CRTM/CRTM_User_Guide.pdf)
        profiles.DateTimes[i,0] = #Year
        profiles.DateTimes[i,1] = #Month
        profiles.DateTimes[i,2] = #Day
        profiles.Pi[i,:] =        #Pressure levels/interfaces in hPa
        profiles.P[i,:] =         #Pressure layers 
        profiles.T[i,:] =         #Temperature Layers
        profiles.Q[i,:] =         #Specific humidity  relative to dry air in g/kg
        profiles.O3[i,:] =        #Ozone concentration ppmv (dry air)
        profiles.clouds[i,:,0,0] = #cloud concentration kg/m**2 (optional)
        profiles.clouds[i,:,0,1] = #cloud effective radius microns (optional)
        profiles.aerosols[i,:,0,0] = #aerosol concentration kg/m**2 (optional)
        profiles.aerosols[i,:,0,1] = #aerosol effective radius microns (optional)
        profiles.aerosolType[i] =    #integer representing aerosol type (optional !  1 = Dust 2 = Sea salt-SSAM  
                                     # 3 = Sea salt-SSCM 4 = Sea salt-SSCM2 5 = Sea salt-SSCM3  6 = Organic carbon 7 = Black carbon 8 = Sulfate)
        profiles.cloudType[i] =      #integer representing cloud type (optional)
        profiles.cloudFraction[i,:] = #cloud fraction (optional)
        profiles.climatology[i] =    #integer representing climatology 1-6 modtran style climatological profiles (optional)
        profiles.surfaceFractions[i,:] = # fractions of land, water, snow, ice
        profiles.surfaceTemperatures[i,:] = #surface temperatures (K) of land, water, snow, ice
        profiles.Salinity[i] =              #salinity in PSU
        profiles.windSpeed10m[i] =          #10m windspeed m/s
        profiles.LAI[i] =                   #Leaf area index (optional)
        profiles.windDirection10m[i] =      #10m wind direction (note opposite of atmospheric convention, uses oceanongrapher convention) deg E from N
        # land, soil, veg, water, snow, ice
        profiles.surfaceTypes[i,0] = #land classification type index refer to CRTM documentation (optional)
        profiles.surfaceTypes[i,1] = #soil classification type index refer to CRTM documentation (optional)
        profiles.surfaceTypes[i,2] = #vegetation type classification type index refer to CRTM documentation (optional)
        profiles.surfaceTypes[i,3] = # water type (1=Sea water) check CRTM documentation for others.
        profiles.surfaceTypes[i,4] = # snow type 1=old snow 2=new snow (optional).
        profiles.surfaceTypes[i,5] = # ice type 1= new ice (optional) (optional)

```


---------------------------------------------------------------------------------------- 

