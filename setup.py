import os, sys, tarfile, configparser, ssl, shutil, socket
from contextlib import closing
import urllib.request
from skbuild import setup
def main():
    #Completely remove previous _skbuild, because cache will remember previous interation and ignore you if you change something.
    try:shutil.rmtree('_skbuild')
    except:pass 
    
    #path of this file.
    scriptDir = os.path.split(os.path.abspath(__file__))[0]
    #read configuration
    download_coef, with_install, coef_path, crtm_install= readSetup('setup.cfg',scriptDir)
    #make sure we convert string to bool.
    if('T' in download_coef.upper()): download_coef= True
    else:download_coef = False
    if('T' in with_install.upper()): with_install = True
    else: with_install = False
    os.environ['CRTM_INSTALL'] = crtm_install
    shutil.copy(os.path.join(scriptDir,'setup.cfg'),os.path.join(scriptDir,'pyCRTM','pycrtm_setup.txt'))
    #If the user selects download_coef
    if(download_coef and not with_install): downloadAndMoveCoef(coef_path)
    elif(with_install): downloadAndMoveCoef(os.path.join(scriptDir,'coefficients'))
    else:
        print("Skipping Download of Coefficients. Hopefully, you know what this means, and downloaded the coefficients somewhere.")
    if(with_install):
        f = open(os.path.join(scriptDir,'MANIFEST.in'),'w')
        f.write('include pyCRTM/pycrtm_setup.txt crtm_io.py coefficients/*.* testCases/*.* testCases/data/*.*')
        f.close()
    else:
        f = open(os.path.join(scriptDir,'MANIFEST.in'),'w')
        f.write('include pyCRTM/pycrtm_setup.txt testCases/*.* testCases/data/*.*')
        f.close()
    requires=['numpy']
    setup(
        name="pyCRTM_JCSDA",
        version='1.0.0',
        description='Python wrapper for the CRTM.',
        author='Bryan Karpowicz',
        requires=requires,
        include_package_data=True,
        packages=['pycrtm'],
        py_modules=['crtm_io', 'pyCRTM'],
        package_data={'pyCRTM':['pyCRTM/setup.txt']})
    if(download_coef):
        shutil.rmtree('fix_crtm-internal_develop')
    if(with_install):
        shutil.rmtree('coefficients')
    os.remove('MANIFEST.in')
def readSetup(setup_file, scriptDir):
    cfg = configparser.ConfigParser()
    cfg.read( os.path.join(scriptDir,'setup.cfg') )
    crtm_install = cfg['Setup']['crtm_install']
    download_coef = cfg['Setup']['download']
    with_install = cfg['Setup']['coef_with_install']
    if('path' in list(cfg['Coefficients'].keys())):
        coef_path = cfg['Coefficients']['path']
    else:
        coef_path = os.path.join(scritpDir,'coefficients')
    return download_coef, with_install, coef_path, crtm_install

def downloadCoef(url):
    print("Downlading coefficients from {}".format(url))
    context = ssl._create_unverified_context()
    with closing(urllib.request.urlopen(url,context=context)) as r:
        with open('fix_REL-2.4.0.tgz', 'wb') as f:
            shutil.copyfileobj(r, f)
    print("Downloaded from {}".format(url))
def extractCoef():
    print("Untarring CRTM Tarball {}".format('fix_REL-2.4.0.tgz') )
    t = tarfile.open( 'fix_REL-2.4.0.tgz'  )
    t.extractall()
    t.close()
    print("Done Untarring.")
def moveCrtmCoefficients(installLocation):
    print("Moving Coefficients.") 
    if( not os.path.isdir(installLocation  ) ):
        os.makedirs( installLocation   )
    cwd = os.getcwd()
    p = os.path.join(cwd,'fix_crtm-internal_develop','SpcCoeff','Little_Endian')
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), installLocation)

    p = os.path.join(cwd,'fix_crtm-internal_develop','TauCoeff','ODPS','Little_Endian')
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), installLocation)

    p = os.path.join(cwd,'fix_crtm-internal_develop','CloudCoeff','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), installLocation)
    
    p = os.path.join(cwd,'fix_crtm-internal_develop','AerosolCoeff','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), installLocation)

    p = os.path.join(cwd,'fix_crtm-internal_develop','EmisCoeff','IR_Ice','SEcategory','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), installLocation)

    p = os.path.join(cwd,'fix_crtm-internal_develop','EmisCoeff','IR_Land','SEcategory','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), installLocation)

    p = os.path.join(cwd,'fix_crtm-internal_develop','EmisCoeff','IR_Snow','SEcategory','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), installLocation)

    p = os.path.join(cwd,'fix_crtm-internal_develop','EmisCoeff','IR_Water','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), installLocation)

    p = os.path.join(cwd,'fix_crtm-internal_develop','EmisCoeff','MW_Water','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), installLocation)

    p = os.path.join(cwd,'fix_crtm-internal_develop','EmisCoeff','VIS_Ice','SEcategory','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), installLocation)

    p = os.path.join(cwd,'fix_crtm-internal_develop','EmisCoeff','VIS_Land','SEcategory','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), installLocation)

    p = os.path.join(cwd,'fix_crtm-internal_develop','EmisCoeff','VIS_Snow','SEcategory','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), installLocation)

    p = os.path.join(cwd,'fix_crtm-internal_develop','EmisCoeff','VIS_Water','SEcategory','Little_Endian') 
    for f in os.listdir(p):
        shutil.copy(os.path.join(p,f), installLocation)
def downloadAndMoveCoef(coef_path):
    #if file is downloaded don't do it again.
    if(not os.path.exists('fix_REL-2.4.0.tgz') ):
        # if you're on discover, don't use ftp, because, well, you can't. Anyone else can talk to UCAR.
        #if('discover' in socket.gethostname()):
        url = 'https://gmao.gsfc.nasa.gov/gmaoftp/bkarpowi/crtm/fix_REL-2.4.0_le.tgz'
        #else:
        #    url = 'ftp://ftp.ucar.edu/pub/cpaess/bjohns/fix_REL-2.4.0.tgz'
        downloadCoef(url)
    extractCoef()
    moveCrtmCoefficients(coef_path)

if __name__ == "__main__":
    main()

