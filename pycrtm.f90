MODULE pycrtm 

REAL(KIND=8), ALLOCATABLE :: outTransmission(:, :, :) ! outTransmission(N_profiles, nChan, N_Layers)

REAL(KIND=8), ALLOCATABLE :: aerosolEffectiveRadius(:,:,:) !(N_Profiles,N_layers, N_aerosols)
REAL(KIND=8), ALLOCATABLE :: aerosolConcentration(:,:,:)   !(N_profiles,N_layers, N_aerosols)
INTEGER,      ALLOCATABLE :: aerosolType(:,:)                   !(N_Profiles, N_aerosols)


REAL(KIND=8), ALLOCATABLE :: cloudEffectiveRadius(:,:,:) !(N_Profiles,N_layers, N_clouds)
REAL(KIND=8), ALLOCATABLE :: cloudConcentration(:,:,:)   !(N_profiles,N_layers, N_clouds)
REAL(KIND=8), ALLOCATABLE :: cloudFraction(:,:)          !(N_profiles,N_layers)
INTEGER,      ALLOCATABLE :: cloudType(:,:)                   !(N_Profiles, N_clouds)

REAL(KIND=8), ALLOCATABLE :: emissivityReflectivity(:,:,:) ! 2,N_profiles, nChan

CONTAINS

SUBROUTINE wrap_forward( coefficientPath, sensor_id_in, IRwaterCoeff_File, MWwaterCoeff_File, & 
                        output_tb_flag, output_transmission_flag,  zenithAngle, scanAngle, azimuthAngle, solarAngle, &
                        output_emissivity_flag, use_passed_emissivity, & 
                        year, month, day, & 
                        nChan, N_Profiles, N_LAYERS, N_trace, &
                        pressureLevels, pressureLayers, temperatureLayers, & 
                        traceConcLayers, trace_IDs, & 
                        climatology, & 
                        surfaceTemperatures, surfaceFractions, LAI, salinity,  windSpeed10m, windDirection10m, & 
                        landType, soilType, vegType, waterType, snowType, iceType, nthreads, &  
                        outTb )      

  ! ============================================================================
  ! STEP 1. **** ENVIRONMENT SETUP FOR CRTM USAGE ****
  !
  ! MODULE usage
  USE CRTM_MODULE
  ! Disable all implicit typing
  IMPLICIT NONE
  ! ============================================================================
  ! variables for interface
  CHARACTER(len=*), INTENT(IN) :: coefficientPath
  CHARACTER(len=*), INTENT(IN) :: sensor_id_in
  CHARACTER(len=*), INTENT(IN) :: IRwaterCoeff_File
  CHARACTER(len=*), INTENT(IN) :: MWwaterCoeff_File
  LOGICAL,          INTENT(IN) :: output_tb_flag, output_transmission_flag, output_emissivity_flag 
  LOGICAL,          INTENT(IN) :: use_passed_emissivity
  ! The scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  INTEGER,      INTENT(IN) :: nChan, N_Profiles, N_Layers, N_trace
  REAL(KIND=8), INTENT(IN) :: zenithAngle(n_profiles), scanAngle(n_profiles) 
  REAL(KIND=8), INTENT(IN) :: azimuthAngle(n_profiles), solarAngle(n_profiles,2)
  INTEGER,      INTENT(IN) :: year(n_profiles), month(n_profiles), day(n_profiles) 
  REAL(KIND=8), INTENT(IN) :: pressureLevels(N_profiles, N_LAYERS+1)
  REAL(KIND=8), INTENT(IN) :: pressureLayers(N_profiles, N_LAYERS), temperatureLayers(N_Profiles,N_Layers)
  REAL(KIND=8), INTENT(IN) :: traceConcLayers(N_Profiles,N_layers,N_trace)
  INTEGER,      INTENT(IN) :: trace_IDs(N_trace)
  INTEGER,      INTENT(IN) :: climatology(N_profiles)
  REAL(KIND=8), INTENT(IN) :: surfaceTemperatures(N_Profiles,4), surfaceFractions(N_profiles, 4)
  REAL(KIND=8), INTENT(IN) :: LAI(N_Profiles), salinity(N_Profiles),  windSpeed10m(N_Profiles), windDirection10m(N_Profiles)
  INTEGER,      INTENT(IN) :: landType(N_Profiles), soilType(N_Profiles), vegType(N_Profiles), waterType(N_Profiles)
  INTEGER,      INTENT(IN) :: snowType(N_Profiles), iceType(N_Profiles) 
  INTEGER,      INTENT(IN) :: nthreads
  REAL(KIND=8), INTENT(OUT) :: outTb(N_Profiles,nChan) 
  CHARACTER(len=256), DIMENSION(1) :: sensor_id
  ! --------------------------
  ! Some non-CRTM-y Parameters
  ! --------------------------
  CHARACTER(*), PARAMETER :: SUBROUTINE_NAME   = 'wrap_forward'
  CHARACTER(*), PARAMETER :: PROGRAM_VERSION_ID = '0.01'



  ! ============================================================================
  ! STEP 2. **** SET UP SOME PARAMETERS FOR THE CRTM RUN ****
  !
  ! Sensor information
  INTEGER     , PARAMETER :: N_SENSORS = 1
  ! ============================================================================
  


  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: version
  INTEGER :: err_stat, alloc_stat
  INTEGER :: n_channels, N_clouds_crtm, N_aerosols_crtm
  INTEGER :: l, n
  LOGICAL :: cloudsOn, aerosolsOn

  ! ============================================================================
  ! STEP 3. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  ! 3a. Define the "non-demoninational" arguments
  ! ---------------------------------------------
  TYPE(CRTM_ChannelInfo_type)             :: chinfo(1)
  TYPE(CRTM_Geometry_type)                :: geo(1)


  ! 3b. Define the FORWARD variables
  ! --------------------------------
  TYPE(CRTM_Atmosphere_type),ALLOCATABLE  :: atm(:)
  TYPE(CRTM_Surface_type), ALLOCATABLE    :: sfc(:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)
  type(crtm_options_type)                 :: options(1)
  sensor_id(1) = sensor_id_in
  ! Program header
  ! --------------
  CALL CRTM_Version( Version )
  CALL Program_Message( SUBROUTINE_NAME, &
    'Running simulation.', &
    'CRTM Version: '//TRIM(Version) )

  IF (.not. allocated(emissivityReflectivity)) THEN
    IF ( output_emissivity_flag ) THEN
      allocate( emissivityReflectivity(2,N_profiles, nChan) )
    END IF
  END IF
  ! ============================================================================
  ! STEP 4. **** INITIALIZE THE CRTM ****
  !
  ! 4a. Initialise all the sensors at once
  ! --------------------------------------
  ! allocate globals in the MODULE based upon user selection through interface.
  CALL check_and_allocate_globals(output_transmission_flag, N_Profiles, nChan, N_layers)
  ! Figure out what needs allocating for Clouds and Aerosols. Are they on?
  CALL aerosols_and_clouds_on(N_aerosols_crtm, N_clouds_crtm, aerosolsOn, cloudsOn)
  
  err_stat = CRTM_Init( sensor_id,  chinfo, &
                        File_Path=coefficientPath, &
                        Load_CloudCoeff = cloudsOn, &  
                        Load_AerosolCoeff = aerosolsOn, &
                        IRwaterCoeff_File = IRwaterCoeff_File, & 
                        MWwaterCoeff_File = MWwaterCoeff_File, & 
                        Quiet=.True. )
  CALL check_allocate_status(err_stat,'Error Initializing CRTM')

  WRITE( *,'(/5x,"Initializing the CRTM...")' )

  ! 4b. Output some channel information
  ! -----------------------------------
  n_channels = CRTM_ChannelInfo_n_Channels(chinfo(1))
  WRITE( *,'(/5x,"Processing a total of ",i0," channels...")' ) n_channels


  ! ============================================================================

  ! ==========================================================================
  ! STEP 5. **** ALLOCATE STRUCTURE ARRAYS ****
  !
  ! 5a. Determine the number of channels
  !     for the current sensor
  ! ------------------------------------
  n_channels = CRTM_ChannelInfo_n_Channels(chinfo(1))

  
  ! 5b. Allocate the ARRAYS
  ! -----------------------
  ! Begin loop over profile
  ! ----------------------
  !$ CALL omp_set_num_threads(nthreads)
  !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(chinfo,emissivityReflectivity,outTb,outTransmission)& 
  !$OMP& shared(zenithAngle, scanAngle, azimuthAngle, solarAngle, output_tb_flag)&
  !$OMP& shared(nChan, N_Profiles, N_LAYERS, N_Clouds_crtm, N_aerosols_crtm, N_trace)&
  !$OMP& shared(use_passed_emissivity, output_emissivity_flag)&
  !$OMP& shared(pressureLevels, pressureLayers, temperatureLayers)& 
  !$OMP& shared(traceConcLayers,trace_IDs, output_transmission_flag)& 
  !$OMP& shared(aerosolEffectiveRadius, aerosolConcentration, aerosolType, cloudsOn, aerosolsOn)& 
  !$OMP& shared(cloudEffectiveRadius, cloudConcentration, cloudType, cloudFraction, climatology)& 
  !$OMP& shared(surfaceTemperatures, surfaceFractions, LAI, salinity,  windSpeed10m, windDirection10m)& 
  !$OMP& shared(landType, soilType, vegType, waterType, snowType, iceType, year, month, day)&
  !$OMP& NUM_THREADS(nthreads) 
  Profile_Loop: DO n = 1, N_Profiles
    n_channels = CRTM_ChannelInfo_n_Channels(chinfo(1))

    ! 5c. Allocate the STRUCTURE INTERNALS
    !     NOTE: Only the Atmosphere structures
    !           are allocated in this example
    ! ----------------------------------------
    ! The input FORWARD structure
    ALLOCATE( rts( n_channels, 1), STAT = alloc_stat )
    CALL check_allocate_status(alloc_stat, "Error allocating Solution rts(n_channels,1).")
    ! allocate 1 profile
    ALLOCATE( atm(1), STAT = alloc_stat )
    CALL check_allocate_status(alloc_stat, "Error allocating atm.")
    ! allocate 1 surface
    ALLOCATE( sfc(1), STAT = alloc_stat )
    CALL check_allocate_status(alloc_stat, "Error allocating sfc.")

    CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_trace, N_CLOUDS_crtm, N_AEROSOLS_crtm )
    CALL check_LOGICAL_status(ANY(.not. CRTM_Atmosphere_Associated(atm) ), "Failed in CRTM_Atmopsphere_Create")

    ! ==========================================================================
    ! STEP 6. **** ASSIGN INPUT DATA ****
    !
    ! 6a. Atmosphere and Surface input
    !     NOTE: that this is the hard part (in my opinion :o). The mechanism by
    !     by which the atmosphere and surface data are loaded in to their
    !     respective structures below was done purely to keep the step-by-step
    !     instructions in this program relatively "clean".
    ! ------------------------------------------------------------------------
    ! ...Profile data
    CALL set_profile(atm, n, climatology(n), pressureLevels(n,:), pressureLayers(n,:), temperatureLayers(n,:),&
                         traceConcLayers(n,:,:), trace_IDs(:), &
                         N_trace, N_aerosols_crtm, N_clouds_crtm, aerosolsOn, cloudsOn)

    ! 6b. Geometry input
    ! ------------------
    CALL CRTM_Geometry_SetValue( geo, &
                                 year = year(n), & 
                                 month = month(n), & 
                                 day = day(n), & 
                                 Sensor_Zenith_Angle  = zenithAngle(n),   &
                                 Sensor_Scan_Angle    = scanAngle(n),     & 
                                 Sensor_Azimuth_Angle = azimuthAngle(n),  &  
                                 Source_Zenith_Angle  = solarAngle(n,1),  & 
                                 Source_Azimuth_Angle = solarAngle(n,2) )
 
    ! ==========================================================================
    ! 4a.1 Profile #1
    ! ---------------
    ! set the surface properties for the profile.
    CALL set_surface(sfc, surfaceFractions(n,:), landType(n), surfaceTemperatures(n,:), LAI(n), & 
                           soilType(n), vegType(n), waterType(n), snowType(n), iceType(n), &
                           windSpeed10m(n), windDirection10m(n), salinity(n))
    ! ==========================================================================
    ! STEP 8. **** CALL THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
    !
    ! 8a. The forward model
    ! ---------------------

    ! Need this to get transmission out of solution, otherwise won't be allocated !!!
    CALL crtm_rtsolution_create( rts, n_layers )
    CALL check_LOGICAL_status( any(.not. crtm_rtsolution_associated( rts ) ),'rts failed to create.') 

    CALL crtm_options_create( options, nChan )
    CALL check_LOGICAL_status( any(.not. crtm_options_associated( options ) ),'options failed to create.' )
    CALL set_emissivity(options, n, use_passed_emissivity)

    !$ CALL omp_set_num_threads(1)
    err_stat = CRTM_Forward( atm        , &  ! Input
                             sfc        , &  ! Input
                             geo        , &  ! Input
                             chinfo     , &  ! Input
                             rts        , &  ! Output
                             options = options ) 

    CALL check_allocate_status(err_stat, "Error CALLing CRTM_Forward.")

    ! ============================================================================
    ! 8c. **** OUTPUT THE RESULTS TO SCREEN **** (Or transfer it into a series of arrays out of this thing!)
    !
    ! User should read the user guide or the source code of the routine
    ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
    ! select the needed variables for outputs.  These variables are contained
    ! in the structure RTSolution.
    IF (output_transmission_flag) THEN
        DO l=1,nChan
            outTransmission(n, l,1:n_layers) = & 
             dexp(-1.0*cumsum( rts(l,1)%Layer_Optical_Depth ) )
        END DO
    END IF

    IF( output_emissivity_flag ) THEN
        emissivityReflectivity(1,n,:) = rts(:,1)%Surface_Emissivity 
        emissivityReflectivity(2,n,:) = rts(:,1)%Surface_Reflectivity
    END IF 

    IF (output_tb_flag) THEN
        outTb(n,:) = rts(:,1)%Brightness_Temperature
    ELSE
        outTb(n,:) = rts(:,1)%Radiance
    END IF 
 
    
    ! ==========================================================================
    ! STEP 9. **** CLEAN UP FOR NEXT PROFILE ****
    !
    ! 9a. Deallocate the structures
    ! -----------------------------


    ! 9b. Deallocate the arrays
    ! -------------------------
    ! ==========================================================================
    CALL CRTM_Atmosphere_Destroy(atm)
    CALL crtm_rtsolution_destroy(rts)
    CALL crtm_options_destroy(options)
    deallocate(atm,stat=alloc_stat)
    CALL check_allocate_status(alloc_stat,"Atm failed to deallocate.")
    deallocate(sfc, stat=alloc_stat)
    CALL check_allocate_status(alloc_stat,"Sfc failed to deallocate.")
    DEALLOCATE(rts, STAT = alloc_stat)
    CALL check_allocate_status(alloc_stat,"Rts failed to deallocate.")
  END DO Profile_Loop
  !$OMP end PARALLEL DO
  
  ! ==========================================================================
  ! 10. **** DESTROY THE CRTM ****
  !
  WRITE( *, '( /5x, "Destroying the CRTM..." )' )
  err_stat = CRTM_Destroy( chinfo )
  CALL check_allocate_status(err_stat, 'Error Destroying the CRTM')
  write(*,*)'wrap_forward done!'

end SUBROUTINE wrap_forward

SUBROUTINE wrap_k_matrix( coefficientPath, sensor_id_in, IRwaterCoeff_File, MWwaterCoeff_File, & 
                        output_tb_flag, output_transmission_flag, zenithAngle, scanAngle, azimuthAngle, solarAngle, &  
                        output_emissivity_flag, use_passed_emissivity, & 
                        year, month, day, & 
                        nChan, N_profiles, N_LAYERS, N_trace, & 
                        pressureLevels, pressureLayers, temperatureLayers, & 
                        traceConcLayers, trace_IDs, & 
                        climatology, & 
                        surfaceTemperatures, surfaceFractions, LAI, salinity, windSpeed10m, windDirection10m, & 
                        landType, soilType, vegType, waterType, snowType, iceType, &  
                        nthreads, outTb, & 
                        temperatureJacobian, traceJacobian, skinK, emisK, reflK, &
                        windSpeedK, windDirectionK )      

  ! ============================================================================
  ! STEP 1. **** ENVIRONMENT SETUP FOR CRTM USAGE ****
  !
  ! MODULE usage
  USE CRTM_MODULE
  ! Disable all implicit typing
  IMPLICIT NONE
  ! ============================================================================
  


  ! --------------------------
  ! Some non-CRTM-y Parameters
  ! --------------------------
  CHARACTER(*), PARAMETER :: SUBROUTINE_NAME   = 'wrap_k_matrix'
  CHARACTER(*), PARAMETER :: PROGRAM_VERSION_ID = '0.01'

  ! variables for interface
  CHARACTER(len=*), INTENT(IN) :: coefficientPath
  CHARACTER(len=*), INTENT(IN) :: sensor_id_in
  CHARACTER(len=*), INTENT(IN) :: IRwaterCoeff_File
  CHARACTER(len=*), INTENT(IN) :: MWwaterCoeff_File
  LOGICAL, INTENT(IN) :: output_tb_flag, output_transmission_flag, output_emissivity_flag
  LOGICAL, INTENT(IN) :: use_passed_emissivity 
  INTEGER, INTENT(IN) :: nChan, N_profiles, N_Layers, N_trace 
  ! The scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  REAL(KIND=8), INTENT(IN) :: zenithAngle(N_profiles), scanAngle(N_profiles)
  REAL(KIND=8), INTENT(IN) :: azimuthAngle(N_profiles), solarAngle(N_profiles,2)
  INTEGER,      INTENT(IN) :: year(n_profiles), month(n_profiles), day(n_profiles)
  REAL(KIND=8), INTENT(IN) :: pressureLevels(N_profiles, N_Layers+1)
  REAL(KIND=8), INTENT(IN) :: pressureLayers(N_profiles, N_layers), temperatureLayers(N_profiles, N_layers)
  REAL(KIND=8), INTENT(IN) :: traceConcLayers(N_profiles, N_layers, N_trace)
  INTEGER,      INTENT(IN) :: trace_IDs(N_trace)
  INTEGER,      INTENT(IN) ::  climatology(N_profiles)
  REAL(KIND=8), INTENT(IN) :: surfaceTemperatures(N_profiles,4), surfaceFractions(N_profiles,4), LAI(N_profiles) 
  REAL(KIND=8), INTENT(IN) :: salinity(N_profiles), windSpeed10m(N_profiles), windDirection10m(N_profiles)
  INTEGER,      INTENT(IN) :: landType(N_profiles), soilType(N_profiles), vegType(N_profiles), waterType(N_profiles) 
  INTEGER,      INTENT(IN) :: snowType(N_profiles), iceType(N_profiles) 
  REAL(KIND=8), INTENT(OUT) :: outTb(N_profiles,nChan)
  REAL(KIND=8), INTENT(OUT) :: skinK(N_profiles,nChan,4), emisK(N_profiles,nChan), reflK(N_profiles,nChan)
  REAL(KIND=8), INTENT(OUT) :: windSpeedK(N_profiles,nChan), windDirectionK(N_profiles,nChan)
  REAL(KIND=8), INTENT(OUT) :: temperatureJacobian(N_profiles, nChan, N_LAYERS)
  REAL(KIND=8), INTENT(OUT) :: traceJacobian(N_profiles, nChan, N_LAYERS, N_trace)
  INTEGER,      INTENT(IN) :: nthreads

  CHARACTER(len=256) :: sensor_id(1)
  ! ============================================================================
  ! STEP 2. **** SET UP SOME PARAMETERS FOR THE CRTM RUN ****
  !
  
  ! Sensor information
  INTEGER, PARAMETER :: N_SENSORS = 1
  ! ============================================================================
  


  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: version
  INTEGER :: err_stat, alloc_stat
  INTEGER :: n_channels, N_aerosols_crtm, N_clouds_crtm
  INTEGER :: l, n, i_abs
  LOGICAL :: cloudsOn, aerosolsOn


  ! ============================================================================
  ! STEP 3. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  ! 3a. Define the "non-demoninational" arguments
  ! ---------------------------------------------
  TYPE(CRTM_ChannelInfo_type)             :: chinfo(1)
  TYPE(CRTM_Geometry_type)                :: geo(1)
  TYPE(crtm_options_type)                 :: options(1)
  ! 3b. Define the FORWARD variables
  ! --------------------------------
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm(:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: sfc(:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)
 
  ! 3c. Define the K-MATRIX variables
  ! ---------------------------------
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_K(:,:)
  TYPE(CRTM_Surface_type)   , ALLOCATABLE :: sfc_K(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_K(:,:)
  ! ============================================================================


  sensor_id(1) = sensor_id_in
  ! Program header
  ! --------------
  CALL CRTM_Version( Version )
  CALL Program_Message( SUBROUTINE_NAME, &
    'Running simulation.', &
    'CRTM Version: '//TRIM(Version) )

  IF (.not. allocated(emissivityReflectivity)) THEN
    IF ( output_emissivity_flag ) THEN
      allocate( emissivityReflectivity( 2,N_profiles, nChan ) )
    END IF
  END IF
 

  ! ============================================================================
  ! STEP 4. **** INITIALIZE THE CRTM ****
  !
  ! 4a. Initialise all the sensors at once
  ! --------------------------------------
  ! allocate globals in the MODULE based upon user selection through interface.
  CALL check_and_allocate_globals(output_transmission_flag, N_Profiles, nChan, N_layers)

  ! figure out how to allocate aerosols/clouds and are the even turned on by the user?
  CALL aerosols_and_clouds_on( N_aerosols_crtm, N_clouds_crtm, aerosolsOn, cloudsOn)

  WRITE( *,'(/5x,"Initializing the CRTM...")' )

  err_stat = CRTM_Init( sensor_id,  chinfo, &
                        File_Path=coefficientPath, &
                        Load_CloudCoeff = cloudsOn, &  
                        Load_AerosolCoeff = aerosolsOn, &
                        IRwaterCoeff_File = IRwaterCoeff_File, & 
                        MWwaterCoeff_File = MWwaterCoeff_File, &   
                        Quiet=.True. )
  CALL check_allocate_status(err_stat, 'Error initializing CRTM')

  ! 4b. Output some channel information
  ! -----------------------------------
  n_channels = CRTM_ChannelInfo_n_Channels(chinfo(1))
  WRITE( *,'(/5x,"Processing a total of ",i0," channels...")' ) n_channels
  WRITE( *,'(7x,i0," from ",a)' )  CRTM_ChannelInfo_n_Channels(chinfo(1)), TRIM(sensor_id(1))

  ! Begin loop over sensors
  ! ----------------------
  !$ CALL omp_set_num_threads(nthreads)
  !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(emissivityReflectivity,outTb)&
  !$OMP& shared(temperatureJacobian, output_tb_flag, output_transmission_flag)& 
  !$OMP& shared(traceJacobian,outTransmission)&
  !$OMP& shared(use_passed_emissivity, output_emissivity_flag)&
  !$OMP& shared(nChan, N_Layers,N_CLOUDS_crtm, N_AEROSOLS_crtm, N_trace)&
  !$OMP& shared(pressureLevels, pressureLayers, temperatureLayers)& 
  !$OMP& shared(traceConcLayers, trace_IDs, cloudsOn, aerosolsOn, zenithAngle,scanAngle,azimuthAngle,solarAngle)& 
  !$OMP& shared(aerosolEffectiveRadius, aerosolConcentration, aerosolType)& 
  !$OMP& shared(cloudEffectiveRadius, cloudConcentration, cloudType, cloudFraction, climatology)& 
  !$OMP& shared(surfaceTemperatures, surfaceFractions, LAI, salinity,  windSpeed10m, windDirection10m)& 
  !$OMP& shared(skinK, emisK, reflK, windSpeedK, windDirectionK)& 
  !$OMP& shared(landType, soilType, vegType, waterType, snowType, iceType)&
  !$OMP& shared(sensor_id,coefficientPath, chinfo, year, month, day, N_profiles)&
 
  !$OMP& NUM_THREADS(nthreads) 
  Profile_Loop: DO n = 1, N_profiles
  
    ! ==========================================================================
    ! STEP 5. **** ALLOCATE STRUCTURE ARRAYS ****
    !
    ! 5a. Determine the number of channels
    !     for the current sensor
    ! ------------------------------------
    n_channels = CRTM_ChannelInfo_n_Channels(chinfo(1))

    
    ! 5b. Allocate the ARRAYS
    ! -----------------------
    allocate(atm(1), STAT = alloc_stat)
    CALL check_allocate_status(alloc_stat,'Error allocating atm(1)')

    allocate(sfc(1), STAT = alloc_stat)
    CALL check_allocate_status(alloc_stat,'Error allocating sfc(1)')

    ALLOCATE( rts( n_channels, 1 ), STAT = alloc_stat )
    CALL check_allocate_status(alloc_stat,'Error allocating rts(n_channels,1)')

    ALLOCATE( atm_K( n_channels, 1 ), STAT = alloc_stat )
    CALL check_allocate_status(alloc_stat,'Error allocating atm_k')

    ALLOCATE( sfc_K( n_channels, 1 ), STAT = alloc_stat )
    CALL check_allocate_status(alloc_stat,'Error allocating sfc_k')

    ALLOCATE( rts_K( n_channels, 1 ), STAT = alloc_stat )
    CALL check_allocate_status(alloc_stat,'Error allocating rts_k')

    ! 5c. Allocate the STRUCTURE INTERNALS
    !     NOTE: Only the Atmosphere structures
    !           are allocated in this example
    ! ----------------------------------------
    ! The input FORWARD structure
    CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_trace, N_CLOUDS_crtm, N_AEROSOLS_crtm )
    CALL check_LOGICAL_status(ANY(.NOT. CRTM_Atmosphere_Associated(atm)), 'Error in CRTM_Atmosphere_Create Atm()')

    ! The output K-MATRIX structure
    CALL CRTM_Atmosphere_Create( atm_K, N_LAYERS, N_trace, N_CLOUDS_crtm, N_AEROSOLS_crtm )
    CALL check_LOGICAL_status(ANY(.NOT. CRTM_Atmosphere_Associated(atm_K)),  'Error in CRTM_Atmosphere_Create Atm_k()')

    CALL crtm_rtsolution_create( rts, n_layers )
    CALL check_LOGICAL_status(any(.not. crtm_rtsolution_associated( rts )),  'Error in crtm_rtsolution_create rts()')
    
    CALL crtm_rtsolution_create( rts_k, n_layers )
    CALL check_LOGICAL_status(any(.not. crtm_rtsolution_associated( rts_k )),  'Error in crtm_rtsolution_create rts_k()')

    ! ==========================================================================
    ! STEP 6. **** ASSIGN INPUT DATA ****
    !
    ! 6a. Atmosphere and Surface input
    !     NOTE: that this is the hard part (in my opinion :o). The mechanism by
    !     by which the atmosphere and surface data are loaded in to their
    !     respective structures below was done purely to keep the step-by-step
    !     instructions in this program relatively "clean".
    ! ------------------------------------------------------------------------
    ! ...Profile data
    CALL set_profile(atm, n, climatology(n), pressureLevels(n,:), pressureLayers(n,:), temperatureLayers(n,:),&
                         traceConcLayers(n,:,:), trace_IDs(:), &
                         N_trace, N_aerosols_crtm, N_clouds_crtm, aerosolsOn, cloudsOn)
 

    ! 6b. Geometry input
    ! ------------------
    ! All profiles are given the same value
    CALL CRTM_Geometry_SetValue( geo, &
                                 year = year(n), & 
                                 month = month(n), & 
                                 day = day(n), & 
                                 Sensor_Zenith_Angle  = zenithAngle(n),   &
                                 Sensor_Scan_Angle    = scanAngle(n),     & 
                                 Sensor_Azimuth_Angle = azimuthAngle(n),  &  
                                 Source_Zenith_Angle  = solarAngle(n,1),  & 
                                 Source_Azimuth_Angle = solarAngle(n,2) )
    ! ==========================================================================

    ! ==========================================================================
    ! STEP 7. **** INITIALIZE THE K-MATRIX ARGUMENTS ****
    !
    ! 7a. Zero the K-matrix OUTPUT structures
    ! ---------------------------------------
    
    CALL CRTM_Atmosphere_Zero( atm_K )
    CALL CRTM_Surface_Zero( sfc_K )

    ! 7b. Inintialize the K-matrix INPUT so
    !     that the results are dTb/dx
    ! -------------------------------------
    IF (output_tb_flag) THEN
        rts_K%Radiance               = ZERO
        rts_K%Brightness_Temperature = ONE
    ELSE 
        rts_K%Radiance               = ONE
        rts_K%Brightness_Temperature = ZERO
    END IF
    ! ==========================================================================
    
    ! 4a.1 Profile #1
    ! ---------------
    ! ...Land surface CHARACTERistics
     CALL set_surface(sfc, surfaceFractions(n,:), landType(n), surfaceTemperatures(n,:), LAI(n), & 
                           soilType(n), vegType(n), waterType(n), snowType(n), iceType(n), &
                           windSpeed10m(n), windDirection10m(n), salinity(n))
    ! ==========================================================================
    ! STEP 8. **** CALL THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
    !
    ! 8b. The K-matrix model
    ! ----------------------
    CALL crtm_options_create( options, nChan )
    CALL check_LOGICAL_status( any(.not. crtm_options_associated( options ) ),'options failed to create' )
    CALL set_emissivity(options, n,  use_passed_emissivity)

    !$ CALL omp_set_num_threads(1)
    err_stat = CRTM_K_Matrix( atm        , &  ! FORWARD  Input
                              sfc        , &  ! FORWARD  Input
                              rts_K      , &  ! K-MATRIX Input
                              geo        , &  ! Input
                              chinfo     , &  ! Input
                              atm_K      , &  ! K-MATRIX Output
                              sfc_K      , &  ! K-MATRIX Output
                              rts        , &  ! FORWARD  Output
                              Options=options)
    CALL check_allocate_status(err_stat,'Error CALLing the CRTM K-Matrix Model')    

    ! ==========================================================================
    ! STEP 9. **** CLEAN UP FOR NEXT SENSOR ****
    !
    ! 9a. Deallocate the structures
    ! -----------------------------
    

    ! 9b. Deallocate the arrays
    ! -------------------------
    ! transfer jacobians out
    DO l=1,nChan
        temperatureJacobian(n, l, 1:n_layers) = atm_k(l, 1)%Temperature(1:n_layers)
        !jacobians of H2O, O3, etc... will be determined by the order in which they were assigned in atm. 
        DO i_abs=1,N_trace
            traceJacobian(n,l, 1:n_layers,i_abs) = atm_k(l,1)%Absorber(1:n_layers,i_abs)
        END DO
        skinK(n,l,1) = sfc_K(l,1)%Land_Temperature
        skinK(n,l,2) = sfc_K(l,1)%Water_Temperature
        skinK(n,l,3) = sfc_K(l,1)%Ice_Temperature
        skinK(n,l,4) = sfc_K(l,1)%Snow_Temperature
        windSpeedK(n,l) = sfc_K(l,1)%Wind_Speed
        windDirectionK(n,l) = sfc_K(l,1)%Wind_Direction
        emisK(n,l) = RTS_K(l,1)%Surface_Emissivity
        reflK(n,l) = RTS_K(l,1)%Surface_Reflectivity
        IF (output_transmission_flag) then 
             outTransmission(n, l,1:n_layers) = & 
             dexp(-1.0* cumsum( rts(l,1)%Layer_Optical_Depth ) ) 
        END IF
    END DO

    IF (output_tb_flag) THEN
        outTb(n,:) = rts(:,1)%Brightness_Temperature
    ELSE
        outTb(n,:) = rts(:,1)%Radiance
    END IF

    IF (output_emissivity_flag) THEN  
        emissivityReflectivity(1,n,:) = rts(:,1)%Surface_Emissivity
        emissivityReflectivity(2,n,:) = rts(:,1)%Surface_Reflectivity
    END IF

    CALL CRTM_Atmosphere_Destroy(atm)
    CALL CRTM_Atmosphere_Destroy(atm_k)
    CALL crtm_options_destroy(options)
    DEALLOCATE(atm_k, STAT = alloc_stat)
    CALL check_allocate_status(alloc_stat, 'Atm_k deallocate failed')

    DEALLOCATE(rts_K, STAT = alloc_stat)
    CALL check_allocate_status(alloc_stat, 'rts_k deallocate failed')

    DEALLOCATE(sfc_k, STAT = alloc_stat)
    CALL check_allocate_status(alloc_stat, 'sfc_k deallocate failed')

    DEALLOCATE(rts, STAT = alloc_stat)
    CALL check_allocate_status(alloc_stat, 'rts deallocate failed')

    DEALLOCATE(atm, STAT = alloc_stat)
    CALL check_allocate_status(alloc_stat, 'atm deallocate failed')

    DEALLOCATE(sfc, STAT = alloc_stat)
    CALL check_allocate_status(alloc_stat, 'sfc deallocate failed')
    ! ==========================================================================

  END DO Profile_Loop
  !$OMP END PARALLEL DO
  ! ==========================================================================
  ! 10. **** DESTROY THE CRTM ****
  !
  WRITE( *, '( /5x, "Destroying the CRTM..." )' )
  err_stat = CRTM_Destroy( chinfo )
  CALL check_allocate_status(err_stat, 'Error destroying the CRTM.')
  ! ==========================================================================
END SUBROUTINE wrap_k_matrix

  SUBROUTINE check_and_allocate_globals(output_transmission_flag, N_Profiles, nChan, N_layers)
  LOGICAL, INTENT(IN) :: output_transmission_flag
  INTEGER, INTENT(IN) :: N_profiles, nChan, N_layers
  IF(output_transmission_flag) then
    IF ( allocated(outTransmission) ) deallocate(outTransmission)
    allocate( outTransmission(N_Profiles, nChan, N_layers) ) 
  END IF
  end SUBROUTINE check_and_allocate_globals

  SUBROUTINE aerosols_and_clouds_on(N_aerosols_crtm, N_clouds_crtm, aerosolsOn, cloudsOn)

  INTEGER, INTENT(OUT) :: N_AEROSOLS_crtm, N_CLOUDS_crtm
  LOGICAL, INTENT(OUT) :: aerosolsOn, cloudsOn
  INTEGER :: shp(2)
  IF( .not. allocated(aerosolType) ) THEN
    N_AEROSOLS_crtm = 0
    aerosolsOn = .False.
  ELSE
    shp = shape(aerosolType)
    N_AEROSOLS_crtm = shp(2)
    aerosolsOn = .True. 
  END IF

  IF( .not. allocated(cloudType) ) THEN
    N_CLOUDS_crtm = 0
    cloudsOn = .False.
  ELSE
    shp = shape(cloudType)
    N_CLOUDS_crtm = shp(2)  
    cloudsOn = .True. 
  END IF
 
  END SUBROUTINE aerosols_and_clouds_on

  SUBROUTINE set_emissivity(options, n, use_passed_emissivity )
    USE crtm_MODULE
    IMPLICIT NONE
    type(crtm_options_type), INTENT(INOUT) :: options(1)
    INTEGER :: n
    LOGICAL:: use_passed_emissivity

    IF ( .not. use_passed_emissivity ) THEN
        Options(1)%Use_Emissivity = .false.   ! compute it
    ELSE
        Options(1)%Use_Emissivity = .true.    ! user supplied
        Options(1)%Emissivity(:) = emissivityReflectivity(1,n,:)
    END IF 

    IF ( .not. use_passed_emissivity ) THEN
        Options(1)%Use_Direct_Reflectivity = .false.
    ELSE
        Options(1)%Use_Direct_Reflectivity = .true.  ! 1: User-supplied
        Options(1)%Direct_Reflectivity(:) = emissivityReflectivity(2,n,:) 
    END IF
  END SUBROUTINE set_emissivity

  SUBROUTINE check_allocate_status(alloc_stat,message)
    INTEGER :: alloc_stat
    CHARACTER(len=*) :: message
    IF ( alloc_stat /= 0 ) THEN
      write(*,*) message
      STOP
    END IF
  end SUBROUTINE check_allocate_status

  SUBROUTINE check_LOGICAL_status(stat,message)
    LOGICAL :: stat
    CHARACTER(len=*) :: message
    IF ( stat ) THEN
      write(*,*) message
      STOP
    END IF
  end SUBROUTINE check_LOGICAL_status

  SUBROUTINE set_profile(atm, n, climatology, pressureLevels, pressureLayers, temperatureLayers,&
                         traceConcLayers, trace_IDs, & 
                         N_trace, N_aerosols_crtm, N_clouds_crtm, aerosolsOn, cloudsOn)
  USE CRTM_MODULE
  TYPE(CRTM_Atmosphere_type), INTENT(INOUT) :: atm(:)
  INTEGER :: n,climatology
  REAL(KIND=8) :: pressureLevels(:), pressureLayers(:), temperatureLayers(:), traceConcLayers(:,:)
  INTEGER :: trace_IDs(:) 
  INTEGER :: N_trace, N_aerosols_crtm, N_clouds_crtm
  LOGICAL :: aerosolsOn, cloudsOn
  INTEGER :: i_abs,species  
 
    atm(1)%Climatology = climatology
    atm(1)%Level_Pressure = pressureLevels(:)
    atm(1)%Pressure = pressureLayers(:)
    atm(1)%Temperature = temperatureLayers(:)
   
    DO i_abs = 1,N_trace 
      atm(1)%Absorber(:,i_abs)      = traceConcLayers(:,i_abs)
      atm(1)%Absorber_Id(i_abs)     = trace_IDs(i_abs)
      IF( trace_IDs(i_abs) == H2O_ID ) THEN 
        atm(1)%absorber_units(i_abs) = MASS_MIXING_RATIO_UNITS
      ELSE 
        atm(1)%absorber_units(i_abs)  = VOLUME_MIXING_RATIO_UNITS
      END IF
    END DO

    IF( aerosolsOn )  THEN
      DO species = 1, N_aerosols_crtm
        atm(1)%Aerosol(species)%Type                = aerosolType(n, species)
        atm(1)%Aerosol(species)%Effective_Radius(:) = aerosolEffectiveRadius(n, :, species)
        atm(1)%Aerosol(species)%Concentration(:)    = aerosolConcentration(n, :, species)
      END DO
    END IF
    IF( cloudsOn ) THEN
      DO species = 1, N_clouds_crtm
        atm(1)%Cloud(species)%Type                = cloudType(n, species)
        atm(1)%Cloud(species)%Effective_Radius(:) = cloudEffectiveRadius(n, :, species)
        atm(1)%Cloud(species)%Water_Content(:)    = cloudConcentration(n, :, species)
      END DO
      atm(1)%Cloud_Fraction(:)            = cloudFraction(n,:)
    END IF

  END SUBROUTINE set_profile


  SUBROUTINE set_surface(sfc, surfaceFractions, landType, surfaceTemperatures, LAI, soilType, & 
                         vegType, waterType, snowType, iceType, windSpeed10m, windDirection10m, & 
                         salinity)
    USE CRTM_MODULE  
    TYPE(CRTM_Surface_type) :: sfc(:)
    REAL(KIND=8) :: surfaceFractions(:), surfaceTemperatures(:), LAI, windSpeed10m, windDirection10m, salinity
    INTEGER :: landType, soilType, vegType, waterType, snowType, iceType

    sfc%Land_Coverage     = surfaceFractions(1)
    sfc%Land_Type         = landType 
    sfc%Land_Temperature  = surfaceTemperatures(1)
    sfc%Lai               = LAI
    sfc%Soil_Type         = soilType 
    sfc%Vegetation_Type   = vegType 
    ! ...Water surface CHARACTERistics
    sfc%Water_Coverage    = surfaceFractions(2)
    sfc%Water_Type        = waterType 
    sfc%Water_Temperature = surfaceTemperatures(2)

    ! ...Snow coverage CHARACTERistics
    sfc%Snow_Coverage    = surfaceFractions(3)
    sfc%Snow_Type        = snowType 
    sfc%Snow_Temperature = surfaceTemperatures(3)
    ! ...Ice surface CHARACTERistics
    sfc%Ice_Coverage    = surfaceFractions(4)
    sfc%Ice_Type        = iceType 
    sfc%Ice_Temperature = surfaceTemperatures(4)

    sfc%Wind_Speed = windSpeed10m
    sfc%Wind_Direction = windDirection10m
    sfc%Salinity = salinity
  END SUBROUTINE set_surface


  FUNCTION cumsum(x) RESULT(xout)
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x
    REAL(KIND=8), DIMENSION(size(x)) :: xout
    INTEGER :: n,j
    n = size(x)
    xout(1) = x(1)
    DO j=2,n
        xout(j) = xout(j-1) + x(j)
    END DO
  END FUNCTION cumsum


END MODULE pycrtm
