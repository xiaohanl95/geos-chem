!-------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!--------------------------------------------------------------------------


! This script can be added to the GEOS-Chem model to build a Plume-in-Grid
! model. The script mainly includes a trajectory track module (SUBROUTINE:
! lagrange_run) and a plume model (SUBROUTINE: plume_run). For any question
! about the code, please contact: Hongwei Sun (hongwei_sun@g.harvard.edu)


! Here are some tips to embed this script to the GEOS-Chem model
!
! First: add a new module (lagrange_mod.f90) into GEOS-Chem:
!(1) put the new module in the folder: GCClassic.v13/src/GEOS-Chem/GeosCore
!(2) edit CMakeLists.txt (in the GoesCore/): add lagrange_mod.F90 into the txt file
!(3) edit main.F90 (in GEOS-Chem/Interfaces/GCClassic/):
!	USE Lagrange_Mod
!	CALL lagrange_init(am_I_root, Input_Opt, State_Chm, State_Grid, State_Met, RC) ! before the time step loop
!	CALL plume_inject(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)  ! inside the time step loop
!	CALL lagrange_cleanup() ! After the time step loop
!
!
! Second: add new variable (PASV_LA) to the GEOS-Chem model: 
! (we added this new variable, but this is not required, you can also use some existed variables)
!(1) modify the input.geos:
!	------------------------+------------------------------------------------------
!	%%% ADVECTED SPECIES MENU %%%:
!	Species name            : PASV_LA3
!	------------------------+------------------------------------------------------
!	%%% PASSIVE SPECIES MENU %%%:
!	Passive species #1    : PASV_LA         98.0      -1       1.0e-999
!(2) modify the species_database.yml:
!	PASV_LA:
! 	  FullName: PASV_LA
!  	  Is_Advected: true
!  	  Is_Gas: true
!  	  MW_g: 98.0


! -----Notes by Bingqing Zhang ---------------
! Dec 16, 2025, Bingqing Zhang
! This code was updated to read point-source emission information and to
! control activation of the Lagrangian model through the geoschem_config.yml
! file, removing the need for hard-coded emission settings.
!
! at a start, I assume all injection occur at the same location, only species are different, so I don't handle plume seg interaction
! Plume location will only initilize based on the first source location
!
! The desired structure
! for every injection
!   CALL lagrange_run(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
!   CALL plume_run(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
!   CALL lagrange_write_std( am_I_Root, RC )
!But the current structure, lagrange_run contains loop for every injection
!
! Known issue: in order to call do_chemistry module, I stored concentration of all chemical species in plume seg, which might require extremely large memory



MODULE Plume_list

  USE precision_mod
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Plume2d_list, Plume1d_list

  TYPE :: Plume2d_list
    INTEGER :: IsNew     = MISSING_INT ! 1: the plume is new injected
    INTEGER :: label     = MISSING_INT! injected rank

    REAL(fp) :: LON = MISSING, LAT = MISSING, LEV = MISSING
    REAL(fp) :: LENGTH = MISSING, ALPHA = MISSING
    REAL(fp) :: LIFE = MISSING
    REAL(fp) :: PDX = MISSING, PDY = MISSING
    REAL(fp), DIMENSION(:,:,:), POINTER :: CONCNT2d => NULL() ! [n_x_max, n_y_max, n_species]

    TYPE(Plume2d_list), POINTER :: next => NULL()
  END TYPE


  TYPE :: Plume1d_list
    INTEGER :: label         = MISSING_INT ! injected rank
    INTEGER :: Is_transfer   = MISSING_INT ! if =1, transfer to the host Euletian model

    REAL(fp) :: LON = MISSING, LAT = MISSING, LEV = MISSING
    REAL(fp) :: LENGTH = MISSING, ALPHA = MISSING
    REAL(fp) :: LIFE = MISSING
    REAL(fp) :: RA = MISSING, RB = MISSING, THETA = MISSING
    REAL(fp), DIMENSION(:,:), POINTER :: CONCNT1d => NULL() ! [n_slab_max,n_species]

    TYPE(Plume1d_list), POINTER :: next => NULL()
  END TYPE

END MODULE Plume_list



MODULE Lagrange_Mod

  USE Plume_list
  USE PRECISION_MOD
  USE ERROR_MOD
  USE ERRCODE_MOD
  USE PhysConstants,   ONLY : PI, Re, g0, AIRMW, AVO, BOLTZ
  USE TIME_MOD,        ONLY : GET_YEAR, GET_MONTH, GET_DAY, GET_HOUR, GET_MINUTE, GET_SECOND
  USE UNITCONV_MOD    
  USE INPUT_OPT_MOD,   ONLY : PlumeSource_t

  IMPLICIT NONE
  PRIVATE

  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: lagrange_init
  PUBLIC :: plume_inject
  PUBLIC :: lagrange_cleanup

  ! PUBLIC VARIABLES:
  PUBLIC :: n_x_max, n_y_max
  PUBLIC :: n_x_mid, n_y_mid
  PUBLIC :: n_slab_max, n_slab_25, n_slab_50, n_slab_75
  PUBLIC :: use_lagrange
  ! variables read from geoschem_config.yml
  LOGICAL                               :: use_lagrange
  LOGICAL                               :: plume_inject_on 
  LOGICAL                               :: plume_diag
  LOGICAL                               :: TROPP_sink
  INTEGER                               :: Num_of_sources
  INTEGER                               :: n_x_max            !number of x grids in 2D, should be (9 x odd)
  INTEGER                               :: n_y_max            !number of y grids in 2D, should be (9 x odd)
  INTEGER                               :: n_species          
  TYPE(PlumeSource_t), ALLOCATABLE      :: Plume_sources(:)
  REAL(fp)                              :: Dx_init
  REAL(fp)                              :: Dy_init
  REAL(fp)                              :: Length_init    ! m
  REAL(fp)                              :: Aircraft_speed ! m/s
  REAL(fp)                              :: plume_interval ! degree

  ! Other variables
  INTEGER               :: Volume_Sort  = 1 ! 1 = use SortList() function, transfer largest (not oldest) plume segment for volume criterion
  INTEGER               :: Calc_entropy = 1 ! 1 = turn on entropy calculation
  REAL(fp)		          :: Entropy0	 ! perfect entropy without diffusion
  INTEGER               :: IIPAR, JJPAR, LLPAR
  INTEGER               :: n_x_mid, n_y_mid, n_x_mid2, n_y_mid2 
  INTEGER               :: n_x_max2, n_y_max2 
  INTEGER               :: n_slab_25, n_slab_50, n_slab_75
  INTEGER               :: n_slab_max, n_slab_max2
  INTEGER               :: N_parcel   ! 131        
  INTEGER               :: Num_inject, Num_Plume2d, Num_Plume1d, Num_dissolve     
  INTEGER               :: tt 
  INTEGER               :: N_total, N_stop_inject
  INTEGER               :: Stop_inject ! 1: stop injecting; 0: keep injecting
  INTEGER               :: nspec

  REAL(fp)              :: DX, DY
  REAL(fp) 		          :: mass_eu, mass_la, mass_la2
  REAL(fp)              :: Length_lat
  REAL(fp), POINTER     :: X_mid(:), Y_mid(:), P_mid(:)
  REAL(fp), POINTER     :: P_edge(:)
  REAL(fp), POINTER     :: X_edge(:), Y_edge(:)
  REAL(fp)              :: X_edge2, Y_edge2

  ! some parameter for sensitive test
  INTEGER, PARAMETER        :: N1_split = 5     ! Cross-section splitting
  INTEGER, PARAMETER        :: N2_split = 5     ! length splitting
  INTEGER, PARAMETER        :: Split_length = 1 ! how many times of DX
  REAL(fp), PARAMETER       :: Dissolve_critiria = 10*0.01
  REAL(fp), PARAMETER       :: Volume_percent    = 30*0.01
  REAL(fp), PARAMETER       :: Critical_day      = 28.0         ! [day]

  ! Species ID flags
  INTEGER :: id_SO2,  id_SO4,  id_OH,   id_O3,   id_NH3,   id_NH4,   id_H2O
  INTEGER :: id_NK01, id_SF01, id_AW01, id_H2SO4

  TYPE(Plume2d_list), POINTER :: Plume2d_tail, Plume2d_head
  TYPE(Plume1d_list), POINTER :: Plume1d_tail, Plume1d_head

  !-------------------------------------------------------------
  ! Some parameters to be retired 
  integer               :: id_PASV_LA3, id_PASV_LA2, id_PASV_LA 
  integer               :: id_PASV_EU2, id_PASV_EU
  integer, parameter    :: i_tracer  = 1
  integer, parameter    :: i_product = 2
  real(fp), parameter       :: Kchem = 1.0e-20_fp ! chemical reaction rate
  ! use Kchem = 1.0e-20_fp for >1 year simulation
  !-------------------------------------------------------------
CONTAINS

!-----------------------------------------------------------------

  SUBROUTINE lagrange_init(am_I_root, Input_Opt, State_Chm, State_Grid, State_Met, RC)

    USE Plume_list

    USE Input_Opt_Mod,   ONLY : OptInput, PlumeSource_t
    USE State_Met_Mod,   ONLY : MetState
    USE State_Chm_Mod,   ONLY : ChmState, Ind_
    USE State_Grid_Mod,  ONLY : GrdState
    USE Species_Mod,     ONLY : SpcConc
    USE TIME_MOD,        ONLY : GET_TS_DYN
    USE TIME_MOD,        ONLY : GET_YEAR, GET_MONTH, GET_DAY, GET_HOUR, GET_MINUTE, GET_SECOND

    LOGICAL,        INTENT(IN)            :: am_I_Root   ! Are we on the root CPU
    TYPE(MetState), intent(in)            :: State_Met
    TYPE(ChmState), intent(inout)         :: State_Chm
    TYPE(OptInput), intent(in)            :: Input_Opt
    TYPE(GrdState), INTENT(IN)            :: State_Grid  ! Grid State objectgg
    INTEGER,        INTENT(OUT)           :: RC         ! Success or failure
    ! Pointers
    TYPE(SpcConc), POINTER                :: Spc(:)
    !TYPE(PlumeSource_t), ALLOCATABLE      :: Plume_sources(:)

    INTEGER                       :: i_box, i_slab
    INTEGER                       :: ii, jj, kk
    INTEGER                       :: id_tracer, N
    INTEGER                       :: i_lon, i_lat, i_lev            !1:IIPAR
    CHARACTER(LEN=255)            :: spc_name
    CHARACTER(LEN=255)            :: FILENAME, FileEntropy, File996
    CHARACTER(LEN=255)            :: FILENAME2, FILENAME3
    CHARACTER(LEN=255)            :: ErrMsg, ThisLoc, LOC
    
    REAL(fp)                      :: lon1, lat1, lon2, lat2
    REAL(fp)                      :: box_lon_edge, box_lat_edge
    REAL(fp)                      :: curr_lon, curr_lat, curr_lev
            
    !REAL(fp), dimension(:,:,:), allocatable :: box_concnt_2D
    !REAL(fp), dimension(:,:), allocatable   :: box_concnt_1D

    !TYPE(Plume2d_list), POINTER :: Plume2d_new, PLume2d, Plume2d_prev
    !TYPE(Plume1d_list), POINTER :: Plume1d_new, Plume1d, Plume1d_prev
    REAL(fp)                      :: Dt

    RC                 =   GC_SUCCESS
    ErrMsg             =   ''
    ThisLoc            =   ' -> at Lagrangian Init (in module GeosCore/lagrange_mod.F90)'

    Num_inject         =    0
    Num_Plume2d        =    0
    Num_Plume1d        =    0
    Num_dissolve       =    0

    n_species          =              State_Chm%nSpecies
    ! Point to chemical species array. 
    ! The unit has been convereted to molec/cm3 before the function and will be convereted back to v/v dry air after the function finish
    Spc                =>             State_Chm%Species
                     
    ! Copy all neccessary variable from geoschem.yml here
    use_lagrange        =             Input_Opt%LagrangianModel_Activate
    plume_inject_on     =             Input_Opt%PlumeInjection_Activate
    plume_diag          =             Input_Opt%PlumeInjection_Diag
    TROPP_sink          =             Input_Opt%TropSink_Activate
    Num_of_sources      =             Input_Opt%Plume_sources_num
    Length_init         =             Input_Opt%Initial_length*1000.0   ! km to m
    Aircraft_speed      =             Input_Opt%Aircraft_speed  !m/s
    plume_interval      =             Input_Opt%plume_interval ! degree

    IF (Num_of_sources > 0) THEN
      ALLOCATE( Plume_sources(Num_of_sources), STAT=RC )
      IF (RC /= 0) THEN
          errMsg = 'Error allocating Plume_sources'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       END IF
       DO N = 1, Num_of_sources
          Plume_sources(N)%lat1        =  Input_Opt%Plume_sources(N)%lat1
          Plume_sources(N)%lat2        =  Input_Opt%Plume_sources(N)%lat2
          Plume_sources(N)%lon         =  Input_Opt%Plume_sources(N)%lon
          Plume_sources(N)%lev         =  Input_Opt%Plume_sources(N)%lev
          Plume_sources(N)%rate        =  Input_Opt%Plume_sources(N)%rate
          Plume_sources(N)%species     =  Input_Opt%Plume_sources(N)%species
       ENDDO
    ENDIF

    n_x_max                    =      Input_Opt%PlumeGrid2d_nx ! odd number to make sure n_x_mid to be integer
    n_y_max                    =      Input_Opt%PlumeGrid2d_ny ! odd number to make sure n_y_mid to be integer
    Dx_init                    =      Input_Opt%PlumeGrid2d_dx
    Dy_init                    =      Input_Opt%PlumeGrid2d_dy
    ! the odd number of n_x_max can ensure a center grid
    n_x_mid                    =      (n_x_max+1)/2 
    n_y_mid                    =      (n_y_max+1)/2  
    n_x_max2                   =      n_x_max+2
    n_y_max2                   =      n_y_max+2
    n_x_mid2                   =      (n_x_max2+1)/2 
    n_y_mid2                   =      (n_y_max2+1)/2  
    ! n_slab_max should be divided by 4, to ensure n_slab_25 is an integer.
    n_slab_max                 =      (INT(n_y_max/4)+1)*4 ! close to n_y_max, number of slabs in 1D
    n_slab_max2                =      n_slab_max+2
  ! add 2 more slab grid to containing background concentration

    ! Debug output, to make sure input is properly read
    ! WRITE( 6, 90  ) 'debug: PLUME MODEL SETTINGS, READ from geoschem_config.yml'
    ! WRITE( 6, 95  ) 'debug: ----------------'
    !WRITE( 6, 100 ) 'debug: Turn on plume injection? : ',                        &
    !                plume_inject_on
    !WRITE( 6, 100 ) 'debug: Turn on plume injection diagnostic? : ',             &
    !                plume_diag
    !WRITE( 6, 110 ) 'Plume injection diagnostic path?    : ',             &
    !                TRIM( Input_Opt%Plume_sources_diag_dir  )
    !WRITE( 6, 100 ) 'debug: Turn on lagrangian model? : ',                       &
    !                use_lagrange
    !WRITE( 6, 100 ) 'debug: Turn on tropospheric sink? : ',                      &
    !                TROPP_sink
    !WRITE(6,*) 'debug: Plume segmentation 2D grid:',                             &
    !    ' nx=', n_x_max,                                 &
    !    ' ny=', n_y_max,                                 &
    !    ' dx=', Dx_init,                                 &
    !    ' dy=', Dy_init
    !WRITE( 6, 105 ) 'debug: Number of sources    : ', Num_of_sources
    ! print information of each sources one by one
    !WRITE( 6, 95  ) 'debug: ------------------------------------------------'
    !WRITE(6, 120) 'debug: No.', 'lat1', 'lat2', 'lon', 'lev', 'rate', 'species'
    !DO N = 1, Num_of_sources
    !  WRITE(6,121) N,                                                    &
    !  Plume_sources(N)%lat1,                                   &
    !  Plume_sources(N)%lat2,                                   &
    !  Plume_sources(N)%lon,                                    &
    !  Plume_sources(N)%lev,                                    &
    !  Plume_sources(N)%rate,                                   &
    !  TRIM(Plume_sources(N)%species)
    !END DO
    !WRITE( 6, 95  ) '------------------------------------------------'
   
	
    Dt = GET_TS_DYN()
    ! N_parcel = NINT(132 *1000/Length_init /600*Dt)
    N_parcel = NINT(Aircraft_speed * Dt / Length_init)
    IF (N_parcel<1) N_parcel=1
    ! use for plume injection
    N_total    = 60 * 1.0e+5 / Length_init
    Length_lat = Length_init / 1.0e+5

    ! define how many plume segments will be injected
    ! -1 means keep inecting in the whole simulation
    N_stop_inject = -1 ! N_total ! -1
    Stop_inject = 0

    ! Check 2D domain grid
    IF(MOD(n_x_max,9).ne.0) WRITE(6,*)"*** ERROR,  n_x_max should be divisible by 9 ***,"
    IF(MOD(n_y_max,9).ne.0) WRITE(6,*)"*** ERROR,  n_y_max should be divisible by 9 ***"


    IIPAR = State_Grid%NX
    JJPAR = State_Grid%NY
    LLPAR = State_Grid%NZ
    DX = State_Grid%DX
    DY = State_Grid%DY
    X_edge => State_Grid%XEdge(:,1) 
    Y_edge => State_Grid%YEdge(1,:)
    P_edge => State_Met%PEDGE(1,1,:) 
    X_edge2       = X_edge(2)
    Y_edge2       = Y_edge(2) 
!----------Output Files need to be organized--------------
    FILENAME2   = 'Plume_lifetime_seconds.txt'
    OPEN( 484,      FILE=TRIM( FILENAME2   ), STATUS='REPLACE',  &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
    CLOSE(484)

    FILENAME3 = 'Plume_number.txt'
    OPEN( 261,      FILE=TRIM( FILENAME3   ), STATUS='REPLACE',  &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
    CLOSE(261)

    FileEntropy   = 'Plume_entropy.txt'
    OPEN( 487,      FILE=TRIM( FileEntropy   ), STATUS='REPLACE',  &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
    CLOSE(487)

    File996   = 'Plume_entropy_check.txt'
    OPEN( 996,      FILE=TRIM( File996  ), STATUS='REPLACE',  &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
    CLOSE(996)
!----------Output Files need to be organized--------------

    !allocate(box_concnt_2D(n_x_max, n_y_max, n_species))
    !allocate(box_concnt_1D(n_slab_max,n_species))


    !WRITE(6,'(a)') '--------------------------------------------------------'
    !WRITE(6,'(a)') '--------------------------------------------------------'
    !WRITE(6,'(a)') ' Initial Lagrnage Module (Using Dynamic time step)'
    !WRITE(6,'(a)') '--------------------------------------------------------'
    !WRITE(6,*) 'time step=', Dt
    !WRITE(6,*) 'Injected plume every time step: ', N_parcel
    !WRITE(6,'(a)') '--------------------------------------------------------'
    !WRITE(6,'(a)') '--------------------------------------------------------'

    ! Note (BZ): As a start, assume single location, but allows multiple types of species
    ! -----------------------------------------------------------
    ! create first node (head) for 2d linked list:
    ! 
    ! At the beginning, there are only one node in the 2d list, 
    ! so Plume2d_tail and Plume2d_head refer to the same node
    ! -----------------------------------------------------------

    IF (.not.plume_inject_on) THEN
      WRITE(6,'(a)') ' No plume injection, no lagrangian module configuration, skip initialization'
      RETURN
    ENDIF
    IF (plume_inject_on .AND. num_of_sources.lt.1) THEN
      ErrMsg = 'Plume injection is turned on but no source specified'
      CALL GC_Error( ErrMsg, RC, ThisLoc )
      RETURN
    ENDIF
    curr_lon    =  Plume_sources(1)%lon
    curr_lat    =  Plume_sources(1)%lat1
    curr_lev    =  Plume_sources(1)%lev
    i_lon = Find_iLonLat(curr_lon, DX, X_edge2)
    if(i_lon>IIPAR) i_lon=i_lon-IIPAR
    if(i_lon<1) i_lon=i_lon+IIPAR
    !WRITE(6,*) 'debug: ilon=', i_lon
    i_lat = Find_iLonLat(curr_lat, DY, Y_edge2)
    if(i_lat>JJPAR) i_lat=JJPAR
    if(i_lat<1) i_lat=1
    !WRITE(6,*) 'debug: ilat=', i_lat
    i_lev = Find_iPLev(curr_lev,P_edge)
    !WRITE(6,*) 'debug: curr_lev=', curr_lev
    !WRITE(6,*) 'debug: P_edge(1)=', curr_lev
    !WRITE(6,*) 'debug: ilev1=', i_lev
    if(i_lev>LLPAR) i_lev=LLPAR
    !WRITE(6,*) 'debug: ilev2=', i_lev

    ! Check that species units are in [kg]
    id_SO2= Ind_('SO2')
    id_SO4= Ind_('SO4')
    WRITE(6,'(a)') 'debug: Unit for SO2 is: ' // TRIM(UNIT_STR(Spc(id_SO2)%Units))
    WRITE(6,'(a)') 'debug: Unit for SO4 is: ' // TRIM(UNIT_STR(Spc(id_SO4)%Units))
    ! Initial Unit for all species is: molec/cm3
    IF ( Spc(id_SO2)%Units /= MOLECULES_SPECIES_PER_CM3 ) THEN
      ErrMsg = 'Incorrect species units: ' // TRIM(UNIT_STR(Spc(id_SO2)%Units))
      CALL GC_Error( ErrMsg, RC, ThisLoc )
    ENDIF


    IF (use_lagrange) THEN
      
      !WRITE(6,*) 'debug, BZ: Tracer id_PASV_EU2  = ', id_PASV_EU2, &
      !        '  Tracer Name: ', TRIM(State_Chm%SpcData(id_PASV_EU2)%Info%Name)
      !id_PASV_EU  = State_Chm%nAdvect-1
      !id_PASV_EU2 = State_Chm%nAdvect
      ! set initial background concentration of injected aerosol as 0 in GCM
      ! output the apecies' name for double check ???
      !State_Chm%Species(id_PASV_EU2)%Conc(:,:,:) = 0.0e+0_fp  ! [kg/kg]
      !State_Chm%Species(id_PASV_EU)%Conc(:,:,:)  = 0.0e+0_fp  ! [kg/kg]
      ! debug purpose, BZ, does species database and SpecData contains the same species and orders?
      ! WRITE(6,*) 'debug, BZ: Tracer id_PASV_EU  = ', id_PASV_EU, &
      !         '  Tracer Name: ', TRIM(State_Chm%SpcData(id_PASV_EU)%Info%Name)
        
      !WRITE(6,*) 'debug, BZ: Tracer id_PASV_EU2 = ', id_PASV_EU2, &
      !        '  Tracer Name: ', TRIM(State_Chm%SpcData(id_PASV_EU2)%Info%Name)

      ! injection is on, source larger than 0, add background to the plume grid
      ! Add ijected species


      ! Initialize the first plume segment
      ALLOCATE(Plume2d_tail)
      Plume2d_tail%IsNew = 1
      Plume2d_tail%label = 1

      !Plume2d_tail%LON = Inject_lon
      !Plume2d_tail%LAT = -29.95e+0_fp
      !Plume2d_tail%LEV = Inject_hPa
      Plume2d_tail%LON   = Plume_sources(1)%lon
      Plume2d_tail%LAT   = Plume_sources(1)%lat1
      Plume2d_tail%LEV   = Plume_sources(1)%lev

      Plume2d_tail%ALPHA  = 0.0e+0_fp
      Plume2d_tail%LIFE = 0.0e+0_fp

      Plume2d_tail%LENGTH = Length_init ! m 
      Plume2d_tail%PDX = Dx_init
      Plume2d_tail%PDY = Dy_init
      
      ALLOCATE(Plume2d_tail%CONCNT2d(n_x_max, n_y_max, n_species+num_of_sources)) ! add additional species here to track species that are newly added
      Plume2d_tail%CONCNT2d = 0.0e+0_fp ! The final unit looks like (g/m3)
      
      ! pass all species as from background as initial concentration
      DO N = 1, n_species
        Plume2d_tail%CONCNT2d(:,:,N) = Spc(N)%Conc(i_lon, i_lat, i_lev)  ! molec/cm3
      ENDDO

      ! add injected species to the center of plume grid
      DO N = 1, num_of_sources
        spc_name = Plume_sources(N)%species
        id_tracer   = Ind_(spc_name)
        Plume2d_tail%CONCNT2d(n_x_mid,n_y_mid,id_tracer) = Plume2d_tail%CONCNT2d(n_x_mid,n_y_mid,id_tracer)   &
                    + (Plume2d_tail%LENGTH * Plume_sources(N)%rate / State_Chm%SpcData(id_tracer)%Info%MW_g * Avo) &
                    /(Plume2d_tail%PDX * Plume2d_tail%PDY *Plume2d_tail%LENGTH*1.E6_fp ) ! molec/cm3
                    
         Plume2d_tail%CONCNT2d(n_x_mid,n_y_mid,n_species+N) = Plume2d_tail%CONCNT2d(n_x_mid,n_y_mid,n_species+N)  &
                    + (Plume2d_tail%LENGTH * Plume_sources(N)%rate / State_Chm%SpcData(id_tracer)%Info%MW_g * Avo) &
                    /(Plume2d_tail%PDX * Plume2d_tail%PDY *Plume2d_tail%LENGTH*1.E6_fp ) ! molec/cm3
      ENDDO
      NULLIFY(Plume2d_tail%next)
      Plume2d_head => Plume2d_tail
      Num_Plume2d = 1
    ELSE
      ! add injected species into eulerian grid
      WRITE(6,'(a)') ' lagrange module was turned off, add injected species to eulerian grid'
      DO N = 1, num_of_sources
        spc_name = Plume_sources(N)%species
        id_tracer   = Ind_(spc_name)
        Spc(id_tracer)%Conc(i_lon, i_lat, i_lev) = Spc(id_tracer)%Conc(i_lon, i_lat, i_lev)  &
                + (Length_init * Plume_sources(N)%rate / State_Chm%SpcData(id_tracer)%Info%MW_g * Avo) &
                    /(State_Met%AIRVOL(i_lon, i_lat, i_lev)*1.E6_fp ) ! molec/cm3
      ENDDO
    ENDIF
    Num_inject = 1
    tt = 0
    !  WRITE(6,'(a)') ' '
    !  WRITE(6,'(a)') '********************************************************'
    !  WRITE(6,'(a)') ' You are using lagrange_mod now'
    !  WRITE(6,'(a)') ' set variable use_lagrange = 0 to turn off lagrnage_mod '
    !  WRITE(6,'(a)') ' '
    !  WRITE(6,'(a)') ' WARNING: all the initial concentration is set as 0 now'
    !  WRITE(6,'(a)') '********************************************************'
    !  WRITE(6,'(a)') ' '

    !  id_PASV_LA  = State_Chm%nAdvect-2
    !  id_PASV_LA2 = State_Chm%nAdvect-1
    !  id_PASV_LA3 = State_Chm%nAdvect

    !  State_Chm%Species(id_PASV_LA)%Conc(:,:,:) = 0.0e+0_fp  ! [kg/kg]
    !  State_Chm%Species(id_PASV_LA2)%Conc(:,:,:) = 0.0e+0_fp  ! [kg/kg]
    !  State_Chm%Species(id_PASV_LA3)%Conc(:,:,:) = 0.0e+0_fp  ! [kg/kg]

	  ! debug purpose, BZ, does species database and SpecData contains the same species and order?
	  !WRITE(6,*) 'debug, BZ: Tracer id_PASV_LA  = ', id_PASV_LA, &
    !        '  Tracer Name: ', TRIM(State_Chm%SpcData(id_PASV_LA)%Info%Name)
			
	  !WRITE(6,*) 'debug, BZ: Tracer id_PASV_LA2  = ', id_PASV_LA2, &
    !        '  Tracer Name: ', TRIM(State_Chm%SpcData(id_PASV_LA2)%Info%Name)
			
	  !WRITE(6,*) 'debug, BZ: Tracer id_PASV_LA3  = ', id_PASV_LA3, &
    !        '  Tracer Name: ', TRIM(State_Chm%SpcData(id_PASV_LA3)%Info%Name)
    !ENDIF


    !deallocate(box_concnt_2D)
    !deallocate(box_concnt_1D)

    !nullify(Plume2d_new)
    !nullify(PLume2d)
    !nullify(Plume2d_prev)
    !nullify(Plume1d_new)
    !nullify(Plume1d)
    !nullify(Plume1d_prev)

90  FORMAT( /, A                 )
95  FORMAT( A                    )
100 FORMAT( A, L5                )
105 FORMAT( A, I0                )
110 FORMAT( A, A                 )
120 FORMAT(A8,1X,"|", A8,1X,"|",A8,1X,"|",1X,A8,1X,"|",1X,A8,1X,"|",1X,A10,1X,"|",1X,A)
121 FORMAT(I8,1X,"|",F8.3,1X,"|",F8.3,1X,"|",1X,F8.3,1X,"|",1X,F8.3,1X,"|",1X,F10.3,1X,"|",1X,A)
  END SUBROUTINE lagrange_init

!=================================================================

  SUBROUTINE plume_inject(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)

    USE Plume_list

    USE Input_Opt_Mod,   ONLY : OptInput, PlumeSource_t
    USE State_Met_Mod,   ONLY : MetState
    USE State_Chm_Mod,   ONLY : ChmState, Ind_
    USE State_Grid_Mod,  ONLY : GrdState
    USE Species_Mod,     ONLY : SpcConc
    USE TIME_MOD,        ONLY : GET_TS_DYN
    USE UnitConv_Mod !,    ONLY : Convert_Spc_Units, MOLECULES_SPECIES_PER_CM3


    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
    TYPE(MetState), intent(in)    :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(OptInput), intent(in)    :: Input_Opt
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure
    
    TYPE(SpcConc), POINTER :: Spc(:)

    !REAL(fp), POINTER :: X_edge(:), Y_edge(:)
    !REAL(fp)          :: X_edge2, Y_edge2


    INTEGER                :: i_box, i_lon, i_lat, i_lev 
    INTEGER                :: nAdv, i_advect1
    INTEGER                :: id_tracer, N
    INTEGER                :: OrigUnit
    REAL(fp)               :: box_lon, box_lat, box_lev
    REAL(fp), POINTER      :: PASV_EU
    REAL(fp)               :: MW_g
    REAL(fp)               :: Dt
    REAL(fp)               :: Entropy2_Concnt, Entropy2_V, Entropy2
    REAL(fp)               :: tracer2_mol, air2_mol, mix2_ratio
    REAL(fp)               :: tracer0_mol, air0_mol, mix0_ratio

    !CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: spc_name
    CHARACTER(LEN=255)     :: ErrMsg
    CHARACTER(LEN=255)     :: FileEntropy
    CHARACTER(LEN=255)     :: ThisLoc
    
    
    
    
    !X_edge => State_Grid%XEdge(:,1) 
    !Y_edge => State_Grid%YEdge(1,:) 
    !X_edge2       = X_edge(2)
    !Y_edge2       = Y_edge(2)
    ErrMsg                 =    ''
    Dt                     =    GET_TS_DYN()
    Spc                    =>   State_Chm%Species
    ThisLoc                =   ' -> at plume_inject (in module GeosCore/lagrange_mod.F90)'

    IF (.not.plume_inject_on) THEN
      !WRITE(6,'(a)') ' No plume injection, no lagrangian module configuration, skip '
      RETURN
    ENDIF
    IF (plume_inject_on .AND. num_of_sources.lt.1) THEN
      ErrMsg = 'Plume injection is turned on but no source specified'
      CALL GC_Error( ErrMsg, RC, ThisLoc )
      RETURN
    ENDIF

    ! -----------------------------------------------------------
    ! add new box every time step
    ! -----------------------------------------------------------
    ! N_stop_inject<0: keep injecting plume in the whole simulation;
    ! Num_inject<N_stop_inject: only inject N_stop_inject plumes at
    !                           beginning;
    IF(N_stop_inject<0.OR.Num_inject<N_stop_inject)THEN
      
      IF (.NOT. use_lagrange) THEN
        DO i_box = Num_inject+1, Num_inject+N_parcel, 1
          box_lon    = Plume_sources(1)%lon
          !box_lat    = ( Plume_sources(1)%lat1 + Length_lat * MOD(i_box,N_total) ) &
          !             * (-1.0)**FLOOR(1.0*i_box/N_total) ! -29.995S:29.995N:0.01
          box_lat    = GetInjectionLat(Plume_sources(1)%lat1,Plume_sources(1)%lat2,plume_interval,i_box)
          write(6,*) 'debug (BZ): injection lat: ', box_lat
          box_lev    = Plume_sources(1)%lev

          i_lon = Find_iLonLat(box_lon, DX, X_edge2)
          if(i_lon>IIPAR) i_lon=i_lon-IIPAR
          if(i_lon<1) i_lon=i_lon+IIPAR
          !WRITE(6,*) 'debug: ilon=', i_lon
          i_lat = Find_iLonLat(box_lat, DY, Y_edge2)
          if(i_lat>JJPAR) i_lat=JJPAR
          if(i_lat<1) i_lat=1
          !WRITE(6,*) 'debug: ilat=', i_lat
          i_lev = Find_iPLev(box_lev,P_edge)
          !WRITE(6,*) 'debug: curr_lev=', curr_lev
          !WRITE(6,*) 'debug: P_edge(1)=', curr_lev
          !WRITE(6,*) 'debug: ilev1=', i_lev
          if(i_lev>LLPAR) i_lev=LLPAR

          ! instantly add injected species into Eulerian grid 
          DO N = 1, num_of_sources
            spc_name = Plume_sources(N)%species
            id_tracer   = Ind_(spc_name)
            Spc(id_tracer)%Conc(i_lon, i_lat, i_lev) = Spc(id_tracer)%Conc(i_lon, i_lat, i_lev)  &
                    + (Length_init * Plume_sources(N)%rate / State_Chm%SpcData(id_tracer)%Info%MW_g * Avo) &
                    /(State_Met%AIRVOL(i_lon, i_lat, i_lev) * 1.0e+6_fp) ! molec/cm3
                    !(DX * DY * State_Met%BXHEIGHT(i_lon, i_lat, i_lev)*1.E6_fp ) ! molec/cm3
          ENDDO

          ! No plume seg created, do I need to calculate entropy?
          !-----------------------------------------------------------
          ! calculte the perfect entropy in Plume-in-Grid model without any
	        ! diffusion, all the injected tracer stay in the center 2-D grid 
          ! cell forever.
          ! State_Met%AIRDEN: [kg m-3]
          ! State_Met%AIRVOL: [m3]

          !tracer0_mol = Length_init*30.0/98.0 ! [mol] ( inject: 30.0 g/m )
          !air0_mol    = Dx_init*Dy_init*Length_init & ! [m]*[m]*[m]
          !*State_Met%AIRDEN(i_lon,i_lat,i_lev)*1000/AIRMW ! [mol]
          !mix0_ratio  = tracer0_mol/air0_mol
          !IF(mix0_ratio>0) Entropy0 = Entropy0 - BOLTZ*tracer0_mol*log(mix0_ratio)
          !-----------------------------------------------------------

          ! ====================================================================
          ! Add concentraion of PASV into conventional Eulerian GEOS-Chem in 
          ! corresponding with injected parcels in Lagrangian model
          ! For conventional GEOS-Chem for comparison with Lagrangian Model:
          !
          !  AD(I,J,L) = grid box dry air mass [kg]
          !  AIRMW     = dry air molecular wt [g/mol]
          !  MW_G(N)   = species molecular wt [g/mol]
          !     
          ! the conversion is:
          ! 
          !====================================================================
          !PASV_EU => State_Chm%Species(id_PASV_EU)%Conc(i_lon,i_lat,i_lev)  ! original unit is v/v, need to convert to [kg/kg]? BZ

          !MW_g = State_Chm%SpcData(nAdv)%Info%MW_g
          ! Here assume the injection rate is 30 kg/km for H2SO4: 
          !PASV_EU = PASV_EU + Length_init*1.0e-3_fp*30.0 &
          !                            /State_Met%AD(i_lon,i_lat,i_lev)

        ENDDO 
        Num_inject = Num_inject + N_parcel
      ELSE ! use_lagrange = 1
        ! State_Chm%Species(id_PASV_LA3)%Conc(:,:,:) = 0.0e+0_fp  ! [kg/kg]

        CALL lagrange_run(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
        CALL plume_run(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
        CALL lagrange_write_std( am_I_Root, RC )

      ENDIF ! IF(N_stop_inject<0.or.Num_inject<N_stop_inject)THEN
      !Num_inject = Num_inject + N_parcel
      tt = tt + 1
    ENDIF
      !======================================================================
      ! Convert species to [molec/cm3] (ewl, 8/16/16)
      !======================================================================
      !CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
      !                        new_units=MOLECULES_SPECIES_PER_CM3, RC = RC, previous_units=OrigUnit )
      !IF ( RC /= GC_SUCCESS ) THEN
      !   ErrMsg = 'Unit conversion error!'
      !   CALL GC_Error( ErrMsg, RC, 'plume_run in lagrange_mod.F90')
      !   RETURN
      !ENDIF


      !======================================================================
      ! run the fake 2nd order chemical reaction
      !======================================================================

      !!$OMP PARALLEL DO           &
      !!$OMP DEFAULT( SHARED     ) &
      !!$OMP PRIVATE( i_lev, i_lat, i_lon )
      !DO i_lev = 1, LLPAR
      !DO i_lat = 1, JJPAR
      !DO i_lon = 1, IIPAR

      ! from EU to EU2
      !State_Chm%Species(id_PASV_EU2)%Conc(i_lon,i_lat,i_lev) = &
      !        State_Chm%Species(id_PASV_EU2)%Conc(i_lon,i_lat,i_lev) &
      !      + Dt* Kchem*State_Chm%Species(id_PASV_EU)%Conc(i_lon,i_lat,i_lev)**2


      !  For all the grid cells in the troposphere, let concentration to be zero
      !IF(TROPP_sink == 1)THEN
      !IF ( State_Met%PEDGE(i_lon,i_lat,i_lev) > State_Met%TROPP(i_lon,i_lat) ) THEN
        ! State_Chm%Species(:)%Conc(i_lon,i_lat,i_lev) = 0.0
		!DO nspec = 1, State_Chm%nSpecies
	!		State_Chm%Species(nspec)%Conc(i_lon,i_lat,i_lev) = 0.0
		!END DO
     ! ENDIF
      !ENDIF


      !ENDDO
      !ENDDO
      !ENDDO
      !!$OMP END PARALLEL DO

    ! Calculate Entropy for Eulerian model
    !IF(mod(tt*NINT(Dt),24*60*60)==0)THEN   ! output once every day (24 hours)

    !  DO i_lon = 1, IIPAR
    !  DO i_lat = 1, JJPAR
    !  DO i_lev = 1, LLPAR
	! AIRDEN [kg m-3]
	! AIRVOL [m3]

    !    Entropy2_Concnt = State_Chm%Species(id_PASV_EU)%Conc(i_lon,i_lat,i_lev)
    !    Entropy2_V = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp

    !    tracer2_mol = Entropy2_Concnt*Entropy2_V/AVO
    !    air2_mol    = Entropy2_V/1e+6_fp*State_Met%AIRDEN(i_lon,i_lat,i_lev)*1000/AIRMW
    !    mix2_ratio  = tracer2_mol/air2_mol

    !    IF(mix2_ratio>0) Entropy2 = Entropy2 - BOLTZ*tracer2_mol*log(mix2_ratio)

    !  ENDDO
    !  ENDDO
    !  ENDDO



    !  FileEntropy   = 'Plume_entropy.txt'
    !  OPEN( 487,      FILE=TRIM( FileEntropy   ), STATUS='OLD',  &
    !       POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

    !  WRITE(*,*) "Entropy:", Entropy2, Entropy0
    !  WRITE(487,*) Entropy2, Entropy0

    !  CLOSE(487)


    !ENDIF


      !=======================================================================
      ! Convert species back to original units (ewl, 8/16/16)
      !=======================================================================
    !  CALL Convert_Spc_Units( Input_Opt, State_Chm,  State_Grid, State_Met, &
    !                          OrigUnit,  RC )

    !  IF ( RC /= GC_SUCCESS ) THEN
    !     ErrMsg = 'Unit conversion error!'
    !     CALL GC_Error( ErrMsg, RC, 'plume_mod.F90' )
    !     RETURN
    !  ENDIF

    !ENDIF ! IF(use_lagrange==0)THEN


    !=======================================================================
    !
    ! use plume model
    !
    !=======================================================================
    !IF(use_lagrange)THEN

    ! call the lagrnage_run() and plume_run() to calculate injected plume



    !ENDIF

    

  END SUBROUTINE plume_inject

!=================================================================

  SUBROUTINE lagrange_run(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)

    USE Plume_list

    USE Input_Opt_Mod,   ONLY : OptInput, PlumeSource_t
    USE State_Chm_Mod,   ONLY : ChmState, Ind_
    USE State_Met_Mod,   ONLY : MetState
    USE Species_Mod,     ONLY : SpcConc
    USE TIME_MOD,        ONLY : GET_TS_DYN
    USE State_Grid_Mod,  ONLY : GrdState
    USE UnitConv_Mod

!    USE GC_GRID_MOD,   ONLY : XEDGE, YEDGE
!    USE CMN_SIZE_Mod,  ONLY : DLAT, DLON !new

    LOGICAL, INTENT(IN)           :: am_I_Root
    TYPE(MetState), INTENT(IN)    :: State_Met
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State objectgg
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    TYPE(SpcConc), POINTER :: Spc(:)
    
           
    INTEGER                :: nAdv
    INTEGER                :: i_box
    INTEGER                :: i_lon, i_lat, i_lev, i_lon_1, i_lat_1, i_lev_1 
    INTEGER                :: next_i_lon, next_i_lat, next_i_lev  
    INTEGER                :: ii, jj, kk
    INTEGER                :: ki
    INTEGER                :: i_species, id_tracer
    INTEGER                :: i_x, i_y
    INTEGER                :: i_slab

    REAL(fp)               :: MW_g
    REAL(fp)               :: Dt, RK_Dt(5)                  ! = 600.0e+0_fp          
    REAL(fp), POINTER      :: PASV_EU    
    REAL(fp)               :: curr_lon, curr_lat, curr_pressure
    REAL(fp)               :: RK_lon, RK_lat
    !REAL(fp)               :: X_edge2, Y_edge2
    REAL(fp)               :: dbox_lon, dbox_lat, dbox_lev
    REAL(fp)               :: dbox_x_PS, dbox_y_PS
    REAL(fp)               :: box_x_PS, box_y_PS
    REAL(fp)               :: RK_Dlon(4), RK_Dlat(4), RK_Dlev(4)
    REAL(fp)               :: RK_x_PS, RK_y_PS
    REAL(fp)               :: RK_Dx_PS, RK_Dy_PS
    REAL(fp)               :: RK_u(4), RK_v(4), RK_omeg(4)

    REAL(fp), POINTER      :: u(:,:,:)
    REAL(fp), POINTER      :: v(:,:,:)
    REAL(fp), POINTER      :: omeg(:,:,:)
    REAL(fp), POINTER      :: Ptemp(:,:,:)
    REAL(fp), POINTER      :: T1(:,:,:)
    REAL(fp), POINTER      :: T2(:,:,:)

    real(fp) :: curr_T1, next_T2, ratio
    real(fp) :: V_prev, V_new

    real(fp) :: box_u, box_v, box_omeg
    real(fp) :: curr_u, curr_v, curr_omeg
    real(fp) :: curr_u_PS, curr_v_PS

    !real(fp), pointer :: X_edge(:), Y_edge(:)


    real(fp) :: length0
    real(fp) :: lon1, lon2, lat1, lat2


    real(fp) :: box_lon_edge, box_lat_edge

!    real(fp) :: box_alpha
    real(fp) :: D_wind, D_x, D_y
    real(fp) :: Ly ! Lyapunov exponent [s-1]


    real(fp) :: box_lon, box_lat, box_lev
    real(fp) :: box_length, box_alpha, box_theta
    real(fp) :: Pdx, Pdy
    real(fp) :: box_Ra, box_Rb
    real(fp) :: box_extra, box_life, box_label

    real(fp), dimension(:,:,:), allocatable :: box_concnt_2D
    real(fp), dimension(:,:), allocatable :: box_concnt_1D

    real(fp)  :: start, finish

    TYPE(Plume2d_list), POINTER :: Plume2d_new, Plume2d_curr, Plume2d_prev
    TYPE(Plume1d_list), POINTER :: Plume1d_new, Plume1d_curr, Plume1d_prev

    CHARACTER(LEN=255)            :: spc_name

    !CHARACTER(LEN=63)      :: OrigUnit
    INTEGER :: OrigUnit
    CHARACTER(LEN=255)     :: ErrMsg
    CHARACTER(LEN=255)     :: ThisLoc

    Dt = GET_TS_DYN()
    ThisLoc                =   ' -> at lagrange_run (in module GeosCore/lagrange_mod.F90)'
    RC     =  GC_SUCCESS
    ErrMsg = ''

    IF(Stop_inject==1) GOTO 400 ! deallocate and nullify -> exit

    ALLOCATE(box_concnt_2D(n_x_max, n_y_max, n_species + Num_of_sources)) ! add additional species here to track species that are newly added
    ALLOCATE(box_concnt_1D(n_slab_max, n_species + Num_of_sources))

    RK_Dt(1) = 0.0
    RK_Dt(2) = 0.5*Dt
    RK_Dt(3) = 0.5*Dt
    RK_Dt(4) = Dt
    RK_Dt(5) = 0.0

    u     => State_Met%U   ! m/s
    v     => State_Met%V   ! m/s
    omeg  => State_Met%OMEGA     ! Updraft velocity [Pa/s]
    Ptemp => State_Met%THETA     ! Potential temperature [K]
    T1    => State_Met%TMPU1     ! Temperature at start of timestep [K]
    T2    => State_Met%TMPU2     ! Temperature at end of timestep [K]

    !X_edge => State_Grid%XEdge(:,1) !XEDGE(:,1,1)   ! IIPAR+1 ! new
    !Y_edge => State_Grid%YEdge(1,:) !YEDGE(1,:,1)  
    ! Use second YEDGE, because sometimes YMID(2)-YMID(1) is not DLAT
    !X_edge2       = X_edge(2)
    !Y_edge2       = Y_edge(2)
     
    


    ! -----------------------------------------------------------
    ! add new box every time step, 
    ! -----------------------------------------------------------
    !IF(N_stop_inject<0.or.Num_inject<N_stop_inject)THEN
    DO i_box = Num_inject+1, Num_inject+N_parcel, 1 !BZ: consider putting this loop in subroutine: plume_inject

      ALLOCATE(Plume2d_new)

      Plume2d_new%IsNew = 1
      Plume2d_new%label = i_box

      Plume2d_new%LON = Plume_sources(1)%lon
      Plume2d_new%LAT = GetInjectionLat(Plume_sources(1)%lat1,Plume_sources(1)%lat2,plume_interval,i_box)
      write(6,*) 'debug (BZ): injection lat: ', box_lat
      Plume2d_new%LEV = Plume_sources(1)%lev
      !Plume2d_new%LAT = ( -30.005e+0_fp + Length_lat * MOD(i_box,N_total) ) &
      !               * (-1.0)**FLOOR(1.0*i_box/N_total) 
!      Plume2d_new%LAT = ( -30.005e+0_fp + 0.01e+0_fp * MOD(i_box,6000) ) &
!                     * (-1.0)**FLOOR(i_box/6000.0) ! -29.995S:29.995N:0.01

      i_lon = Find_iLonLat(Plume2d_new%LON, DX, X_edge2)
      if(i_lon>IIPAR) i_lon=i_lon-IIPAR
      if(i_lon<1) i_lon=i_lon+IIPAR
      i_lat = Find_iLonLat(Plume2d_new%LAT, DY, Y_edge2)
      if(i_lat>JJPAR) i_lat=JJPAR
      if(i_lat<1) i_lat=1
      i_lev = Find_iPLev(Plume2d_new%LEV,P_edge)
      if(i_lev>LLPAR) i_lev=LLPAR

      Plume2d_new%LENGTH = Length_init ! 1000m 
      Plume2d_new%ALPHA  = 0.0e+0_fp
      Plume2d_new%LIFE   = 0.0e+0_fp
      Plume2d_new%PDX    = Dx_init
      Plume2d_new%PDY    = Dy_init


      ALLOCATE(Plume2d_new%CONCNT2d(n_x_max, n_y_max, n_species+num_of_sources))

      ! pass background concentration
      DO i_species = 1, n_species
        Plume2d_new%CONCNT2d(:,:,i_species) = Spc(i_species)%Conc(i_lon, i_lat, i_lev)  ! molec/cm3
      ENDDO
      ! add tracer concentration
      DO i_species = 1, num_of_sources
        spc_name = Plume_sources(i_species)%species
        id_tracer   = Ind_(spc_name)
        Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,id_tracer) = Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,id_tracer)   &
                    + (Plume2d_new%LENGTH * Plume_sources(i_species)%rate / State_Chm%SpcData(id_tracer)%Info%MW_g * Avo) &
                    /(Plume2d_new%PDX * Plume2d_new%PDY *Plume2d_new%LENGTH*1.E6_fp ) ! molec/cm3
                    
         Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,n_species+i_species) = Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,n_species+i_species)  &
                    + (Plume2d_new%LENGTH * Plume_sources(i_species)%rate / State_Chm%SpcData(id_tracer)%Info%MW_g * Avo) &
                    /(Plume2d_new%PDX * Plume2d_new%PDY *Plume2d_new%LENGTH*1.E6_fp ) ! molec/cm3
      ENDDO
      NULLIFY(Plume2d_new%next)
      Plume2d_tail%next => Plume2d_new
      Plume2d_tail => Plume2d_tail%next
      
      ! track total mass injected

      
      !Plume2d_new%CONCNT2d = 0.0e+0_fp
      !Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,i_tracer) = Plume2d_new%LENGTH*30.0 &
      !             / (Plume2d_new%PDX*Plume2d_new%Pdy*Plume2d_new%LENGTH)
      ! From [g/m3] to [molec/cm3], 98.0 g/mol for H2SO4
      !Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,i_tracer) = &
      !          Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,i_tracer)/1.0e+6_fp/98.0*AVO


      !mass_la = mass_la + Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,i_tracer)* &
      !       (Plume2d_new%PDX*Plume2d_new%Pdy*Plume2d_new%LENGTH*1.0e+6_fp)



    ENDDO ! i_box = Num_inject+1, Num_inject+N_parcel, 1

    Num_Plume2d = Num_Plume2d + N_parcel
    Num_inject = Num_inject + N_parcel


    !ENDIF ! IF(N_stop_inject<0.or.Num_inject<N_stop_inject)THEN

    !=======================================================================
    ! for 2D plume: Run Lagrangian trajectory-track HERE
    !=======================================================================

    IF(.NOT.ASSOCIATED(Plume2d_head)) GOTO 401

    Plume2d_curr => Plume2d_head
    i_box = 0

    DO WHILE(ASSOCIATED(Plume2d_curr))

      i_box = i_box+1

      box_lon    = Plume2d_curr%LON
      box_lat    = Plume2d_curr%LAT
      box_lev    = Plume2d_curr%LEV
      box_length = Plume2d_curr%LENGTH
      box_alpha  = Plume2d_curr%ALPHA

      box_label  = Plume2d_curr%label
      box_life   = Plume2d_curr%LIFE

      Pdx    = Plume2d_curr%PDX
      Pdy    = Plume2d_curr%PDY

      box_concnt_2D = Plume2d_curr%CONCNT2d


      box_life = box_life + Dt

      !-----------------------------------------------------------------------
      ! For the center of the plume
      !-----------------------------------------------------------------------

      ! make sure the location is not out of range
      do while (box_lat > Y_edge(JJPAR+1))
         box_lat = Y_edge(JJPAR+1) - ( box_lat-Y_edge(JJPAR+1) )
      end do
      do while (box_lat < Y_edge(1))
         box_lat = Y_edge(1) + ( box_lat-Y_edge(1) )
      end do
      do while (box_lon > X_edge(IIPAR+1))
         box_lon = box_lon - 360.0
      end do
      do while (box_lon < X_edge(1))
         box_lon = box_lon + 360.0
      end do

      i_lon = Find_iLonLat(box_lon, DX, X_edge2)
      if(i_lon>IIPAR) i_lon=i_lon-IIPAR
      if(i_lon<1) i_lon=i_lon+IIPAR
      i_lat = Find_iLonLat(box_lat, DY, Y_edge2)
      if(i_lat>JJPAR) i_lat=JJPAR
      if(i_lat<1) i_lat=1
      i_lev = Find_iPLev(box_lev,P_edge)
      if(i_lev>LLPAR) i_lev=LLPAR

      curr_lon      = box_lon
      curr_lat      = box_lat
      curr_pressure = box_lev

      DO Ki = 1,4,1

       !------------------------------------------------------------------
       ! For vertical wind speed:
       ! pay attention for the polar region * * *
       !------------------------------------------------------------------
       if(abs(curr_lat)>Y_mid(JJPAR))then
          curr_omeg = Interplt_wind_RLL_polar(omeg, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       else
          curr_omeg = Interplt_wind_RLL(omeg, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       endif

       RK_omeg(Ki) = curr_omeg
       RK_Dlev(Ki)   = Dt * curr_omeg / 100.0     ! Pa => hPa

       curr_pressure = box_lev + RK_Dt(Ki+1) * curr_omeg / 100.0

       if(curr_pressure<P_mid(LLPAR)) &
             curr_pressure = P_mid(LLPAR) !+ ( P_mid(LLPAR) - curr_pressure )
       if(curr_pressure>P_mid(1)) &
             curr_pressure = P_mid(1) !- ( curr_pressure - P_mid(1) )


       !------------------------------------------------------------------
       ! For the region where lat<72, use Regualr Longitude-Latitude Mesh:
       !------------------------------------------------------------------
       if(abs(curr_lat)<=72.0)then

         curr_u = Interplt_wind_RLL(u, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
         curr_v = Interplt_wind_RLL(v, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

         RK_u(Ki) = curr_u
         RK_v(Ki) = curr_v

         dbox_lon  = (RK_Dt(Ki+1)*curr_u) &
                        / (2.0*PI*Re*COS(box_lat*PI/180.0)) * 360.0
         dbox_lat  = (RK_Dt(Ki+1)*curr_v) / (PI*Re) * 180.0

         curr_lon  = box_lon + dbox_lon
         curr_lat  = box_lat + dbox_lat
 

         RK_Dlon(Ki) = (Dt*curr_u) &
                      / (2.0*PI*Re*COS(box_lat*PI/180.0)) * 360.0
         RK_Dlat(Ki) = (Dt*curr_v) / (PI*Re) * 180.0

       endif

       !------------------------------------------------------------------
       ! For the polar region (lat>=72), use polar sterographic
       !------------------------------------------------------------------
       if(abs(curr_lat)>72.0)then

         if(abs(curr_lat)>Y_mid(JJPAR))then 
         curr_u_PS = Interplt_uv_PS_polar(1, u, v, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
         curr_v_PS = Interplt_uv_PS_polar(0, u, v, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
         else
         curr_u_PS = Interplt_uv_PS(1, u, v, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)    
         curr_v_PS = Interplt_uv_PS(0, u, v, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)    
         endif

         RK_u(Ki) = curr_u_PS
         RK_v(Ki) = curr_v_PS

         dbox_x_PS = RK_Dt(Ki+1)*curr_u_PS
         dbox_y_PS = RK_Dt(Ki+1)*curr_v_PS


         RK_Dx_PS = Dt*curr_u_PS
         RK_Dy_PS = Dt*curr_v_PS


         !------------------------------------------------------------------
         ! change from (lon,lat) in RLL to (x,y) in PS: 
         !------------------------------------------------------------------
         if(box_lat<0)then
           box_x_PS = -1.0* Re* COS(box_lon*PI/180.0) &
                                / TAN(box_lat*PI/180.0)
           box_y_PS = -1.0* Re* SIN(box_lon*PI/180.0) &
                                / TAN(box_lat*PI/180.0)
         else
           box_x_PS = Re* COS(box_lon*PI/180.0) &
                        / TAN(box_lat*PI/180.0)
           box_y_PS = Re* SIN(box_lon*PI/180.0) &
                        / TAN(box_lat*PI/180.0)
         endif

         RK_x_PS  = box_x_PS + RK_Dx_PS
         RK_y_PS  = box_y_PS + RK_Dy_PS

         box_x_PS  = box_x_PS + dbox_x_PS
         box_y_PS  = box_y_PS + dbox_y_PS


         !------------------------------------------------------------------
         ! change from (x,y) in PS to (lon,lat) in RLL
         !------------------------------------------------------------------
         if(box_x_PS>0.0)then
           curr_lon = ATAN( box_y_PS / box_x_PS )*180.0/PI 
         endif
         if(box_x_PS<0.0 .and. box_y_PS<=0.0)then
           curr_lon = ATAN( box_y_PS / box_x_PS )*180.0/PI -180.0
         endif
         if(box_x_PS<0.0 .and. box_y_PS>0.0)then
           curr_lon = ATAN( box_y_PS / box_x_PS )*180.0/PI +180.0
         endif
           
         if(curr_lat<0.0)then
           curr_lat= -1* ATAN( Re/SQRT(box_x_PS**2+box_y_PS**2) ) *180.0/PI
         else
           curr_lat= ATAN( Re / SQRT(box_x_PS**2+box_y_PS**2) ) *180.0/PI
         endif

         !------------------------------------------------------------------
         ! For 4th order Runge Kutta
         !------------------------------------------------------------------
         if(RK_x_PS>0.0)then
           RK_lon = ATAN( RK_y_PS / RK_x_PS )*180.0/PI
         endif
         if(RK_x_PS<0.0 .and. RK_y_PS<=0.0)then
           RK_lon = ATAN( RK_y_PS / RK_x_PS )*180.0/PI -180.0
         endif
         if(RK_x_PS<0.0 .and. RK_y_PS>0.0)then
           RK_lon = ATAN( RK_y_PS / RK_x_PS )*180.0/PI +180.0
         endif

         if(box_lat<0.0)then
           RK_lat = -1 * ATAN( Re / SQRT(RK_x_PS**2+RK_y_PS**2) ) *180.0/PI
         else
           RK_lat = ATAN( Re / SQRT(RK_x_PS**2+RK_y_PS**2) ) *180.0/PI
         endif

         RK_Dlon(Ki) = RK_lon - box_lon
         RK_Dlat(Ki) = RK_lat - box_lat

       endif ! if(abs(curr_lat)>72.0)then

      ENDDO ! Ki = 1,4,1

      box_lon = box_lon + &
                (RK_Dlon(1)+2.0*RK_Dlon(2)+2.0*RK_Dlon(3)+RK_Dlon(4))/6.0
      box_lat = box_lat + &
                (RK_Dlat(1)+2.0*RK_Dlat(2)+2.0*RK_Dlat(3)+RK_Dlat(4))/6.0
      box_lev = box_lev + &
                (RK_Dlev(1)+2.0*RK_Dlev(2)+2.0*RK_Dlev(3)+RK_Dlev(4))/6.0


      box_u    = ( RK_u(1) + 2.0*RK_u(2) + 2.0*RK_u(3) + RK_u(4) ) / 6.0
      box_v    = ( RK_v(1) + 2.0*RK_v(2) + 2.0*RK_v(3) + RK_v(4) ) / 6.0
      box_omeg = ( RK_omeg(1) + 2.0*RK_omeg(2) + 2.0*RK_omeg(3) + RK_omeg(4) ) / 6.0


      !--------------------------------------------------------------------
      ! interpolate temperature for plume volumn change (PV=nRT):
      !--------------------------------------------------------------------
      if(abs(curr_lat)>Y_mid(JJPAR))then
         curr_T1 = Interplt_wind_RLL_polar(T1, i_lon, i_lat, &
                                 i_lev, curr_lon, curr_lat, curr_pressure)
      else
         curr_T1 = Interplt_wind_RLL(T1, i_lon, i_lat, i_lev, &
                                        curr_lon, curr_lat, curr_pressure)
      endif

      i_lon = Find_iLonLat(curr_lon, DX, X_edge2)
      if(i_lon>IIPAR) i_lon=i_lon-IIPAR
      if(i_lon<1) i_lon=i_lon+IIPAR

      i_lat = Find_iLonLat(curr_lat, DY, Y_edge2)
      if(i_lat>JJPAR) i_lat=JJPAR
      if(i_lat<1) i_lat=1

      i_lev = Find_iPLev(box_lev,P_edge)
      if(i_lev>LLPAR) i_lev=LLPAR


      if(abs(box_lat)>Y_mid(JJPAR))then
         next_T2 = Interplt_wind_RLL_polar(T2, i_lon_1, i_lat_1, i_lev_1, box_lon, box_lat, box_lev)
      else
         next_T2 = Interplt_wind_RLL(T2, i_lon_1, i_lat_1, i_lev_1, box_lon, box_lat, box_lev)
      endif

      ! PV=nRT, V2 = T2/P2 : T1/P1 * V1
      ratio = ( (next_T2/box_lev)/(curr_T1/curr_pressure) )**(1/2)

      ! assume the volume change mainly apply to the cross-section, 
      ! the box_length would not change 

      V_prev = Pdx*Pdy*box_length*1.0e+6_fp

      Pdx = Pdx *ratio
      Pdy = Pdy *ratio

      V_new = Pdx*Pdy*box_length*1.0e+6_fp


!      call cpu_time(start)
!      box_concnt_2D(:,:,:) = box_concnt_2D(:,:,:)*V_prev/V_new
!      call cpu_time(finish)
!      WRITE(6,*)'Time1 (finish-start) for 2D:', i_box, finish-start


      DO i_species = 1, n_species + num_of_sources, 1
      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( i_y, i_x )
      DO i_y = 1,n_y_max,1
      DO i_x = 1,n_x_max,1
         box_concnt_2D(i_x,i_y,i_species) = &
                        box_concnt_2D(i_x,i_y,i_species)*V_prev/V_new
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO
      ENDDO


      !------------------------------------------------------------------
      ! calcualte the box_alpha [0,2*PI)
      ! angle between plume length and eastward (from west to east)
      !------------------------------------------------------------------
      if((box_u**2+box_v**2)==0)then
          box_alpha = 0.0
      else
        IF(box_v>=0)THEN
          box_alpha = ACOS( box_u/SQRT(box_u**2+box_v**2) ) 
        ELSE
          box_alpha = 2*PI - ACOS( box_u/SQRT(box_u**2+box_v**2) )
        ENDIF
      endif

      !------------------------------------------------------------------
      ! calcualte the Lyaponov exponent (Ly), unit: s-1
      !------------------------------------------------------------------
      Ly = Calc_Ly(u, v, i_lon, i_lat, i_lev, box_alpha, box_lon, box_lat)


      !------------------------------------------------------------------
      ! Horizontal stretch:
      ! Adjust the length/radius of box based on Lyaponov exponent (Ly)
      !------------------------------------------------------------------

      length0           = box_length
      box_length = EXP(Ly*Dt) * length0


      Pdx = Pdx*SQRT(length0/box_length)
      Pdy = Pdy*SQRT(length0/box_length)

      ! ----------------------------------------------------------------
      ! Entrainment/detrainment should be considered here
      ! ----------------------------------------------------------------

      ! ----------------------------------------------------------------
      ! update     
      ! ----------------------------------------------------------------
      Plume2d_curr%LON = box_lon
      Plume2d_curr%LAT = box_lat
      Plume2d_curr%LEV = box_lev
      Plume2d_curr%LENGTH = box_length
      Plume2d_curr%ALPHA  = box_alpha
      Plume2d_curr%label   = box_label
      Plume2d_curr%LIFE   = box_life


      Plume2d_curr%PDX = Pdx
      Plume2d_curr%PDY = Pdy

      Plume2d_curr%CONCNT2d    = box_concnt_2D

      Plume2d_curr => Plume2d_curr%next

    ENDDO  ! DO WHILE(ASSOCIATED(Plume2d))

401 CONTINUE


    !=======================================================================
    ! For 1d plume: Run Lagrangian trajectory-track HERE
    !=======================================================================

    IF(.NOT.ASSOCIATED(Plume1d_head)) GOTO 400


    Plume1d_curr => Plume1d_head
    i_box = 0

    DO WHILE(ASSOCIATED(Plume1d_curr))

      i_box = i_box+1


      box_lon    = Plume1d_curr%LON
      box_lat    = Plume1d_curr%LAT
      box_lev    = Plume1d_curr%LEV
      box_length = Plume1d_curr%LENGTH
      box_alpha  = Plume1d_curr%ALPHA
      box_label  = Plume1d_curr%label
      box_life   = Plume1d_curr%LIFE

      box_Ra    = Plume1d_curr%RA
      box_Rb    = Plume1d_curr%RB


      box_concnt_1D = Plume1d_curr%CONCNT1d


      box_life = box_life + Dt

      !-----------------------------------------------------------------------
      ! For the center of the plume
      !-----------------------------------------------------------------------

      ! make sure the location is not out of range
      do while (box_lat > Y_edge(JJPAR+1))
         box_lat = Y_edge(JJPAR+1) - ( box_lat-Y_edge(JJPAR+1) )
      end do

      do while (box_lat < Y_edge(1))
         box_lat = Y_edge(1) + ( box_lat-Y_edge(1) )
      end do

      do while (box_lon > X_edge(IIPAR+1))
         box_lon = box_lon - 360.0
      end do

      do while (box_lon < X_edge(1))
         box_lon = box_lon + 360.0
      end do


      curr_lon      = box_lon
      curr_lat      = box_lat
      curr_pressure = box_lev      ! hPa


      i_lon = Find_iLonLat(curr_lon, DX, X_edge2)
      if(i_lon>IIPAR) i_lon=i_lon-IIPAR
      if(i_lon<1) i_lon=i_lon+IIPAR

      i_lat = Find_iLonLat(curr_lat, DY, Y_edge2)
      if(i_lat>JJPAR) i_lat=JJPAR
      if(i_lat<1) i_lat=1

      i_lev = Find_iPLev(curr_pressure,P_edge)
      if(i_lev>LLPAR) i_lev=LLPAR

      DO Ki = 1,4,1

       !------------------------------------------------------------------
       ! For vertical wind speed:
       ! pay attention for the polar region * * *
       !------------------------------------------------------------------
       if(abs(curr_lat)>Y_mid(JJPAR))then
          curr_omeg = Interplt_wind_RLL_polar(omeg, i_lon, i_lat, i_lev, &
                                        curr_lon, curr_lat, curr_pressure)
       else
          curr_omeg = Interplt_wind_RLL(omeg, i_lon, i_lat, i_lev, &
                                        curr_lon, curr_lat, curr_pressure)
       endif

       RK_omeg(Ki) = curr_omeg

       RK_Dlev(Ki)   = Dt * curr_omeg / 100.0     ! Pa => hPa
       curr_pressure = box_lev + RK_Dt(Ki+1) * curr_omeg / 100.0

       if(curr_pressure<P_mid(LLPAR)) &
             curr_pressure = P_mid(LLPAR) !+ ( P_mid(LLPAR) - curr_pressure )
       if(curr_pressure>P_mid(1)) &
             curr_pressure = P_mid(1) !- ( curr_pressure - P_mid(1) )


       !------------------------------------------------------------------
       ! For the region where lat<72, use Regualr Longitude-Latitude Mesh:
       !------------------------------------------------------------------
       if(abs(curr_lat)<=72.0)then

         curr_u = Interplt_wind_RLL(u, i_lon, i_lat, i_lev, curr_lon, &
                                                  curr_lat, curr_pressure)
         curr_v = Interplt_wind_RLL(v, i_lon, i_lat, i_lev, curr_lon, &
                                                  curr_lat, curr_pressure)

         RK_u(Ki) = curr_u
         RK_v(Ki) = curr_v

         dbox_lon  = (RK_Dt(Ki+1)*curr_u) &
                        / (2.0*PI*Re*COS(box_lat*PI/180.0)) * 360.0
         dbox_lat  = (RK_Dt(Ki+1)*curr_v) / (PI*Re) * 180.0

         curr_lon  = box_lon + dbox_lon
         curr_lat  = box_lat + dbox_lat


         RK_Dlon(Ki) = (Dt*curr_u) &
                      / (2.0*PI*Re*COS(box_lat*PI/180.0)) * 360.0
         RK_Dlat(Ki) = (Dt*curr_v) / (PI*Re) * 180.0

       endif


       !------------------------------------------------------------------
       ! For the polar region (lat>=72), use polar sterographic
       !------------------------------------------------------------------
       if(abs(curr_lat)>72.0)then

         if(abs(curr_lat)>Y_mid(JJPAR))then
         curr_u_PS = Interplt_uv_PS_polar(1, u, v, i_lon, i_lat, i_lev, &
                                        curr_lon, curr_lat, curr_pressure)
         curr_v_PS = Interplt_uv_PS_polar(0, u, v, i_lon, i_lat, i_lev, &
                                        curr_lon, curr_lat, curr_pressure)
         else
         curr_u_PS = Interplt_uv_PS(1, u, v, i_lon, i_lat, i_lev, &
                                        curr_lon, curr_lat, curr_pressure) 
         curr_v_PS = Interplt_uv_PS(0, u, v, i_lon, i_lat, i_lev, &
                                        curr_lon, curr_lat, curr_pressure) 
         endif

         RK_u(Ki) = curr_u_PS
         RK_v(Ki) = curr_v_PS

         dbox_x_PS = RK_Dt(Ki+1)*curr_u_PS
         dbox_y_PS = RK_Dt(Ki+1)*curr_v_PS


         RK_Dx_PS = Dt*curr_u_PS
         RK_Dy_PS = Dt*curr_v_PS

         !------------------------------------------------------------------
         ! change from (lon,lat) in RLL to (x,y) in PS: 
         !------------------------------------------------------------------
         if(box_lat<0)then
           box_x_PS = -1.0* Re* COS(box_lon*PI/180.0) &
                                / TAN(box_lat*PI/180.0)
           box_y_PS = -1.0* Re* SIN(box_lon*PI/180.0) &
                                / TAN(box_lat*PI/180.0)
         else
           box_x_PS = Re* COS(box_lon*PI/180.0) &
                        / TAN(box_lat*PI/180.0)
           box_y_PS = Re* SIN(box_lon*PI/180.0) &
                        / TAN(box_lat*PI/180.0)
         endif

         RK_x_PS  = box_x_PS + RK_Dx_PS
         RK_y_PS  = box_y_PS + RK_Dy_PS

         box_x_PS  = box_x_PS + dbox_x_PS
         box_y_PS  = box_y_PS + dbox_y_PS


         !------------------------------------------------------------------
         ! change from (x,y) in PS to (lon,lat) in RLL
         !------------------------------------------------------------------
         if(box_x_PS>0.0)then
           curr_lon = ATAN( box_y_PS / box_x_PS )*180.0/PI
         endif
         if(box_x_PS<0.0 .and. box_y_PS<=0.0)then
           curr_lon = ATAN( box_y_PS / box_x_PS )*180.0/PI -180.0
         endif
         if(box_x_PS<0.0 .and. box_y_PS>0.0)then
           curr_lon = ATAN( box_y_PS / box_x_PS )*180.0/PI +180.0
         endif

         if(curr_lat<0.0)then
           curr_lat= -1* ATAN( Re/SQRT(box_x_PS**2+box_y_PS**2) ) *180.0/PI
         else
           curr_lat= ATAN( Re / SQRT(box_x_PS**2+box_y_PS**2) ) *180.0/PI
         endif



         !------------------------------------------------------------------
         ! For 4th order Runge Kutta
         !------------------------------------------------------------------
         if(RK_x_PS>0.0)then
           RK_lon = ATAN( RK_y_PS / RK_x_PS )*180.0/PI
         endif
         if(RK_x_PS<0.0 .and. RK_y_PS<=0.0)then
           RK_lon = ATAN( RK_y_PS / RK_x_PS )*180.0/PI -180.0
         endif
         if(RK_x_PS<0.0 .and. RK_y_PS>0.0)then
           RK_lon = ATAN( RK_y_PS / RK_x_PS )*180.0/PI +180.0
         endif

         if(box_lat<0.0)then
           RK_lat = -1 * ATAN( Re / SQRT(RK_x_PS**2+RK_y_PS**2) ) *180.0/PI
         else
           RK_lat = ATAN( Re / SQRT(RK_x_PS**2+RK_y_PS**2) ) *180.0/PI
         endif

         RK_Dlon(Ki) = RK_lon - box_lon
         RK_Dlat(Ki) = RK_lat - box_lat

       endif ! if(abs(curr_lat)>72.0)then

      ENDDO ! Ki = 1,4,1

      box_lon = box_lon + &
                (RK_Dlon(1)+2.0*RK_Dlon(2)+2.0*RK_Dlon(3)+RK_Dlon(4))/6.0
      box_lat = box_lat + &
                (RK_Dlat(1)+2.0*RK_Dlat(2)+2.0*RK_Dlat(3)+RK_Dlat(4))/6.0
      box_lev = box_lev + &
                (RK_Dlev(1)+2.0*RK_Dlev(2)+2.0*RK_Dlev(3)+RK_Dlev(4))/6.0



      box_u    = ( RK_u(1) + 2.0*RK_u(2) &
                         + 2.0*RK_u(3) + RK_u(4) ) / 6.0
      box_v    = ( RK_v(1) + 2.0*RK_v(2) &
                         + 2.0*RK_v(3) + RK_v(4) ) / 6.0
      box_omeg = ( RK_omeg(1) + 2.0*RK_omeg(2) &
                         + 2.0*RK_omeg(3) + RK_omeg(4) ) / 6.0



      !--------------------------------------------------------------------
      ! interpolate temperature for plume volumn change 
      ! (PV=nRT):
      !--------------------------------------------------------------------
      if(abs(curr_lat)>Y_mid(JJPAR))then
         curr_T1 = Interplt_wind_RLL_polar(T1, i_lon, i_lat, &
                                 i_lev, curr_lon, curr_lat, curr_pressure)
      else
         curr_T1 = Interplt_wind_RLL(T1, i_lon, i_lat, i_lev, &
                                        curr_lon, curr_lat, curr_pressure)
      endif



      i_lon = Find_iLonLat(curr_lon, DX, X_edge2)
      if(i_lon>IIPAR) i_lon=i_lon-IIPAR
      if(i_lon<1) i_lon=i_lon+IIPAR

      i_lat = Find_iLonLat(curr_lat, DY, Y_edge2)
      if(i_lat>JJPAR) i_lat=JJPAR
      if(i_lat<1) i_lat=1

      i_lev = Find_iPLev(box_lev,P_edge)
      if(i_lev>LLPAR) i_lev=LLPAR


      if(abs(box_lat)>Y_mid(JJPAR))then
         next_T2 = Interplt_wind_RLL_polar(T2, i_lon, i_lat, i_lev, &
                           box_lon, box_lat, box_lev)
      else
         next_T2 = Interplt_wind_RLL(T2, i_lon, i_lat, i_lev, &
                           box_lon, box_lat, box_lev)
      endif


      ! PV=nRT, V2 = T2/P2 : T1/P1 * V1
      ratio = ( (next_T2/box_lev)/(curr_T1/curr_pressure) )**(1/2)


      ! assume the volume change mainly apply to the cross-section, 
      ! the box_length would not change 

      V_prev = box_Ra*box_Rb*box_length*1.0e+6_fp

      box_Ra = box_Ra *ratio
      box_Rb = box_Rb *ratio

      V_new = box_Ra*box_Rb*box_length*1.0e+6_fp

      box_concnt_1D(:,:) = box_concnt_1D(:,:)*V_prev/V_new



      !------------------------------------------------------------------
      ! calcualte the box_alpha [0,2*PI)
      !------------------------------------------------------------------
      if((box_u**2+box_v**2)==0)then
          box_alpha = 0.0
      else
        IF(box_v>=0)THEN
          box_alpha = ACOS( box_u/SQRT(box_u**2+box_v**2) )
        ELSE
          box_alpha = 2*PI - ACOS( box_u/SQRT(box_u**2+box_v**2) )
        ENDIF
      endif

      !------------------------------------------------------------------
      ! calcualte the Lyaponov exponent (Ly), unit: s-1
      !------------------------------------------------------------------
      Ly = Calc_Ly(u, v, i_lon, i_lat, i_lev, box_alpha, box_lon, box_lat)



      !------------------------------------------------------------------
      ! Horizontal stretch:
      ! Adjust the length/radius of box based on Lyaponov exponent (Ly)
      !------------------------------------------------------------------
      length0           = box_length
      box_length = EXP(Ly*Dt) * length0


      box_Ra = box_Ra*SQRT(length0/box_length)
      box_Rb = box_Rb*SQRT(length0/box_length)



      Plume1d_curr%LON    = box_lon
      Plume1d_curr%LAT    = box_lat
      Plume1d_curr%LEV    = box_lev
      Plume1d_curr%LENGTH = box_length
      Plume1d_curr%ALPHA  = box_alpha
      Plume1d_curr%label  = box_label
      Plume1d_curr%LIFE   = box_life
      Plume1d_curr%RA     = box_Ra
      Plume1d_curr%RB     = box_Rb

      Plume1d_curr%CONCNT1d = box_concnt_1D


      Plume1d_curr => Plume1d_curr%next

    ENDDO  ! DO WHILE(ASSOCIATED(Plume1d))



400 CONTINUE


    !------------------------------------------------------------------
    ! Everything is done, clean up pointers
    !------------------------------------------------------------------

    IF(allocated(box_concnt_2D)) deallocate(box_concnt_2D)
    IF(allocated(box_concnt_1D)) deallocate(box_concnt_1D)

    IF(ASSOCIATED(u)) nullify(u)
    IF(ASSOCIATED(v)) nullify(v)
    IF(ASSOCIATED(omeg)) nullify(omeg)

    IF(ASSOCIATED(Ptemp)) nullify(Ptemp)
    IF(ASSOCIATED(T1)) nullify(T1)
    IF(ASSOCIATED(T2)) nullify(T2)

    IF(ASSOCIATED(X_edge)) nullify(X_edge)
    IF(ASSOCIATED(Y_edge)) nullify(Y_edge)

    IF(ASSOCIATED(PASV_EU)) nullify(PASV_EU)

    IF(ASSOCIATED(Plume2d_new)) nullify(Plume2d_new)
    IF(ASSOCIATED(Plume2d_curr)) nullify(Plume2d_curr)
    IF(ASSOCIATED(Plume2d_prev)) nullify(Plume2d_prev)
    IF(ASSOCIATED(Plume1d_new)) nullify(Plume1d_new)
    IF(ASSOCIATED(Plume1d_curr)) nullify(Plume1d_curr)
    IF(ASSOCIATED(Plume1d_prev)) nullify(Plume1d_prev)

  END SUBROUTINE lagrange_run

!------------------------------------------------------------------
! Get lattitude of injection at specific time step
! This make sure plume seg created from lat1:dlat:lat2 -> lat2:-dlat,lat1 -> ...
!------------------------------------------------------------------
Real(fp) FUNCTION GetInjectionLat(lat1, lat2, dlat, global_id) 
    implicit none
    REAL(fp), INTENT(IN) :: lat1, lat2, dlat
    INTEGER,  INTENT(IN) :: global_id
    INTEGER              :: n_total, cycle_length
    INTEGER              :: pos, idx

    n_total = INT((lat2 - lat1)/dlat) + 1
    cycle_length = 2*n_total - 2

    pos = MOD(global_id-1, cycle_length)

    IF (pos < n_total) THEN
        idx = pos
    ELSE
        idx = cycle_length - pos
    END IF

    GetInjectionLat = lat1 + idx*dlat

END FUNCTION GetInjectionLat
!------------------------------------------------------------------
! calcualte the Lyaponov exponent (Ly), unit: s-1
!------------------------------------------------------------------
  REAL(fp) FUNCTION Calc_Ly(u, v, i_lon, i_lat, i_lev, box_alpha, box_lon, box_lat)

    implicit none

    real(fp)    :: box_alpha, box_lon, box_lat
    real(fp), pointer :: u(:,:,:), v(:,:,:)

    real(fp)    :: D_wind, DX_m, DY_m, Ly

    integer     :: i_lon, i_lat, i_lev
    integer     :: next_i_lon, next_i_lat


      IF(box_alpha>=1.75*PI)THEN
        if(box_lon>=X_mid(i_lon))then
          next_i_lon = i_lon + 1
          DX_m   = DX/360.0 * 2*PI *Re*COS(box_lat/180*PI)
        else
          next_i_lon = i_lon - 1
          DX_m   = -1*DX/360.0 * 2*PI *Re*COS(box_lat/180*PI)
        endif

        next_i_lat = i_lat

        if(next_i_lon>IIPAR) next_i_lon=next_i_lon-IIPAR
        if(next_i_lon<1)     next_i_lon=next_i_lon+IIPAR
        if(next_i_lat>JJPAR) next_i_lat=JJPAR
        if(next_i_lat<1)     next_i_lat=1

        D_wind = u(next_i_lon, next_i_lat, i_lev)-u(i_lon, i_lat, i_lev)
        Ly     = D_wind/DX_m

      ELSEIF(box_alpha>=1.25*PI)THEN
        if(box_lat>=Y_mid(i_lat))then
          next_i_lat = i_lat + 1
          DY_m   = DY/360.0 * 2*PI*Re
        else
          next_i_lat = i_lat - 1
          DY_m   = -1*DY/360.0 * 2*PI*Re
        endif

        next_i_lon = i_lon

        if(next_i_lon>IIPAR) next_i_lon=next_i_lon-IIPAR
        if(next_i_lon<1)     next_i_lon=next_i_lon+IIPAR
        if(next_i_lat>JJPAR) next_i_lat=JJPAR
        if(next_i_lat<1)     next_i_lat=1

        D_wind = v(next_i_lon, next_i_lat, i_lev)-v(i_lon, i_lat, i_lev)
        Ly     = D_wind/DY_m

      ELSEIF(box_alpha>=0.75*PI)THEN
        if(box_lon>=X_mid(i_lon))then
          next_i_lon = i_lon + 1
          DX_m   = DX/360.0 * 2*PI *Re*COS(box_lat/180*PI)
        else
          next_i_lon = i_lon - 1
          DX_m   = -1*DX/360.0 * 2*PI *Re*COS(box_lat/180*PI)
        endif

        next_i_lat = i_lat

        if(next_i_lon>IIPAR) next_i_lon=next_i_lon-IIPAR
        if(next_i_lon<1)     next_i_lon=next_i_lon+IIPAR
        if(next_i_lat>JJPAR) next_i_lat=JJPAR
        if(next_i_lat<1)     next_i_lat=1

        D_wind = u(next_i_lon, next_i_lat, i_lev)-u(i_lon, i_lat, i_lev)
        Ly     = D_wind/DX_m

      ELSEIF(box_alpha>=0.25*PI)THEN
        if(box_lat>=Y_mid(i_lat))then
          next_i_lat = i_lat + 1
          DY_m   = DY/360.0 * 2*PI*Re
        else
          next_i_lat = i_lat - 1
          DY_m   = -1*DY/360.0 * 2*PI*Re
        endif

        next_i_lon = i_lon

        if(next_i_lon>IIPAR) next_i_lon=next_i_lon-IIPAR
        if(next_i_lon<1)     next_i_lon=next_i_lon+IIPAR
        if(next_i_lat>JJPAR) next_i_lat=JJPAR
        if(next_i_lat<1)     next_i_lat=1

        D_wind = v(next_i_lon, next_i_lat, i_lev)-v(i_lon, i_lat, i_lev)
        Ly     = D_wind/DY_m

      ELSE
        if(box_lon>=X_mid(i_lon))then
          next_i_lon = i_lon + 1
          DX_m   = DX/360.0 * 2*PI *Re*COS(box_lat/180*PI)
        else
          next_i_lon = i_lon - 1
          DX_m   = -1*DX/360.0 * 2*PI *Re*COS(box_lat/180*PI)
        endif

        next_i_lat = i_lat

        if(next_i_lon>IIPAR) next_i_lon=next_i_lon-IIPAR
        if(next_i_lon<1)     next_i_lon=next_i_lon+IIPAR
        if(next_i_lat>JJPAR) next_i_lat=JJPAR
        if(next_i_lat<1)     next_i_lat=1

        D_wind = u(next_i_lon, next_i_lat, i_lev)-u(i_lon, i_lat, i_lev)
        Ly     = D_wind/DX_m

      ENDIF


      Calc_Ly = Ly

    return
  END FUNCTION

!------------------------------------------------------------------
! functions to interpolate wind speed (u,v,omeg) 
! based on the surrounding 4 points.

  real(fp) function Interplt_wind_RLL(wind, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: wind(:,:,:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          :: distance(2,2), Weight(2,2)
    real(fp)          :: wind_lonlat(2), wind_lonlat_lev

    ! Identify wether particle is exactly located on the grid point
    if(curr_pressure==P_mid(i_lev))then
    if(curr_lon==X_mid(i_lon))then
    if(curr_lat==Y_mid(i_lat))then

          Interplt_wind_RLL = wind(i_lon, i_lat, i_lev)
          return

    endif
    endif
    endif


    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southwest of the particle or under
    ! the particle
    if(curr_lon>=X_mid(i_lon))then
      init_lon = i_lon
    else
      init_lon = i_lon - 1
    endif

    if(curr_lat>=Y_mid(i_lat))then
      init_lat = i_lat
    else
      init_lat = i_lat - 1
    endif

    ! For pressure level, P_mid(1) is about surface pressure
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif

    if(init_lev==0) init_lev = 1
    if(init_lev==LLPAR) init_lev = LLPAR-1


    ! calculate the distance between particle and grid point
    do i = 1,2
    do j = 1,2
      ii = i + init_lon - 1
      jj = j + init_lat - 1

      ! For some special circumstance:
      if(ii==0)then
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), Y_mid(jj))
      else if(ii==(IIPAR+1))then
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii-IIPAR), Y_mid(jj))
      else
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
      endif

    enddo
    enddo


    ! Calculate the inverse distance weight
    do i = 1,2
    do j = 1,2

       if(distance(i,j)==0)then
          Weight(:,:) = 0
          Weight(i,j) = 1
          GOTO 100
       endif

       Weight(i,j) = 1.0/distance(i,j) / SUM( 1.0/distance(:,:) )

    enddo
    enddo

 100 CONTINUE


    do k = 1,2     
        kk = k + init_lev - 1 
        if(init_lon==0)then
            wind_lonlat(k) =  Weight(1,1) * wind(IIPAR,init_lat,kk) &
                          + Weight(1,2) * wind(IIPAR,init_lat+1,kk) &
                          + Weight(2,1) * wind(init_lon+1,init_lat,kk) &
                          + Weight(2,2) * wind(init_lon+1,init_lat+1,kk)
        else if(init_lon==IIPAR)then
            wind_lonlat(k) =  Weight(1,1) * wind(init_lon,init_lat,kk) &
                          + Weight(1,2) * wind(init_lon,init_lat+1,kk) &
                          + Weight(2,1) * wind(1,init_lat,kk)   &
                          + Weight(2,2) * wind(1,init_lat+1,kk)
        else
            wind_lonlat(k) =  Weight(1,1) * wind(init_lon,init_lat,kk) &
                          + Weight(1,2) * wind(init_lon,init_lat+1,kk) &
                          + Weight(2,1) * wind(init_lon+1,init_lat,kk) &
                          + Weight(2,2) * wind(init_lon+1,init_lat+1,kk)
        endif
    enddo


    ! second interpolate vertically (Linear)
    IF(P_mid(init_lev+1)==P_mid(init_lev))THEN
      WRITE(6,*)"*** WARNING: two same pressure level happens! ***"
      WRITE(6,*)"init_lev, P_mid(init_lev), P_mid(init_lev+1):"
      WRITE(6,*)init_lev, P_mid(init_lev), P_mid(init_lev+1)
      wind_lonlat_lev = wind_lonlat(1)
    ELSE
      wind_lonlat_lev = wind_lonlat(1) + (wind_lonlat(2)-wind_lonlat(1)) &
                                 / (P_mid(init_lev+1)-P_mid(init_lev)) &
                                     * (curr_pressure-P_mid(init_lev))
    ENDIF

    !Line_Interplt( wind_lonlat(1), wind_lonlat(2), P_mid(i_lev), P_mid(i_lev+1), curr_pressure )

    Interplt_wind_RLL = wind_lonlat_lev


    return
  end function


! functions to interpolate vertical wind speed (w)
! based on the surrounding 3 points, one of the points is the north/south polar
! point. The w value at polar point is the average of all surrounding grid points.

  real(fp) function Interplt_wind_RLL_polar(wind, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: wind(:,:,:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          :: distance(3), Weight(3)
    real(fp)          :: wind_lonlat(2), wind_lonlat_lev
    real(fp)          :: wind_polar

    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southwest of the particle or under
    ! the particle
    if(curr_lon>=X_mid(i_lon))then
      init_lon = i_lon
    else
      init_lon = i_lon - 1
    endif

    if(curr_lat>=Y_mid(i_lat))then
      init_lat = i_lat
    else
      init_lat = i_lat - 1
    endif
    
    if(init_lat==0)then
      init_lat = 1
    endif

    ! For pressure level, P_mid(1) is about surface pressure
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif

    if(init_lev==0) init_lev = 1
    if(init_lev==LLPAR) init_lev = LLPAR-1


    ! calculate the distance between particle and grid point
    j = 1
    do i = 1,2
      ii = i + init_lon - 1
      jj = j + init_lat - 1

      ! For some special circumstance:
      if(ii==0)then
      distance(i)= Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), Y_mid(jj))
      else if(ii==(IIPAR+1))then
      distance(i)= Distance_Circle(curr_lon, curr_lat, X_mid(1), Y_mid(jj))
      else
      distance(i)= Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
      endif

    enddo


   if(ii==0)then
    distance(3)= Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), 90.0e+0_fp)
   else if(ii==(IIPAR+1))then
    distance(3)= Distance_Circle(curr_lon, curr_lat, X_mid(1), 90.0e+0_fp)
   else
    distance(3)= Distance_Circle(curr_lon, curr_lat, X_mid(ii), 90.0e+0_fp)
   endif


   IF(distance(3)==0.0)THEN
       do k=1,2
         kk = k + init_lev - 1
         wind_lonlat(k) = SUM(wind(:,init_lat,kk))/IIPAR
       enddo
   ELSE
       ! Calculate the inverse distance weight
       do i=1,3
          Weight(i) = 1.0/distance(i) / sum( 1.0/distance(:) )
       enddo

       do k=1,2
         kk = k + init_lev - 1

         wind_polar = SUM(wind(:,init_lat,kk))/IIPAR      

         if(init_lon==0)then    
             wind_lonlat(k) =  Weight(1)*wind(IIPAR,init_lat,kk)   &
                             + Weight(2)*wind(init_lon+1,init_lat,kk)   &
                             + Weight(3)*wind_polar
         else if(init_lon==IIPAR)then
             wind_lonlat(k) =  Weight(1)*wind(init_lon,init_lat,kk)   &
                             + Weight(2)*wind(1,init_lat,kk)   &
                             + Weight(3)*wind_polar
         else
             wind_lonlat(k) =  Weight(1)*wind(init_lon,init_lat,kk)   &
                             + Weight(2)*wind(init_lon+1,init_lat,kk)   &
                             + Weight(3)*wind_polar
         endif
     enddo
   ENDIF

    ! second interpolate vertically (Linear)
    IF(P_mid(init_lev+1)==P_mid(init_lev))THEN
      WRITE(6,*)"*** WARNING: two same pressure level happens! ***"
      WRITE(6,*)"init_lev, P_mid(init_lev), P_mid(init_lev+1):"
      WRITE(6,*)init_lev, P_mid(init_lev), P_mid(init_lev+1)
      wind_lonlat_lev = wind_lonlat(1)
    ELSE
      wind_lonlat_lev =   wind_lonlat(1) &
                      + (wind_lonlat(2)-wind_lonlat(1)) &
                      / (P_mid(init_lev+1)-P_mid(init_lev)) &
                      * (curr_pressure-P_mid(init_lev))
    ENDIF

    !Line_Interplt( wind_lonlat(1), wind_lonlat(2), P_mid(i_lev),
    !P_mid(i_lev+1), curr_pressure )

    Interplt_wind_RLL_polar = wind_lonlat_lev

    return
  end function

!------------------------------------------------------------------
! functions to interpolate wind speed (u,v) 
! based on the surrounding 3 points, one of the points is the north/south polar
! point. The u/v value at polar point is the average of all surrounding grid points.

  real(fp) function Interplt_uv_PS_polar(i_uv, u_RLL, v_RLL, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

    implicit none

    real(fp)          :: curr_lon, curr_lat, curr_pressure
    real(fp)          :: curr_x, curr_y
    real(fp), pointer :: u_RLL(:,:,:), v_RLL(:,:,:)

    integer           :: i_uv
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk

    real(fp)          :: x_PS(3), y_PS(3)  ! the third value x_PS(3) is the polar point
    real(fp)          :: uv_PS(3,2)
    real(fp)          :: uv_polars(IIPAR)
    real(fp)          :: distance_PS(3), Weight_PS(3)
    real(fp)          :: uv_xy(2), uv_xy_lev

    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southwest of the particle or under
    ! the particle
    if(curr_lon>=X_mid(i_lon))then
      init_lon = i_lon
    else
      init_lon = i_lon - 1
    endif

    if(curr_lat>=Y_mid(i_lat))then
      init_lat = i_lat
    else
      init_lat = i_lat - 1
    endif

    ! For pressure level, P_mid(1) is about surface pressure, has biggerst
    ! value.
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif

    if(init_lev==0) init_lev = 1
    if(init_lev==LLPAR) init_lev = LLPAR-1


    ! change from (lon,lat) in RLL to (x,y) in PS: 
    if(curr_lat<0)then
      curr_x = -1.0* Re* COS(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
      curr_y = -1.0* Re* SIN(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
    else
      curr_x = Re* COS(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
      curr_y = Re* SIN(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
    endif


    if(init_lat==0)then
           jj = init_lat + 1    ! south polar
       else
           jj = init_lat
    endif
    
    do i=1,2  ! Interpolate location and wind of grid points into Polar Stereographic Plane
   
       ii = i + init_lon - 1

       ! For lon=180 deg:
       if(ii==IIPAR+1)then
          ii = 1
       endif
       if(ii==0)then
          ii = IIPAR
       endif

      ! Interpolate location and wind into Polar Stereographic Plane
      if(Y_mid(jj)>0)then
        x_PS(i)= Re* COS(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
        y_PS(i)= Re* SIN(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
      else
        x_PS(i)= -1.0* Re* COS(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
        y_PS(i)= -1.0* Re* SIN(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
      endif


       do k=1,2
          kk = k + init_lev - 1
          IF(i_uv==1)THEN ! i_ux==1 for u
          if(Y_mid(jj)>0)then
            uv_PS(i,k) = -1.0* ( u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                          + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2) )
          else
            uv_PS(i,k) = u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                        + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
          endif
          ENDIF

          IF(i_uv==0)THEN ! for v
          if(Y_mid(jj)>0)then
            uv_PS(i,k) = u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                       - v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
          else
            uv_PS(i,k) = -1* u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                           + v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
          endif
          ENDIF
       enddo
    enddo

    ! Third grid point: south/north polar point
    x_PS(3) = 0.0
    y_PS(3) = 0.0
    
    do k=1,2
      kk = k + init_lev - 1
      IF(i_uv==1)THEN ! i_ux==1 for u
      ! interpolate all the grid points surrounding the polar point:
      do ii = 1,IIPAR
        if(Y_mid(jj)>0)then
          uv_polars(ii) = -1.0* ( u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                         / SIN(Y_mid(jj)*PI/180.0) &
                         + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                         / (SIN(Y_mid(jj)*PI/180.0)**2) )
          else
            uv_polars(ii) = u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                          + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
          endif
        enddo

            uv_PS(3,k) = SUM(uv_polars)/IIPAR
          ENDIF

          IF(i_uv==0)THEN ! for v
          do ii = 1,IIPAR
             if(Y_mid(jj)>0)then
               uv_polars(ii) = u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                             - v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
             else
               uv_polars(ii) = -1* u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                                 + v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
             endif
          enddo
             uv_PS(3,k) = SUM(uv_polars)/IIPAR
          ENDIF
    enddo



    ! calculate the distance between particle and grid point
    do i = 1,3
       distance_PS(i)= SQRT( (x_PS(i)-curr_x)**2.0 + (y_PS(i)-curr_y)**2.0 )
    enddo

    ! Calculate the inverse distance weight
    do i = 1,3
        Weight_PS(i) = 1.0/distance_PS(i) / sum( 1.0/distance_PS(:) )
    enddo


    do k = 1,2
      uv_xy(k) = Weight_PS(1) * uv_PS(1,k) &
               + Weight_PS(2) * uv_PS(2,k) &
               + Weight_PS(3) * uv_PS(3,k)
    enddo


    ! second interpolate vertically (Linear)

    uv_xy_lev = uv_xy(1)+ (uv_xy(2)-uv_xy(1)) &
     / (P_mid(init_lev+1)-P_mid(init_lev)) * (curr_pressure-P_mid(init_lev))

    Interplt_uv_PS_polar = uv_xy_lev

    return
  end function


!------------------------------------------------------------------
! functions to interpolate wind speed (u,v,omeg) 
! based on the surrounding 4 points.

  real(fp) function Interplt_uv_PS(i_uv, u_RLL, v_RLL, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

    implicit none

    real(fp)          :: curr_lon, curr_lat, curr_pressure
    real(fp)          :: curr_x, curr_y
    real(fp), pointer :: u_RLL(:,:,:), v_RLL(:,:,:)

    integer           :: i_uv
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk

    real(fp)          :: x_PS(2,2), y_PS(2,2)
    real(fp)          :: uv_PS(2,2,2)
    real(fp)          :: distance_PS(2,2), Weight_PS(2,2)
    real(fp)          :: uv_xy(2), uv_xy_lev


    ! first interpolate horizontally (Inverse Distance Weighting)

 ! identify the grid point located in the southwest of the particle or under
    ! the particle
    if(curr_lon>=X_mid(i_lon))then
      init_lon = i_lon
    else
      init_lon = i_lon - 1
    endif

    if(curr_lat>=Y_mid(i_lat))then
      init_lat = i_lat
    else
      init_lat = i_lat - 1
    endif

!For pressure level, P_mid(1) is about surface pressure, has biggerst value.
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif

    if(init_lev==0) init_lev = 1
    if(init_lev==LLPAR) init_lev = LLPAR-1

    
    ! change from (lon,lat) in RLL to (x,y) in PS: 
    if(curr_lat<0)then
      curr_x = -1.0* Re* COS(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
      curr_y = -1.0* Re* SIN(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
    else
      curr_x = Re* COS(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
      curr_y = Re* SIN(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
    endif

    ! i,j means the four grid point value that surround around the particle
    do i=1,2
    do j=1,2

      ii = i + init_lon - 1
      jj = j + init_lat - 1

      ! Get the right ii and jj for interpolation at polar point:
      ! For South Polar Point:
      if(jj==0)then
        jj = jj+1
      endif
      ! For North Polar Point:
      if(jj==JJPAR+1)then
        jj = jj-1
      endif

    
      ! For lon=180 deg:
      if(ii==IIPAR+1)then
        ii = 1
      endif
      if(ii==0)then
        ii = IIPAR
      endif


      ! Interpolate location and wind into Polar Stereographic Plane
      if(Y_mid(jj)>0)then

        x_PS(i,j) = Re* COS(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)  
        y_PS(i,j) = Re* SIN(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
          
        do k=1,2
           kk = k + init_lev - 1

           if(i_uv==1)then ! i_ux==1 for u
             uv_PS(i,j,k)= -1.0* ( u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                         / SIN(Y_mid(jj)*PI/180.0) &
                                + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2) )
           endif

           if(i_uv==0)then ! for v
             uv_PS(i,j,k) = u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                          - v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
           endif
        enddo

      else

      x_PS(i,j)= -1.0* Re* COS(X_mid(ii)*PI/180.0) /TAN(Y_mid(jj)*PI/180.0)
      y_PS(i,j)= -1.0* Re* SIN(X_mid(ii)*PI/180.0) /TAN(Y_mid(jj)*PI/180.0)

        do k=1,2
           kk = k + init_lev - 1

           if(i_uv==1)then
             uv_PS(i,j,k) = u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                          + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
           endif

           if(i_uv==0)then
             uv_PS(i,j,k) = -1* u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                              + v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
           endif
        enddo

      endif

    enddo
    enddo

    ! calculate the distance between particle and grid point
    do i = 1,2
    do j = 1,2
       distance_PS(i,j) = SQRT( (x_PS(i,j)-curr_x)**2.0 &
                                                + (y_PS(i,j)-curr_y)**2.0 )
    enddo
    enddo

    ! Calculate the inverse distance weight
    do i = 1,2
    do j = 1,2
        Weight_PS(i,j) = 1.0/distance_PS(i,j) / sum( 1.0/distance_PS(:,:) )
    enddo
    enddo


    do k = 1,2
      uv_xy(k) =  Weight_PS(1,1) * uv_PS(1,1,k) &
                 + Weight_PS(1,2) * uv_PS(1,2,k) &
                 + Weight_PS(2,1) * uv_PS(2,1,k) &
                 + Weight_PS(2,2) * uv_PS(2,2,k)
    enddo


    ! second interpolate vertically (Linear)

    uv_xy_lev = uv_xy(1)+ &
              (uv_xy(2)-uv_xy(1)) / (P_mid(init_lev+1)-P_mid(init_lev)) &
                                        * (curr_pressure-P_mid(init_lev))

    Interplt_uv_PS = uv_xy_lev

    return
  end function


!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------

  SUBROUTINE plume_run(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)

    USE Plume_list

    USE Input_Opt_Mod,   ONLY : OptInput

    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Met_Mod,   ONLY : MetState

    USE TIME_MOD              ! For computing date & time 
    USE TIME_MOD,        ONLY : GET_TS_DYN
    USE TIME_MOD,        ONLY : ITS_TIME_FOR_EXIT

    USE TIME_MOD,        ONLY : GET_YEAR
    USE TIME_MOD,        ONLY : GET_MONTH
    USE TIME_MOD,        ONLY : GET_DAY
    USE TIME_MOD,        ONLY : GET_HOUR
    USE TIME_MOD,        ONLY : GET_MINUTE

    USE State_Grid_Mod,  ONLY : GrdState
!    USE GC_GRID_MOD,     ONLY : XEDGE, YEDGE, XMID, YMID
!    USE CMN_SIZE_Mod,    ONLY : IIPAR, JJPAR, LLPAR, DLAT, DLON
    ! DLAT( IIPAR, JJPAR, LLPAR ), DLON( IIPAR, JJPAR, LLPAR )
    ! XEDGE( IM+1, JM,   L ), YEDGE( IM,   JM+1, L ), IM=IIPAR, JM=JJPAR

    USE UnitConv_Mod,    ONLY : Convert_Spc_Units


    INTERFACE
      RECURSIVE FUNCTION SortList(head0) RESULT (new_head)
        USE Plume_list
        TYPE(Plume1d_list), POINTER :: head0, new_head
      END FUNCTION
    END INTERFACE



    logical, intent(in)           :: am_I_Root
    TYPE(MetState), intent(in)    :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State objectgg
    TYPE(OptInput), intent(in)    :: Input_Opt

    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    REAL(fp) :: Dt          ! = 600.0e+0_fp          
    REAL(fp) :: Pdt
    REAL(fp) :: Dt2


    integer  :: i_box   !, i_species
    integer  :: ii_box
    integer  :: i_x, i_y, i, j
    integer  :: i_slab
    integer  :: i_cell

    integer  :: i_lon, i_lat, i_lev
    integer  :: i_species, i_advect
    integer  :: i_advect1, i_advect2

    integer  :: Nt

    integer  :: t1s

    integer :: alloc_stat

    integer :: Stop_loop 

    integer :: Type_transfer

    real(fp) :: curr_lon, curr_lat, curr_pressure
!    real(fp) :: next_lon, next_lat, next_pressure

    real(fp) :: X_edge2, Y_edge2

    real(fp), pointer :: u(:,:,:)
    real(fp), pointer :: v(:,:,:)
    real(fp), pointer :: omeg(:,:,:)
    real(fp), pointer :: Ptemp(:,:,:)
    real(fp), pointer :: P_BXHEIGHT(:,:,:)

    ! cumulative volume of all plume in a Eulerian grid cell
    real(fp), dimension(:), allocatable  :: SumV_Plume, Entropy_V
!    real(fp)          :: curr_u, curr_v, curr_omeg


    real(fp)          :: curr_Ptemp, Ptemp_shear  ! potential temperature
    real(fp)          :: U_shear, V_shear, UV_shear

    real(fp)          :: wind_s_shear
    real(fp)          :: theta_previous


    real(fp), pointer :: X_edge(:)
    real(fp), pointer :: Y_edge(:)


    real(fp)  :: Outer(n_slab_max), Inner(n_slab_max)


    real(fp)  :: eddy_v, eddy_h
    real(fp)  :: eddy_A, eddy_B
    real(fp)  :: Cv, Ch, Omega_N, N_BV

    real(fp)  :: grid_volumn
    real(fp)  :: V_grid_1D, V_grid_2D
    real(fp)  :: mass_plume, mass_plume_new
    real(fp)  :: D_mass_plume !, D_mass_extra ! mass change in plume and extra
                                            ! model every time step [molec]

    real(fp)  :: Ratio_radius, Ra_Rb

    real(fp)  :: L_b, L_O, Ee

    real(fp)  :: CFL
    real(fp)  :: CFL_2d(n_x_max2,n_y_max2) ! for 2D advection

    real(fp)  :: Pu(n_x_max2,n_y_max2) ! for 2D advection


    real(fp)  :: Pc_middle
    real(fp)  :: Pc_bottom, Pc_top
    real(fp)  :: Pc_left, Pc_right


    real(fp)  :: Pc(n_x_max,n_y_max),Pc2(n_x_max,n_y_max) !, Ec(n_x_max,n_y_max)
    real(fp)  :: Pc_bdy(n_x_max2,n_y_max2)
!    real(fp)  :: D_concnt(n_x_max,n_y_max)
    real(fp)  :: Xscale, Yscale, frac_mass
    real(fp)  :: C2d_prev(n_x_max,n_y_max) !, C2d_prev_extra(n_x_max,n_y_max)

    real(fp)  :: Cslab(n_slab_max) !, Extra_Cslab(n_slab_max)

    real(fp)  :: Concnt2D_bdy(n_x_max2, n_y_max2)
    real(fp)  :: Concnt1D_bdy(n_slab_max2)

    real(fp)  :: C1d_prev(n_slab_max) !, C1d_prev_extra(n_slab_max)

    real(fp)  :: DlonN, DlatN ! split 1D from 1 to 5 segments

    real(fp)  :: backgrd_concnt

    real(fp)  :: Prod_before, Prod_after, Eul_concnt

    real(fp), allocatable :: box_lon_new(:), box_lat_new(:)

    real(fp)  :: start, finish
    real(fp)  :: start2, finish2


    real(fp) :: box_lon, box_lat, box_lev
    real(fp) :: box_length, box_alpha, box_theta
    real(fp) :: Pdx, Pdy
    real(fp) :: box_Ra, box_Rb
    real(fp) :: box_extra, box_life, box_label

    real(fp) :: box_concnt_2D(n_x_max, n_y_max, n_species)
    real(fp) :: box_concnt_1D(n_slab_max, n_species)

    real(fp) :: Max_Concnt_Grid
    real(fp) :: Entropy_Concnt, Entropy
    real(fp) :: tracer_mol, air_mol, mix_ratio


    TYPE(Plume2d_list), POINTER :: Plume2d_new, PLume2d, Plume2d_prev
    TYPE(Plume1d_list), POINTER :: Plume1d_new, Plume1d, Plume1d_prev


    CHARACTER(LEN=255)     :: FILENAME2
    !CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: ThisLoc
    CHARACTER(LEN=255)     :: ErrMsg

    INTEGER :: YEAR
    INTEGER :: MONTH
    INTEGER :: DAY
    INTEGER :: HOUR
    INTEGER :: MINUTE
    INTEGER :: SECOND
    INTEGER :: OrigUnit

    CHARACTER(LEN=255)     :: FileConcnt, FileXYZ, FileEntropy, File996

    CHARACTER(LEN=25) :: YEAR_C
    CHARACTER(LEN=25) :: MONTH_C
    CHARACTER(LEN=25) :: DAY_C
    CHARACTER(LEN=25) :: HOUR_C
    CHARACTER(LEN=25) :: MINUTE_C
    CHARACTER(LEN=25) :: SECOND_C


    Dt = GET_TS_DYN()
    ThisLoc                =   ' -> at plume_run (in module GeosCore/lagrange_mod.F90)'
    Entropy = 0.0


    YEAR        = GET_YEAR()
    MONTH       = GET_MONTH()
    DAY         = GET_DAY()
    HOUR        = GET_HOUR()
    MINUTE      = GET_MINUTE()
    SECOND      = GET_SECOND()

    WRITE(YEAR_C,*) YEAR
    WRITE(MONTH_C,*) MONTH
    WRITE(DAY_C,*) DAY
    WRITE(HOUR_C,*) HOUR
    WRITE(MINUTE_C,*) MINUTE
    WRITE(SECOND_C,*) SECOND



    IF(Stop_inject==1) GOTO 999

    Stop_loop = 0


    RC        =  GC_SUCCESS
    ErrMsg    =  ''


    FILENAME2   = 'Plume_lifetime_seconds.txt'

    OPEN( 484,      FILE=TRIM( FILENAME2   ), STATUS='OLD',  &
          POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )


    Dt = GET_TS_DYN()

    u => State_Met%U ! [m/s]
    v => State_Met%V ! V [m s-1]
    omeg => State_Met%OMEGA  ! Updraft velocity [Pa/s]

    Ptemp => State_Met%THETA ! ! Potential temperature [K]

    P_BXHEIGHT => State_Met%BXHEIGHT  ![IIPAR,JJPAR,KKPAR]



    allocate(SumV_Plume(LLPAR*JJPAR*IIPAR))
    allocate(Entropy_V(LLPAR*JJPAR*IIPAR))
    SumV_Plume = 0.0
    Entropy_V  = 0.0



    X_edge => State_Grid%XEdge(:,1) !XEDGE(:,1,1)   ! IIPAR+1 ! new
    Y_edge => State_Grid%YEdge(1,:) !YEDGE(1,:,1) 
    ! Use second YEDGE, because sometimes YMID(2)-YMID(1) is not DLAT

    X_edge2       = X_edge(2)
    Y_edge2       = Y_edge(2)

    !======================================================================
    ! Convert species to [molec/cm3] (ewl, 8/16/16)
    !======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                              new_units=MOLECULES_SPECIES_PER_CM3, RC = RC, previous_units=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'plume_run in lagrange_mod.F90')
       RETURN
    ENDIF


    !======================================================================
    ! run the fake 2nd order chemical reaction
    !======================================================================


    !-----------------------------------------------------------------------

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( i_lev, i_lat, i_lon, i_cell )
    DO i_lev = 1, LLPAR
    DO i_lat = 1, JJPAR
    DO i_lon = 1, IIPAR

      ! set initial value for Entropy_V
      i_cell = (i_lev-1)*IIPAR*JJPAR+(i_lat-1)*IIPAR+i_lon
      Entropy_V(i_cell) = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp


      ! from LA to LA2
      State_Chm%Species(id_PASV_LA2)%Conc(i_lon,i_lat,i_lev) = &
                 State_Chm%Species(id_PASV_LA2)%Conc(i_lon,i_lat,i_lev) &
               + Dt*Kchem*State_Chm%Species(id_PASV_LA)%Conc(i_lon,i_lat,i_lev)**2


      !  For all the grid cells in the troposphere, let concentration to be zero
      IF(TROPP_sink == 1)THEN
      IF ( State_Met%PEDGE(i_lon,i_lat,i_lev)>State_Met%TROPP(i_lon,i_lat) ) THEN
        !State_Chm%Species(:)%Conc(i_lon,i_lat,i_lev) = 0.0
		DO nspec = 1, State_Chm%nSpecies
			State_Chm%Species(nspec)%Conc(i_lon,i_lat,i_lev) = 0.0
		END DO
      ENDIF
      ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !=====================================================================
    ! For 2D plume: Run distortion & dilution HERE
    !=====================================================================

    IF(ASSOCIATED(Plume2d_prev)) NULLIFY(Plume2d_prev)
    IF(.NOT.ASSOCIATED(Plume2d_head)) GOTO 900



    Plume2d => Plume2d_head
    i_box = 0

    DO WHILE(ASSOCIATED(Plume2d))

      i_box = i_box+1

      box_lon    = Plume2d%LON
      box_lat    = Plume2d%LAT
      box_lev    = Plume2d%LEV
      box_length = Plume2d%LENGTH
      box_alpha  = Plume2d%ALPHA

      box_label  = Plume2d%label
      box_life   = Plume2d%LIFE

      Pdx    = Plume2d%PDX
      Pdy    = Plume2d%PDY

      box_concnt_2D = Plume2d%CONCNT2d



      ! make sure the location is not out of range
      do while (box_lat > Y_edge(JJPAR+1))
         box_lat = Y_edge(JJPAR+1) &
                        - ( box_lat-Y_edge(JJPAR+1) )
      end do

      do while (box_lat < Y_edge(1))
         box_lat = Y_edge(1) + ( box_lat-Y_edge(1) )
      end do

      do while (box_lon > X_edge(IIPAR+1))
         box_lon = box_lon - 360.0
      end do

      do while (box_lon < X_edge(1))
         box_lon = box_lon + 360.0
      end do


      curr_lon      = box_lon
      curr_lat      = box_lat
      curr_pressure = box_lev      ! hPa



       i_lon = Find_iLonLat(curr_lon, DX, X_edge2)
       if(i_lon>IIPAR) i_lon=i_lon-IIPAR
       if(i_lon<1) i_lon=i_lon+IIPAR

       i_lat = Find_iLonLat(curr_lat, DY, Y_edge2)
       if(i_lat>JJPAR) i_lat=JJPAR
       if(i_lat<1) i_lat=1

       i_lev = Find_iPLev(curr_pressure,P_edge)


!       backgrd_concnt(i_box,:) = State_Chm%Species(i_lon,i_lat,i_lev,:)
!       backgrd_concnt = State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1)
       ! [molec/cm3]

!       Plume_I(i_box) = i_lon
!       Plume_J(i_box) = i_lat
!       Plume_L(i_box) = i_lev




      backgrd_concnt = State_Chm%Species(id_PASV_LA)%Conc(i_lon,i_lat,i_lev)
      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( i_y, i_x )
      DO i_y = 1,n_y_max,1
      DO i_x = 1,n_x_max,1

      ! run the fake 2nd order chemical reaction
      box_concnt_2D(i_x,i_y,2) = box_concnt_2D(i_x,i_y,2) &
         + Dt* Kchem*(box_concnt_2D(i_x,i_y,1)+backgrd_concnt)**2 &
         - Dt* Kchem*backgrd_concnt**2

      ENDDO
      ENDDO
      !$OMP END PARALLEL DO




       !====================================================================
       ! calculate the wind shear along plume corss-section
       ! clock-wise 90 degree from the plume length direction (box_alpha)
       ! calculate the diffusivity in horizontal and vertical direction
       !====================================================================

       ! calculate the wind_s shear along pressure direction
       wind_s_shear = Wind_shear_s(u, v, P_BXHEIGHT, box_alpha, i_lon, &
                        i_lat, i_lev,curr_lon, curr_lat, curr_pressure)


       ! Calculate vertical eddy diffusivity (U.Schumann, 2012) :
       Cv = 0.2
       Omega_N = 0.1
       Ptemp_shear = Vertical_shear(Ptemp, P_BXHEIGHT, i_lon, i_lat, &
                             i_lev,curr_lon, curr_lat, curr_pressure)


       !--------------------------------------------------------------------
       ! interpolate potential temperature for Plume module:
       !--------------------------------------------------------------------
       if(abs(curr_lat)>Y_mid(JJPAR))then
         curr_Ptemp = Interplt_wind_RLL_polar(Ptemp, i_lon, i_lat, &
                                 i_lev, curr_lon, curr_lat, curr_pressure)
       else
         curr_Ptemp = Interplt_wind_RLL(Ptemp, i_lon, i_lat, i_lev, &
                                        curr_lon, curr_lat, curr_pressure)
       endif


       N_BV = SQRT(Ptemp_shear*g0/curr_Ptemp)
       IF(N_BV<=0.001) N_BV = 0.001

       ! diffusivity unit: [m2/s]
       eddy_v = Cv * Omega_N**2 / N_BV
       eddy_h = 10.0


       ! Define the wind field based on wind shear
       DO i_y = 1, n_y_max2
         Pu(:,i_y) = (i_y-n_y_mid2)*Pdy * wind_s_shear
         ! [m s-1]
       ENDDO



         !--------------------------------------------------------------
         ! if Pdx become smaller than half Dx_initm,
         ! combine 9 grids into 1 grids. 
         ! Update: box_concnt_2D(), Pdx(), Pdy(),  Extra_mass_2D()
         !--------------------------------------------------------------
         IF(Pdx<=0.5*Dx_init)THEN

600     CONTINUE

           DO i_species = 1, n_species
         
           Pc = box_concnt_2D(:,:,i_species) ! [molec cm-3]
           box_concnt_2D(:,:,i_species) = 0.0

           DO i = 1, n_x_max/3, 1
           DO j = 1, n_y_max/3, 1

             i_x = (i-1)*3+1
             i_y = (j-1)*3+1

             box_concnt_2D(i+n_x_max/3,j+n_y_max/3,i_species) = &
                                      SUM(Pc(i_x:i_x+2,i_y:i_y+2))/9

           ENDDO
           ENDDO

           ENDDO ! DO i_species = 1, n_species


           Pdx = Pdx*3
           Pdy = Pdy*3



           IF(Pdx<=0.5*Dx_init) WRITE(6,*) &
               "ERROR 0: more combination: ", i_box, Plume2d%label, Pdx, Pdy

           IF(Pdx<=0.5*Dx_init) GOTO 600


         ENDIF ! IF(Pdx(i_box)<0.5*Dx_init)THEN



       !-------------------------------------------------------------------
       ! Calculate the advection-diffusion in 2D grids
       !-------------------------------------------------------------------
       V_grid_2D       = Pdx*Pdy*box_length*1.0e+6_fp


       DO i_species=1,n_species,1

       C2d_prev(1:n_x_max,1:n_y_max) = &
                              box_concnt_2D(1:n_x_max,1:n_y_max,i_species)

       Nt = CEILING(Dt/120)
       Pdt = Dt/Nt ! Dt=600 ! FLOOR(Pdx/Pu(1,1)/10)*10

       ! Find the best Pdt to meet CFL condition:
 700     CONTINUE


       CFL = Pdt*Pu(1,1)/Pdx
       IF(MAX( ABS(CFL), ABS(2*eddy_h*Pdt/(Pdx**2)), &
                         ABS(2*eddy_v*Pdt/(Pdy**2)) ) > 0.8)THEN

          Nt = Nt+1
          Pdt = Dt/Nt
          GOTO 700

        ENDIF


       Concnt2D_bdy(:,:) = 0.0
       Concnt2D_bdy(2:n_x_max2-1,2:n_y_max2-1) = box_concnt_2D(:,:,i_species)

        
       IF( abs(Dt/Pdt-Nt) > 0.00001 ) THEN
         WRITE(6,*) "*** ERROR: Check Pdt ***"
         WRITE(6,*) Pdt, Nt, Dt
         WRITE(6,*) Pdx/Pu(1,1), Pdy**2/(2*eddy_v), Pdx**2/(2*eddy_h)
       ENDIF
       



       DO t1s = 1, NINT(Dt/Pdt)

         ! advection ----------------------------------------------------
         ! diffusion ----------------------------------------------------

         Pc_bdy(:,:) = 0.0
         Pc_bdy(2:n_x_max2-1,2:n_y_max2-1) = Concnt2D_bdy(2:n_x_max2-1,2:n_y_max2-1)


         ! Only calculate the vertical half 2D domain         

         !$OMP PARALLEL DO           &
         !$OMP DEFAULT( SHARED     ) &
         !$OMP PRIVATE(i_y,i_x,CFL,Pc_middle,Pc_top,Pc_bottom,Pc_right,Pc_left)
         DO i_y = 1, n_y_mid2, 1
         DO i_x = 2, n_x_max2-1, 1

           Pc_middle = Pc_bdy( i_x,   i_y  )
           Pc_top    = Pc_bdy( i_x,   i_y+1)
           Pc_bottom = Pc_bdy( i_x,   i_y-1)
           Pc_right  = Pc_bdy( i_x+1, i_y  )
           Pc_left   = Pc_bdy( i_x-1, i_y  )

           CFL       = Pdt*Pu(i_x,i_y)/Pdx

           Concnt2D_bdy(i_x, i_y) = Pc_middle           &
              - 0.5 * CFL    * ( Pc_right - Pc_left )   &
              + 0.5 * CFL**2 * ( Pc_right - 2*Pc_middle + Pc_left )         &
              + Pdt*( eddy_h*( Pc_right -2*Pc_middle +Pc_left   ) /(Pdx**2) &
                     +eddy_v*( Pc_top   -2*Pc_middle +Pc_bottom ) /(Pdy**2) )

         ENDDO
         ENDDO
         !$OMP END PARALLEL DO

         ! update the other half based on vertical symmetry 
         DO i_y = n_y_mid2+1, n_y_max2-1, 1
           Concnt2D_bdy(2:n_x_max2-1:1, i_y) = Concnt2D_bdy(n_x_max2-1:2:-1,n_y_max2+1-i_y)
         ENDDO

       ENDDO ! DO t1s = 1, NINT(Dt/Pdt)


         box_concnt_2D(:,:,i_species) = Concnt2D_bdy(2:n_x_max2-1,2:n_y_max2-1)

         !================================================================
         ! Update the concentration in the background grid cell
         ! after the interaction with 2D plume
         !================================================================

         grid_volumn     = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ! [cm3]
         V_grid_2D       = Pdx*Pdy*box_length*1.0e+6_fp


         ! the boundary always represents the background concentration
         ! [molec]
         mass_plume      = V_grid_2D*SUM(C2d_prev(:,:))

         mass_plume_new = V_grid_2D * SUM(box_concnt_2D(:,:,i_species))


         D_mass_plume = mass_plume_new - mass_plume


         i_advect = id_PASV_LA +i_species -1

         backgrd_concnt = State_Chm%Species(i_advect)%Conc(i_lon,i_lat,i_lev)

         backgrd_concnt = ( backgrd_concnt*grid_volumn &
                                           - D_mass_plume) /grid_volumn

         State_Chm%Species(i_advect)%Conc(i_lon,i_lat,i_lev) = backgrd_concnt


       ENDDO ! DO i_species=1,n_species,1



         !====================================================================
         ! Change from 2D to 1D, 
         ! once the tilting degree is bigger than 88 deg (88/180*3.14)
         !====================================================================
         

         IF(Plume2d%LIFE>4.0*3600.0)THEN

         frac_mass = 0.95


         Pc2 = box_concnt_2D(:,:,1)
         Xscale = Get_XYscale(Pc2, Pdx, Pdy, frac_mass, 2)
         Yscale = Get_XYscale(Pc2, Pdx, Pdy, frac_mass, 1)

         box_theta = ATAN( Xscale/Yscale )



         IF( Xscale/Yscale>25.0 &
                        .or. Plume2d%LIFE>12.0*3600.0 ) THEN !!! shw ???

           DO i_species=1,n_species,1

           Pc = box_concnt_2D(:,:,i_species)


           IF(box_theta>1.569) WRITE(6,*)&
            '*** ERROR *** i_box, box_theta:', i_box, box_theta,&
                                                     Plume2d%LIFE

	   !-------------------------------------------------------------------
           ! assign length/width and concentration to 1D slab model 
           ! return initial box_Ra(i_box), box_Rb(i_box), box_concnt_1D(i_box)
	   ! 
	   ! There are two interpolation options for: bilnear and conservative
	   !-------------------------------------------------------------------

           Cslab = box_concnt_1D(:,i_species)
	
	   ! bilinear:
           CALL Slab_init_bilinear(Pdx, Pdy, box_theta, Pc, Yscale, box_Ra, box_Rb, Cslab, n_slab_max)

	   ! conservative:
!           CALL Slab_init_conservative(Pdx, Pdy, box_theta, Pc, Yscale, box_Ra, box_Rb, Cslab, n_slab_max, n_x_max, n_y_max)

           box_concnt_1D(:,i_species) = Cslab


           V_grid_1D = box_Ra*box_Rb*box_length*1.0e+6_fp

           mass_plume_new = V_grid_1D &
                            * SUM(box_concnt_1D(1:n_slab_max,i_species))

           mass_plume = V_grid_2D * SUM(box_concnt_2D(:,:,i_species))



           !!! update the background concentration:
           D_mass_plume = mass_plume_new - mass_plume



           IF(ABS(D_mass_plume/mass_plume)>0.01)THEN
             WRITE(6,*) '*********************************************'
             WRITE(6,*) '    more than 1% mass lost from 2-D to 1-D'
             WRITE(6,*) ABS(D_mass_plume/mass_plume)
             WRITE(6,*) '*********************************************'
           ENDIF



           i_advect = id_PASV_LA +i_species -1

           backgrd_concnt = State_Chm%Species(i_advect)%Conc(i_lon,i_lat,i_lev)

           backgrd_concnt = ( backgrd_concnt*grid_volumn &
                                            - D_mass_plume ) /grid_volumn

           State_Chm%Species(i_advect)%Conc(i_lon,i_lat,i_lev) = backgrd_concnt

           ENDDO ! DO i_species=1,n_species,1


           ! ----------------------------------
           ! delete the node for 2D plume:
           ! ----------------------------------

           IF(.NOT.ASSOCIATED(Plume2d%next) .and. &
                          .NOT.ASSOCIATED(Plume2d_prev))THEN

             WRITE(6,*)'                '
             WRITE(6,*)'*** Only left the last Plume2d', Num_Plume2d
             WRITE(6,*)'                '
             DEALLOCATE(Plume2d_head)

             i_box = i_box-1
             Num_Plume2d = Num_Plume2d - 1
 
             Stop_loop = 1
             GOTO 901

           ENDIF


           IF(.NOT.ASSOCIATED(Plume2d%next))THEN ! delete last node
              Plume2d_tail%next => Plume2d_prev
              Plume2d_tail => Plume2d_tail%next

              IF(ASSOCIATED(Plume2d_prev%next)) NULLIFY(Plume2d_prev%next)
              DEALLOCATE(Plume2d%CONCNT2d)
              DEALLOCATE(Plume2d)
              i_box = i_box-1
              Num_Plume2d = Num_Plume2d - 1

              Stop_loop = 1
              GOTO 901

           ELSEIF(ASSOCIATED(Plume2d_prev))THEN ! delete node, not head/tail
              Plume2d => Plume2d_prev%next
              Plume2d_prev%next => Plume2d%next
              DEALLOCATE(Plume2d%CONCNT2d)
              DEALLOCATE(Plume2d)
              Plume2d => Plume2d_prev%next
           ELSE ! delete head node
              Plume2d => Plume2d_head
              Plume2d_head => Plume2d%next
              DEALLOCATE(Plume2d%CONCNT2d)
              DEALLOCATE(Plume2d)
              Plume2d => Plume2d_head
           ENDIF

           i_box = i_box-1
           Num_Plume2d = Num_Plume2d - 1


           ! ----------------------------------
           ! add new node for 1D plume:
           ! ----------------------------------
901   CONTINUE

           IF(.NOT.ASSOCIATED(Plume1d_head))THEN ! create first node for 1d plume

             ALLOCATE(Plume1d_tail)
             Plume1d_tail%LON = box_lon
             Plume1d_tail%LAT = box_lat
             Plume1d_tail%LEV = box_lev

             Plume1d_tail%LENGTH = box_length
             Plume1d_tail%ALPHA  = box_alpha

             Plume1d_tail%Is_transfer    = 0

             Plume1d_tail%label = box_label
             Plume1d_tail%LIFE = box_life

             Plume1d_tail%RA = box_Ra
             Plume1d_tail%RB = box_Rb

             ALLOCATE(Plume1d_tail%CONCNT1d(n_slab_max,n_species))
             Plume1d_tail%CONCNT1d = box_concnt_1D

             NULLIFY(Plume1d_tail%next)
             Plume1d_head => Plume1d_tail


             WRITE(6,*) '*** CREATE first Plume1d node ***', box_label

           ELSE

             ALLOCATE(Plume1d_new, stat=alloc_stat)
             IF(alloc_stat/=0) WRITE(6,*)'ERROR 5:', i_box, alloc_stat

             Plume1d_new%label = box_label

             Plume1d_new%LON = box_lon
             Plume1d_new%LAT = box_lat
             Plume1d_new%LEV = box_lev

             Plume1d_new%LENGTH = box_length
             Plume1d_new%ALPHA = box_alpha

             Plume1d_new%Is_transfer = 0

             Plume1d_new%label  = box_label
             Plume1d_new%LIFE = box_life

             Plume1d_new%RA = box_Ra
             Plume1d_new%RB = box_Rb

             ALLOCATE(Plume1d_new%CONCNT1d(n_slab_max,n_species))
             Plume1d_new%CONCNT1d = box_concnt_1D

             NULLIFY(Plume1d_new%next)
             Plume1d_tail%next => Plume1d_new
             Plume1d_tail => Plume1d_tail%next


           ENDIF

           Num_Plume1d = Num_Plume1d + 1

           IF(Stop_loop==1)THEN
             GOTO 900 ! skip the whole loop, go to Plume1d
           ELSE
             GOTO 800 ! skip this box, go to next box
           ENDIF

         ENDIF ! IF(box_theta(i_box)>98/180*3.14)
         ENDIF ! IF(Lifetime(i_box)>3*3600.0)THEN



         !===================================================================
         ! Fourth judge:
         ! If the plume touch the tropopause, dissolve the plume
         !===================================================================
         IF ( box_lev>State_Met%TROPP(i_lon,i_lat) ) THEN

           Num_dissolve = Num_dissolve+1


           DO i_species=1,n_species

           i_advect = id_PASV_LA +i_species -1

           backgrd_concnt = State_Chm%Species(i_advect)%Conc(i_lon,i_lat,i_lev)

           backgrd_concnt = ( SUM( box_concnt_2D(:,:,i_species) )*V_grid_2D &
                             + backgrd_concnt*grid_volumn ) / grid_volumn


           State_Chm%Species(i_advect)%Conc(i_lon,i_lat,i_lev) = backgrd_concnt

           ENDDO

           ! ----------------------------------
           ! delete the node for 2D plume:
           ! ----------------------------------
           IF(.NOT.ASSOCIATED(Plume2d%next) .and. &
                          .NOT.ASSOCIATED(Plume2d_prev))THEN
              GOTO 900
           ENDIF


           IF(.NOT.ASSOCIATED(Plume2d%next))THEN ! delete last node
              Plume2d_tail%next => Plume2d_prev
              Plume2d_tail => Plume2d_tail%next

!            Plume2d => Plume2d_prev%next !!! ???
              IF(ASSOCIATED(Plume2d_prev%next)) NULLIFY(Plume2d_prev%next)
              DEALLOCATE(Plume2d%CONCNT2d)
              DEALLOCATE(Plume2d)
              i_box = i_box-1
              Num_Plume2d = Num_Plume2d - 1
              GOTO 900
           ELSEIF(ASSOCIATED(Plume2d_prev))THEN ! delete node, not head/tail
              Plume2d => Plume2d_prev%next
              Plume2d_prev%next => Plume2d%next
              DEALLOCATE(Plume2d%CONCNT2d)
              DEALLOCATE(Plume2d)
              Plume2d => Plume2d_prev%next
           ELSE ! delete head node
              Plume2d => Plume2d_head
              Plume2d_head => Plume2d%next
              DEALLOCATE(Plume2d%CONCNT2d)
              DEALLOCATE(Plume2d)
              Plume2d => Plume2d_head
           ENDIF

           i_box = i_box-1
           Num_Plume2d = Num_Plume2d - 1


           WRITE(6,*)'333'

           GOTO 800 ! skip this box, go to next box

         ENDIF ! IF ( box_lev>State_Met%TROPP(i_lon,i_lat) ) THEN




       ! add plume volume to SumV_Plume
       i_cell = (i_lev-1)*IIPAR*JJPAR+(i_lat-1)*IIPAR+i_lon
       SumV_Plume(i_cell) = SumV_Plume(i_cell) + V_grid_2D*n_x_max*n_y_max
       


       Plume2d%LON = box_lon
       Plume2d%LAT = box_lat
       Plume2d%LEV = box_lev
       Plume2d%LENGTH = box_length
       Plume2d%ALPHA  = box_alpha
 
       Plume2d%PDX = Pdx
       Plume2d%PDY = Pdy
 
       Plume2d%CONCNT2d    = box_concnt_2D
 

       Plume2d_prev => Plume2d
       Plume2d => Plume2d%next
    

       !------------------------------------------------------------
       ! put the un-dissolved plume concentration into PASV_LA3
       !------------------------------------------------------------
       backgrd_concnt = State_Chm%Species(id_PASV_LA3)%Conc(i_lon,i_lat,i_lev)

       backgrd_concnt = &
                   ( SUM( box_concnt_2D(:,:,1) )*V_grid_2D  &
                        + backgrd_concnt*grid_volumn ) / grid_volumn

       State_Chm%Species(id_PASV_LA3)%Conc(i_lon,i_lat,i_lev) = backgrd_concnt
        

800     CONTINUE

    ENDDO ! DO WHILE(ASSOCIATED(Plume2d))

900     CONTINUE



       !====================================================================
       ! sort the Plume1D linked list based on the grid volume, which will be
       ! used for the volume criterion.
       ! WARNING: the sort list after SortList() will change both Plume1d_head
       ! and Plume1d_tail, make sure the Plume1d_tail is also updated !!! 
       !====================================================================
       IF(Volume_Sort==1) Plume1d_head => SortList(Plume1d_head)


       !====================================================================
       ! begin 1D slab model only for sum volume of plume
       !====================================================================

       IF(.NOT.ASSOCIATED(Plume1d_head)) GOTO 500 ! no plume1d, skip loop


       IF(ASSOCIATED(Plume1d_prev)) NULLIFY(Plume1d_prev)
       Plume1d => Plume1d_head
       i_box = 0

       DO WHILE(ASSOCIATED(Plume1d)) ! how to delete last plume???

         i_box = i_box+1


         box_lat    = Plume1d%LAT
         box_lon    = Plume1d%LON
         box_lev    = Plume1d%LEV
         box_length = Plume1d%LENGTH
         box_alpha  = Plume1d%ALPHA

         box_label  = Plume1d%label
         box_life  = Plume1d%LIFE

         box_Ra    = Plume1d%RA
         box_Rb    = Plume1d%RB


         ! make sure the location is not out of range
         do while (box_lat > Y_edge(JJPAR+1))
            box_lat = Y_edge(JJPAR+1) &
                           - ( box_lat-Y_edge(JJPAR+1) )
         end do

         do while (box_lat < Y_edge(1))
            box_lat = Y_edge(1) + ( box_lat-Y_edge(1) )
         end do

         do while (box_lon > X_edge(IIPAR+1))
            box_lon = box_lon - 360.0
         end do

         do while (box_lon < X_edge(1))
            box_lon = box_lon + 360.0
         end do


         curr_lon      = box_lon
         curr_lat      = box_lat
         curr_pressure = box_lev      ! hPa



         i_lon = Find_iLonLat(curr_lon, DX, X_edge2)
         if(i_lon>IIPAR) i_lon=i_lon-IIPAR
         if(i_lon<1) i_lon=i_lon+IIPAR

         i_lat = Find_iLonLat(curr_lat, DY, Y_edge2)
         if(i_lat>JJPAR) i_lat=JJPAR
         if(i_lat<1) i_lat=1

         i_lev = Find_iPLev(curr_pressure,P_edge)
         if(i_lev>LLPAR) i_lev=LLPAR


         V_grid_1D        = box_Ra*box_Rb*box_length*1.0e+6_fp

         ! add plume volume to SumV_Plume
         i_cell = (i_lev-1)*IIPAR*JJPAR+(i_lat-1)*IIPAR+i_lon
         SumV_Plume(i_cell) = SumV_Plume(i_cell) + V_grid_1D*n_slab_max


         Plume1d_prev => Plume1d
         Plume1d      => Plume1d%next



     ENDDO ! DO WHILE(ASSOCIATED(Plume1d))

     
     !====================================================================
     ! Always make sure the Plume1d_tail refer to the last node of the list!
     !
     ! WARNING: the sort list after SortList() will change both Plume1d_head
     ! and Plume1d_tail, make sure the Plume1d_tail is also updated !!! 
     !==================================================================== 
     Plume1d_tail => Plume1d_prev
    




       !====================================================================
       ! begin 1D slab model
       !====================================================================
       IF(.NOT.ASSOCIATED(Plume1d_head)) GOTO 500 ! no plume1d, skip loop


       IF(ASSOCIATED(Plume1d_prev)) NULLIFY(Plume1d_prev)
       Plume1d => Plume1d_head
       i_box = 0

       DO WHILE(ASSOCIATED(Plume1d)) ! how to delete last plume???

         i_box = i_box+1


         box_lat    = Plume1d%LAT
         box_lon    = Plume1d%LON
         box_lev    = Plume1d%LEV
         box_length = Plume1d%LENGTH
         box_alpha  = Plume1d%ALPHA

         box_label  = Plume1d%label
         box_life  = Plume1d%LIFE

         box_Ra    = Plume1d%RA
         box_Rb    = Plume1d%RB
 
         box_concnt_1D = Plume1d%CONCNT1d




         ! make sure the location is not out of range
         do while (box_lat > Y_edge(JJPAR+1))
            box_lat = Y_edge(JJPAR+1) &
                           - ( box_lat-Y_edge(JJPAR+1) )
         end do

         do while (box_lat < Y_edge(1))
            box_lat = Y_edge(1) + ( box_lat-Y_edge(1) )
         end do

         do while (box_lon > X_edge(IIPAR+1))
            box_lon = box_lon - 360.0
         end do

         do while (box_lon < X_edge(1))
            box_lon = box_lon + 360.0
         end do


         curr_lon      = box_lon
         curr_lat      = box_lat
         curr_pressure = box_lev      ! hPa



         i_lon = Find_iLonLat(curr_lon, DX, X_edge2)
         if(i_lon>IIPAR) i_lon=i_lon-IIPAR
         if(i_lon<1) i_lon=i_lon+IIPAR

         i_lat = Find_iLonLat(curr_lat, DY, Y_edge2)
         if(i_lat>JJPAR) i_lat=JJPAR
         if(i_lat<1) i_lat=1

         i_lev = Find_iPLev(curr_pressure,P_edge)
         if(i_lev>LLPAR) i_lev=LLPAR






         ! run the fake 2nd order chemical reaction
         backgrd_concnt = State_Chm%Species(id_PASV_LA)%Conc(i_lon,i_lat,i_lev)

         box_concnt_1D(:,2) = box_concnt_1D(:,2) &
             + Dt* Kchem* (box_concnt_1D(:,1)+backgrd_concnt)**2 &
             - Dt* Kchem* backgrd_concnt**2



        !====================================================================
        ! calculate the wind shear along plume corss-section
        ! calculate the diffusivity in horizontal and vertical direction
        !====================================================================
        ! calculate the wind_s shear along pressure direction
        wind_s_shear = Wind_shear_s(u, v, P_BXHEIGHT, box_alpha, i_lon, &
                         i_lat, i_lev, curr_lon, curr_lat, curr_pressure)


        ! Calculate vertical eddy diffusivity (U.Schumann, 2012) :
        Cv = 0.2
        Omega_N = 0.1
        Ptemp_shear = Vertical_shear(Ptemp, P_BXHEIGHT, i_lon, i_lat, &
                              i_lev, curr_lon, curr_lat, curr_pressure)


        !--------------------------------------------------------------------
        ! interpolate potential temperature for Plume module:
        !--------------------------------------------------------------------
        if(abs(curr_lat)>Y_mid(JJPAR))then
          curr_Ptemp = Interplt_wind_RLL_polar(Ptemp, i_lon, i_lat, &
                                  i_lev, curr_lon, curr_lat, curr_pressure)
        else
          curr_Ptemp = Interplt_wind_RLL(Ptemp, i_lon, i_lat, i_lev, &
                                         curr_lon, curr_lat, curr_pressure)
        endif


        N_BV = SQRT(Ptemp_shear*g0/curr_Ptemp)
        IF(N_BV<=0.001) N_BV = 0.001

        ! diffusivity unit: [m2/s]
        eddy_v = Cv * Omega_N**2 / N_BV
        eddy_h = 10.0  

       ! ============================================
       ! Decide the time step based on CFL condition
       ! 2*k*Dt/Dr<1 (or Dr-2*k*Dt>0) for diffusion
       !==================================================================

       ! ignore the diffusivity in long radius direction
       eddy_B = eddy_v*SIN(abs(box_theta))
       eddy_B = eddy_B + eddy_h*COS(box_theta) ! b


       Dt2 = Dt

       IF((2*eddy_B*Dt2/(box_Rb**2))>0.5) THEN
          Dt2 = Dt*0.1
       ENDIF


       !--------------------------------------------------------------------
       ! Begin the loop to calculate the advection & diffusion in slab model
       !--------------------------------------------------------------------
!       DO i_species = 1,N_species
!       DO i_species = 1, 1
       do t1s=1,NINT(Dt/Dt2)

       !--------------------------------------------------------------------
       ! For deformation of cross-section caused by wind shear 
       ! (A.D.Naiman et al., 2010):
       ! This should be moved to the outer loop !!! shw
       !--------------------------------------------------------------------

       V_grid_1D        = box_Ra*box_Rb*box_length*1.0e+6_fp

       theta_previous   = box_theta
       box_theta = ATAN( TAN(box_theta) + wind_s_shear*Dt2 )

       ! make sure use box_theta or TAN(box_theta)  ??? 
       box_Ra = box_Ra  * (TAN(box_theta)**2+1)**0.5 &
                          / (TAN(theta_previous)**2+1)**0.5
       box_Rb = box_Rb * (TAN(box_theta)**2+1)**(-0.5) &
                         / (TAN(theta_previous)**2+1)**(-0.5)

       box_Rb = V_grid_1D/(box_Ra*box_length*1.0e+6_fp) !!! shw ???

       !--------------------------------------------------------------------
       ! For the concentration change caused by eddy diffusion:
       ! box_concnt(n_box,N_slab,N_species)
       !--------------------------------------------------------------------

       ! ignore the diffusivity in long radius direction
       eddy_B = eddy_v * SIN(abs(box_theta))
       eddy_B = eddy_B + eddy_h*COS(box_theta) ! b



       IF((2*eddy_B*Dt2/(box_Rb**2))>0.5) THEN

         300     CONTINUE

         !!! shw put this as a function

         !-------------------------------------------------------------------
         ! Decrease half of slab by combining 2 slab into 1 slab
         !-------------------------------------------------------------------
         DO i_species = 1, n_species, 1

         Cslab       = box_concnt_1D(1:n_slab_max,i_species)

         DO i_slab = 1, n_slab_50
           box_concnt_1D(i_slab+N_slab_25,i_species) = &
              (  Cslab(i_slab*2-1) + Cslab(i_slab*2) ) /2
         ENDDO

         !-------------------------------------------------------------------
         ! Meanwhile, add new slab into plume to compensate lost slab
         !-------------------------------------------------------------------
         box_concnt_1D(1:n_slab_25,i_species) = 0.0
         box_concnt_1D(n_slab_75+1:n_slab_max,i_species) = 0.0

         ENDDO ! DO i_species = 1, n_species, 1


         box_Rb = box_Rb*2


       ENDIF ! IF( 2*eddy_B*Dt2/(box_Rb(i_box)**2)>1 )

       IF(2*eddy_B*Dt2/(box_Rb**2)>0.5 )THEN
          GOTO 300
       ENDIF

       V_grid_1D        = box_Ra*box_Rb*box_length*1.0e+6_fp



       DO i_species=1,n_species,1

       C1d_prev       = box_concnt_1D(1:n_slab_max,i_species)

       !---------------------------------------------------------------------
       ! Calculate the diffusion in slab grids
       !---------------------------------------------------------------------
       Cslab       = box_concnt_1D(1:n_slab_max,i_species)

       Concnt1D_bdy(1)          = 0.0
       Concnt1D_bdy(n_slab_max2) = 0.0
       Concnt1D_bdy(2:n_slab_max2-1) = box_concnt_1D(1:n_slab_max,i_species)


       IF(2.0*eddy_B*Dt2/(box_Rb**2)>1) WRITE(6,*)'*** CFL ERROR ***'


       Concnt1D_bdy(2:n_slab_max2-1) = Concnt1D_bdy(2:n_slab_max2-1)     &
        + Dt2*eddy_B*(   Concnt1D_bdy(3:n_slab_max2)   &
                      -2*Concnt1D_bdy(2:n_slab_max2-1) &
                      +  Concnt1D_bdy(1:n_slab_max2-2) ) /(box_Rb**2)


       box_concnt_1D(1:n_slab_max,i_species) = Concnt1D_bdy(2:n_slab_max2-1)


       !====================================================================
       ! Update the concentration in the background grid cell
       ! after the interaction with 1D plume
       !====================================================================

       grid_volumn     = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ![cm3]

       mass_plume = V_grid_1D * SUM( C1d_prev(1:n_slab_max) )
       mass_plume_new = V_grid_1D * SUM(box_concnt_1D(1:n_slab_max,i_species))

       ! [molec]
       D_mass_plume = mass_plume_new - mass_plume

       i_advect = id_PASV_LA +i_species -1

       backgrd_concnt = State_Chm%Species(i_advect)%Conc(i_lon,i_lat,i_lev)

       backgrd_concnt = ( backgrd_concnt*grid_volumn &
                                              - D_mass_plume) / grid_volumn

       State_Chm%Species(i_advect)%Conc(i_lon,i_lat,i_lev) = backgrd_concnt

       ENDDO ! DO i_species=1,n_species,1



       !===================================================================
       ! Second judge:
       ! According to a 2-order process rate, if there is no big difference
       ! between PIG and Eulerian model, then dissolve the PIG into
       ! background grid cell.                                          shw
       ! test this once big time step
       !===================================================================

!       ! (1) Only consider the product producted by the tracer in Plume model
!       ! (2) Consider the interaction between Plume model and Eulerian model
!       ! Reference: Karamchandani et al., (2000)
!       !            Korsakissok and Mallet, (2010)

       backgrd_concnt = State_Chm%Species(id_PASV_LA)%Conc(i_lon,i_lat,i_lev)

       Prod_before = SUM( ( box_concnt_1D(1:n_slab_max,1)+backgrd_concnt )**2 &
                                               - backgrd_concnt**2 ) *V_grid_1D

       Eul_concnt = ( SUM(box_concnt_1D(1:n_slab_max,1)) *V_grid_1D &
                    + backgrd_concnt*grid_volumn ) / grid_volumn

       Prod_after = (Eul_concnt**2 - backgrd_concnt**2) *grid_volumn



       ! ------------------------------------------------------------------
       ! Volume criterion
       ! ------------------------------------------------------------------
       i_cell = (i_lev-1)*IIPAR*JJPAR +(i_lat-1)*IIPAR +i_lon


       ! Compare SumV_Plume with Eulerian grid cell volume
       ! determine the Eulerian grid cell containing too much plume segment


       ! After sorting the plume from largest to smallest by SortList()
       ! [cm3] always delete the largest plume first
       IF(SumV_Plume(i_cell) &
           >= Volume_percent*State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp)THEN

         Plume1d%Is_transfer = 1


         SumV_Plume(i_cell) = SumV_Plume(i_cell) - V_grid_1D*n_slab_max

       ENDIF


       ! ------------------------------------------------------------------
       ! 1. Product criterion
       ! 2. Volume criterion
       ! 3. Location critetion
       ! 4. Time criterion
       ! ------------------------------------------------------------------
       IF(      ABS(Prod_before-Prod_after)/Prod_after<Dissolve_critiria &
           .OR. Plume1d%Is_transfer==1 &
           .OR. box_lev>State_Met%TROPP(i_lon,i_lat) &
           .OR. box_life/3600/24>Critical_day  )THEN 


         IF( ABS(Prod_before-Prod_after)/Prod_after<Dissolve_critiria ) &
                                                      Type_transfer = 11
         IF( Plume1d%Is_transfer==1 )                 Type_transfer = 22
         IF( box_lev>State_Met%TROPP(i_lon,i_lat) )   Type_transfer = 33
         IF( box_life/3600/24>Critical_day )          Type_transfer = 44



         Num_dissolve = Num_dissolve+1

         WRITE(484,*) box_length/1000, Type_transfer, box_life, NINT(box_label)


         mass_la2 = mass_la2 + SUM(box_concnt_1D(1:n_slab_max,1)) * V_grid_1D

         DO i_species=1,n_species,1
         i_advect = id_PASV_LA +i_species -1

         backgrd_concnt = &
                  State_Chm%Species(i_advect)%Conc(i_lon,i_lat,i_lev)

         backgrd_concnt = &
           ( SUM(box_concnt_1D(1:n_slab_max,i_species)) * V_grid_1D   &
                + backgrd_concnt*grid_volumn ) / grid_volumn

         State_Chm%Species(i_advect)%Conc(i_lon,i_lat,i_lev) = backgrd_concnt
         ENDDO

         ! ----------------------------------
         ! delete the node for 1D plume:
         ! ----------------------------------

         IF(.NOT.ASSOCIATED(Plume1d%next) .and. &
                        .NOT.ASSOCIATED(Plume1d_prev))THEN

            WRITE(6,*)'                 '
            WRITE(6,*)'*** Only left the last Plume1d', i_box, box_label, Num_Plume1d
            WRITE(6,*)'                 '
            DEALLOCATE(Plume1d_head)

            i_box = i_box-1
            Num_Plume1d = Num_Plume1d - 1

            GOTO 500 ! Only 1 plume1d left, skip loop
         ENDIF

         IF(.NOT.ASSOCIATED(Plume1d%next))THEN ! delete last node
            Plume1d_tail%next => Plume1d_prev
            Plume1d_tail => Plume1d_tail%next

!            Plume1d => Plume1d_prev%next !!! ???
            IF(ASSOCIATED(Plume1d_prev%next)) NULLIFY(Plume1d_prev%next)
            DEALLOCATE(Plume1d%CONCNT1d)
            DEALLOCATE(Plume1d)
            i_box = i_box-1
            Num_Plume1d = Num_Plume1d - 1
            GOTO 500 ! delete last node, skip the loop
         ELSEIF(ASSOCIATED(Plume1d_prev))THEN ! delete node, not head/tail
            Plume1d => Plume1d_prev%next
            Plume1d_prev%next => Plume1d%next
            DEALLOCATE(Plume1d%CONCNT1d)
            DEALLOCATE(Plume1d)
            Plume1d => Plume1d_prev%next
         ELSE ! delete head node
            Plume1d => Plume1d_head
            Plume1d_head => Plume1d%next
            DEALLOCATE(Plume1d%CONCNT1d)
            DEALLOCATE(Plume1d)
            Plume1d => Plume1d_head
         ENDIF

         i_box = i_box-1
         Num_Plume1d = Num_Plume1d - 1


         GOTO 200 ! skip this box, go to next box

       ENDIF ! IF(ABS(Prod_before-Prod_after)/Prod_after<0.01)THEN



       enddo ! do t1s=1,NINT(Dt/Dt2)
!       ENDDO ! do i_Species = 1,nSpecies



       !===================================================================
       ! For filament structure:
       !
       ! 1. splitting for plume segment crtoss-section
       ! If slab length(Ra) is bigger than 2* horizontal resolution (2*DX),
       ! split slab into five smaller segment (Ra/5)
       !
       ! 2. splitting for plume segment length
       !===================================================================

       ! ----------------------------------------------
       ! 1. Splitting for plume segment cross-section
       ! ----------------------------------------------
       IF(N1_split==1) GOTO 111


       IF(.NOT.ALLOCATED(box_lon_new))THEN
         ALLOCATE(box_lon_new(N1_split-1))
         ALLOCATE(box_lat_new(N1_split-1))
       ENDIF


       IF( box_Ra > Split_length* DX*1.0e+5 ) THEN
!         WRITE(6,*)"SPLIT:", i_box, box_label

         !-----------------------------------------------------------------------
         ! set initial location for new added boxes from splitting
         !-----------------------------------------------------------------------
         ! Based on Great_circle distance
!         DlonN = 2.0 * 180/PI *ASIN( ABS( SIN(box_Ra(i_box)*COS(box_alpha)/N1_split /Re *0.5) &
!                                     / COS(box_lat(i_box)/180*PI)**2 ) )
        
         DlonN = ABS(box_Ra*SIN(box_alpha)/N1_split) / &
                        (2*PI *(Re*COS(box_lat/180*PI)) ) * 360

         DlatN = ABS(box_Ra*COS(box_alpha)/N1_split) / (2*PI*Re) * 360


         DO i = 1, (N1_split-1)/2, 1
           box_lon_new(i) = box_lon - i*DlonN
           box_lat_new(i) = box_lat - i*DlatN
         ENDDO

         DO i = 1, (N1_split-1)/2, 1
           box_lon_new((N1_split-1)/2+i) = box_lon + i*DlonN
           box_lat_new((N1_split-1)/2+i) = box_lat + i*DlatN
         ENDDO 

         box_Ra    = box_Ra/N1_split
!         box_extra = box_extra/N1_split

         do ii_box = 1, N1_split-1, 1

           ALLOCATE(Plume1d_new)
           Plume1d_new%label = Num_inject+1

           Plume1d_new%LON = box_lon_new(ii_box)
           Plume1d_new%LAT = box_lat_new(ii_box)
           Plume1d_new%LEV = box_lev

           Plume1d_new%LENGTH = box_length
           Plume1d_new%ALPHA  = box_alpha

           Plume1d_new%LIFE = box_life

           Plume1d_new%RA = box_Ra
           Plume1d_new%RB = box_Rb

           ALLOCATE(Plume1d_new%CONCNT1d(n_slab_max,n_species))
           Plume1d_new%CONCNT1d = box_concnt_1D

           ! add new splitting plume into the end of linked list
           NULLIFY(Plume1d_new%next)
           Plume1d_tail%next => Plume1d_new
           Plume1d_tail => Plume1d_tail%next

           Num_Plume1d = Num_Plume1d + 1
           Num_inject  = Num_inject + 1

         enddo ! do ii_box = 1, N1_split-1, 1


       ENDIF ! IF( box_Ra(i_box) > 2*DX )

111     CONTINUE


       ! -------------------------------------
       ! 2. Splitting for plume segment length
       ! -------------------------------------
       IF(N2_split==1) GOTO 222


       IF(.NOT.ALLOCATED(box_lon_new))THEN
         ALLOCATE(box_lon_new(N2_split-1))
         ALLOCATE(box_lat_new(N2_split-1))
       ENDIF


       IF( box_length > Split_length* DX*1.0e+5 ) THEN

         DlonN = ABS(box_length*COS(box_alpha)/N2_split) / &
                        (2*PI *(Re*COS(box_lat/180*PI)) ) * 360

         DlatN = ABS(box_length*SIN(box_alpha)/N2_split) / (2*PI*Re) * 360


         DO i = 1, (N2_split-1)/2, 1
           box_lon_new(i) = box_lon - i*DlonN
           box_lat_new(i) = box_lat - i*DlatN
         ENDDO

         DO i = 1, (N2_split-1)/2, 1
           box_lon_new((N2_split-1)/2+i) = box_lon + i*DlonN
           box_lat_new((N2_split-1)/2+i) = box_lat + i*DlatN
         ENDDO

         box_length    = box_length/N2_split
!         box_extra = box_extra/N2_split

         do ii_box = 1, N2_split-1, 1

           ALLOCATE(Plume1d_new)
           Plume1d_new%label = Num_inject+1

           Plume1d_new%LON = box_lon_new(ii_box)
           Plume1d_new%LAT = box_lat_new(ii_box)
           Plume1d_new%LEV = box_lev

           Plume1d_new%LENGTH = box_length
           Plume1d_new%ALPHA  = box_alpha

           Plume1d_new%LIFE = box_life

           Plume1d_new%RA = box_Ra
           Plume1d_new%RB = box_Rb

           ALLOCATE(Plume1d_new%CONCNT1d(n_slab_max,n_species))
           Plume1d_new%CONCNT1d = box_concnt_1D

           ! add new splitting plume into the end of linked list
           NULLIFY(Plume1d_new%next)
           Plume1d_tail%next => Plume1d_new
           Plume1d_tail => Plume1d_tail%next

           Num_Plume1d = Num_Plume1d + 1
           Num_inject  = Num_inject + 1

         enddo ! do ii_box = 1, N2_split-1, 1


       ENDIF ! IF( box_Ra(i_box) > 2*DX )

222     CONTINUE


       ! update the 1d plume for current node
       Plume1d%LON = box_lon
       Plume1d%LAT = box_lat
       Plume1d%LEV = box_lev
       Plume1d%LENGTH = box_length
       Plume1d%ALPHA = box_alpha

       Plume1d%LIFE = box_life

       Plume1d%RA = box_Ra
       Plume1d%RB = box_Rb

       Plume1d%CONCNT1d = box_concnt_1D


       !WRITE(6,*)"Plume run for 1D: assign Plume1d_prev"
       Plume1d_prev => Plume1d
       Plume1d      => Plume1d%next


       !---------------------------------------------------------
       ! put the un-dissolved plume concentration into PASV_LA3
       !---------------------------------------------------------
       backgrd_concnt = State_Chm%Species(id_PASV_LA3)%Conc(i_lon,i_lat,i_lev)

       backgrd_concnt = &
                     ( SUM(box_concnt_1D(1:n_slab_max,1)) * V_grid_1D  &
                            + backgrd_concnt*grid_volumn ) / grid_volumn

       State_Chm%Species(id_PASV_LA3)%Conc(i_lon,i_lat,i_lev) = backgrd_concnt


200     CONTINUE


     ENDDO ! DO WHILE(ASSOCIATED(Plume1d))


500     CONTINUE



    !===============================================================================
    ! The interaction between plume model and GEOS-Chem has finished
    ! Begin to calculate the entropy here
    !===============================================================================

    ! output once every day (24 hours)
    IF(mod(tt*NINT(Dt),30*24*60*60)==0 .AND. Calc_entropy==1)THEN 



    File996   = 'Plume_entropy_check.txt'
    OPEN( 996,      FILE=TRIM( File996   ), STATUS='OLD',  &
           POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )




    ! ------------------------------------------------------------------------------
    ! First, for 2-D plume model

    IF(ASSOCIATED(Plume2d_prev)) NULLIFY(Plume2d_prev)

    IF(ASSOCIATED(Plume2d_head))THEN

    Plume2d => Plume2d_head

    DO WHILE(ASSOCIATED(Plume2d))


      box_lon    = Plume2d%LON
      box_lat    = Plume2d%LAT
      box_lev    = Plume2d%LEV
      box_length = Plume2d%LENGTH
      box_alpha  = Plume2d%ALPHA

      box_label  = Plume2d%label
      box_life   = Plume2d%LIFE

      Pdx    = Plume2d%PDX
      Pdy    = Plume2d%PDY

      box_concnt_2D = Plume2d%CONCNT2d

      ! make sure the location is not out of range
      do while (box_lat > Y_edge(JJPAR+1))
         box_lat = Y_edge(JJPAR+1) &
                        - ( box_lat-Y_edge(JJPAR+1) )
      end do

      do while (box_lat < Y_edge(1))
         box_lat = Y_edge(1) + ( box_lat-Y_edge(1) )
      end do

      do while (box_lon > X_edge(IIPAR+1))
         box_lon = box_lon - 360.0
      end do

      do while (box_lon < X_edge(1))
         box_lon = box_lon + 360.0
      end do


      curr_lon      = box_lon
      curr_lat      = box_lat
      curr_pressure = box_lev      ! hPa



      i_lon = Find_iLonLat(curr_lon, DX, X_edge2)
      if(i_lon>IIPAR) i_lon=i_lon-IIPAR
      if(i_lon<1) i_lon=i_lon+IIPAR

      i_lat = Find_iLonLat(curr_lat, DY, Y_edge2)
      if(i_lat>JJPAR) i_lat=JJPAR
      if(i_lat<1) i_lat=1

      i_lev = Find_iPLev(curr_pressure,P_edge)

      V_grid_2D       = Pdx*Pdy*box_length*1.0e+6_fp


      DO i_y = 1,n_y_max,1
      DO i_x = 1,n_x_max,1

        Entropy_Concnt = box_concnt_2D(i_x,i_y,i_tracer) + &
                        State_Chm%Species(id_PASV_LA)%Conc(i_lon,i_lat,i_lev)

        tracer_mol = Entropy_Concnt * V_grid_2D / AVO
        air_mol    = V_grid_2D/1e+6_fp*State_Met%AIRDEN(i_lon,i_lat,i_lev)*1000/AIRMW
        mix_ratio  = tracer_mol/air_mol

        IF(mix_ratio>0) Entropy = Entropy - BOLTZ*tracer_mol*log(mix_ratio)

        i_cell = (i_lev-1)*IIPAR*JJPAR+(i_lat-1)*IIPAR+i_lon
        Entropy_V(i_cell) = Entropy_V(i_cell) - V_grid_2D


        ! Only output Day 1:
        IF(tt*NINT(Dt)==30*24*60*60 .AND. Calc_entropy==1)THEN
          WRITE(996,*)  2, tracer_mol, mix_ratio
        ENDIF



      ENDDO
      ENDDO


       Plume2d_prev => Plume2d
       Plume2d => Plume2d%next


    ENDDO ! DO WHILE(ASSOCIATED(Plume2d))

    ENDIF ! IF(ASSOCIATED(Plume2d_head))THEN



    ! ------------------------------------------------------------------------------
    ! Second, for 1-D plume model
    IF(ASSOCIATED(Plume1d_head))THEN


    IF(ASSOCIATED(Plume1d_prev)) NULLIFY(Plume1d_prev)
    Plume1d => Plume1d_head

    DO WHILE(ASSOCIATED(Plume1d)) ! how to delete last plume???


      box_lat    = Plume1d%LAT
      box_lon    = Plume1d%LON
      box_lev    = Plume1d%LEV
      box_length = Plume1d%LENGTH
      box_alpha  = Plume1d%ALPHA

      box_label  = Plume1d%label
      box_life  = Plume1d%LIFE

      box_Ra    = Plume1d%RA
      box_Rb    = Plume1d%RB

      box_concnt_1D = Plume1d%CONCNT1d

      ! make sure the location is not out of range
      do while (box_lat > Y_edge(JJPAR+1))
        box_lat = Y_edge(JJPAR+1) &
                      - ( box_lat-Y_edge(JJPAR+1) )
      end do

      do while (box_lat < Y_edge(1))
        box_lat = Y_edge(1) + ( box_lat-Y_edge(1) )
      end do

      do while (box_lon > X_edge(IIPAR+1))
        box_lon = box_lon - 360.0
      end do

      do while (box_lon < X_edge(1))
        box_lon = box_lon + 360.0
      end do


      curr_lon      = box_lon
      curr_lat      = box_lat
      curr_pressure = box_lev      ! hPa



      i_lon = Find_iLonLat(curr_lon, DX, X_edge2)
      if(i_lon>IIPAR) i_lon=i_lon-IIPAR
      if(i_lon<1) i_lon=i_lon+IIPAR

      i_lat = Find_iLonLat(curr_lat, DY, Y_edge2)
      if(i_lat>JJPAR) i_lat=JJPAR
      if(i_lat<1) i_lat=1

      i_lev = Find_iPLev(curr_pressure,P_edge)
      if(i_lev>LLPAR) i_lev=LLPAR


      ! Entropy
      V_grid_1D     = box_Ra*box_Rb*box_length*1.0e+6_fp

      DO i_slab = 1, n_slab_max, 1

        Entropy_Concnt = box_concnt_1D(i_slab,i_tracer) + &
                        State_Chm%Species(id_PASV_LA)%Conc(i_lon,i_lat,i_lev)

        tracer_mol = Entropy_Concnt * V_grid_1D / AVO
        air_mol    = V_grid_1D/1e+6_fp*State_Met%AIRDEN(i_lon,i_lat,i_lev)*1000/AIRMW
        mix_ratio  = tracer_mol/air_mol

        IF(mix_ratio>0) Entropy = Entropy - BOLTZ*tracer_mol*log(mix_ratio)


        i_cell = (i_lev-1)*IIPAR*JJPAR+(i_lat-1)*IIPAR+i_lon
        Entropy_V(i_cell) = Entropy_V(i_cell) - V_grid_1D



        ! Only output Day 1:
        IF(tt*NINT(Dt)==30*24*60*60 .AND. Calc_entropy==1)THEN
          WRITE(996,*)  1, tracer_mol, mix_ratio
        ENDIF




      ENDDO


      Plume1d_prev => Plume1d
      Plume1d      => Plume1d%next


    ENDDO ! DO WHILE(ASSOCIATED(Plume1d))

    ENDIF ! IF(ASSOCIATED(Plume1d_head))THEN



    ! ------------------------------------------------------------------------------
    ! Third, for modified GEOS-Chem grid
    DO i_lon = 1, IIPAR
    DO i_lat = 1, JJPAR
    DO i_lev = 1, LLPAR

       i_cell = (i_lev-1)*IIPAR*JJPAR +(i_lat-1)*IIPAR +i_lon

       Entropy_Concnt = State_Chm%Species(id_PASV_LA)%Conc(i_lon,i_lat,i_lev)

       tracer_mol = Entropy_Concnt*Entropy_V(i_cell)/AVO
       air_mol    = Entropy_V(i_cell)/1e+6_fp*State_Met%AIRDEN(i_lon,i_lat,i_lev)*1000/AIRMW
       mix_ratio  = tracer_mol/air_mol
              
       IF(mix_ratio>0) Entropy = Entropy - BOLTZ*tracer_mol*log(mix_ratio)


        ! Only output Day 1:
        IF(tt*NINT(Dt)==30*24*60*60 .AND. Calc_entropy==1)THEN
          WRITE(996,*) 0, tracer_mol, mix_ratio
        ENDIF


    ENDDO
    ENDDO
    ENDDO

!    WRITE(6,*)'Entropy3:', Entropy

    FileEntropy   = 'Plume_entropy.txt'
    OPEN( 487,      FILE=TRIM( FileEntropy   ), STATUS='OLD',  &
           POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

    WRITE(*,*) "Entropy:", Entropy
    WRITE(487,*) MONTH, DAY, Entropy

    CLOSE(487)



    CLOSE(996)


    ENDIF ! IF(mod(tt*NINT(Dt),24*60*60)==0 .AND. Calc_entropy==1)THEN

    !=======================================================================
    ! Convert species back to original units (ewl, 8/16/16)
    !=======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm,  State_Grid, State_Met, &
                            OrigUnit,  RC )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'plume_mod.F90' )
       RETURN
    ENDIF

    ! Everything is done, clean up pointers
    CLOSE(484)

999      CONTINUE

!    WRITE(6,*) 'Num test:  ', Num_Plume1d, Num_Plume2d, Num_inject, Num_dissolve
!    WRITE(6,*) 'Total mass:', mass_eu, mass_la, mass_la2


    ! if there is no plume left, stop using Lagrangian plume module
    !IF(Num_Plume1d+Num_Plume2d==0) Input_Opt%LagrangianModel_Activate = .False.


    ! Everything is done, clean up pointers


    IF(ALLOCATED(box_lon_new)) deallocate(box_lon_new)
    IF(ALLOCATED(box_lat_new)) deallocate(box_lat_new)

    IF(ASSOCIATED(u)) nullify(u)
    IF(ASSOCIATED(v)) nullify(v)
    IF(ASSOCIATED(omeg)) nullify(omeg)

    IF(ASSOCIATED(Ptemp)) nullify(Ptemp)
    IF(ASSOCIATED(P_BXHEIGHT)) nullify(P_BXHEIGHT)

    IF(ASSOCIATED(X_edge)) nullify(X_edge)
    IF(ASSOCIATED(Y_edge)) nullify(Y_edge)


    IF(ASSOCIATED(Plume2d_new)) nullify(Plume2d_new)
    IF(ASSOCIATED(PLume2d)) nullify(PLume2d)
    IF(ASSOCIATED(Plume2d_prev)) nullify(Plume2d_prev)
    IF(ASSOCIATED(Plume1d_new)) nullify(Plume1d_new)
    IF(ASSOCIATED(Plume1d)) nullify(Plume1d)
    IF(ASSOCIATED(Plume1d_prev)) nullify(Plume1d_prev)

  END SUBROUTINE plume_run

!--------------------------------------------------------------------
! functions to which grid cell (i,j,k) contains the box
!--------------------------------------------------------------------

  integer function Find_iLonLat(curr_xy,  Dxy,  XY_edge2)
    implicit none
    real(fp) :: curr_xy, Dxy, XY_edge2
    Find_iLonLat = INT( (curr_xy - (XY_edge2 - Dxy)) / Dxy )+1
    ! Notice the difference between INT(), FLOOR(), AINT()
    ! for lon: Xedge_Sec - DX = Xedge_first
    return
  end function


  integer function Find_iPLev(curr_pressure,P_edge)
    implicit none
    real(fp) :: curr_pressure
    real(fp), pointer :: P_edge(:)
    integer :: i_lon, i_lat
    integer :: locate(1)
    locate = MINLOC(abs( P_edge(:)-curr_pressure ))

    if(P_edge(locate(1))-curr_pressure >= 0 )then
       Find_iPLev = locate(1)
    else
       Find_iPLev = locate(1) - 1
    endif

    return
  end function


  !-------------------------------------------------------------------
  ! calculation the great-circle distance between two points on the earth
  ! surface
  !-------------------------------------------------------------------

  real(fp) function Distance_Circle(x1, y1, x2, y2)
    implicit none
    real(fp)     :: x1, y1, x2, y2      ! unit is degree
    real(fp)     :: xx1, yy1, xx2, yy2  ! unit is radian
    !real(fp) :: PI, Re

    xx1 = x1/180.0*PI
    yy1 = y1/180.0*PI
    xx2 = x2/180.0*PI
    yy2 = y2/180.0*PI

!     Distance_Circle = Re * 
!       ATAN( SQRT( (COS(y2)*SIN(x2-x1))**2 + &
!              ( COS(y1)*SIN(y2)-SIN(y1)*OCS(y2)COS((x2-x1)) )**2 ) &
!            / (SIN(y1)*SIN(y2)+COS(y1)*COS(y2)*COS(x2-x1)) )

    ! output distance is in unit of [m]
    Distance_Circle = Re * 2.0 * ASIN(SQRT( (SIN((yy1-yy2)*0.5))**2.0 &
                     +COS(yy1)*COS(yy2)*(SIN((xx1-xx2)*0.5))**2.0 ))
    return
  end function


  !-------------------------------------------------------------------  
  ! calculate the wind_s (inside a plume sross-section) shear along pressure
  ! direction
  !-------------------------------------------------------------------  

  real(fp) function Wind_shear_s(u, v, P_BXHEIGHT, plume_alpha, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    real(fp)          :: plume_alpha
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: u(:,:,:), v(:,:,:)
    real(fp), pointer :: P_BXHEIGHT(:,:,:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          :: distance(2,2), Weight(2,2)
    real(fp)          :: u_lonlat(2), v_lonlat(2)
    real(fp)          :: wind_s(2),  Delt_height

    ! first interpolate horizontally (Inverse Distance Weighting)

    if(curr_lon>=X_mid(i_lon))then
      init_lon = i_lon
    else
      init_lon = i_lon - 1
    endif

    if(curr_lat>=Y_mid(i_lat))then
      init_lat = i_lat
    else
      init_lat = i_lat - 1
    endif

    ! For pressure level, P_mid(1) is about surface pressure
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif


    do i = 1,2
    do j = 1,2

      ii = i + init_lon - 1
      jj = j + init_lat - 1

      ! For some special circumstance:
      if(ii==0)then
        distance(i,j) = &
           Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), Y_mid(jj))
      else if(ii==(IIPAR+1))then
        distance(i,j) = &
           Distance_Circle(curr_lon, curr_lat, X_mid(ii-IIPAR), Y_mid(jj))
      else
        distance(i,j) = &
           Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
      endif

    enddo
    enddo

    do i = 1,2
      do j = 1,2
        Weight(i,j) = 1.0/distance(i,j) / sum( 1.0/distance(:,:) )
      enddo
    enddo


    do k = 1,2

      kk = k + init_lev - 1

      if(init_lon==0)then
          u_lonlat(k) =  Weight(1,1) * u(IIPAR,init_lat,kk) &
                       + Weight(1,2) * u(IIPAR,init_lat+1,kk) &
                       + Weight(2,1) * u(init_lon+1,init_lat,kk) &
                       + Weight(2,2) * u(init_lon+1,init_lat+1,kk)
          v_lonlat(k) =  Weight(1,1) * v(IIPAR,init_lat,kk) &
                       + Weight(1,2) * v(IIPAR,init_lat+1,kk) &
                       + Weight(2,1) * v(init_lon+1,init_lat,kk) &
                       + Weight(2,2) * v(init_lon+1,init_lat+1,kk)
      else if(init_lon==IIPAR)then
          u_lonlat(k) =  Weight(1,1) * u(init_lon,init_lat,kk) &
                       + Weight(1,2) * u(init_lon,init_lat+1,kk) &
                       + Weight(2,1) * u(1,init_lat,kk)   &
                       + Weight(2,2) * u(1,init_lat+1,kk)
          v_lonlat(k) =  Weight(1,1) * v(init_lon,init_lat,kk) &
                       + Weight(1,2) * v(init_lon,init_lat+1,kk) &
                       + Weight(2,1) * v(1,init_lat,kk)   &
                       + Weight(2,2) * v(1,init_lat+1,kk)
      else
          u_lonlat(k) =  Weight(1,1) * u(init_lon,init_lat,kk) &
                       + Weight(1,2) * u(init_lon,init_lat+1,kk) &
                       + Weight(2,1) * u(init_lon+1,init_lat,kk) &
                       + Weight(2,2) * u(init_lon+1,init_lat+1,kk)
          v_lonlat(k) =  Weight(1,1) * v(init_lon,init_lat,kk) &
                       + Weight(1,2) * v(init_lon,init_lat+1,kk) &
                       + Weight(2,1) * v(init_lon+1,init_lat,kk) &
                       + Weight(2,2) * v(init_lon+1,init_lat+1,kk)
      endif

      wind_s(k) = u_lonlat(k)*COS(plume_alpha-0.5*PI) + v_lonlat(k)*COS(plume_alpha-PI)

    enddo


    ! second vertical shear of wind_s

    ! This code should be changed !!!
   ! Because it is the pressure center in [hPa] instead of height center in
    ! [m]
    ! Delt_height    = 0.5 * ( P_BXHEIGHT(init_lon,init_lat,init_lev) +
    ! P_BXHEIGHT(init_lon,init_lat,init_lev+1) )
    if(init_lon==0) init_lon=IIPAR
    if(init_lon==IIPAR+1) init_lon=1

    Delt_height = Pa2meter( P_BXHEIGHT(IIPAR,init_lat,init_lev),    &
                          P_edge(init_lev), P_edge(init_lev+1), 1 ) &   
                + Pa2meter( P_BXHEIGHT(IIPAR,init_lat,init_lev+1),   &
                          P_edge(init_lev), P_edge(init_lev+1), 0 )


    ! find the z height of each pressure level in GEOS-Chem

    Wind_shear_s = ( wind_s(2) - wind_s(1) ) / Delt_height

    return
  end function

!-------------------------------------------------------------------
! transform the pressure level [Pa] to the height level [m]
!-------------------------------------------------------------------

  real(fp) function Pa2meter(Box_height, P1, P2, Judge)
  ! Judge: 0 for bottom, 1 for top
  ! 
    implicit none
    real(fp) :: Box_height, P1, P2
    integer  :: Judge

      if(Judge==1)then
        ! calculate the height of the top half of the grid box
        Pa2meter = Box_height * ( DLOG(P2) - DLOG(0.5*(P2+P1)) ) &
                        / ( DLOG(P2) - DLOG(P1) )
      else
        ! calculate the height of the bottom half of the grid box
        Pa2meter = Box_height * ( DLOG(0.5*(P2+P1)) - DLOG(P1) ) &
                        / ( DLOG(P2) - DLOG(P1) )
      endif

    return
  end function


!-------------------------------------------------------------------
! Calculate eddy diffusivity in stratisphere. (U.Schumann, 2012)


  !-------------------------------------------------------------------  
  ! calculate the vertical shear for calculating eddy difussivity

  real(fp) function Vertical_shear(var, P_BXHEIGHT, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: var(:,:,:)
    real(fp), pointer :: P_BXHEIGHT(:,:,:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          :: distance(2,2), Weight(2,2)
    real(fp)          :: var_lonlat(2)
    real(fp)          :: Delt_height

    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southeast of the particle or under
    ! the particle
    if(curr_lon>=X_mid(i_lon))then
      init_lon = i_lon
    else
      init_lon = i_lon - 1
    endif

    if(curr_lat>=Y_mid(i_lat))then
      init_lat = i_lat
    else
      init_lat = i_lat - 1
    endif

    ! For pressure level, P_mid(1) is about surface pressure
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif


    do i = 1,2
    do j = 1,2

      ii = i + init_lon - 1
      jj = j + init_lat - 1

      ! For some special circumstance:
      if(ii==0)then
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), Y_mid(jj))
      else if(ii==(IIPAR+1))then
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii-IIPAR), Y_mid(jj))
      else
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
      endif

    enddo
    enddo

    do i = 1,2
      do j = 1,2
        Weight(i,j) = 1.0/distance(i,j) / sum( 1.0/distance(:,:) )
      enddo
    enddo


    do k = 1,2
      kk            = k + init_lev - 1
      
      IF(init_lon==0)THEN
        var_lonlat(k) =  Weight(1,1) *var(IIPAR,i_lat,kk)   &
                       + Weight(1,2) *var(IIPAR,i_lat+1,kk)   &
                       + Weight(2,1) *var(1,i_lat,kk) &
                       + Weight(2,2) *var(1,i_lat+1,kk)
      ELSE IF(init_lon==IIPAR)THEN
        var_lonlat(k) =  Weight(1,1) *var(IIPAR,i_lat,kk)   &
                       + Weight(1,2) *var(IIPAR,i_lat+1,kk)   &
                       + Weight(2,1) *var(1,i_lat,kk) &
                       + Weight(2,2) *var(1,i_lat+1,kk)
      ELSE
        var_lonlat(k) =  Weight(1,1) *var(init_lon,i_lat,kk)   &
                       + Weight(1,2) *var(init_lon,i_lat+1,kk)   &
                       + Weight(2,1) *var(init_lon+1,i_lat,kk) &
                       + Weight(2,2) *var(init_lon+1,i_lat+1,kk)
      ENDIF
    enddo


    ! second vertical shear of wind
    if(init_lon==0)then
      Delt_height = Pa2meter( P_BXHEIGHT(IIPAR,init_lat,init_lev),    &
                            P_edge(init_lev), P_edge(init_lev+1), 1 ) &
                 + Pa2meter( P_BXHEIGHT(IIPAR,init_lat,init_lev+1),   &
                            P_edge(init_lev), P_edge(init_lev+1), 0 )
    else if(init_lon==IIPAR)then
      Delt_height = Pa2meter( P_BXHEIGHT(1,init_lat,init_lev),        &
                            P_edge(init_lev), P_edge(init_lev+1), 1 ) &
                 + Pa2meter( P_BXHEIGHT(1,init_lat,init_lev+1),       &
                            P_edge(init_lev), P_edge(init_lev+1), 0 )
    else
      Delt_height = Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev), &
                            P_edge(init_lev), P_edge(init_lev+1), 1 ) &   
                 + Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev+1),&
                            P_edge(init_lev), P_edge(init_lev+1), 0 )
    endif
    ! find the z height of each pressure level in GEOS-Chem

    Vertical_shear = ( var_lonlat(2) - var_lonlat(1) ) / Delt_height

    return
  end function


!===================================================================
!
!===================================================================
  REAL(fp) FUNCTION Get_XYscale(concnt1_2D, Pdx, Pdy, frac, axis)

    IMPLICIT NONE

    REAL(fp) :: Pdx, Pdy, concnt1_2D(n_x_max, n_y_max), frac
    INTEGER  :: axis

    INTEGER  :: i, j
    REAL(fp) :: temp, C_sum, mass_total

    REAL(fp) :: D_len
    INTEGER  :: N_max, N_frac

    REAL(fp), allocatable :: concnt1_2D_sum(:)
    integer :: alloc_stat


      IF(axis==2)THEN ! sum y
        N_max = n_x_max
        D_len = Pdx
      ELSE IF(axis==1)THEN ! sum x
        N_max = n_y_max
        D_len = Pdy
      ENDIF

      ALLOCATE(concnt1_2D_sum(N_max), stat=alloc_stat)
      IF(alloc_stat/=0) WRITE(6,*)'ERROR 6:', alloc_stat

      concnt1_2D_sum = SUM(concnt1_2D, DIM=axis)


      ! sort the array from high value to low:
      DO i = N_max-1, 1, -1
      DO j = 1, i
        IF(concnt1_2D_sum(j)<concnt1_2D_sum(j+1))THEN
          temp = concnt1_2D_sum(j)
          concnt1_2D_sum(j) = concnt1_2D_sum(j+1)
          concnt1_2D_sum(j+1) = temp
        ENDIF
      ENDDO
      ENDDO


      mass_total = SUM(concnt1_2D_sum)
      C_sum = 0.0

      ! find the XYscale containing frac of total mass
      DO i = 1, N_max, 1
        C_sum = C_sum + concnt1_2D_sum(i)
        if(C_sum>mass_total*frac)then
          N_frac = i
          EXIT
        endif
      ENDDO


      Get_XYscale = N_frac * D_len

      DEALLOCATE(concnt1_2D_sum)

    RETURN

  END FUNCTION



!======================================================================
! Conservatice interpolation
!======================================================================

  SUBROUTINE Slab_init_conservative(Pdx, Pdy, theta1, Pc_2D, Height1, &
						box_Ra, box_Rb, Cslab, n_slab_max, n_x_max, n_y_max)

    IMPLICIT NONE

    REAL(fp)    :: Pdx, Pdy, theta1, Height1
    REAL(fp), INTENT(INOUT)  :: box_Ra, box_Rb
    REAL(fp), INTENT(INOUT)  :: Cslab(n_slab_max)
    REAL(fp)    :: Pc_2D(n_x_max,n_y_max) !, Ec_2D(n_x_max,n_y_max)


    !INTEGER, parameter     :: Nb = n_slab_max +2
    !INTEGER, parameter     :: Na = (INT(n_x_max/4)+1)*4 +2
    INTEGER, INTENT(IN)     :: n_slab_max, n_x_max, n_y_max            
    INTEGER                 :: Nb, Na  
    
    INTEGER     :: Nb_mid, Na_mid

    REAL(fp), ALLOCATABLE    :: X2d(:,:), Y2d(:,:), C2d(:,:) !, Extra_C2d(Na,Nb)

    REAL(fp)    :: LenB, LenA
    REAL(fp)    :: Adx, Ady, Bdx, Bdy
    REAL(fp)    :: Prod, M, Lb, La

    REAL(fp), ALLOCATABLE    :: C_slab(:) !, Extra_slab(Nb)

    REAL(fp)    :: start, finish

    INTEGER     :: i, j


    REAL(fp)	:: X2d_edge(4), Y2d_edge(4)
    REAL(fp)	:: Area_target
    REAL(fp)	:: X_max, X_min, Y_max, Y_min
    REAL(fp)	:: Xx(n_x_max), Yy(n_y_max)
    REAL(fp)	:: source_Xedge(4), source_Yedge(4)
    REAL(fp)	:: source_X1, source_X2, source_Y1, source_Y2
    REAL(fp)	:: target_X1, target_X2, target_Y1, target_Y2
    REAL(fp)	:: X_intersect, Y_intersect
    REAL(fp)	:: Points_X(8), Points_Y(8)
    REAL(fp)	:: mass_total, area_total, area_overlap

    real(fp), dimension(:), allocatable   :: Pts_X, Pts_Y

    INTEGER	:: ix_max, ix_min, iy_max, iy_min
    INTEGER	:: ix, iy
    INTEGER	:: Nx_dim, Ny_dim
    INTEGER	:: ip, ip_S, i_point
    INTEGER	:: N_diff

    LOGICAL	:: Is_inside

      Nb = n_slab_max +2
      Na = (INT(n_x_max/4)+1)*4 +2
      ALLOCATE(X2d(Na,Nb))
      ALLOCATE(Y2d(Na,Nb))
      ALLOCATE(C2d(Na,Nb))
      ALLOCATE(C_slab(Nb))
      ! define the coordinate system
      DO i=1, n_x_max
        Xx(i) = Pdx*(i-n_x_mid)
      ENDDO
      DO j=1, n_y_max
        Yy(j) = Pdy*(j-n_y_mid)
      ENDDO



      Nb_mid = INT(Nb/2)
      Na_mid = INT(Na/2)


      LenB = Pdy ! 8.0 *Height1 / Nb /SIN(theta1)
      LenA = Pdx ! 1.3 *Height1 *TAN(theta1) /Na /SIN(theta1)


      ! interval in long radius
      Adx = LenA*SIN(theta1)
      Ady = LenA*COS(theta1)

      ! interval in short radius
      Bdy = LenB*SIN(theta1)
      Bdx = LenB*COS(theta1)


      ! find the location of 1D grid in 2D XY grids
      DO i=1,Na,1
        X2d(i,Nb_mid) = -Adx*Na_mid + Adx*(i-0.5)
        Y2d(i,Nb_mid) = -Ady*Na_mid + Ady*(i-0.5)
      ENDDO


      X2d(:,Nb_mid) = X2d(:,Nb_mid) + 0.5*Bdx
      Y2d(:,Nb_mid) = Y2d(:,Nb_mid) - 0.5*Bdy

      DO j=Nb_mid+1, Nb, 1
        X2d(:,j) = X2d(:,j-1) - Bdx
        Y2d(:,j) = Y2d(:,j-1) + Bdy
      ENDDO


      DO j=Nb_mid-1, 1, -1
        X2d(:,j) = X2d(:,j+1) + Bdx
        Y2d(:,j) = Y2d(:,j+1) - Bdy
      ENDDO



      ! begin the conservative interpolation
      DO i=2,Na-1,1
      DO j=2,Nb-1,1


	! four edge points of selected target grid cell
	X2d_edge(1) = ( X2d(i,j) + X2d(i-1,j-1) )/2
	Y2d_edge(1) = ( Y2d(i,j) + Y2d(i-1,j-1) )/2

        X2d_edge(2) = ( X2d(i,j) + X2d(i+1,j-1) )/2
        Y2d_edge(2) = ( Y2d(i,j) + Y2d(i+1,j-1) )/2

        X2d_edge(3) = ( X2d(i,j) + X2d(i+1,j+1) )/2
        Y2d_edge(3) = ( Y2d(i,j) + Y2d(i+1,j+1) )/2

        X2d_edge(4) = ( X2d(i,j) + X2d(i-1,j+1) )/2
        Y2d_edge(4) = ( Y2d(i,j) + Y2d(i-1,j+1) )/2


	! Check whether target grid cell area is equal to Dx*Dy
	Area_target = 							&
	    Heron_Formula( X2d_edge(1), Y2d_edge(1), 			&
		X2d_edge(3), Y2d_edge(3), X2d_edge(2), Y2d_edge(2) ) 	&
	  + Heron_Formula( X2d_edge(1), Y2d_edge(1), 			&
                X2d_edge(3), Y2d_edge(3), X2d_edge(4), Y2d_edge(4) )


	! Find all the source grid cells that may overlap with the target grid
	! cell, based on the four edge points of the target grid cell
	X_max = MAXVAL( X2d_edge(:) )
	X_min = MINVAL( X2d_edge(:) )
  
	Y_max = MAXVAL( Y2d_edge(:) )
	Y_min = MINVAL( Y2d_edge(:) )


	ix_max = CEILING( (X_max - ( Xx(1)-Pdx/2 )) /Pdx )
	ix_min = CEILING( (X_min - ( Xx(1)-Pdx/2 )) /Pdx )
 
        iy_max = CEILING( (Y_max - ( Yy(1)-Pdy/2 )) /Pdy )
        iy_min = CEILING( (Y_min - ( Yy(1)-Pdy/2 )) /Pdy )


	! for the target grid cell located in the edge of the source grid domain
	if(ix_min<1) ix_min = 1
	if(iy_min<1) iy_min = 1

	if(ix_max>n_x_max) ix_max = n_x_max
	if(iy_max>n_y_max) iy_max = n_y_max

	
	! when the target grid cell is totally outside the source grid domain
	if(ix_min>n_x_max)then
	  ix_min = n_x_max
	  ix_max = n_x_max
	endif

	if(iy_min>n_y_max)then
	  iy_min = n_y_max
	  iy_max = n_y_max
	endif

	if(ix_max<1)then
	  ix_max = 1
	  ix_min = 1
	endif

	if(iy_max<1)then
	  iy_max=1
	  iy_min=1
	endif

	! loop for all selected source grid cell, calculate the overlap area
	! between the selected source grid cell and the target grid cell
        Nx_dim = ix_max-ix_min+1
	Ny_dim = iy_max-iy_min+1


	area_overlap 	 = 0.0

	mass_total	 = 0.0
	area_total	 = 0.0


	do ix = ix_min, ix_max, 1
	do iy = iy_min, iy_max, 1
	 
 
          area_overlap = 0.0

          Points_X(:) = 0.0
          Points_Y(:) = 0.0


	  i_point = 0
	  N_diff  = 0

	  ! four edge points of the selected source grid cell
	  source_Xedge(1) = Xx(ix) - Pdx/2
          source_Yedge(1) = Yy(iy) - Pdy/2
	  
          source_Xedge(2) = Xx(ix) - Pdx/2
          source_Yedge(2) = Yy(iy) + Pdy/2

          source_Xedge(3) = Xx(ix) + Pdx/2
          source_Yedge(3) = Yy(iy) + Pdy/2

          source_Xedge(4) = Xx(ix) + Pdx/2
          source_Yedge(4) = Yy(iy) - Pdy/2
	  

	  ! (1) check whether target grid cell edge points are inside the source
	  ! grid cell
	  do ip=1,4,1
	  if(      ABS(X2d_edge(ip)-Xx(ix))<Pdx/2 &
	     .and. ABS(Y2d_edge(ip)-Yy(iy))<Pdy/2 )then
		i_point = i_point+1
		Points_X(i_point) = X2d_edge(ip)
		Points_Y(i_point) = Y2d_edge(ip)
	  endif
	  enddo


	  ! (2) Check whether source grid cell edge points are inside the target
	  ! grid cell
	  do ip_S=1,4,1
	    Is_inside = In_box(source_Xedge(ip_S), source_Yedge(ip_S), &
							X2d_edge, Y2d_edge)
	    if(Is_inside)then
		i_point = i_point+1
                Points_X(i_point) = source_Xedge(ip_S)
                Points_Y(i_point) = source_Yedge(ip_S)
	    endif
	  enddo


	  ! (3) Check the intersection point between source grid cell side and
	  ! target grid cell side
	  do ip=1,4,1
	    if(ip==4)then
		target_X1 = X2d_edge(ip)
		target_Y1 = Y2d_edge(ip)
		target_X2 = X2d_edge(1)
		target_Y2 = Y2d_edge(1)
	    else
                target_X1 = X2d_edge(ip)
                target_Y1 = Y2d_edge(ip)
                target_X2 = X2d_edge(ip+1)
                target_Y2 = Y2d_edge(ip+1)
	    endif

	    do ip_S=1,4,1
	      if(ip_S==4)then
                  source_X1 = source_Xedge(ip_S)
                  source_Y1 = source_Yedge(ip_S)
                  source_X2 = source_Xedge(1)
                  source_Y2 = source_Yedge(1)
	      else
		  source_X1 = source_Xedge(ip_S)
		  source_Y1 = source_Yedge(ip_S)
		  source_X2 = source_Xedge(ip_S+1)
		  source_Y2 = source_Yedge(ip_S+1)
	      endif



	      ! calculate the intersect point for the two given lines.
	      !
	      ! In this special case: (1) the side of the source grid cell must 
	      ! be parallel to either x or y axis; (2) the side slope of the 
	      ! target grid cell must not be parallel to both x and y axis.

	      if(source_X1==source_X2)then
		X_intersect = source_X1
		Y_intersect = (target_Y2-target_Y1) / (target_X2-target_X1) &
					* (X_intersect-target_X1) + target_Y1
	      endif


	      if(source_Y1==source_Y2)then
		Y_intersect = source_Y1
		X_intersect = (target_X2-target_X1) / (target_Y2-target_Y1) &
					* (Y_intersect-target_Y1) + target_X1
	      endif


	      ! check whether this intersect point is located in the both line
	      ! segments
	      if( (X_intersect-target_X1)*(X_intersect-target_X2)<=0 .and. &
		  (Y_intersect-target_Y1)*(Y_intersect-target_Y2)<=0 .and. &
		  (X_intersect-source_X1)*(X_intersect-source_X2)<=0 .and. &
		  (Y_intersect-source_Y1)*(Y_intersect-source_Y2)<=0 )then

                i_point = i_point+1
                Points_X(i_point) = X_intersect
                Points_Y(i_point) = Y_intersect

	      endif

            enddo ! do ip_S=1,4,1
	  enddo ! do ip=1,4,1



          ! tempprary put all the last elements to be same as the i_point
	  ! all the repeated elements will be deleted later
          if(i_point==0)then
            Points_X(:) = 0
            Points_Y(:) = 0
	  else
	    Points_X(i_point:8) = Points_X(i_point)
	    Points_Y(i_point:8) = Points_Y(i_point)
	  endif



	  ! Sort all the intersect points in one direction (clock/anti-clock)
	  CALL Sort_points(Points_X, Points_Y, N_diff)


	  if(N_diff>=3)then
	    ! delete all the repeate points

	    allocate( Pts_X(N_diff) )
	    allocate( Pts_Y(N_diff) )
	  
	    Pts_X(:) = Points_X(1:N_diff)
	    Pts_Y(:) = Points_Y(1:N_diff)

	    ! Calculate the total overlapping area based on the final sorted
	    ! points
	    do ip=2,N_diff-1,1
	      area_overlap = area_overlap + 		&
		Heron_Formula(Pts_X(1), Pts_Y(1), Pts_X(ip), Pts_Y(ip), &
					      Pts_X(ip+1), Pts_Y(ip+1) )
	    enddo


	    area_total = area_total + area_overlap
	    mass_total = mass_total + area_overlap * Pc_2d(ix,iy)

            deallocate(Pts_X)
            deallocate(Pts_Y)


          endif ! if(N_diff>=3)then


	enddo ! iy = iy_min, iy_max, 1
	enddo ! ix = ix_min, ix_max, 1



	if(area_total<=0)then
	  C2d(i, j) = 0.0
	else
	  C2d(i, j) = mass_total/area_total
	endif


      ENDDO ! DO i=1,Na,1
      ENDDO ! DO j=1,Nb,1


      ! define the length/width of slab based on 95% total mass
      Lb = LenB
      La = Height1 *TAN(theta1)


      ! Nb = n_slab_max +2
      ! C2d(Na,Nb)

      DO j=1,n_slab_max,1
        C_slab(j)     = SUM(C2d(:,j+1)) *(LenA*LenB)/ (La*Lb)
      ENDDO


      ! Return:
      box_Ra = La
      box_Rb = Lb
      Cslab(1:n_slab_max) = C_slab(1:n_slab_max)



  END SUBROUTINE Slab_init_conservative




  ! Calculate triangle's area based on its three sides
  REAL(fp) FUNCTION Heron_Formula(x1, y1, x2, y2, x3, y3)
    IMPLICIT NONE

    real(fp)	:: x1, y1, x2, y2, x3, y3
    real(fp)	:: L1, L2, L3, Ss
    real(fp)	:: Area

    L1 = SQRT( (x1-x2)**2 + (y1-y2)**2 )
    L2 = SQRT( (x2-x3)**2 + (y2-y3)**2 )
    L3 = SQRT( (x1-x3)**2 + (y1-y3)**2 )

    Ss = (L1+L2+L3)/2.0


    if( (Ss-L1)*(Ss-L2)*(Ss-L3)<=0 )then
!	WRITE(6,*)"*** ERROR in Heron_Formula ***"
!	WRITE(6,*)"*** three points may in one line ***"
!	WRITE(6,*)x1, x2, x3
!	WRITE(6,*)y1, y2, y3
	Area = 0.0
    else
	Area = SQRT( Ss*(Ss-L1)*(Ss-L2)*(Ss-L3) )
    endif

    Heron_Formula = Area

    RETURN
  END FUNCTION




  ! check whether source point is always in the same side of four line of the
  ! target grid cell, if so, the source point is inside the target grid cell
  LOGICAL FUNCTION In_box(X_source, Y_source, Xs_target, Ys_target)
    IMPLICIT NONE

    real(fp)	:: X_source, Y_source
    real(fp)	:: Xs_target(4), Ys_target(4)
    real(fp)	:: Vec_x1, Vec_x2, Vec_y1, Vec_y2
    real(fp)	:: dot, det, angle

    integer	:: i, signal

    signal = 0

    ! for the first 3 lines of the target grid cell
    do i=1,3,1
	Vec_x1 = X_source - Xs_target(i)
	Vec_y1 = Y_source - Ys_target(i)

	Vec_x2 = Xs_target(i+1) - Xs_target(i)
	Vec_y2 = Ys_target(i+1) - Ys_target(i)

	dot = Vec_x1*Vec_x2 + Vec_y1*Vec_y2
	det = Vec_x1*Vec_y2 - Vec_y1*Vec_x2

	angle  = ATAN2(det,dot)
	signal = signal + NINT(angle/ABS(angle))
    enddo

    ! for the last line of the target grid cell
    Vec_x1 = X_source - Xs_target(4)
    Vec_y1 = Y_source - Ys_target(4)

    Vec_x2 = Xs_target(4) - Xs_target(4)
    Vec_y2 = Y_source - Ys_target(4)

    dot = Vec_x1*Vec_x2 + Vec_y1*Vec_y2
    det = Vec_x1*Vec_y2 - Vec_y1*Vec_x2

    angle  = ATAN2(det,dot)
    signal = signal + NINT(angle/ABS(angle))


    ! check whether source point is always in the same side of four line of the
    ! target grid cell, if so, the source point is inside the target grid cell
    if(signal==4 .or. signal==-4)then
	In_box = 1
    else
	In_box = 0
    endif


    RETURN
  END FUNCTION



  SUBROUTINE Sort_points(Points_X, Points_Y, N_diff)
    IMPLICIT NONE

    real(fp), intent(inout)	:: Points_X(8),  Points_Y(8)
    real(fp), dimension(8)	:: Xs_ascend,    Ys_ascend
    real(fp)			:: Angles(8), Angles_ascend(8)
    real(fp)			:: Vec_X1, Vec_Y1, Vec_X2, Vec_Y2
    real(fp)			:: dot, det
    real(fp)			:: X_center, Y_center


    integer, intent(inout)	:: N_diff
    integer			:: ip
    integer			:: idx(1)

    logical, dimension(8)	:: mask
    
    X_center = SUM(Points_X)/8
    Y_center = SUM(Points_Y)/8



    Vec_X1 = Points_X(1) - X_center
    Vec_Y1 = Points_Y(1) - Y_center

    Angles(1) = 0.0

    do ip=2,8,1
	Vec_X2 = Points_X(ip) - X_center
	Vec_Y2 = Points_Y(ip) - Y_center

	dot = Vec_X1*Vec_X2 + Vec_Y1*Vec_Y2	! dot product
	det = Vec_X1*Vec_Y2 - Vec_Y1*Vec_X2	! determinant

	Angles(ip) = ATAN2(det, dot)
	if(Angles(ip)<0) Angles(ip) = 2*PI + Angles(ip)
    enddo



    mask(:) = .TRUE.

    DO ip = 1, 8, 1
	idx = MINLOC(Angles,mask)
	Angles_ascend(ip) = Angles( idx(1) )
	Xs_ascend(ip) 	  = Points_X( idx(1) )
	Ys_ascend(ip) 	  = Points_Y( idx(1) )
	mask(idx(1)) 	  = .FALSE.
    END DO



    ! select all the identical elements
    N_diff = 1
    Points_X(N_diff) = Xs_ascend(1)
    Points_Y(N_diff) = Ys_ascend(1)

    do ip = 1, 8-1, 1
    if(    Xs_ascend(ip)/=Xs_ascend(ip+1) &
       .or.Ys_ascend(ip)/=Ys_ascend(ip+1))then

	N_diff = N_diff + 1
	Points_X(N_diff) = Xs_ascend(ip+1)
	Points_Y(N_diff) = Ys_ascend(ip+1)

    endif
    enddo


  END SUBROUTINE Sort_points



!======================================================================
! Bilinear interpolation
!======================================================================

  SUBROUTINE Slab_init_bilinear(Pdx, Pdy, theta1, Pc_2D, Height1, &
						box_Ra, box_Rb, Cslab, n_slab_max)
 
    IMPLICIT NONE

    REAL(fp)    :: Pdx, Pdy, theta1, Height1
    REAL(fp), INTENT(INOUT)  :: box_Ra, box_Rb
    REAL(fp), INTENT(INOUT)  :: Cslab(n_slab_max)
    REAL(fp)    :: Pc_2D(n_x_max,n_y_max) !, Ec_2D(n_x_max,n_y_max)
    INTEGER, INTENT(IN)      :: n_slab_max

    !INTEGER, parameter     :: Nb = n_slab_max
    INTEGER, parameter     :: Na = 128
    !INTEGER                :: Nb  
    
    INTEGER     :: Nb_mid, Na_mid

    REAL(fp)    :: X2d(Na,n_slab_max), Y2d(Na,n_slab_max), C2d(Na,n_slab_max) !, Extra_C2d(Na,Nb)

    REAL(fp)    :: LenB, LenA
    REAL(fp)    :: Adx, Ady, Bdx, Bdy
    REAL(fp)    :: Prod, M, Lb, La

    REAL(fp)    :: C_slab(n_slab_max) !, Extra_slab(Nb)

    real(fp)  :: start, finish

    INTEGER     :: i, j

      Nb_mid = INT(n_slab_max/2)
      Na_mid = INT(Na/2)


      LenB = Pdy ! 8.0 *Height1 / Nb /SIN(theta1)
      LenA = Pdx ! 1.3 *Height1 *TAN(theta1) /Na /SIN(theta1)


      ! interval in long radius
      Adx = LenA*SIN(theta1)
      Ady = LenA*COS(theta1)

      ! interval in short radius
      Bdy = LenB*SIN(theta1)
      Bdx = LenB*COS(theta1)


      ! find the location of 1D grid in 2D XY grids
      DO i=1,Na,1
        X2d(i,Nb_mid) = -Adx*Na_mid + Adx*(i-0.5)
        Y2d(i,Nb_mid) = -Ady*Na_mid + Ady*(i-0.5)
      ENDDO


      X2d(:,Nb_mid) = X2d(:,Nb_mid) + 0.5*Bdx
      Y2d(:,Nb_mid) = Y2d(:,Nb_mid) - 0.5*Bdy

      DO j=Nb_mid+1, n_slab_max, 1
        X2d(:,j) = X2d(:,j-1) - Bdx
        Y2d(:,j) = Y2d(:,j-1) + Bdy
      ENDDO


      DO j=Nb_mid-1, 1, -1
        X2d(:,j) = X2d(:,j+1) + Bdx
        Y2d(:,j) = Y2d(:,j+1) - Bdy
      ENDDO


!      call cpu_time(start)

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( i, j )
      DO i=1,Na,1
      DO j=1,n_slab_max,1
        C2d(i,j)       = Interplt_2D(Pdx, Pdy, X2d(i,j), Y2d(i,j), Pc_2D, 2)
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO


        ! define the length/width of slab based on 95% total mass
        Lb = LenB
        La = Height1 *TAN(theta1)

        DO i=1,n_slab_max,1
          C_slab(i)     = SUM(C2d(:,i)) *(LenA*LenB)/ (La*Lb)
        ENDDO


      ! Return:
      box_Ra = La
      box_Rb = Lb
      Cslab(1:n_slab_max) = C_slab(1:n_slab_max)

  END SUBROUTINE Slab_init_bilinear

  REAL(fp) FUNCTION Interplt_2D(Pdx, Pdy, x0, y0, C_2D, ids)
    ! n_x_max, n_y_max, Pdx, Pdy are global variables
    IMPLICIT NONE


    REAL(fp)    :: Pdx, Pdy
    REAL(fp)    :: x0, y0, C_2D(n_x_max,n_y_max)
    REAL(fp)    :: X1d(n_x_max), Y1d(n_y_max)
    REAL(fp)    :: C1, C2

    integer     :: ids ! 1 for Find_theta(); 2 for Slab_init()
    integer     :: Ix0, Iy0
    integer     :: i, j

      ! define the coordinate system
      DO i=1, n_x_max
        X1d(i) = Pdx*(i-n_x_mid)
      ENDDO
      DO j=1, n_y_max
        Y1d(j) = Pdy*(j-n_y_mid)
      ENDDO


      Ix0 = floor( (x0-X1d(1)) / Pdx ) + 1
      Iy0 = floor( (y0-Y1d(1)) / Pdy ) + 1

      IF(Ix0<1 .or. Ix0>=n_x_max)THEN
        Interplt_2D = 0.0
      ELSE IF((Iy0<1 .or. Iy0>=n_y_max))THEN
        Interplt_2D = 0.0
      ELSE
        C1 = Interplt_linear(x0, X1d(Ix0), X1d(Ix0+1), &
                C_2D(Ix0,Iy0),C_2D(Ix0+1,Iy0))

        C2 = Interplt_linear(x0, X1d(Ix0), X1d(Ix0+1), &
                C_2D(Ix0,Iy0+1), C_2D(Ix0+1,Iy0+1))

        Interplt_2D = Interplt_linear(y0, Y1d(Iy0), &
                                          Y1d(Iy0+1), C1, C2)
      ENDIF

    return

  END FUNCTION Interplt_2D


  REAL(fp) FUNCTION Interplt_linear(xx, x1, x2, Cx1, Cx2)

    IMPLICIT NONE

    REAL(fp)    :: xx, x1, x2, Cx1, Cx2


    Interplt_linear = ( Cx1*(x2-xx) + Cx2*(xx-x1) ) /(x2-x1)

    return

  END FUNCTION Interplt_linear



  SUBROUTINE lagrange_write_std( am_I_Root, RC )
     
    USE Plume_list
   
    USE m_netCDF_io_define
    USE m_netcdf_io_read
    USE m_netcdf_io_open
    USE Ncdf_Mod,            ONLY : NC_Open
    USE Ncdf_Mod,            ONLY : NC_Read_Time
    USE Ncdf_Mod,            ONLY : NC_Read_Arr
    USE Ncdf_Mod,            ONLY : NC_Create
    USE Ncdf_Mod,            ONLY : NC_Close
    USE Ncdf_Mod,            ONLY : NC_Var_Def
    USE Ncdf_Mod,            ONLY : NC_Var_Write
    USE Ncdf_Mod,            ONLY : NC_Get_RefDateTime
    USE CHARPAK_Mod,         ONLY : TRANLC
    USE JulDay_Mod,          ONLY : JulDay

    USE Input_Opt_Mod, ONLY : OptInput

    USE TIME_MOD,        ONLY : GET_TS_DYN

    USE TIME_MOD,        ONLY : GET_YEAR
    USE TIME_MOD,        ONLY : GET_MONTH
    USE TIME_MOD,        ONLY : GET_DAY
    USE TIME_MOD,        ONLY : GET_HOUR
    USE TIME_MOD,        ONLY : GET_MINUTE
    USE TIME_MOD,        ONLY : GET_SECOND

    ! Parameters for netCDF routines
    include "netcdf.inc"


    LOGICAL,          INTENT(IN   ) :: am_I_Root   ! root CPU?
    INTEGER,          INTENT(INOUT) :: RC          ! Failure or success

    INTEGER :: i_box

    INTEGER :: YEAR
    INTEGER :: MONTH
    INTEGER :: DAY
    INTEGER :: HOUR
    INTEGER :: MINUTE
    INTEGER :: SECOND

    CHARACTER(LEN=255)            :: FILENAME2, FILENAME3

    CHARACTER(LEN=25) :: YEAR_C
    CHARACTER(LEN=25) :: MONTH_C
    CHARACTER(LEN=25) :: DAY_C
    CHARACTER(LEN=25) :: HOUR_C
    CHARACTER(LEN=25) :: MINUTE_C
    CHARACTER(LEN=25) :: SECOND_C


    TYPE(Plume2d_list), POINTER :: Plume2d_new, PLume2d, Plume2d_prev
    TYPE(Plume1d_list), POINTER :: Plume1d_new, Plume1d, Plume1d_prev

    REAL(fp) :: Dt


    Dt = GET_TS_DYN()

    YEAR        = GET_YEAR()
    MONTH       = GET_MONTH()
    DAY         = GET_DAY()
    HOUR        = GET_HOUR()
    MINUTE      = GET_MINUTE()
    SECOND      = GET_SECOND()

    WRITE(YEAR_C,*) YEAR
    WRITE(MONTH_C,*) MONTH
    WRITE(DAY_C,*) DAY
    WRITE(HOUR_C,*) HOUR
    WRITE(MINUTE_C,*) MINUTE
    WRITE(SECOND_C,*) SECOND

    !-------------------------------------------------------------------
    ! Output box location
    !-------------------------------------------------------------------
    IF(mod(tt*NINT(Dt),24*3600)==0)THEN   ! output once every day (24 hours)

    FILENAME3 = 'Plume_number.txt'

    OPEN( 261,      FILE=TRIM( FILENAME3   ), STATUS='OLD',  &
          POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )



    WRITE(261,*) Num_Plume1d, Num_Plume2d, Num_dissolve, Num_inject


    ENDIF ! IF(mod(tt,1440)==0)THEN



    CLOSE(261)

    IF(ASSOCIATED(Plume2d_new)) nullify(Plume2d_new)
    IF(ASSOCIATED(PLume2d)) nullify(PLume2d)
    IF(ASSOCIATED(Plume2d_prev)) nullify(Plume2d_prev)
    IF(ASSOCIATED(Plume1d_new)) nullify(Plume1d_new)
    IF(ASSOCIATED(Plume1d)) nullify(Plume1d)
    IF(ASSOCIATED(Plume1d_prev)) nullify(Plume1d_prev)


  END SUBROUTINE lagrange_write_std

!-------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------

  subroutine lagrange_cleanup()


    IF (.not.plume_inject_on) THEN
      WRITE(6,'(a)') ' No plume injection, no lagrangian module configuration, skip cleanup'
      RETURN
    ENDIF


!
!    if (allocated(box_lon))      deallocate(box_lon)
!    if (allocated(box_lat))      deallocate(box_lat)
!    if (allocated(box_lev))      deallocate(box_lev)
!    if (allocated(box_length))   deallocate(box_length)
!    if (allocated(box_Ra))       deallocate(box_Ra)
!    if (allocated(box_Rb))       deallocate(box_Rb)
!    if (allocated(box_theta))    deallocate(box_theta)
!
    WRITE(6,'(a)') '--> Lagrange and Plume Module Cleanup <--'
!
!
  end subroutine lagrange_cleanup


END MODULE Lagrange_Mod



RECURSIVE FUNCTION SortList(head0) RESULT (new_head)

  USE Plume_list

  INTERFACE
      FUNCTION MergeSort(head1, head2)
       USE Plume_list
       TYPE(Plume1d_list), POINTER :: head1, head2
       TYPE(Plume1d_list), POINTER :: MergeSort
      END FUNCTION
  END INTERFACE


    TYPE(Plume1d_list), POINTER :: head0, new_head
    TYPE(Plume1d_list), POINTER :: fast, slow


    IF(.NOT.ASSOCIATED(head0) .OR. .NOT.ASSOCIATED(head0%next))THEN
        new_head => head0
    ELSE

        ! use fast/slow search for the middle point
        fast => head0
        slow => head0

        DO WHILE(ASSOCIATED(fast%next).AND.ASSOCIATED(fast%next%next))
          fast => fast%next%next
          slow => slow%next
        ENDDO

        ! cut the whole list into two parts:
        ! fast: left half
        ! slow: right half
        fast => slow
        slow => slow%next
        NULLIFY(fast%next)

        ! recursion
        fast => SortList(head0)
        slow => SortList(slow)

        new_head => MergeSort( fast, slow )
    ENDIF

END FUNCTION SortList




FUNCTION MergeSort(head1, head2)
   USE Plume_list

   TYPE(Plume1d_list), POINTER :: head1, head2, head3
   TYPE(Plume1d_list), POINTER :: p1, p2, p
   TYPE(Plume1d_list), POINTER :: MergeSort

   p1 => head1
   p2 => head2

   ! First determine the new head
   IF( head1%RA*head1%RB*head1%LENGTH &
                > head2%RA*head2%RB*head2%LENGTH )THEN
        head3 =>head1
        p1 => p1%next
   ELSE
        head3 => head2
        p2 => p2%next
   ENDIF


   ! Second, sort between two half lists
   p => head3
   DO WHILE(ASSOCIATED(p1).AND.ASSOCIATED(p2))
     IF( p1%RA*p1%RB*p1%LENGTH &
                > p2%RA*p2%RB*p2%LENGTH )THEN
        p%next => p1
        p1 => p1%next
        p => p%next
     ELSE
        p%next => p2
        p2 => p2%next
        p => p%next
     ENDIF
   ENDDO

   IF(ASSOCIATED(p1)) p%next => p1
   IF(ASSOCIATED(p2)) p%next => p2

   MergeSort => head3

   RETURN

END FUNCTION MergeSort
