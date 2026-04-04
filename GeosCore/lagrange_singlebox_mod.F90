! Note (BZ):
! - Currently not interpolate temperature and other meteorology in solving KPP-related chemistry
! - Diagnostic file (to be added)
!     1. Plume lifetime
!         Plume_id, Plume_life
!     2. Plume_number
!         Dt, Num_Plume_inject, Num_Plume_2d, Num_Plume_1d, Num_Plume_dissolve
!     3. Sulfate Mass
!         Dt, Mass_S_injected, Mass_S_released
! - Occasional KPP integration error:
!     Forced exit from Rosenbrock due to the following error:
!     --> Step size too small: T + 10*H = T or H < Roundoff
!     T=  1.565122350165089E-016 and H=  1.565122350165089E-016
!     ### INTEGRATE RETURNED ERROR AT:           43          16          39
!- CFL condition and plume box resize issue
! - Better treat H2O for KPP-related Chemistry
!     Currently read Eulerian background in Set_inPlume_2d_Kpp_GridBox_Values
! Sulfate-SS and sulfate-Cloud heteorogeneous reactions are included in KPP, 
! Not read in SS to supress SS-Sulfate heteorogeneous reaction
! In cloud sulfate oxidation? 
MODULE Lagrange_singlebox_Mod
  USE Plume_list_Mod
  USE PRECISION_MOD
  USE ERROR_MOD
  USE ERRCODE_MOD
  USE PhysConstants,   ONLY : PI, Re, g0, AIRMW, AVO, BOLTZ
  USE TIME_MOD,        ONLY : GET_YEAR, GET_MONTH, GET_DAY, GET_HOUR, GET_MINUTE, GET_SECOND
  USE TIME_MOD,        ONLY : ITS_TIME_FOR_EXIT
  USE UNITCONV_MOD    
  USE INPUT_OPT_MOD,   ONLY : PlumeSource_t
  USE gckpp_Parameters, ONLY: NREACT

  IMPLICIT NONE
  PRIVATE

  !PUBLIC MEMBER FUNCTIONS
  PUBLIC :: lagrange_init_box
  PUBLIC :: plume_inject_box
  PUBLIC :: plume_model_box
  PUBLIC :: plume_mod_cleanup_box

  ! PUBLIC VARIABLES:
  PUBLIC :: RXNRATE_CONST_KPP
  !PUBLIC :: n_x_max, n_y_max

  ! variables read from geoschem_config.yml
  LOGICAL                               :: use_lagrange
  LOGICAL                               :: plume_inject_on 
  LOGICAL                               :: plume_diag
  LOGICAL                               :: TROPP_sink
  INTEGER                               :: Num_of_sources
  INTEGER                               :: n_x_max            !number of x grids in 2D, should be (9 x odd)
  INTEGER                               :: n_y_max            !number of y grids in 2D, should be (9 x odd)
  INTEGER                               :: N_stop_inject      ! How many parcel to be injected, -1 means continuous injection           
  TYPE(PlumeSource_t), ALLOCATABLE      :: Plume_sources(:)
  REAL(fp)                              :: Dx_init
  REAL(fp)                              :: Dy_init
  REAL(fp)                              :: Length_init    ! m
  REAL(fp)                              :: Aircraft_speed ! m/s
  REAL(fp), ALLOCATABLE                 :: RXNRATE_CONST_KPP(:,:,:,:)
  ! Other variables
  INTEGER               :: n_species
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
  INTEGER               :: N_total
  INTEGER               :: Stop_inject ! 1: stop injecting; 0: keep injecting


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
  SUBROUTINE lagrange_init_box(am_I_root, Input_Opt, State_Chm, State_Grid, State_Met, RC)

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
    
    INTEGER                       :: i_box, i_slab
    INTEGER                       :: ii, jj, kk
    INTEGER                       :: id_tracer, N
    INTEGER                       :: i_lon, i_lat, i_lev            !1:IIPAR
    INTEGER                       :: previous_units, previous_units_temp
    CHARACTER(LEN=255)            :: spc_name
    CHARACTER(LEN=255)            :: FILENAME, FileEntropy, File996
    CHARACTER(LEN=255)            :: FILENAME2, FILENAME3
    CHARACTER(LEN=255)            :: ErrMsg, ThisLoc
    
    REAL(fp)                      :: lon1, lat1, lon2, lat2
    REAL(fp)                      :: box_lon_edge, box_lat_edge
    REAL(fp)                      :: curr_lon, curr_lat, curr_lev
    REAL(fp)                      :: Dt
    REAL(fp)                      :: plume_len_deg, plume_rem_deg, add_deg 
    REAL(fp), PARAMETER           :: eps = 1.0e-10  

    RC                 =   GC_SUCCESS
    ErrMsg             =   ''
    ThisLoc            =   ' -> at lagrange_init_box (in module GeosCore/lagrange_singlebox_mod.F90)'
    Spc                =>   State_Chm%Species

    Num_inject         =    0
    Num_Plume2d        =    0
    Num_Plume1d        =    0
    Num_dissolve       =    0

    n_species          =              State_Chm%nSpecies
    IIPAR              =              State_Grid%NX
    JJPAR              =              State_Grid%NY
    LLPAR              =              State_Grid%NZ
    DX                 =              State_Grid%DX
    DY                 =              State_Grid%DY
    X_edge             =>             State_Grid%XEdge(:,1) 
    Y_edge             =>             State_Grid%YEdge(1,:)
    P_edge             =>             State_Met%PEDGE(1,1,:) 
    X_edge2            =              X_edge(2)
    Y_edge2            =              Y_edge(2)
    X_mid              =>             State_Grid%XMid(:,1) ! Grid box longitude [degrees] ! XMID(:,1,1)   ! IIPAR ! new
    Y_mid              =>             State_Grid%YMid(1,:) ! Grid box latitude center [degree] ! YMID(1,:,1)
    P_mid              =>             State_Met%PMID(1,1,:)  ! Pressure at level centers (hPa)
    
    Write (6, *) "Debug: (BZ) Num of reactions from KPP is: ", NREACT

    ALLOCATE( RXNRATE_CONST_KPP(IIPAR, JJPAR, LLPAR, NREACT ), STAT=RC )
    IF (RC /= 0) THEN
        errMsg = 'Error allocating RXNRATE_CONST_KPP'
        CALL ERROR_STOP( errMsg, thisLoc)
        RETURN
    ENDIF

    ! Copy all neccessary variable from geoschem.yml here
    use_lagrange        =             Input_Opt%LagrangianModel_Activate
    plume_inject_on     =             Input_Opt%PlumeInjection_Activate
    plume_diag          =             Input_Opt%PlumeInjection_Diag
    TROPP_sink          =             Input_Opt%TropSink_Activate
    Num_of_sources      =             Input_Opt%Plume_sources_num
    Length_init         =             Input_Opt%Initial_length*1000.0   ! km to m
    Aircraft_speed      =             Input_Opt%Aircraft_speed  !m/s
    !plume_interval      =             Input_Opt%plume_interval ! degree
    N_stop_inject       =             Input_Opt%Num_of_injection
    
    IF (Num_of_sources > 0) THEN
      ALLOCATE( Plume_sources(Num_of_sources), STAT=RC )
      IF (RC /= 0) THEN
          errMsg = 'Error allocating Plume_sources'
          !CALL GC_Error( errMsg, RC, thisLoc )
          CALL ERROR_STOP( errMsg, thisLoc)
          RETURN
       END IF
       DO N = 1, Num_of_sources
          Plume_sources(N)%lat1        =  MIN(Input_Opt%Plume_sources(N)%lat1, Input_Opt%Plume_sources(N)%lat2)
          Plume_sources(N)%lat2        =  MAX(Input_Opt%Plume_sources(N)%lat1, Input_Opt%Plume_sources(N)%lat2)
          ! The total meridional distance is approximated as (lat2 - lat1) * 111 km 
          ! If the total distance is not divisible by plume_length, the final plume
          ! segment would be shorter than plume_length and would partially overlap
          ! the previous segment near the boundary. To avoid handling partial plume
          ! segments, lat2 is slightly adjusted so that the total distance is an
          ! integer multiple of plume_length.
          plume_len_deg                =  Length_init / 1000.0 / 111.0  ! m -> km -> degrees
          plume_rem_deg                =  mod((Plume_sources(N)%lat2-Plume_sources(N)%lat1),plume_len_deg)
          IF (plume_rem_deg> eps) THEN
            add_deg                    =  plume_len_deg - plume_rem_deg
            Plume_sources(N)%lat2      =  Plume_sources(N)%lat2 + add_deg
          ENDIF
          Plume_sources(N)%lon         =  Input_Opt%Plume_sources(N)%lon
          Plume_sources(N)%lev         =  Input_Opt%Plume_sources(N)%lev
          Plume_sources(N)%rate        =  Input_Opt%Plume_sources(N)%rate
          Plume_sources(N)%species     =  Input_Opt%Plume_sources(N)%species
       ENDDO
    ENDIF

    n_x_max                    =      Input_Opt%PlumeGrid2d_nx ! odd number to make sure n_x_mid to be integer and divisible by 9
    n_y_max                    =      Input_Opt%PlumeGrid2d_ny ! odd number to make sure n_y_mid to be integer and divisible by 9
    Dx_init                    =      Input_Opt%PlumeGrid2d_dx
    Dy_init                    =      Input_Opt%PlumeGrid2d_dy
    ! Check n_x_max and n_y_max to be odd number and divisible by 9
    IF((MOD(n_x_max,9).ne.0).or.(MOD(n_x_max,2).eq.0) ) THEN
      ErrMsg = '*** ERROR,  n_x_max should be an odd number and divisible by 9 ***'
      !CALL GC_Error( ErrMsg, RC, ThisLoc )
      CALL ERROR_STOP (ErrMsg, ThisLoc)
    ENDIF
    IF((MOD(n_y_max,9).ne.0).or.(MOD(n_y_max,2).eq.0) ) THEN
      ErrMsg = '*** ERROR,  n_y_max should be an odd number and divisible by 9 ***'
      !CALL GC_Error( ErrMsg, RC, ThisLoc )
      CALL ERROR_STOP (ErrMsg, ThisLoc)
    ENDIF
    ! the odd number of n_x_max can ensure a center grid
    n_x_mid                    =      (n_x_max+1)/2 
    n_y_mid                    =      (n_y_max+1)/2  
    n_x_max2                   =      n_x_max+2
    n_y_max2                   =      n_y_max+2
    n_x_mid2                   =      (n_x_max2+1)/2 
    n_y_mid2                   =      (n_y_max2+1)/2

    
    Dt = GET_TS_DYN()
    N_parcel = NINT(Aircraft_speed * Dt / Length_init)
    IF (N_parcel<1) N_parcel=1
    WRITE(6,*) 'Debug (BZ): Injected plume every time step: ', N_parcel
  
    Stop_inject = 0

    NULLIFY(Plume2d_head, Plume2d_tail)
    NULLIFY(Plume1d_head, Plume1d_tail)
    
    IF (.not.plume_inject_on) THEN
      WRITE(6,'(a)') ' No plume injection, no lagrangian module configuration, skip initialization'
      RETURN
    ENDIF

    IF (plume_inject_on .AND. num_of_sources.lt.1) THEN
      ErrMsg = 'Plume injection is turned on but no source specified, skip initialization'
      CALL GC_Warning( ErrMsg, RC, ThisLoc )
      !CALL ERROR_STOP (ErrMsg, ThisLoc)
      RETURN
    ENDIF

    ! Convert unit from mol/mol dry to kg/kg dry to molec/cm3
    ! unit is mol/mol dry, first convert mol/mol dry to kg/kg dry
    !CALL Convert_Spc_Units(                                            &
    !            Input_Opt      = Input_Opt,                                   &
    !            State_Chm      = State_Chm,                                   &
    !            State_Grid     = State_Grid,                                  &
    !            State_Met      = State_Met,                                   &
    !            new_units      = KG_SPECIES_PER_KG_DRY_AIR,                   &
    !            previous_units = previous_units,                              &
    !            RC             = RC                                          )
    ! then convert kg/kg dry to molec/cm3, because no direct conversion from mol/mol to molec/cm3
    !CALL Convert_Spc_Units(                                            &
    !            Input_Opt      = Input_Opt,                                   &
    !            State_Chm      = State_Chm,                                   &
    !            State_Grid     = State_Grid,                                  &
    !            State_Met      = State_Met,                                   &
    !            new_units      = MOLECULES_SPECIES_PER_CM3,                   &
    !            previous_units = previous_units_temp,                              &
    !            RC             = RC                                          )
    !WRITE(6,'(a)') 'debug in plume_initialization: before initialization, the new unit is ' // TRIM(UNIT_STR(previous_units))
    !WRITE(6,'(a)') 'debug in plume_initialization: during initialization, the new unit is ' // TRIM(UNIT_STR(State_Chm%Species(1)%Units))
    ! Check that species units are in  molec/cm3
    !id_SO4= Ind_('SO4')
    !IF ( Spc(id_SO4)%Units /= MOLECULES_SPECIES_PER_CM3 ) THEN
    !  ErrMsg = 'Incorrect species units: ' // TRIM(UNIT_STR(Spc(id_SO4)%Units))
      !CALL GC_Error( ErrMsg, RC, ThisLoc )
    !   CALL ERROR_STOP (ErrMsg, ThisLoc)
    !ENDIF
    ! Read injection location at the first time step, always start from lat1
    !curr_lon    =  Plume_sources(1)%lon
    !curr_lat    =  Plume_sources(1)%lat1
    !curr_lev    =  Plume_sources(1)%lev
    !i_lon = Find_iLonLat(curr_lon, DX, X_edge2)
    !if(i_lon>IIPAR) i_lon=i_lon-IIPAR
    !if(i_lon<1) i_lon=i_lon+IIPAR
    !WRITE(6,*) 'debug: ilon=', i_lon
    !i_lat = Find_iLonLat(curr_lat, DY, Y_edge2)
    !if(i_lat>JJPAR) i_lat=JJPAR
    !if(i_lat<1) i_lat=1
    !WRITE(6,*) 'debug: ilat=', i_lat
    !i_lev = Find_iPLev(curr_lev,P_edge)
    !WRITE(6,*) 'debug: ilev1=', i_lev
    !if(i_lev>LLPAR) i_lev=LLPAR

    !IF (use_lagrange) THEN
      
      !! Initialize the first plume segment
      !ALLOCATE(Plume2d_tail)
      !Plume2d_tail%IsNew = 1
      !Plume2d_tail%label = 1
      !Plume2d_tail%active_x_min = 1
      !Plume2d_tail%active_x_max = n_x_max
      !Plume2d_tail%active_y_min = 1
      !Plume2d_tail%active_y_max = n_y_max
      !Plume2d_tail%lat_ind = i_lat
      !Plume2d_tail%lon_ind = i_lon
      !Plume2d_tail%lev_ind = i_lev

      !Plume2d_tail%LON = Inject_lon
      !Plume2d_tail%LAT = -29.95e+0_fp
      !Plume2d_tail%LEV = Inject_hPa
      !Plume2d_tail%LON   = Plume_sources(1)%lon
      !Plume2d_tail%LAT   = Plume_sources(1)%lat1
      !Plume2d_tail%LEV   = Plume_sources(1)%lev

      !Plume2d_tail%ALPHA  = 0.0e+0_fp
      !Plume2d_tail%LIFE = 0.0e+0_fp

      !Plume2d_tail%LENGTH = Length_init ! m 
      !Plume2d_tail%PDX = Dx_init
      !Plume2d_tail%PDY = Dy_init
      
      !ALLOCATE(Plume2d_tail%CONCNT2d(n_x_max, n_y_max, n_species)) 
      !Plume2d_tail%CONCNT2d = 0.0e+0_fp ! molec/cm3
      
      ! pass all species as from background as initial concentration
      !DO N = 1, n_species
      !  Plume2d_tail%CONCNT2d(:,:,N) = Spc(N)%Conc(i_lon, i_lat, i_lev)  ! molec/cm3
      !ENDDO

      ! add injected species to the center of plume grid
      !DO N = 1, num_of_sources
      !  spc_name = Plume_sources(N)%species
      !  id_tracer   = Ind_(spc_name)
      !  Plume2d_tail%CONCNT2d(n_x_mid,n_y_mid,id_tracer) = Plume2d_tail%CONCNT2d(n_x_mid,n_y_mid,id_tracer)   &
      !              + (Plume2d_tail%LENGTH * Plume_sources(N)%rate / State_Chm%SpcData(id_tracer)%Info%MW_g * Avo) &
      !              /(Plume2d_tail%PDX * Plume2d_tail%PDY *Plume2d_tail%LENGTH*1.E6_fp ) ! molec/cm3
                    
         !Plume2d_tail%CONCNT2d(n_x_mid,n_y_mid,n_species+N) = Plume2d_tail%CONCNT2d(n_x_mid,n_y_mid,n_species+N)  &
         !           + (Plume2d_tail%LENGTH * Plume_sources(N)%rate / State_Chm%SpcData(id_tracer)%Info%MW_g * Avo) &
         !           /(Plume2d_tail%PDX * Plume2d_tail%PDY *Plume2d_tail%LENGTH*1.E6_fp ) ! molec/cm3
      !ENDDO
      !NULLIFY(Plume2d_tail%next)
      !Plume2d_head => Plume2d_tail
      !Num_Plume2d = 1
    !ELSE
      ! add injected species into eulerian grid
      !WRITE(6,'(a)') ' lagrange module was turned off, add injected species to eulerian grid'
      !DO N = 1, num_of_sources
      !  spc_name = Plume_sources(N)%species
      !  id_tracer   = Ind_(spc_name)
      !  Spc(id_tracer)%Conc(i_lon, i_lat, i_lev) = Spc(id_tracer)%Conc(i_lon, i_lat, i_lev)  &
      !          + (Length_init * Plume_sources(N)%rate / State_Chm%SpcData(id_tracer)%Info%MW_g * Avo) &
      !              /(State_Met%AIRVOL(i_lon, i_lat, i_lev)*1.E6_fp ) ! molec/cm3
      !ENDDO
    !ENDIF
    !Num_inject = 1
    !tt = 0
    ! convert back from molec/cm3 to kg/kg dry
    !CALL Convert_Spc_Units(                                            &
    !            Input_Opt  = Input_Opt,                                       &
    !            State_Chm  = State_Chm,                                       &
    !            State_Grid = State_Grid,                                      &
    !            State_Met  = State_Met,                                       &
    !            new_units  = previous_units_temp,                                  &
    !            RC         = RC                                              )
    ! convert back from kg/kg dry to v/ v dry
    !CALL Convert_Spc_Units(                                            &
    !            Input_Opt  = Input_Opt,                                       &
    !            State_Chm  = State_Chm,                                       &
    !            State_Grid = State_Grid,                                      &
    !            State_Met  = State_Met,                                       &
    !            new_units  = previous_units,                                  &
    !            RC         = RC                                              )
    !WRITE(6,'(a)') 'debug in plume_initialize: after initialization, the unit is ' // TRIM(UNIT_STR(State_Chm%Species(1)%Units))
90  FORMAT( /, A                 )
95  FORMAT( A                    )
100 FORMAT( A, L5                )
105 FORMAT( A, I0                )
110 FORMAT( A, A                 )
120 FORMAT(A8,1X,"|", A8,1X,"|",A8,1X,"|",1X,A8,1X,"|",1X,A8,1X,"|",1X,A10,1X,"|",1X,A)
121 FORMAT(I8,1X,"|",F8.3,1X,"|",F8.3,1X,"|",1X,F8.3,1X,"|",1X,F8.3,1X,"|",1X,F10.3,1X,"|",1X,A)

  END SUBROUTINE lagrange_init_box
  
  SUBROUTINE plume_inject_box(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)

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
    
    TYPE(SpcConc), POINTER        :: Spc(:)
    TYPE(Plume2d_list), POINTER   :: Plume2d_new

    !REAL(fp), POINTER :: X_edge(:), Y_edge(:)
    !REAL(fp)          :: X_edge2, Y_edge2


    INTEGER                :: i_box, i_lon, i_lat, i_lev, i_species
    INTEGER                :: nAdv, i_advect1
    INTEGER                :: id_tracer
    INTEGER                :: previous_units, previous_units_temp
    INTEGER                :: Num_Stop
    REAL(fp)               :: box_lon, box_lat, box_lev
    REAL(fp), POINTER      :: PASV_EU
    REAL(fp)               :: MW_g
    REAL(fp)               :: Dt
    REAL(fp)               :: Vgrid_2D, Vgrid_1D
    REAL(fp)               :: Entropy2_Concnt, Entropy2_V, Entropy2
    REAL(fp)               :: tracer2_mol, air2_mol, mix2_ratio
    REAL(fp)               :: tracer0_mol, air0_mol, mix0_ratio
    !REAL(fp)               :: plume_length_deg, plane_route_deg, plane_loc_last

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
    ThisLoc                =   ' -> at plume_inject_box (in module GeosCore/lagrange_singlebox_mod.F90)'

    id_SO4= Ind_('SO4')

    IF (.NOT.plume_inject_on) THEN
      !WRITE(6,'(a)') ' No plume injection, no lagrangian module configuration, skip '
      RETURN
    ENDIF
    IF (plume_inject_on .AND. num_of_sources.LT.1) THEN
      !ErrMsg = 'Plume injection is turned on but no source specified'
      !CALL GC_Error( ErrMsg, RC, ThisLoc )
      !CALL ERROR_STOP (ErrMsg, ThisLoc)
      RETURN
    ENDIF

    ! In theory, before plume injection, the unit should be kg/kg
    CALL Convert_Spc_Units(                                            &
            Input_Opt      = Input_Opt,                                   &
            State_Chm      = State_Chm,                                   &
            State_Grid     = State_Grid,                                  &
            State_Met      = State_Met,                                   &
            new_units      = MOLECULES_SPECIES_PER_CM3,                   &
            previous_units = previous_units_temp,                              &
            RC             = RC                                          )
    WRITE(6,'(a)') 'debug : Before plume injection, the unit is ' // TRIM(UNIT_STR(previous_units_temp))
    WRITE(6,'(a)') 'debug : During plume injection, the unit is ' // TRIM(UNIT_STR(State_Chm%Species(id_SO4)%Units))
    
    ! Check that species units are in  molec/cm3
    IF ( Spc(id_SO4)%Units /= MOLECULES_SPECIES_PER_CM3 ) THEN
      ErrMsg = 'Incorrect species units: ' // TRIM(UNIT_STR(Spc(id_SO4)%Units))
      !CALL GC_Error( ErrMsg, RC, ThisLoc )
       CALL ERROR_STOP (ErrMsg, ThisLoc)
    ENDIF
    ! -----------------------------------------------------------
    ! add new box every time step
    ! -----------------------------------------------------------
    IF(N_stop_inject<0.OR.Num_inject<N_stop_inject)THEN
      IF (N_stop_inject < 0) THEN
        Num_Stop = Num_inject + N_parcel
        ELSE
        Num_Stop = MIN(Num_inject + N_parcel, N_stop_inject)
      END IF
        IF (.NOT. use_lagrange) THEN
            DO i_box = Num_inject+1, Num_Stop, 1
                
              !IF (i_box .GT. N_stop_inject) EXIT 
                
                box_lon    = Plume_sources(1)%lon
                box_lat    = GetInjectionLat(Plume_sources(1)%lat1,Plume_sources(1)%lat2, Length_init, i_box - 1)
                write(6,*) 'debug (BZ): injection lat: ', box_lat, 'Num_inject: ', i_box
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
                !WRITE(6,*) 'debug: ilev1=', i_lev
                if(i_lev>LLPAR) i_lev=LLPAR
                ! instantly add injected species into Eulerian grid 
                DO i_species = 1, num_of_sources
                    spc_name = Plume_sources(i_species)%species
                    id_tracer   = Ind_(TRIM(spc_name))
                    write(6,*) 'debug (BZ): species: ', TRIM(spc_name), 'id: ', id_tracer, 'conc before injection: ', Spc(id_tracer)%Conc(i_lon, i_lat, i_lev)
                    Spc(id_tracer)%Conc(i_lon, i_lat, i_lev) = Spc(id_tracer)%Conc(i_lon, i_lat, i_lev)  &
                            + (Length_init * Plume_sources(i_species)%rate / State_Chm%SpcData(id_tracer)%Info%MW_g * Avo) &
                            /(State_Met%AIRVOL(i_lon, i_lat, i_lev) * 1.0e+6_fp) ! molec/cm3
                    write(6,*) 'debug (BZ): species: ', spc_name, 'id: ', id_tracer, 'conc after injection: ', Spc(id_tracer)%Conc(i_lon, i_lat, i_lev)

                ENDDO
            ENDDO
            Num_inject = Num_Stop 

        ELSE
            DO i_box = Num_inject+1, Num_Stop, 1
                ALLOCATE(Plume2d_new)
                Plume2d_new%IsNew           = 1
                Plume2d_new%label           = i_box

                Plume2d_new%LON             = Plume_sources(1)%lon
                Plume2d_new%LAT             = GetInjectionLat(Plume_sources(1)%lat1,Plume_sources(1)%lat2, Length_init, i_box - 1)
                write(6,*) 'debug (BZ): injection lat: ', Plume2d_new%LAT, 'Num_inject: ', i_box
                Plume2d_new%LEV             = Plume_sources(1)%lev
                 
                i_lon = Find_iLonLat(Plume2d_new%LON, DX, X_edge2)
                if(i_lon>IIPAR) i_lon=i_lon-IIPAR
                if(i_lon<1) i_lon=i_lon+IIPAR
                i_lat = Find_iLonLat(Plume2d_new%LAT, DY, Y_edge2)
                if(i_lat>JJPAR) i_lat=JJPAR
                if(i_lat<1) i_lat=1
                i_lev = Find_iPLev(Plume2d_new%LEV,P_edge)
                if(i_lev>LLPAR) i_lev=LLPAR
                Plume2d_new%lat_ind = i_lat
                Plume2d_new%lon_ind = i_lon
                Plume2d_new%lev_ind = i_lev

                Plume2d_new%LENGTH = Length_init ! 1000m 
                Plume2d_new%ALPHA  = 0.0e+0_fp
                Plume2d_new%LIFE   = 0.0e+0_fp
                Plume2d_new%PDX    = Dx_init
                Plume2d_new%PDY    = Dy_init

                Vgrid_2D           = (Plume2d_new%PDX * Plume2d_new%PDY * Plume2d_new%LENGTH ) *1.E6_fp ! cm3
                ALLOCATE(Plume2d_new%CONCNT2d(n_x_max, n_y_max, n_species))
                ALLOCATE(Plume2d_new%MassRef2d(n_species))
                Plume2d_new%CONCNT2d = 0.0e+0_fp
                ! pass background concentration
                DO i_species = 1, n_species
                    Plume2d_new%CONCNT2d(:,:,i_species) = Spc(i_species)%Conc(i_lon, i_lat, i_lev)  ! molec/cm3
                    Plume2d_new%MassRef2d(i_species) = Spc(i_species)%Conc(i_lon, i_lat, i_lev) &
                                                              * Vgrid_2D  *  n_x_max * n_y_max ! molec 
                ENDDO
                ! add tracer concentration
                DO i_species = 1, num_of_sources
                    spc_name = Plume_sources(i_species)%species
                    id_tracer   = Ind_(TRIM(spc_name))
                    Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,id_tracer) = Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,id_tracer)   &
                                + (Plume2d_new%LENGTH * Plume_sources(i_species)%rate / State_Chm%SpcData(id_tracer)%Info%MW_g * Avo) &
                                /(Plume2d_new%PDX * Plume2d_new%PDY *Plume2d_new%LENGTH*1.E6_fp ) ! molec/cm3
                                
                ENDDO
                NULLIFY(Plume2d_new%next)
                IF ( .NOT. ASSOCIATED(Plume2d_head) ) THEN
                  ! First node in the list
                  Plume2d_head => Plume2d_new
                  Plume2d_tail => Plume2d_new
                ELSE
                  ! Append to list
                  Plume2d_tail%next => Plume2d_new
                  Plume2d_tail => Plume2d_tail%next
                ENDIF
                Num_Plume2d = Num_Plume2d + 1
            ENDDO
            Num_inject = Num_Stop 
        ENDIF
    ENDIF

    CALL Convert_Spc_Units(                                            &
               Input_Opt      = Input_Opt,                                   &
               State_Chm      = State_Chm,                                   &
               State_Grid     = State_Grid,                                  &
               State_Met      = State_Met,                                   &
               new_units      = previous_units_temp,                   &
               RC             = RC                                          )
    WRITE(6,'(a)') 'debug: after plume injection, the unit is ' // TRIM(UNIT_STR(State_Chm%Species(id_SO4)%Units))

  END SUBROUTINE plume_inject_box

  SUBROUTINE plume_model_box(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
  
    USE Input_Opt_Mod,   ONLY : OptInput, PlumeSource_t
    USE State_Chm_Mod,   ONLY : ChmState, Ind_
    USE State_Met_Mod,   ONLY : MetState
    USE Species_Mod,     ONLY : SpcConc
    USE TIME_MOD,        ONLY : GET_TS_DYN
    USE State_Grid_Mod,  ONLY : GrdState
    !USE State_Diag_Mod,           ONLY : DgnState
    !USE State_Diag_Mod,           ONLY : DgnMap
    USE UnitConv_Mod

    LOGICAL, INTENT(IN)           :: am_I_Root
    TYPE(MetState), INTENT(IN)    :: State_Met
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State objectgg
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    !TYPE(DgnState), INTENT(INOUT) :: State_Diag ! Diagnostics State object
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    INTEGER                :: previous_units, previous_units_temp

    REAL(fp)               :: Dt

    !CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: spc_name
    CHARACTER(LEN=255)     :: ErrMsg
    CHARACTER(LEN=255)     :: ThisLoc
    
    TYPE(SpcConc), POINTER :: Spc(:)
    !TYPE(Plume2d_list), POINTER :: Plume2d_tail=> NULL(), Plume2d_head=> NULL()

    !X_edge => State_Grid%XEdge(:,1) 
    !Y_edge => State_Grid%YEdge(1,:) 
    !X_edge2       = X_edge(2)
    !Y_edge2       = Y_edge(2)
    ErrMsg                 =    ''
    Dt                     =    GET_TS_DYN()
    Spc                    =>   State_Chm%Species
    ThisLoc                =   ' -> at plume_model_box (in module GeosCore/lagrange_singlebox_mod.F90)'

    id_SO4= Ind_('SO4')
    
    IF (use_lagrange .AND. plume_inject_on) THEN

      ! In theory, before plume injection, the unit should be kg/kg
      CALL Convert_Spc_Units(                                            &
              Input_Opt      = Input_Opt,                                   &
              State_Chm      = State_Chm,                                   &
              State_Grid     = State_Grid,                                  &
              State_Met      = State_Met,                                   &
              new_units      = MOLECULES_SPECIES_PER_CM3,                   &
              previous_units = previous_units_temp,                              &
              RC             = RC                                          )
      WRITE(6,'(a)') 'debug : Before plume box model, the unit is ' // TRIM(UNIT_STR(previous_units_temp))
      WRITE(6,'(a)') 'debug : During plume box model, the unit is ' // TRIM(UNIT_STR(State_Chm%Species(id_SO4)%Units))

      ! Check that species units are in  molec/cm3
      IF ( Spc(id_SO4)%Units /= MOLECULES_SPECIES_PER_CM3 ) THEN
        ErrMsg = 'Incorrect species units: ' // TRIM(UNIT_STR(Spc(id_SO4)%Units))
        !CALL GC_Error( ErrMsg, RC, ThisLoc )
        CALL ERROR_STOP (ErrMsg, ThisLoc)
      ENDIF

      CALL plume_physics(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
      ! Plume physical evolution: 
      ! Update plume lifetime 
      ! movement (advection); plume stretching; adiabatic volume change due to thermodynamics (p-T-V); 
      ! Update in-plume concentration due to volume change 
      ! In-plume species concentration: 
      ! advection, diffusion; entrainment

      ! BZ: See Chemistry_mod.F90
      ! Before Chemistry, do I need to set CO2 to 421ppm?This is set for Eulerian grid before chemistry
      ! This is necessary to reduce the error norm in KPP.
      ! See https://github.com/geoschem/geos-chem/issues/1529.

      CALL plume_chem_microphysics(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
      ! New module update chemistry and microphysics

      
      CALL plume_structure_change(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
      ! Plume structure evolution:
      ! 2D to 1D
      ! Dissolve when meet the criteria 

      ! write diagnostic output
      CALL lagrange_write_std( am_I_Root, RC )
      ! convert unit back
      CALL Convert_Spc_Units(                                            &
               Input_Opt      = Input_Opt,                                   &
               State_Chm      = State_Chm,                                   &
               State_Grid     = State_Grid,                                  &
               State_Met      = State_Met,                                   &
               new_units      = previous_units_temp,                   &
               RC             = RC                                          )
      WRITE(6,'(a)') 'debug: after plume box model, the unit is ' // TRIM(UNIT_STR(State_Chm%Species(id_SO4)%Units))

    ELSE
      RETURN
    ENDIF
    
  END SUBROUTINE plume_model_box
  
  SUBROUTINE plume_mod_cleanup_box( RC)

     USE ErrCode_Mod

    INTEGER, INTENT(OUT) :: RC          ! Success or failure?

     ! Initialize
    RC = GC_SUCCESS

    IF ((.NOT.plume_inject_on).OR.(.NOT.use_lagrange)) THEN
      WRITE(6,'(a)') ' No plume injection / no lagrangian module configuration, skip cleanup'
      RETURN
    ENDIF

    IF ( ALLOCATED(RXNRATE_CONST_KPP))  THEN
          DEALLOCATE(RXNRATE_CONST_KPP, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:RXNRATE_CONST_KPP', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
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
    !WRITE(6,'(a)') '--> Lagrange and Plume Module Cleanup <--'
!
!
  END SUBROUTINE plume_mod_cleanup_box


  SUBROUTINE plume_physics(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
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

  !TYPE(SpcConc), POINTER        :: Spc(:)
  !TYPE(Plume2d_list), POINTER :: Plume2d_head, Plume2d_tail
          
  INTEGER                :: i_box, i_lon, i_lat, i_lev
  INTEGER                :: ii, jj, kk, i, j
  INTEGER                :: ki
  INTEGER                :: i_species, id_tracer
  INTEGER                :: i_x, i_y
  INTEGER                :: i_slab
  INTEGER                :: Nt
  INTEGER                :: t1s
  INTEGER                :: OrigUnit

  REAL(fp)               :: MW_g
  REAL(fp)               :: Dt
  REAL(fp)               :: ratio
  real(fp)               :: V_prev, V_new
  real(fp)               :: grid_volume,  V_grid_2D
  real(fp)               :: Ly ! Lyapunov exponent [s-1]
  real(fp)               :: length0
  real(fp)               :: Pdx, Pdy, Pdt
  REAL(fp)               :: box_lon, box_lat, box_lev                          
  real(fp)               :: box_length, box_alpha, box_theta
  real(fp)               :: box_life, box_label
  REAL(fp)               :: box_x_PS, box_y_PS
  real(fp)               :: box_u, box_v, box_omeg
  REAL(fp)               :: dbox_lon, dbox_lat, dbox_lev
  REAL(fp)               :: dbox_x_PS, dbox_y_PS
  REAL(fp)               :: curr_lon, curr_lat, curr_pressure
  real(fp)               :: curr_u, curr_v, curr_omeg, curr_ptemp
  real(fp)               :: curr_u_PS, curr_v_PS
  real(fp)               :: curr_T1, next_T2 
  REAL(fp)               :: RK_Dx_PS, RK_Dy_PS
  REAL(fp)               :: RK_x_PS, RK_y_PS
  REAL(fp)               :: RK_lon, RK_lat
  real(fp)               :: wind_s_shear, Ptemp_shear
  real(fp)               :: Cv, Ch, Omega_N, N_BV
  real(fp)               :: eddy_v, eddy_h 
  real(fp)               :: CFL
  real(fp)               :: Pc_middle, Pc_bottom, Pc_top, Pc_left, Pc_right
  real(fp)               :: background_conc
  real(fp)               :: mass_plume, mass_plume_new, D_mass_plume, mass_plume_scale
  real(fp)               :: background_mass, background_mass_new
  real(fp)               :: excess_mass
  ! real(fp)               :: Massref_species
  REAL(fp)               :: RK_Dt(5)
  REAL(fp)               :: RK_u(4), RK_v(4), RK_omeg(4)
  REAL(fp)               :: RK_Dlon(4), RK_Dlat(4), RK_Dlev(4)
  real(fp)               :: Pu(n_x_max2,n_y_max2) ! for 2D advection
  !real(fp)               :: Pc(n_x_max,n_y_max), Pc2(n_x_max,n_y_max) !, Ec(n_x_max,n_y_max)
  !real(fp)               :: Pc_bdy(n_x_max2,n_y_max2)
  real(fp)               :: C2d_prev(n_x_max,n_y_max) !, C2d_prev_extra(n_x_max,n_y_max)
  real(fp)               :: C2d_new(n_x_max2,n_y_max2)
  real(fp)               :: C2d_bg(n_x_max2,n_y_max2)
  !eal(fp)               :: Concnt2D_bdy(n_x_max2, n_y_max2)
  

  REAL(fp), dimension(:,:,:), allocatable :: box_concnt_2D
  !real(fp), dimension(:,:), allocatable :: box_concnt_1D
  REAL(fp), POINTER      :: u(:,:,:)
  REAL(fp), POINTER      :: v(:,:,:)
  REAL(fp), POINTER      :: omeg(:,:,:)
  REAL(fp), POINTER      :: Ptemp(:,:,:)
  REAL(fp), POINTER      :: T1(:,:,:)
  REAL(fp), POINTER      :: T2(:,:,:)
  real(fp), pointer      :: P_BXHEIGHT(:,:,:)

  TYPE(SpcConc), POINTER        :: Spc(:)
  TYPE(Plume2d_list), POINTER :: Plume2d_new, Plume2d_curr, Plume2d_prev
  !TYPE(Plume1d_list), POINTER :: Plume1d_new, Plume1d_curr, Plume1d_prev


  
  CHARACTER(LEN=255)     :: ErrMsg
  CHARACTER(LEN=255)     :: ThisLoc
  CHARACTER(LEN=255)     :: spc_name
  

  !CHARACTER(LEN=63)      :: OrigUnit
  
! unused vars -------
  !INTEGER                :: nAdv
  !REAL(fp), POINTER      :: PASV_EU 
  !real(fp)  :: eddy_A, eddy_B 
  !real(fp)  :: CFL_2d(n_x_max2,n_y_max2) ! for 2D advection  
  !real(fp)  :: Xscale, Yscale, frac_mass
  !real(fp)  :: D_concnt(n_x_max,n_y_max)
  !real(fp)  :: Cslab(n_slab_max) !, Extra_Cslab(n_slab_max)
  !real(fp) :: lon1, lon2, lat1, lat2
  !real(fp) :: box_lon_edge, box_lat_edge
  !real(fp) :: D_wind, D_x, D_y
  !real(fp)  :: start, finish


  Spc                    =>   State_Chm%Species
  
  Dt = GET_TS_DYN()
  ThisLoc                =   ' -> at plume_physics (in module GeosCore/lagrange_singlebox_mod.F90)'
  RC     =  GC_SUCCESS
  ErrMsg = ''
  NULLIFY(Plume2d_new, Plume2d_curr, Plume2d_prev)
  !IF(Stop_inject==1) GOTO 400 ! deallocate and nullify -> exit
  
  ALLOCATE(box_concnt_2D(n_x_max, n_y_max, n_species )) !

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
  P_BXHEIGHT => State_Met%BXHEIGHT  ![IIPAR,JJPAR,KKPAR]
  !=======================================================================
  ! for 2D plume: Run Lagrangian trajectory-track HERE
  !=======================================================================

  IF(.NOT.ASSOCIATED(Plume2d_head)) GOTO 401
  Plume2d_curr => Plume2d_head
  i_box = 0
  
  DO WHILE(ASSOCIATED(Plume2d_curr))
   
    !i_box         = Plume2d_curr%label
    box_lon       = Plume2d_curr%LON
    box_lat       = Plume2d_curr%LAT
    box_lev       = Plume2d_curr%LEV
    i_lon         = Plume2d_curr%lon_ind
    i_lat         = Plume2d_curr%lat_ind
    i_lev         = Plume2d_curr%lev_ind

    box_length    = Plume2d_curr%LENGTH
    box_alpha     = Plume2d_curr%ALPHA
    box_label     = Plume2d_curr%label
    box_life      = Plume2d_curr%LIFE

    Pdx           = Plume2d_curr%PDX
    Pdy           = Plume2d_curr%PDY

    box_concnt_2D = Plume2d_curr%CONCNT2d

    write(6,*) 'debug (BZ): solve plume physics in plume box: ', box_label
    
    box_life = box_life + Dt

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

    ! update index based on new location
    i_lon = Find_iLonLat(box_lon, DX, X_edge2)
    if(i_lon>IIPAR) i_lon=i_lon-IIPAR
    if(i_lon<1) i_lon=i_lon+IIPAR

    i_lat = Find_iLonLat(box_lat, DY, Y_edge2)
    if(i_lat>JJPAR) i_lat=JJPAR
    if(i_lat<1) i_lat=1

    i_lev = Find_iPLev(box_lev,P_edge)
    if(i_lev>LLPAR) i_lev=LLPAR


    if(abs(box_lat)>Y_mid(JJPAR))then
        next_T2 = Interplt_wind_RLL_polar(T2, i_lon, i_lat, i_lev, box_lon, box_lat, box_lev)
    else
        next_T2 = Interplt_wind_RLL(T2, i_lon, i_lat, i_lev, box_lon, box_lat, box_lev)
    endif

    ! PV=nRT, V2 = T2/P2 : T1/P1 * V1
    ratio = ( (next_T2/box_lev)/(curr_T1/curr_pressure) )**(1/2)

    ! assume the volume change mainly apply to the cross-section, 
    ! the box_length would not change 

    V_prev = Pdx*Pdy*box_length*1.0e+6_fp

    Pdx = Pdx *ratio
    Pdy = Pdy *ratio

    V_new = Pdx*Pdy*box_length*1.0e+6_fp

    DO i_species = 1, n_species, 1
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
    IF((box_u**2+box_v**2)==0)then
        box_alpha = 0.0
    ELSE
      IF(box_v>=0)THEN
        box_alpha = ACOS( box_u/SQRT(box_u**2+box_v**2) ) 
      ELSE
        box_alpha = 2*PI - ACOS( box_u/SQRT(box_u**2+box_v**2) )
      ENDIF
    ENDIF

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
    ! update     
    ! ----------------------------------------------------------------
    Plume2d_curr%IsNew   = 0
    
    Plume2d_curr%LON     = box_lon
    Plume2d_curr%LAT     = box_lat
    Plume2d_curr%LEV     = box_lev
    Plume2d_curr%lat_ind = i_lat
    Plume2d_curr%lon_ind = i_lon
    Plume2d_curr%lev_ind = i_lev

    Plume2d_curr%LENGTH  = box_length
    Plume2d_curr%ALPHA   = box_alpha
    Plume2d_curr%label   = box_label
    Plume2d_curr%LIFE    = box_life


    Plume2d_curr%PDX = Pdx
    Plume2d_curr%PDY = Pdy

    Plume2d_curr%CONCNT2d    = box_concnt_2D

    !-------- solve in-plume concentration here
    ! advection and diffusion
    
    !box_lon          = Plume2d_curr%LON
    !box_lat          = Plume2d_curr%LAT
    !box_lev          = Plume2d_curr%LEV
    !box_length       = Plume2d_curr%LENGTH
    !box_alpha        = Plume2d_curr%ALPHA
    !box_label        = Plume2d_curr%label
    !box_life         = Plume2d_curr%LIFE
    !Pdx              = Plume2d_curr%PDX
    !Pdy              = Plume2d_curr%PDY
    !box_concnt_2D    = Plume2d_curr%CONCNT2d
    !i_lat            = Plume2d_curr%lat_ind
    !i_lon            = Plume2d_curr%lon_ind
    !i_lev            = Plume2d_curr%lev_ind

    curr_lon         = box_lon
    curr_lat         = box_lat
    curr_pressure    = box_lev      ! hPa
    
    grid_volume     = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ! [cm3]
    V_grid_2D       = Pdx*Pdy*box_length*1.0e+6_fp ! [cm3]
    
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
    IF(abs(curr_lat)>Y_mid(JJPAR))then
      curr_Ptemp = Interplt_wind_RLL_polar(Ptemp, i_lon, i_lat, &
                              i_lev, curr_lon, curr_lat, curr_pressure)
    ELSE
      curr_Ptemp = Interplt_wind_RLL(Ptemp, i_lon, i_lat, i_lev, &
                                      curr_lon, curr_lat, curr_pressure)
    ENDIF


    N_BV = SQRT(Ptemp_shear*g0/curr_Ptemp)
    IF(N_BV<=0.001) N_BV = 0.001

    ! diffusivity unit: [m2/s]
    eddy_v = Cv * Omega_N**2 / N_BV
    eddy_h = 10.0


    ! Define the wind field based on wind shear
    DO i_y = 1, n_y_max2
      Pu(:,i_y) = (i_y-n_y_mid2)*Pdy * wind_s_shear ! [m s-1]
    ENDDO
    !--------------------------------------------------------------
    ! if Pdx become smaller than half Dx_init,
    ! combine 9 grids into 1 grids. 
    ! Update: box_concnt_2D(), Pdx(), Pdy(),  Extra_mass_2D()
    ! (BZ): The logic here need to be clarified, temporiraly disable
    ! Resize and throw an error when not meet CFL condition
    !--------------------------------------------------------------
    IF(Pdx<=0.5*Dx_init) THEN
      errMsg = 'Size are too small to meet the CFL condition! '
          !CALL GC_Error( errMsg, RC, thisLoc )
      CALL ERROR_STOP( errMsg, thisLoc)
    ENDIF
!    IF(Pdx<=0.5*Dx_init)THEN
!
!600     CONTINUE
!      DO i_species = 1, n_species
!        Pc = box_concnt_2D(:,:,i_species) ! [molec cm-3]
!        !box_concnt_2D(:,:,i_species) = 0.0
!        ! Should fill with background concentration
!        box_concnt_2D(:,:,i_species) = Spc(i_species)%Conc(i_lon, i_lat, i_lev)
!
!        DO i = 1, n_x_max/3, 1
!        DO j = 1, n_y_max/3, 1
!
!          i_x = (i-1)*3+1
!          i_y = (j-1)*3+1
!
!          box_concnt_2D(i+n_x_max/3,j+n_y_max/3,i_species) = &
!                                    SUM(Pc(i_x:i_x+2,i_y:i_y+2))/9
!
!        ENDDO
!        ENDDO
!      ENDDO ! DO i_species = 1, n_species
!
!      Pdx = Pdx*3
!      Pdy = Pdy*3
!
!
!
!      IF(Pdx<=0.5*Dx_init) WRITE(6,*) &
!          "ERROR 0: more combination: ", Plume2d_curr%label, Pdx, Pdy
!
!      IF(Pdx<=0.5*Dx_init) GOTO 600
!
!
!    ENDIF ! IF(Pdx(i_box)<0.5*Dx_init)THEN
    !-------------------------------------------------------------------
       ! Calculate the advection-diffusion in 2D grids
       !-------------------------------------------------------------------
       !V_grid_2D       = Pdx*Pdy*box_length*1.0e+6_fp


        DO i_species= 1, n_species, 1
          
          ! Massref_species=Plume2d_curr%MassRef2d(i_species)
          background_conc = Spc(i_species)%Conc(i_lon, i_lat, i_lev)
          C2d_prev(1:n_x_max,1:n_y_max) = &
                                  box_concnt_2D(1:n_x_max,1:n_y_max,i_species)

          Nt = CEILING(Dt/120)
          Pdt = Dt/Nt ! Dt=600 ! FLOOR(Pdx/Pu(1,1)/10)*10

          ! Find the best Pdt to meet CFL condition:
 700    CONTINUE


          CFL = Pdt*Pu(1,1)/Pdx
          IF(MAX( ABS(CFL), ABS(2*eddy_h*Pdt/(Pdx**2)), &
                            ABS(2*eddy_v*Pdt/(Pdy**2)) ) > 0.8)THEN

            Nt = Nt+1
            Pdt = Dt/Nt
            GOTO 700

          ENDIF

          !Concnt2D_bdy(:,:) = 0.0 
          ! Maybe fill the unused outer cells with background concentration
          C2d_bg(:,:)  = background_conc
          C2d_bg(2:n_x_max2-1,2:n_y_max2-1) =C2d_prev(1:n_x_max,1:n_y_max)
          C2d_new (:,:) = C2d_bg (:,:)
        
          IF( abs(Dt/Pdt-Nt) > 0.00001 ) THEN
            WRITE(6,*) "*** ERROR: Check Pdt ***"
            WRITE(6,*) Pdt, Nt, Dt
            WRITE(6,*) Pdx/Pu(1,1), Pdy**2/(2*eddy_v), Pdx**2/(2*eddy_h)
          ENDIF
       

          DO t1s = 1, NINT(Dt/Pdt)

            ! advection ----------------------------------------------------
            ! diffusion ----------------------------------------------------

            !Pc_bdy(:,:) = Spc(i_species)%Conc(i_lon, i_lat, i_lev)
            !Pc_bdy(2:n_x_max2-1,2:n_y_max2-1) = Concnt2D_bdy(2:n_x_max2-1,2:n_y_max2-1)
            C2d_bg(:,:)  = background_conc
            C2d_bg(2:n_x_max2-1,2:n_y_max2-1) =C2d_new(2:n_x_max2-1,2:n_y_max2-1) 
            ! Only calculate the vertical half 2D domain         

            !$OMP PARALLEL DO           &
            !$OMP DEFAULT( SHARED     ) &
            !$OMP PRIVATE(i_y,i_x,CFL,Pc_middle,Pc_top,Pc_bottom,Pc_right,Pc_left)
            DO i_y = 2, n_y_mid2, 1
            DO i_x = 2, n_x_max2-1, 1
              Pc_middle = C2d_bg( i_x,   i_y  )
              Pc_top    = C2d_bg( i_x,   i_y+1)
              Pc_bottom = C2d_bg( i_x,   i_y-1)
              Pc_right  = C2d_bg( i_x+1, i_y  )
              Pc_left   = C2d_bg( i_x-1, i_y  )
           
              CFL       = Pdt*Pu(i_x,i_y)/Pdx

              C2d_new(i_x, i_y) = Pc_middle           &
                - 0.5 * CFL    * ( Pc_right - Pc_left )   &
                + 0.5 * CFL**2 * ( Pc_right - 2*Pc_middle + Pc_left )         &
                + Pdt*( eddy_h*( Pc_right -2*Pc_middle +Pc_left   ) /(Pdx**2) &
                       +eddy_v*( Pc_top   -2*Pc_middle +Pc_bottom ) /(Pdy**2) )
              ! update the other half based on vertical symmetry 
              C2d_new(i_x, n_y_max2+1-i_y) = C2d_new(i_x, i_y) 

            ENDDO
            ENDDO
            !$OMP END PARALLEL DO

          ! update the other half based on vertical symmetry 
            !DO i_y = n_y_mid2+1, n_y_max2-1, 1
            !  Concnt2D_bdy(2:n_x_max2-1:1, i_y) = Concnt2D_bdy(n_x_max2-1:2:-1,n_y_max2+1-i_y)
            !ENDDO

          ENDDO ! DO t1s = 1, NINT(Dt/Pdt)


          

          
         !================================================================
         ! Calculate the mass exchange of plume to background cell
         ! Update the concentration in the background and plume
         ! accordingly
         !================================================================

          ! the boundary always represents the background concentration
          mass_plume      = V_grid_2D * SUM(C2d_prev(:,:))! molec
          mass_plume_new  = V_grid_2D * SUM(C2d_new(2:n_x_max2-1,2:n_y_max2-1))
          !mass_plume_edge = V_grid_2D * ( SUM(C2d_new(2, 2:n_y_max2-1)) + &
          !                          SUM(C2d_new(n_x_max2-1, 2:n_y_max2-1)) +   &
          !                          SUM(C2d_new(3:n_x_max2-2 , 2))   + &
          !                          SUM(C2d_new(3:n_x_max2-2 , n_y_max2-1)) )

          D_mass_plume    = mass_plume_new - mass_plume
          background_mass     = background_conc * grid_volume
          excess_mass = D_mass_plume - background_mass

          ! If mass need to enter the plume significantly larger than background mass, decrease the 
          ! mass entered, by scaling the whole plume conc
          ! If background = 0 but D_MASS_PLUME > (Maybe later change to a lower threshold)
          IF ((background_mass <= 0.0_fp).AND. (D_mass_plume.GT. 0.0_fp)) THEN
            ! just print this value to see how small it would be
            !ErrMsg = "(Debug: BZ) Mass enter the plume but background is 0, Plume num: ", &
            !        Plume2d_curr%label, '; Species: ', i_species, 'D_mass_plume: ', D_mass_plume
            !GC_Warning( ErrMsg, RC, ThisLoc )
            WRITE (6, *) "(Debug: BZ) Mass enter the plume but background is 0, Plume num: ", &
                    Plume2d_curr%label, '; Species: ', i_species, 'D_mass_plume: ', D_mass_plume

          ELSEIF (excess_mass .GT. 1.0e-2_fp * background_mass) THEN
            ! Later might need to adjust the species that enter the plume to ensure mass conservation
            !ErrMsg = "(Debug: BZ) Mass enter the plume larger than mass in background, Plume num: ", &
            !        Plume2d_curr%label, '; Species: ', i_species, 'D_mass_plume: ', D_mass_plume, &
            !        'background_mass: ', background_mass
            !GC_Warning( ErrMsg, RC, ThisLoc )
            WRITE (6, *) "(Debug: BZ) Mass enter the plume larger than mass in background, Plume num: ", &
                    Plume2d_curr%label, '; Species: ', i_species, 'D_mass_plume: ', D_mass_plume, &
                    'background_mass: ', background_mass
            ! mass_plume_new      = mass_plume + background_mass
            !mass_edge_scale     = (mass_plume_edge - excess_mass) / mass_plume_edge 
            !C2d_new(2:n_x_max2-1,2:n_y_max2-1) = C2d_new(2:n_x_max2-1,2:n_y_max2-1) * &
            !        (mass_plume + background_mass) / mass_plume_new
          ENDIF

          background_mass_new = MAX( 0.0_fp, background_mass - D_mass_plume )

          !i_advect = id_PASV_LA +i_species -1

          !backgrd_concnt = State_Chm%Species(i_species)%Conc(i_lon,i_lat,i_lev)

          !backgrd_concnt = ( backgrd_concnt*grid_volume &
                                          !- D_mass_plume) /grid_volume
          box_concnt_2D(:,:,i_species) = C2d_new(2:n_x_max2-1,2:n_y_max2-1)

          
          Spc(i_species)%Conc(i_lon,i_lat,i_lev) = background_mass_new / grid_volume

        ENDDO ! DO i_species=1,n_species,1
    Plume2d_curr%CONCNT2d    = box_concnt_2D
    Plume2d_curr => Plume2d_curr%next
  ENDDO  ! DO WHILE(ASSOCIATED(Plume2d))


401 CONTINUE
  !=======================================================================
  ! For 1d plume: Run Lagrangian trajectory-track HERE
  !=======================================================================
  IF(.NOT.ASSOCIATED(Plume1d_head)) GOTO 400

400 CONTINUE
  
  !------------------------------------------------------------------
  ! Everything is done, clean up pointers
  !------------------------------------------------------------------
  ! deallocate unused space
  IF(allocated(box_concnt_2D)) deallocate(box_concnt_2D)
  !IF(allocated(box_concnt_1D)) deallocate(box_concnt_1D)
  
  ! Nullify pointers
  IF(ASSOCIATED(u)) nullify(u)
  IF(ASSOCIATED(v)) nullify(v)
  IF(ASSOCIATED(omeg)) nullify(omeg)
  IF(ASSOCIATED(Ptemp)) nullify(Ptemp)
  IF(ASSOCIATED(T1)) nullify(T1)
  IF(ASSOCIATED(T2)) nullify(T2)
  IF(ASSOCIATED(P_BXHEIGHT)) nullify(P_BXHEIGHT)
  IF(ASSOCIATED(Spc)) nullify(Spc)
  !IF(ASSOCIATED(X_edge)) nullify(X_edge)
  !IF(ASSOCIATED(Y_edge)) nullify(Y_edge)

  !IF(ASSOCIATED(PASV_EU)) nullify(PASV_EU)

  IF(ASSOCIATED(Plume2d_new)) nullify(Plume2d_new)
  IF(ASSOCIATED(Plume2d_curr)) nullify(Plume2d_curr)
  IF(ASSOCIATED(Plume2d_prev)) nullify(Plume2d_prev)
  !IF(ASSOCIATED(Plume1d_new)) nullify(Plume1d_new)
  !IF(ASSOCIATED(Plume1d_curr)) nullify(Plume1d_curr)
  !IF(ASSOCIATED(Plume1d_prev)) nullify(Plume1d_prev)

  END SUBROUTINE plume_physics

  SUBROUTINE plume_chem_microphysics(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
    USE Input_Opt_Mod,            ONLY : OptInput, PlumeSource_t
    USE State_Chm_Mod,            ONLY : ChmState, Ind_
    USE State_Met_Mod,            ONLY : MetState
    USE Species_Mod,              ONLY : SpcConc, Species
    USE TIME_MOD,                 ONLY : GET_TS_DYN
    USE State_Grid_Mod,           ONLY : GrdState
    USE UnitConv_Mod
    !USE State_Diag_Mod,           ONLY : DgnState
    !USE State_Diag_Mod,           ONLY : DgnMap
    ! KPP-related module
#ifdef KPP_INTEGRATOR_AUTOREDUCE
    USE fullchem_AutoReduceFuncs, ONLY : fullchem_AR_KeepHalogensActive
    USE fullchem_AutoReduceFuncs, ONLY : fullchem_AR_SetKeepActive
    USE fullchem_AutoReduceFuncs, ONLY : fullchem_AR_UpdateKppDiags
    USE fullchem_AutoReduceFuncs, ONLY : fullchem_AR_SetIntegratorOptions
#endif
    USE GcKpp_Global
    USE GcKpp_Parameters
    USE Gckpp_Monitor,            ONLY : SPC_NAMES, Eqn_Names, Fam_Names
    USE GcKpp_Rates,              ONLY : UPDATE_RCONST, RCONST
    USE GcKpp_Integrator,         ONLY : Integrate
    USE GcKpp_Function
    ! GEOS-Chem related module
    USE UCX_MOD,                  ONLY : SO4_PHOTFRAC

    LOGICAL, INTENT(IN)           :: am_I_Root
    TYPE(MetState), INTENT(IN)    :: State_Met
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State objectgg
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    !TYPE(DgnState), INTENT(INOUT) :: State_Diag ! Diagnostics State object
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    
    
    INTEGER                       :: i_box, i_lon, i_lat, i_lev
    INTEGER                       :: i_species, i_phot, i_kpp, i_rxn
    INTEGER                       :: i_x, i_y
    INTEGER                       :: n_species 
    INTEGER                       :: Thread, IERR,  P, F, errorCount
    INTEGER                       :: SpcID, KppID
    INTEGER                       :: RXN_O3_1, RXN_O3_2
    !INTEGER                       :: ind_SO2, ind_SO4, ind_OH
    
    INTEGER                       :: ISTATUS(20)
    INTEGER                       :: ICNTRL (20)


    !REAL(fp)                      :: Dt
    REAL(fp)                      :: SO4_FRAC,   SR,        LWC
    REAL(dp)                      :: RCNTRL (20)
    REAL(dp)                      :: RSTATE (20)
    REAL(dp)                      :: C_before_integrate(NSPEC)
    REAL(dp)                      :: local_RCONST(NREACT)
    REAL(fp)                      :: H2SO4_RATE_2d(n_x_max,n_y_max) ! H2SO4 prod rate [kg s-1]
    REAL(fp)                      :: PSO4AQ_RATE_2d(n_x_max,n_y_max) ! Cld chem sulfate prod rate [kg s-1]
    
    LOGICAL                       :: Failed2x,  Size_Res, doSuppress

    CHARACTER(LEN=255)            :: ErrMsg
    CHARACTER(LEN=255)            :: ThisLoc

    TYPE(SpcConc), POINTER        :: Spc(:)
    TYPE(Species), POINTER        :: SpcInfo

    TYPE(Plume2d_list), POINTER :: Plume2d_next, Plume2d_curr, Plume2d_prev
    
    ! SAVEd scalars
    LOGICAL,  SAVE         :: FIRSTCHEM = .False.
#ifdef MODEL_CLASSIC
#ifndef NO_OMP
    INTEGER, EXTERNAL      :: OMP_GET_THREAD_NUM
#endif
#endif
    ! Parameters below copied from GEOS-Chem fullchem_mod.F90
    ! Defines the slot in which the H-value from the KPP integrator is stored
    ! This should be the same as the value of Nhnew in gckpp_Integrator.F90
    ! Define this locally in order to break a compile-time dependency.
    !    -- Bob Yantosca (05 May 2022)
    INTEGER,     PARAMETER :: Nhnew = 3

    ! Add Nhexit, the last timestep length -- Obin Sturm (30 April 2024)
    INTEGER,     PARAMETER :: Nhexit = 2

    ! Suppress printing out KPP error messages after this many errors occur
    INTEGER,     PARAMETER :: INTEGRATE_FAIL_TOGGLE = 20
    
    !========================================================================
    ! plume_chem_microphysics begins here!
    ! Mainly adapted from GEOS-Chem fullchem_mod.F90
    ! Currently only solve KPP-related gas phase chemistry
    ! Currently not consider other processes listed in chemistry_mod.F90, including:
    ! 1) UCX: Calc_Strat_Aer
    ! No solid particle formation inside, liquid fraction taken from Eulerian grid process
    ! 2) Aerosol_Mod: Aerosol_Conc, RdAer
    ! 3) Dust_Mod: RDust_Online (Dust OD)
    ! aerosol AOD related calculation
    ! aerosol surface area for heteorogenous chem
    ! 4) ChemSulfate -> chem_SO2:
    ! In-cloud chemistry and cloud pH 
    ! e.g., sulfate loss by H2O2, O3, HOBr, HCHO,...
    ! 5) photolysis rate modification 
    ! all taken from Eulerian grid process diagnostic
    !========================================================================

    ! Initialization
    NULLIFY(Plume2d_next, Plume2d_curr, Plume2d_prev)

    Spc                    =>   State_Chm%Species
    SpcInfo                =>   NULL()
    Plume2d_curr           =>   Plume2d_head
    ThisLoc                =    ' -> at plume_chem_microphysics (in module GeosCore/lagrange_singlebox_mod.F90)'
    RC                     =    GC_SUCCESS
    ErrMsg                 =    ''
    i_box                  =    0
    n_species              =    State_Chm%nSpecies
    Thread                 =    1
    errorCount             =    0
    Failed2x               =   .FALSE.
    doSuppress             =   .FALSE.
    id_SO2                 =    Ind_('SO2')
    id_SO4                 =    Ind_('SO4')
    id_OH                  =    Ind_('OH')
    H2SO4_RATE_2d          =    0.0d0
    PSO4AQ_RATE_2d         =    0.0d0
    ! RXN_O3_1 specifies: O3 + hv -> O2 + O
    ! RXN_O3_2 specifies: O3 + hv -> O2 + O(1D)
    ! (BZ): For debug purpose
    RXN_O3_1              = State_Chm%Phot%RXN_O3_1
    RXN_O3_2              = State_Chm%Phot%RXN_O3_2
    ! Noted that in GEOS-Chem, the default chemistry timestep is 20min, dynamic time step is 10min
    ! Here we implement chemsitry timestep using 10min
    Dt                     =    GET_TS_DYN()
    

    ! Set up integration convergence conditions and timesteps
    ! This is defined in gckpp_Global and set to be public
    ATOL = State_Chm%KPP_AbsTol   ! Absolute tolerance
    RTOL = State_Chm%KPP_RelTol   ! Relative tolerance

    ! IF (State_Diag%Archive_RxnConst        ) Write(6, *) "Debug: (BZ) ; Archive_RxnRate", State_Diag%Archive_RxnRate
    ! For debug process, print rate constant 202: SO2 + OH {+M} = SO4 + HO2 + PH2SO4 :
    Write (6, *) "Debug: (BZ): In Plume  (Before Plume Chem): rate constant for RXN 202 = ", &
      RXNRATE_CONST_KPP(23, 40, 39 ,202)
    DO WHILE(ASSOCIATED(Plume2d_curr))
      i_box = i_box+1
      i_lon         = Plume2d_curr%lon_ind
      i_lat         = Plume2d_curr%lat_ind
      i_lev         = Plume2d_curr%lev_ind
      ! Test if we need to do the chemistry for box (I,J,L),
      ! otherwise move onto the next box.
      ! MaxChemLev = MaxStratLev = 59 
      ! Hard coded in GeosUtil/gc_grid_mod.F90
      ! IF ( .not. State_Met%InChemGrid(i_lon,i_lat,i_lev) ) CYCLE
      IF ( .not. State_Met%InChemGrid(i_lon,i_lat,i_lev) ) THEN
        WRITE(6,*) 'Debug (BZ): Outside chem grid: (i_box, i_lon, i_lat, i_lev): ',    &
        Plume2d_curr%label, i_lon, i_lat, i_lev
        plume2d_curr => plume2d_curr%next
        ! BZ: Maybe if plume reach out of chem grid, directly release species to Eulerian grid?
        CYCLE
      ENDIF
      WRITE(6,*) 'Debug (BZ): Euleria grid: Conc of SO2 ', Spc(id_SO2)%Conc(i_lon,i_lat,i_lev)
      WRITE(6,*) 'Debug (BZ): Euleria grid: Conc of SO4 ', Spc(id_SO4)%Conc(i_lon,i_lat,i_lev)
      WRITE(6,*) 'Debug (BZ): Euleria grid: Conc of OH ', Spc(id_OH)%Conc(i_lon,i_lat,i_lev)

      WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc before Chem: Conc of SO2 ', SUM(Plume2d_curr%CONCNT2d(:,:,id_SO2))/(n_x_max*n_y_max)
      WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc before Chem: Conc of SO4 ', SUM(Plume2d_curr%CONCNT2d(:,:,id_SO4))/(n_x_max*n_y_max)
      WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc before Chem: Conc of OH ', SUM(Plume2d_curr%CONCNT2d(:,:,id_OH))/(n_x_max*n_y_max)

      WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc before Chem: Conc of SO2 ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_SO2)
      WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc before Chem: Conc of SO4 ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_SO4)
      WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc before Chem: Conc of OH ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_OH)

      ! Some species (counter or diagnostic species) needed to be zero out before solving KPP chemistry
      ! See details in Do_Chemistry in Fullchem_mod.F90
      DO i_species = 1, n_species
        ! Get info about this species from the species database
        SpcInfo      => State_Chm%SpcData(i_species)%Info

        ! isoprene oxidation counter species
        IF ( TRIM( SpcInfo%Name ) == 'LISOPOH' .or. &
              TRIM( SpcInfo%Name ) == 'LISOPNO3' ) THEN
            Plume2d_curr%CONCNT2d(:,:,i_species) = 0.0_fp
        ENDIF

       ! aromatic oxidation counter species
        IF ( Input_Opt%LSOA .or. Input_Opt%LSVPOA ) THEN
            SELECT CASE ( TRIM( SpcInfo%Name ) )
              CASE ( 'LBRO2H', 'LBRO2N', 'LTRO2H', 'LTRO2N', &
                      'LXRO2H', 'LXRO2N', 'LNRO2H', 'LNRO2N' )
                  Plume2d_curr%CONCNT2d(:,:,i_species) = 0.0_fp
            END SELECT
        ENDIF

        ! Sulfate gas/cloud prod diagnostic species
        IF ( TRIM( SpcInfo%Name ) == 'PH2SO4' .or. &
              TRIM( SpcInfo%Name ) == 'PSO4AQ' ) THEN
            Plume2d_curr%CONCNT2d(:,:,i_species) = 0.0_fp
        ENDIF

        ! Free pointer
        SpcInfo => NULL()
      ENDDO

      !!! (BZ) Test output: Will Photolysis rate array be initalized every dynamic timestep?
      ! State_Chm%Phot%ZPJ(i_lev,n_photoRxn,i_lat,i_lon)
      WRITE(6,*) 'Debug (BZ): Do_Chemistry  (In plume): PhotoRxn: O3 + hv -> O2 + O; rate: ', State_Chm%Phot%ZPJ(39,RXN_O3_1,23,40)
      WRITE(6,*) 'Debug (BZ): Do_Chemistry  (In plume): PhotoRxn:  O3 + hv -> O2 + O(1D); rate: ', State_Chm%Phot%ZPJ(39,RXN_O3_2,23,40)

      !========================================================================
      ! MAIN LOOP: Compute reaction rates and call chemical solver
      !========================================================================
      
      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE(i_y,i_x, Thread)&
      !$OMP PRIVATE( SO4_FRAC, IERR,     RCNTRL,  ISTATUS,   RSTATE     )&
      !$OMP PRIVATE( ICNTRL,   C_before_integrate )&
      !$OMP PRIVATE( SpcID,    KppID,    F,       P            )&
      !$OMP PRIVATE( SR               )&
      !$OMP PRIVATE( SIZE_RES, LWC                                            )&
      !$OMP COLLAPSE( 2                    )&
      !$OMP SCHEDULE( DYNAMIC, 24          )&
      !$OMP REDUCTION( +:errorCount        )
      DO i_y = 1, n_y_max, 1
        DO i_x = 1, n_x_max, 1
          ! Skip to the end of the loop if we have failed integration twice
          IF ( Failed2x ) CYCLE
           ! Initialize private loop variables for each (i_x, i_y)
          IERR      = 0                        ! KPP success or failure flag
          ISTATUS   = 0.0_dp                   ! Rosenbrock output
          ICNTRL    = 0                        ! Rosenbrock input (integer)
          RCNTRL    = 0.0_fp                   ! Rosenbrock input (real)
          RSTATE    = 0.0_dp                   ! Rosenbrock output
          SO4_FRAC  = 0.0_fp                   ! Frac of SO4 avail for photolysis
          P         = 0                        ! GEOS-Chem photolyis species ID
          !LCH4      = 0.0_fp                   ! P/L diag: Methane loss rate
          !PCO_TOT   = 0.0_fp                   ! P/L diag: Total P(CO)
          !PCO_CH4   = 0.0_fp                   ! P/L diag: P(CO) from CH4
          !PCO_NMVOC = 0.0_fp                   ! P/L diag: P(CO) from NMVOC
          SR        = 0.0_fp                   ! Enhancement to O2 catalysis rate
          LWC       = 0.0_fp                   ! Liquid water content
          SIZE_RES  = .FALSE.                  ! Size resolved calculation?
          C         = 0.0_dp                   ! KPP species conc's
          RCONST    = 0.0_dp                   ! KPP rate constants
          PHOTOL    = 0.0_dp                   ! Photolysis array for KPP
          K_CLD     = 0.0_dp                   ! Sulfur in-cloud rxn het rates
          K_MT      = 0.0_dp                   ! Sulfur sea salt rxn het rates
          CFACTOR   = 1.0_dp                   ! KPP conversion factor
          SRO3      = 0.0_dp                   ! Enhanced sulfate production of
          SRHOBr    = 0.0_dp                   !  O3, HOBr, HCl in size-resolved
          SRHOCl    = 0.0_dp                   !  cloud droplets
#ifdef MODEL_CLASSIC
#ifndef NO_OMP
          Thread    = OMP_GET_THREAD_NUM() + 1 ! OpenMP thread number
#endif
#endif
#ifdef KPP_INTEGRATOR_AUTOREDUCE
       ! Per discussions for Lin et al., force keepActive throughout the
       ! atmosphere if keepActive option is enabled. (hplin, 2/9/22)
       CALL fullchem_AR_SetKeepActive( option=.TRUE. )
#endif
          ! Get photolysis rates (daytime only)
          ! Update SUNCOSmid threshold from 0 to cos(98 degrees)
          ! Loop over the FAST-JX photolysis species
          IF ( State_Met%SUNCOSmid(i_lon,i_lat) > -0.1391731e+0_fp ) THEN
            ! Only proceed if doing photolysis
            IF ( Input_Opt%Do_Photolysis ) THEN
              DO i_phot = 1, State_Chm%Phot%nMaxPhotRxns

                ! Copy photolysis rate from FAST_JX into KPP PHOTOL array
                PHOTOL(i_phot) = State_Chm%Phot%ZPJ(i_lev,i_phot,i_lon,i_lat)
                ! In GEOS-Chem, maybe useful to archieve instantaneous photolysis rate [s -1] and noon time photolysis rate [s -1]
                ! The mapping between the GEOS-Chem photolysis species and
                ! the FAST-JX photolysis species is contained in the lookup
                ! table in input file FJX_j2j.dat.
                ! Some GEOS-Chem photolysis species may have multiple
                ! branches for photolysis reactions.  These will be
                ! represented by multiple entries in the FJX_j2j.dat
                ! lookup table.
                !    NOTE: For convenience, we have stored the GEOS-Chem
                !    photolysis species index (range: 1..State_Chm%nPhotol)
                !    for each of the FAST-JX photolysis species (range;
                !    1..State_Chm%Phot%nMaxPhotRxns) in the GC_PHOTO_ID array

                ! GC photolysis species index
                P = State_Chm%Phot%GC_Photo_Id(i_phot)
                ! Below is diagnostic steps and might be ignored for now
                ! If this FAST_JX photolysis species maps to a valid
                ! GEOS-Chem photolysis species (for this simulation)...
                !IF ( P > 0 .and. P <= State_Chm%nPhotol ) THEN
                  !
                !ELSE IF ( P == State_Chm%nPhotol+1 ) THEN
                  ! J(O3_O1D).  This used to be stored as the nPhotol+1st
                  ! diagnostic in Jval, but needed to be broken off
                  ! to facilitate cleaner diagnostic indexing (bmy, 6/3/20)
                !ELSE IF ( P == State_Chm%nPhotol+2 ) THEN
                  ! J(O3_O3P).  This used to be stored as the nPhotol+2nd
                  ! diagnostic in Jval, but needed to be broken off
                  ! to facilitate cleaner diagnostic indexing (bmy, 6/3/20)
                !ENDIF

              ENDDO
            ENDIF
          ENDIF
          
          ! Initialize the KPP "C" vector of species concentrations [molec/cm3]
          DO i_species = 1, NSPEC
            SpcID = State_Chm%Map_KppSpc(i_species)
            C(i_species)  = 0.0_dp
            IF ( SpcId > 0 ) C(i_species) =  Plume2d_curr%CONCNT2d(i_x,i_y,i_species)
          ENDDO

          !=====================================================================
          ! CHEMISTRY MECHANISM INITIALIZATION (#1)
          !
          ! Populate KPP global variables and arrays in gckpp_global.F90
          !
          ! NOTE: This has to be done before Set_Sulfur_Chem_Rates, so that
          ! the NUMDEN and SR_TEMP KPP variables will be populated first.
          ! Otherwise this can lead to differences in output that are evident
          ! when running with different numbers of OpenMP cores.
          ! See https://github.com/geoschem/geos-chem/issues/1157
          !    -- Bob Yantosca (08 Mar 2022)
          !=====================================================================
          ! Copy values into the various KPP global variables
          CALL Set_inPlume_2d_Kpp_GridBox_Values( I_EU          = i_lon,                          &
                                        J_EU          = i_lat,                          &
                                        L_EU         = i_lev,                          &
                                        I_LA         =   i_x,                     &
                                        J_LA        = i_y,                       &
                                        Input_Opt  = Input_Opt,                  &
                                        State_Chm  = State_Chm,                  &
                                        State_Grid = State_Grid,                 &
                                        State_Met  = State_Met,                  &
                                        RC         = RC                         )
          !=====================================================================
          ! CHEMISTRY MECHANISM INITIALIZATION (#2)
          !
          ! Update reaction rates [1/s] for sulfur chemistry in cloud and on
          ! seasalt.  These will be passed to the KPP chemical solver.
          !
          ! NOTE: This has to be done before fullchem_SetStateHet so that
          ! State_Chm%HSO3_aq and State_Chm%SO3_aq will be populated first.
          ! These are copied into State_Het%HSO3_aq and State_Het%SO3_aq.
          ! See https://github.com/geoschem/geos-chem/issues/1157
          !    -- Bob Yantosca (08 Mar 2022)
          !=====================================================================
          ! Compute sulfur chemistry reaction rates [1/s]
          ! If size_res = T, we'll call fullchem_HetDropChem below.
          ! Could remove the use of State_Diag
          ! defines the variables State_Chm%HSO3_aq and State_Chm%SO3aq
          ! Therefore, we must call Set_Sulfur_Chem_Rates after
          !  Set_KPP_GridBox_Values, but before fullchem_SetStateHet.  Otherwise we
          !  will not be able to copy State_Chm%HSO3_aq to State_Het%HSO3_aq and
          !  State_Chm%SO3_aq to State_Het%SO3_aq properly.
          !CALL Set_Sulfur_Chem_Rates( I          = i_lon,                           &
          !                            J          = i_lat,                           &
          !                            L          = i_lev,                           &
          !                            Input_Opt  = Input_Opt,                   &
          !                            State_Chm  = State_Chm,                   &
          !                            State_Diag = State_Diag,                  &
          !                            State_Grid = State_Grid,                  &
          !                            State_Met  = State_Met,                   &
          !                            size_res   = size_res,                    &
          !                            RC         = RC                          )

          !=====================================================================
          ! CHEMISTRY MECHANISM INITIALIZATION (#3)
          !
          ! Populate the various fields of the State_Het object.
          !
          ! NOTE: This has to be done after fullchem_SetStateHet so that
          ! State_Chm%HSO3_aq and State_Chm%SO3_aq will be populated first.
          ! These are copied into State_Het%HSO3_aq and State_Het%SO3_aq.
          ! See https://github.com/geoschem/geos-chem/issues/1157
          !    -- Bob Yantosca (08 Mar 2022)
          !=====================================================================

          ! Populate fields of the State_Het object
          ! These values are used in the heterogenous chemistry reaction rate computations.
          ! Ignore heteogeneous reaction for now and complete later
          !CALL fullchem_SetStateHet( I         = I,                             &
          !                            J         = J,                             &
          !                            L         = L,                             &
          !                            id_SALA   = id_SALA,                       &
          !                            id_SALAAL = id_SALAAL,                     &
          !                            id_SALC   = id_SALC,                       &
          !                            id_SALCAL = id_SALCAL,                     &
          !                            Input_Opt = Input_Opt,                     &
          !                            State_Chm = State_Chm,                     &
          !                            State_Met = State_Met,                     &
          !                            H         = State_Het,                     &
          !                            RC        = RC                            )

          !=====================================================================
          ! (ignore for now) CHEMISTRY MECHANISM INITIALIZATION (#5)
          !
          ! Call Het_Drop_Chem (formerly located in sulfate_mod.F90) to
          ! estimate the in-cloud sulfate production rate in heterogeneous
          ! cloud droplets based on the Yuen et al., 1996 parameterization.
          ! Code by Becky Alexander (2011) with updates by Mike Long and Bob
          ! Yantosca (2021).
          !
          ! We will only call Het_Drop_Chem if:
          ! (1) It is at least 0.01% cloudy
          ! (2) We are doing a size-resolved computation
          ! (3) The grid box is over water
          ! (4) The temperature is above -5C
          ! (5) Liquid water content is nonzero
          !=====================================================================

          !=====================================================================
          ! Prepare arrays
          !=====================================================================

          ! Zero out dummy species index in KPP
          !! Since we did not initialize PL_Kpp_Id, need to define dummy species from another way

          DO i_kpp = 1, NFAM
              KppID = Ind_( TRIM ( Fam_Names(i_kpp) ), 'K' )
              ! Exit if an invalid ID is encountered
              !IF ( KppId <= 0 ) THEN
              !    ErrMsg = 'Invalid KPP ID for prod/loss species: '            // &
              !         TRIM( Fam_Names(i_kpp) )
              !    CALL GC_Error( ErrMsg, RC, ThisLoc )
              !    RETURN
              !ENDIF
              IF ( KppID > 0 ) C(KppID) = 0.0_dp
          ENDDO

          !=====================================================================
          ! Update reaction rates
          !=====================================================================

          ! Update the array of rate constants
          ! Mannually Set K_CLoud and K_MT = 0?
          ! Read from last timestep
          ! CALL Update_RCONST()
          
          !=====================================================================
          ! HISTORY (aka netCDF diagnostics)
          !
          ! Archive KPP reaction rates [molec cm-3 s-1]
          ! See gckpp_Monitor.F90 for a list of chemical reactions
          !=====================================================================

#ifdef KPP_INTEGRATOR_AUTOREDUCE
          !=====================================================================
          ! This is currently an empty function, might be discarded?
          ! Set options for the KPP integrator in vectors ICNTRL and RCNTRL
          ! This now needs to be done within the parallel loop
          !=====================================================================
          CALL fullchem_AR_SetIntegratorOptions( Input_Opt, State_Chm,          &
                                                 State_Met, FirstChem,          &
                                                 i_lon,    i_lat,  i_lev,       &
                                                 ICNTRL,    RCNTRL             )
#endif
          !=====================================================================
          ! Integrate the box forwards
          !=====================================================================
          C_before_integrate = C
          CALL Integrate( 0.0_dp, DT, ICNTRL, RCNTRL, ISTATUS, RSTATE, IERR )
          ! Print grid box indices to screen if integrate failed
          IF ( IERR < 0 ) THEN

              ! Turn off error output after a certain limit is reached
              IF ( .not. doSuppress ) THEN
                WRITE( 6, * ) '### INTEGRATE RETURNED ERROR AT: ', i_x, i_y, Plume2d_curr%label
                errorCount = errorCount + 1
                IF ( errorCount > INTEGRATE_FAIL_TOGGLE ) THEN
                    WRITE( 6, '(a)' ) &
                      '### Further error output has been switched off'
                    doSuppress = .TRUE.
                ENDIF
              ENDIF
            ENDIF
          !=====================================================================
          ! HISTORY: Archive KPP solver diagnostics
          !=====================================================================
          
          !=====================================================================
          ! Try another time if it failed
          !=====================================================================
          IF ( IERR < 0 ) THEN
            ! Zero the first time step (Hstart).  Also reset C with
            ! concentrations prior to the 1st call to "Integrate".
            RCNTRL(3) = 0.0_dp
            C         = C_before_integrate

            ! Disable auto-reduce solver for the second iteration for safety
            IF ( Input_Opt%Use_AutoReduce ) THEN
              RCNTRL(12) = -1.0_dp ! without using ICNTRL
            ENDIF

            ! Update rates again
            CALL Update_RCONST()

            ! Call the integrator
            ! NOTE: Some integrators (like LSODE) will overwrite the TIN value
            ! upon exit.  To prevent this, pass 0.0, DT as the 1st 2 arguments.
            CALL Integrate( 0.0_dp, DT, ICNTRL, RCNTRL, ISTATUS, RSTATE, IERR )
            
            !==================================================================
            ! Exit upon the second failure
            !==================================================================
            IF ( IERR < 0 ) THEN
              ! Print error message
              WRITE(6,     '(a   )' ) '## INTEGRATE FAILED TWICE !!! '
              WRITE(ERRMSG,'(a,i3)' ) 'Integrator error code :', IERR
             
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             ! Make sure only one thread at a time executes this block
             !$OMP CRITICAL
             !
             ! Set a flag to break out of loop gracefully
             ! NOTE: You can set a GDB breakpoint here to examine the error
             Failed2x = .TRUE.

             ! Print concentrations at failure grid box
             PRINT*, REPEAT( '#', 79 )
             PRINT*, '### KPP DEBUG OUTPUT!'
             PRINT*, '### Species concentrations at problem box ',i_lon, i_lat, Plume2d_curr%label
             PRINT*, REPEAT( '#', 79 )
             DO i_species = 1, n_species
                PRINT*, C(i_species), TRIM( ADJUSTL( SPC_NAMES(i_species) ) )
             ENDDO

             ! Print rate constants at failure grid box
             PRINT*, REPEAT( '#', 79 )
             PRINT*, '### KPP DEBUG OUTPUT!'
             PRINT*, '### Reaction rates at problem box ', i_lon, i_lat, Plume2d_curr%label
             PRINT*, REPEAT( '#', 79 )
             DO i_rxn = 1, NREACT
                PRINT*, RCONST(i_rxn), TRIM( ADJUSTL( EQN_NAMES(i_rxn) ) )
             ENDDO
             !
             !$OMP END CRITICAL
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             ! Start skipping to end of loop upon 2 failures in a row
             CYCLE
            ENDIF
          ENDIF

          !=====================================================================
          ! Check we have no negative values and copy the concentrations
          ! calculated from the C array back into Plume conc array
          !=====================================================================

          ! Loop over KPP species
          DO i_species = 1, NSPEC

              ! GEOS-Chem species ID
              SpcID = State_Chm%Map_KppSpc(i_species)

              ! Skip if this is not a GEOS-Chem species
              IF ( SpcID <= 0 ) CYCLE


              ! Set negative concentrations to zero
              C(i_species) = MAX( C(i_species), 0.0_dp )

              ! Copy concentrations back into State_Chm%Species
               Plume2d_curr%CONCNT2d(i_x,i_y,i_species) = REAL( C(i_species), kind=fp )

          ENDDO
#ifdef TOMAS
              !-----------------------------------------------------------------
              ! FOR TOMAS MICROPHYSICS:
              !
              ! Obtain P/L with a unit [kg S] for tracing
              ! gas-phase sulfur species production (SO2, SO4, MSA)
              ! (win, 8/4/09)
              !
              ! TODO: Abstract this to a subroutine, to simplify DO_FULLCHEM
              !-----------------------------------------------------------------
              ! Calculate H2SO4 production rate [kg s-1] in each
              ! time step (win, 8/4/09)
              H2SO4_RATE_2d(i_x, i_y)= C(ind_PH2SO4) / AVO * 98.e-3_fp * &
                           State_Met%AIRVOL(i_lon,i_lat,i_lev)    * &
                           1.0e+6_fp / DT  ! kg s-1 box-1
        
              IF ( H2SO4_RATE_2d(i_x, i_y) < 0.0d0) THEN
                !ErrMsg = "H2SO4_RATE_2D negative in (Plumeid, x, y):", &
                !    Plume2d_curr%label, i_x, i_y, "was:", H2SO4_RATE_2d(i_x, i_y), "  setting to 0.0d0"
                !CALL GC_Warning( ErrMsg, RC, ThisLoc )
                WRITE(6, *) "H2SO4_RATE_2D negative in (Plumeid, x, y):", &
                    Plume2d_curr%label, i_x, i_y, "was:", H2SO4_RATE_2d(i_x, i_y), "  setting to 0.0d0"
                H2SO4_RATE_2d(i_x, i_y) = 0.0d0
              ENDIF

              PSO4AQ_RATE_2d(i_x, i_y) = C(ind_PSO4AQ) / AVO * 98.e-3_fp * &
                            State_Met%AIRVOL(i_lon,i_lat,i_lev)    * &
                            1.0e+6_fp ! kg per timestep box-1

              IF ( PSO4AQ_RATE_2d(i_x, i_y) < 0.0d0) THEN
                !ErrMsg = "PSO4AQ_RATE_2D negative in (Plumeid, x, y):", &
                !    Plume2d_curr%label, i_x, i_y, "was:", PSO4AQ_RATE_2d(i_x, i_y), "  setting to 0.0d0"
                
                !CALL GC_Warning( ErrMsg, RC, ThisLoc )
                WRITE(6, *) "PSO4AQ_RATE_2D negative in (Plumeid, x, y):", &
                      Plume2d_curr%label, i_x, i_y, "was:", PSO4AQ_RATE_2d(i_x, i_y), "  setting to 0.0d0"
                PSO4AQ_RATE_2d(i_x, i_y) = 0.0d0
              ENDIF
#endif
              !====================================================================
              ! Archieve KPP diagnostic output and write diagnostics files
              ! e.g., Chemical production and loss
              ! OH reactivity: inverse of its life-time
              ! GcKpp_Util -> Get_OHreactivity
              !====================================================================
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
      
      WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc after Chem: Conc of SO2 ', SUM(Plume2d_curr%CONCNT2d(:,:,id_SO2))/(n_x_max*n_y_max)
      WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc after Chem: Conc of SO4 ', SUM(Plume2d_curr%CONCNT2d(:,:,id_SO4))/(n_x_max*n_y_max)
      WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc after Chem: Conc of OH ', SUM(Plume2d_curr%CONCNT2d(:,:,id_OH))/(n_x_max*n_y_max)

      WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc after Chem: Conc of SO2 ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_SO2)
      WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc after Chem: Conc of SO4 ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_SO4)
      WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc after Chem: Conc of OH ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_OH)
      !=======================================================================
      ! Return gracefully if integration failed 2x anywhere
      ! (as we cannot break out of a parallel DO loop!)
      !=======================================================================
      IF ( Failed2x ) THEN
        ErrMsg = 'KPP failed to converge after 2 iterations!'
        CALL GC_Error( ErrMsg, RC, ThisLoc )
        RETURN
      ENDIF
#ifdef TOMAS
      !-----------------------------------------------------------------
      ! TOMAS microphysics:
      ! [Currently not consider]: additional process: in-cloud oxidation: 
      ! SO4 production from aqueous chemistry of SO2 >> PSO4AQ_RATE_2d 
      ! Distributed onto size-resolved aerosol num and sulfate mass >> 
      !    fullchem_mod -> TOMAS_SO4_AQ >> TOMAS_MOD -> AQOXID
      !-----------------------------------------------------------------
       
#endif
      !-----------------------------------------------------------------
      ! [Currently not considered]: Sea salt chemistry: ChemSeaSalt <<SEASALT_MOD
      ! - SALA, SALACL, SALAAL wet settling
      ! - SALC, SALCCL, SALCAL wet settling
      ! - Marine organic aerosol MOPO-> MOPI, e-folding time 1.15 days: CHEM_MOPO AND CHEM_MOPI <<SEASALT_MOD
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! [Currently not considered]: Recalculate PSC properties: Calc_Strat_Aer<< UCX_MOD
      ! In non PSC formation regime, calculate liquid phase aerosol (SLA)
      ! Based on partitioning of total SO4 as H2SO4, 
      ! CALL TERNARY( PCENTER,TCENTER,H2OSUM,H2SO4_BOX_L, &
		  !			   0.e+0_fp   ,HClSUM,HOClSUM,HBrSUM,HOBrSUM, &
			!		   W_H2SO4,W_H2O,W_HNO3,W_HCl,W_HOCl,W_HBr,W_HOBr, &
			!		   HNO3GASFRAC,HClGASFRAC,HOClGASFRAC, &
			!		   HBrGASFRAC,HOBrGASFRAC,VOL_SLA,RHO_AER_BOX)
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      ! [Currently not considered]: sulfate chemistry: ChemSulfate << SULFATE_MOD
      ! SO4s [kg] gravitational settling
      ! NITs [kg] gravitational settling
      ! Stratospheric aerosol gravitational settling: SETTLE_STRAT_AER << UCX
      !    - Settling SLAs
      !      - change concentration of SLAs including: SO4, HNO3, HCl, HOCl, HBr, HOBr, H2O, 
      !      - Settle some BCPI
      !    - Settling SPAs
      !      - NIT, H2O, and aerosol mass
      ! SO2 chemistry: CHEM_SO2 << sulfate_mod
      !    - SO2 production:
      !      - DMS + OH, DMS + NO3 (saved in CHEM_DMS)
      !      - HMS -> SO2 + HCHO (aq)
      !    - SO2 loss:
      !      - SO2 + OH  -> SO4 [This is for aerosol-only simulation, otherwise is included in KPP]
      !      - SO2       -> drydep [This is done in mixing_mod.F90]
      !      - Sea salt alkalinity and sea salt-sulfate reaction [Not included in KPP]
      !           condition: alkalinity>0, SO2 present and excess O3 present
      !      - Acid uptake on dust particles: 
      !           condition: Alkalinity >0, SO2 present, and O3 excess
      !      - SO2 in cloud chemistry (with cloud, with SO2, and T>-15 C, with LWC)
      !      - Metal catalyzed oxidation of SO2 pathway (on dust)
      !      - SO2 loss by H2O2 (not included in KPP)
      !      - SO2 loss by O3 (not included in KPP) 
      !           - Pathway 1: L3S incloud liquid phase?
      !           - Pathway 2: L3S_1 in cloud solid phase?          
      !      - SO2 loss by HOBr (not included in KPP)
      !      - SO2 + HCHO (aq)-> HMS
      !      - SO2 + HMS -> 2 SO4
      !    - SO2 = SO2_0 * exp(-bt) +  PSO2_DMS/bt * [1-exp(-bt)]
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! [Currently not considered]: Carbon chemistry: ChemCarbon << CARBON_MOD
      ! BCPO/ECOB -> BCPI/ECIL and OCPO/OCOB -> OCPI/OCIL:   e-folding time 1.15 days
      ! SOAP -> SOAS, for TOMAS: SOA condensation (COACOND << TOMAS_MOD)
      ! SOA chemistry SOA_CHEMISTRY << CARBON_MOD
      Plume2d_curr => Plume2d_curr%next
    ENDDO
    !Write (6, *) "Debug: (BZ): In Plume  (After Plume Chem): rate constant for RXN 202 = ", &
    !    State_Diag%RxnConst(23, 40, 39 ,202)
  END SUBROUTINE plume_chem_microphysics

  SUBROUTINE plume_structure_change(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
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

    TYPE(SpcConc), POINTER        :: Spc(:)

    INTEGER                       :: i_box, i_lon, i_lat, i_lev
    INTEGER                       :: i_species, n_species
    INTEGER                       :: Stop_loop

    real(fp)                      :: box_lon, box_lat, box_lev
    real(fp)                      :: box_length, box_alpha, box_theta
    real(fp)                      :: box_extra, box_life, box_label
    real(fp)                      :: Pdx, Pdy
    real(fp)                      :: grid_volume, V_grid_2D
    real(fp)                      :: conc_background, mass_release
    REAL(fp)                      :: Dt

    CHARACTER(LEN=255)            :: ErrMsg
    CHARACTER(LEN=255)            :: ThisLoc

    real(fp), dimension(:,:,:), allocatable :: box_concnt_2D

    TYPE(Plume2d_list), POINTER :: Plume2d_new, Plume2d_curr, Plume2d_prev

    Spc                    =>  State_Chm%Species
    n_species              =   State_Chm%nSpecies
    Dt                     =   GET_TS_DYN()
    ThisLoc                =   ' -> at plume_structure_change (in module GeosCore/lagrange_singlebox_mod.F90)'
    RC                     =   GC_SUCCESS
    ErrMsg                 =   ''
    NULLIFY(Plume2d_new, Plume2d_curr, Plume2d_prev)

    !ALLOCATE(box_concnt_2D(n_x_max, n_y_max, n_species))
    ! dissolve 2D plume seg 
    ! lifetime larger than 1-Day
    IF(.NOT.ASSOCIATED(Plume2d_head)) GOTO 401
    Plume2d_curr => Plume2d_head
    !NULLIFY(Plume2d_prev)
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

      Pdx        = Plume2d_curr%PDX
      Pdy        = Plume2d_curr%PDY

      !box_concnt_2D = Plume2d_curr%CONCNT2d
      i_lon         = Plume2d_curr%lon_ind
      i_lat         = Plume2d_curr%lat_ind
      i_lev         = Plume2d_curr%lev_ind

      grid_volume     = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ! [cm3]
      V_grid_2D       = Pdx*Pdy*box_length*1.0e+6_fp
      ! --------------------------------------------------------------------
      ! delete the node for 2D plume (current criteria: plume lifetime > 1days)
      ! At the end of simulation, automatically dissolve all plumes
      ! --------------------------------------------------------------------
      IF((box_life .GT. 1.0*24.0*60.0*60).OR. (ITS_TIME_FOR_EXIT())) THEN
    
      ! delete the only node
      IF(.NOT.ASSOCIATED(Plume2d_curr%next) .AND. &
                      .NOT.ASSOCIATED(Plume2d_prev))THEN
        !i_box = i_box-1
        !Stop_loop = 1
        ! Release mass in Plume2d_head in Eulerian grid
        DO i_species = 1, n_species
          !conc_background = Spc(i_species)%Conc(i_lon, i_lat, i_lev)
          !mass_release = (SUM(Plume2d_curr%CONCNT2d(:,:,i_species)) -      &
          !             conc_background *n_x_max * n_y_max ) * V_grid_2D
          mass_release = SUM(Plume2d_curr%CONCNT2d(:,:,i_species))  * V_grid_2D - &
                      Plume2d_curr%MassRef2d(i_species)
          Spc(i_species)%Conc(i_lon, i_lat, i_lev) =       &
                         MAX(0.0_fp, Spc(i_species)%Conc(i_lon, i_lat, i_lev) + &
                         mass_release / grid_volume)
        ENDDO

        WRITE(6,*)'                '
        WRITE(6,*)'debug (BZ): grid_volume=', grid_volume, '(i_lat, i_lon, i_lev) = (', i_lat, i_lon, i_lev,')'
		WRITE(6,*)'*** Deleting the last Plume2d segment: Num_Plume2d = ', Num_Plume2d, 'label = ', Plume2d_curr%label
        WRITE(6,*)'                '
        IF (ASSOCIATED(Plume2d_curr%CONCNT2d)) DEALLOCATE(Plume2d_curr%CONCNT2d)
        DEALLOCATE(Plume2d_curr)
        NULLIFY(Plume2d_prev)
        NULLIFY(Plume2d_head)
        NULLIFY(Plume2d_tail)
        Num_Plume2d = Num_Plume2d - 1

        GOTO  401! add corresponding 1-D seg

      ! delete the tail node 
      ELSEIF(.NOT.ASSOCIATED(Plume2d_curr%next))THEN 
        WRITE(6,*)'debug (BZ): grid_volume=', grid_volume, '(i_lat, i_lon, i_lev) = (', i_lat, i_lon, i_lev,')'
		WRITE(6,*)'*** Deleting tail node: Num_Plume2d = ', Num_Plume2d, 'label = ', Plume2d_curr%label
        ! Release mass in Plume2d_curr in Eulerian grid
        DO i_species = 1, n_species
          !conc_background = Spc(i_species)%Conc(i_lon, i_lat, i_lev)
          !mass_release = (SUM(Plume2d_curr%CONCNT2d(:,:,i_species)) -      &
          !             conc_background *n_x_max * n_y_max ) * V_grid_2D
          mass_release = SUM(Plume2d_curr%CONCNT2d(:,:,i_species))  * V_grid_2D - &
                      Plume2d_curr%MassRef2d(i_species) 
          Spc(i_species)%Conc(i_lon, i_lat, i_lev) =       &
                         MAX(0.0_fp, Spc(i_species)%Conc(i_lon, i_lat, i_lev) + &
                         mass_release / grid_volume )
        ENDDO

          Plume2d_tail => Plume2d_prev
          NULLIFY(Plume2d_tail%next)
          IF (ASSOCIATED(Plume2d_curr%CONCNT2d)) DEALLOCATE(Plume2d_curr%CONCNT2d)
          DEALLOCATE(Plume2d_curr)
          i_box = i_box-1
          Num_Plume2d = Num_Plume2d - 1

          Stop_loop = 1
          GOTO 401
          
        ! delete the head node
      ELSEIF (.NOT. ASSOCIATED(Plume2d_prev)) THEN
        Plume2d_head => Plume2d_curr%next
        !Plume2d_prev => Plume2d_curr%next
        ! Release mass in Plume2d_curr in Eulerian grid
        WRITE(6,*)'debug (BZ): grid_volume=', grid_volume, '(i_lat, i_lon, i_lev) = (', i_lat, i_lon, i_lev,')'
		WRITE(6,*)'*** Deleting head node: Num_Plume2d = ', Num_Plume2d, 'label = ', Plume2d_curr%label
        DO i_species = 1, n_species
          !conc_background = Spc(i_species)%Conc(i_lon, i_lat, i_lev)
          !mass_release = (SUM(Plume2d_curr%CONCNT2d(:,:,i_species)) -      &
          !             conc_background *n_x_max * n_y_max ) * V_grid_2D
          mass_release = SUM(Plume2d_curr%CONCNT2d(:,:,i_species))  * V_grid_2D - &
                      Plume2d_curr%MassRef2d(i_species) 
           Spc(i_species)%Conc(i_lon, i_lat, i_lev) =       &
                         MAX(0.0_fp, Spc(i_species)%Conc(i_lon, i_lat, i_lev) + &
                         mass_release / grid_volume )
        ENDDO
        IF (ASSOCIATED(Plume2d_curr%CONCNT2d)) DEALLOCATE(Plume2d_curr%CONCNT2d)
        DEALLOCATE(Plume2d_curr)
        Plume2d_curr => Plume2d_head
        Num_Plume2d = Num_Plume2d - 1
        ! delete middle node
      ELSE
          Plume2d_prev%next => Plume2d_curr%next
          ! Release mass in Plume2d_curr in Eulerian grid
          WRITE(6,*)'debug (BZ): grid_volume=', grid_volume, '(i_lat, i_lon, i_lev) = (', i_lat, i_lon, i_lev,')'
          WRITE(6,*)'*** Deleting middle node: Num_Plume2d = ', Num_Plume2d, 'label = ', Plume2d_curr%label
          DO i_species = 1, n_species
            !conc_background = Spc(i_species)%Conc(i_lon, i_lat, i_lev)
            !mass_release = (SUM(Plume2d_curr%CONCNT2d(:,:,i_species)) -      &
            !            conc_background *n_x_max * n_y_max ) * V_grid_2D
            mass_release = SUM(Plume2d_curr%CONCNT2d(:,:,i_species))  * V_grid_2D - &
                      Plume2d_curr%MassRef2d(i_species) 
            Spc(i_species)%Conc(i_lon, i_lat, i_lev) =       &
                          MAX(0.0_fp, Spc(i_species)%Conc(i_lon, i_lat, i_lev) + &
                          mass_release / grid_volume)
          ENDDO
          IF (ASSOCIATED(Plume2d_curr%CONCNT2d)) DEALLOCATE(Plume2d_curr%CONCNT2d)
          DEALLOCATE(Plume2d_curr)
          Plume2d_curr => Plume2d_prev%next
          Num_Plume2d = Num_Plume2d - 1
      ENDIF

    ELSE
      ! IF deletion occur, no need to move Plume_prev and Plume_curr
      Plume2d_prev => Plume2d_curr
      Plume2d_curr => Plume2d_curr%next 
    ENDIF
    
    ENDDO
401 CONTINUE
    !IF(.NOT.ASSOCIATED(Plume1d_head)) GOTO 400

400 CONTINUE
    ! Cleanup pointer
    !IF (ALLOCATED(box_concnt_2D)) DEALLOCATE(box_concnt_2D)

  END SUBROUTINE plume_structure_change
    

  SUBROUTINE lagrange_write_std(am_I_Root, RC)

    LOGICAL,          INTENT(IN   ) :: am_I_Root   ! root CPU?
    INTEGER,          INTENT(INOUT) :: RC          ! Failure or success

  END SUBROUTINE lagrange_write_std
!------------------------------------------------------------------
! all functions listed below
!------------------------------------------------------------------
Real(fp) FUNCTION GetInjectionLat(lat1, lat2, plume_length, inject_num) 
    !------------------------------------------------------------------
    ! Get lattitude of injection at specific time step
    ! This make sure plume seg created from lat1:dlat:lat2 -> lat2:-dlat,lat1 -> ...
    !------------------------------------------------------------------

    implicit none
    REAL(fp), INTENT(IN) :: lat1, lat2
    !REAL(fp), INTENT(IN) :: plane_speed  !m/s
    REAL(fp), INTENT(IN) :: plume_length !m
    !REAL(fp), INTENT(IN) :: dt    ! s 
    INTEGER,  INTENT(IN) :: inject_num ! number of injection before this cycle 
    
    REAL(fp)             :: plane_route_deg, plume_length_deg
    REAL(fp)             :: lat_span_deg
    REAL(fp)             :: route_pos_deg
    INTEGER              :: pass_index


    plume_length_deg      =  plume_length / 1000.0 / 111.0
    lat_span_deg   = lat2 - lat1
    plane_route_deg = real(inject_num) * plume_length_deg
    pass_index       = int(plane_route_deg / lat_span_deg)
    route_pos_deg    = mod(plane_route_deg, lat_span_deg)

    if (mod(pass_index,2) == 0) then
      GetInjectionLat = lat1 + route_pos_deg
    else
      GetInjectionLat = lat2 - route_pos_deg
    endif
    
    
    ! old function below
    !INTEGER              :: n_total, cycle_length
    !INTEGER              :: pos, idx

    !n_total = INT((lat2 - lat1)/dlat) + 1
    !cycle_length = 2*n_total - 2

    !pos = MOD(global_id-1, cycle_length)

    !IF (pos < n_total) THEN
    !    idx = pos
    !ELSE
    !    idx = cycle_length - pos
    !END IF

    !GetInjectionLat = lat1 + idx*dlat

END FUNCTION GetInjectionLat

  REAL(fp) FUNCTION Calc_Ly(u, v, i_lon, i_lat, i_lev, box_alpha, box_lon, box_lat)
    !------------------------------------------------------------------
    ! calcualte the Lyaponov exponent (Ly), unit: s-1
    !------------------------------------------------------------------

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

  real(fp) function Interplt_wind_RLL(wind, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    !------------------------------------------------------------------
    ! functions to interpolate wind speed (u,v,omeg) 
    ! based on the surrounding 4 points.
    !------------------------------------------------------------------

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
  end function  Interplt_wind_RLL
  real(fp) function Interplt_wind_RLL_polar(wind, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
  !-------------------------------------------------------------------------------
  ! functions to interpolate vertical wind speed (w)
  ! based on the surrounding 3 points, one of the points is the north/south polar
  ! point. The w value at polar point is the average of all surrounding grid points.
  !-------------------------------------------------------------------------------
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
  end function Interplt_wind_RLL_polar

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
  end function Interplt_uv_PS_polar
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
  end function Interplt_uv_PS

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

    INTEGER               :: axis
    REAL(fp)              :: Pdx, Pdy, frac
    REAL(fp)              :: concnt1_2D(n_x_max, n_y_max)

    INTEGER               :: i, j
    INTEGER               :: N_max, N_frac
    INTEGER               :: alloc_stat
    REAL(fp)              :: temp, C_sum, mass_total
    REAL(fp)              :: D_len
    REAL(fp), allocatable :: concnt1_2D_sum(:)
    
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

SUBROUTINE Set_inPlume_2d_Kpp_GridBox_Values( I_EU,J_EU, L_EU, I_LA, J_LA, Input_Opt, State_Chm, State_Grid, State_Met, RC)

  ! Temporarily set all environmental variables = Euleria grid, later consider variable interpolation within plume grid
  USE ErrCode_Mod
  USE GcKpp_Global
  USE GcKpp_Parameters
  USE Input_Opt_Mod,          ONLY : OptInput
  USE PhysConstants,          ONLY : CONSVAP, RGASLATM, RSTARG
  USE Pressure_Mod,           ONLY : Get_Pcenter
  USE State_Chm_Mod,          ONLY : ChmState
  USE State_Grid_Mod,         ONLY : GrdState
  USE State_Met_Mod,          ONLY : MetState

  INTEGER,        INTENT(IN)  :: I_EU, J_EU, L_EU, I_LA, J_LA
  TYPE(OptInput), INTENT(IN)  :: Input_Opt
  TYPE(ChmState), INTENT(IN)  :: State_Chm
  TYPE(GrdState), INTENT(IN)  :: State_Grid
  TYPE(MetState), INTENT(IN)  :: State_Met

  INTEGER,        INTENT(OUT) :: RC

  ! LOCAL VARIABLES:
  INTEGER            :: F,       N,        NA,    KppId,    SpcId
  REAL(f8)           :: CONSEXP, VPRESH2O
  ! Strings
  CHARACTER(LEN=255) :: ErrMsg, ThisLoc

  ! Initialization
  RC      = GC_SUCCESS
  NA      = State_Chm%nAeroType
  ErrMsg  = ''
  ThisLoc = &
    ' -> at Set_inPlume_2d_Kpp_GridBox_Values (in module GeosCore/lagrange_singlebox_mod.F90)'
  
  !========================================================================
  ! Populate global variables in gckpp_Global.F90
  !========================================================================
  ! Solar quantities
  SUNCOS          = State_Met%SUNCOSmid(I_EU,J_EU)

  ! Pressure and density quantities
  NUMDEN          = State_Met%AIRNUMDEN(I_EU,J_EU,L_EU)
  H2O             = State_Met%AVGW(I_EU,J_EU,L_EU) * NUMDEN
  PRESS           = Get_Pcenter( I_EU, J_EU, L_EU )

  ! Temperature quantities
  TEMP            = State_Met%T(I_EU,J_EU,L_EU)
  INV_TEMP        = 1.0_dp   / TEMP
  TEMP_OVER_K300  = TEMP     / 300.0_dp
  K300_OVER_TEMP  = 300.0_dp / TEMP
  SR_TEMP         = SQRT( TEMP )
  FOUR_R_T        = 4.0_dp * CON_R    * TEMP
  FOUR_RGASLATM_T = 4.0_dp * RGASLATM * TEMP
  EIGHT_RSTARG_T  = 8.0_dp * RSTARG   * TEMP

  ! Relative humidity quantities
  CONSEXP         = 17.2693882_dp * (TEMP - 273.16_dp) / (TEMP - 35.86_dp)
  VPRESH2O        = CONSVAP * EXP( CONSEXP ) / TEMP
  RELHUM          = ( H2O / VPRESH2O ) * 100_dp 
END SUBROUTINE Set_inPlume_2d_Kpp_GridBox_Values
END MODULE Lagrange_singlebox_Mod
