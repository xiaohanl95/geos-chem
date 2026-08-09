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
! April 10, 2026
! Note:
! One possible model logic:
! The PiG state stores only concentration anomaly relative to the current
! Eulerian background. This is intentional. Some species are strongly coupled
! to the full chemical mechanism but USE ERROR_MOare not intended to evolve as an
! independent long-memory plume reservoir. Therefore chemistry is solved on
! reconstructed absolute concentrations, and the result is then converted back
! to anomaly with respect to the updated background state.
! Chemistry update in anomaly form:
!
!   C_plume_abs_before = C_bg_before + C_anom_before
!   C_plume_abs_after  = KPP( C_plume_abs_before )
!   C_bg_after         = KPP( C_bg_before )
!   C_anom_after       = C_plume_abs_after - C_bg_after
!
! Only C_anom is stored in the plume state.
! April 15, 2026
! Another possible model logic:
! 

! August 4, 2026:
! (To do) Could store reaction constant needed only instead of full 4-D array to save space
! RXNRATE_CONST_KPP(NX_GC,NY_GC,NZ_GC,NREACT)
! Accumulator for chemistry, OH, HO2, NH3, NH4, can save space by storing active grid cells
! TOMAS module
! NK_TOT threshold 1e-5 for doing COND_NUC, need to check if this threshold is reasonable
MODULE Lagrange_singlebox_Mod
  USE Plume_list_Mod
  USE PRECISION_MOD
  USE ERROR_MOD
  USE ERRCODE_MOD
  USE PhysConstants,   ONLY : PI, Re, g0, AIRMW, AVO, BOLTZ
  !USE TIME_MOD,        ONLY : GET_YEAR, GET_MONTH, GET_DAY, GET_HOUR, GET_MINUTE, GET_SECOND
  !USE TIME_MOD,        ONLY : ITS_TIME_FOR_EXIT, GET_TS_DYN, GET_TAU, GET_TAUb
  !USE TIME_MOD
  USE UNITCONV_MOD    
  USE INPUT_OPT_MOD,   ONLY : PlumeSource_t
  ! USE gckpp_Parameters, ONLY: NREACT, NSPEC

  IMPLICIT NONE
  PRIVATE

  !PUBLIC MEMBER FUNCTIONS
  PUBLIC :: lagrange_init_box
  PUBLIC :: lagrange_init_box_tomas
  PUBLIC :: plume_inject_box
  PUBLIC :: plume_model_box
  PUBLIC :: plume_mod_cleanup_box

  ! PUBLIC VARIABLES:
  PUBLIC :: RXNRATE_CONST_KPP
  !PUBLIC :: SpcConc_BEFORE_KPP
  !PUBLIC :: SpcConc_AFTER_KPP
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
  REAL(fp)                              :: Length_init             ! m
  REAL(fp)                              :: Aircraft_speed          ! m/s
  REAL(fp)                              :: Critical_day_2D         ! Maximum 2-D plume lifetime allowed, day
  REAL(fp)                              :: Critical_day_1D         ! Maximum total plume lifetime allowed, day
  ! Species ID flags correspond to GEOS-Chem species
  INTEGER                               :: id_SO2,  id_SO4,  id_NH3,   id_NH4,  id_PASVLA
  INTEGER                               :: id_OH,   id_O3,   id_H2O,   id_HO2,  id_PH2SO4
  INTEGER                               :: id_NK01, id_SF01, id_AW01,  id_H2SO4
  INTEGER                               :: id_SO2pl, id_SO4pl

  ! Species ID flags correspond to Plume species
  INTEGER, PARAMETER                    :: nspc_p_bulk           = 6 ! Num of species in Plume, bulk
  INTEGER, PARAMETER                    :: nspc_p_tomas_tracer   = 3 ! Num of tomas tracer type in Plume, currently consider NK, SF, AW
  INTEGER, PARAMETER                    :: id_SO2_p              = 1
  INTEGER, PARAMETER                    :: id_SO4_p              = 2
  INTEGER, PARAMETER                    :: id_OH_p               = 3
  INTEGER, PARAMETER                    :: id_HO2_p              = 4
  INTEGER, PARAMETER                    :: id_PH2SO4_p           = 5
  INTEGER, PARAMETER                    :: id_PASVLA_p           = 6
  INTEGER, PARAMETER                    :: id_NH4_p              = 7
  INTEGER, PARAMETER                    :: id_NH3_p              = 8
  INTEGER, PARAMETER                    :: id_H2SO4_p            = 9
  INTEGER, PARAMETER                    :: id_H2O_p              = 10
  INTEGER, PARAMETER                    :: id_NK01_p             = 11
  INTEGER, PARAMETER                    :: id_SF01_p             = 12
  INTEGER, PARAMETER                    :: id_AW01_p             = 13 ! TOMAS tracer need to in this order: NK, ..., ..., AW
  CHARACTER(LEN=16), PARAMETER          :: spc_names_p(*)        = [ 'SO2', 'SO4', 'OH', &
                                  'HO2', 'PH2SO4','PASVLA','NH4','NH3', 'H2SO4','H2O', 'NK01','SF01','AW01' ] ! Uppercase

  ! hard-coded reaction id
  ! 202 is the reaction ID of SO2 + OH -> H2SO4 in fullchem KPP, need to adjust if using different mechanism or more species 
  INTEGER, PARAMETER    :: SO2_OH_RXN_ID = 202 
  
  
  ! Other variables
  INTEGER               :: nspc_GC ! Num of species in GEOS-Chem
  INTEGER               :: nspc_p ! Num of species in Plume
  INTEGER               :: nBins  ! Number of TOMAS bins
  !INTEGER               :: Volume_Sort  = 1 ! 1 = use SortList() function, transfer largest (not oldest) plume segment for volume criterion
  !INTEGER               :: Calc_entropy = 1 ! 1 = turn on entropy calculation
  !REAL(fp)		          :: Entropy0	 ! perfect entropy without diffusion
  INTEGER               :: NX_GC, NY_GC, NZ_GC
  INTEGER               :: n_x_mid, n_y_mid, n_x_mid2, n_y_mid2 
  INTEGER               :: n_x_max2, n_y_max2 
  INTEGER               :: n_slab_25, n_slab_50, n_slab_75
  INTEGER               :: n_slab_max, n_slab_max2
  INTEGER               :: N_parcel   ! 131        
  INTEGER               :: Num_inject, Num_Plume2d, Num_Plume1d, Num_Plume1d_new, Num_dissolve_2D, Num_dissolve_1D 
  INTEGER               :: Num_transfer_2D, Num_Plume2d_acc, Num_Plume1d_acc    
  INTEGER               :: tt 
  INTEGER               :: N_total
  INTEGER               :: Stop_inject ! 1: stop injecting; 0: keep injecting

  REAL(fp)              :: DX_GC, DY_GC
  REAL(fp)              :: Length_lat
  REAL(fp)              :: time_elapsed
  REAL(fp), POINTER     :: X_mid(:), Y_mid(:), P_mid(:)
  REAL(fp), POINTER     :: P_edge(:)
  REAL(fp), POINTER     :: X_edge(:), Y_edge(:)
  REAL(fp)              :: X_edge2, Y_edge2
  REAL(fp), ALLOCATABLE :: RXNRATE_CONST_KPP(:,:,:,:)
!   REAL(fp), ALLOCATABLE :: SpcConc_BEFORE_KPP(:,:,:,:)
!   REAL(fp), ALLOCATABLE :: SpcConc_AFTER_KPP(:,:,:,:)
  REAL(fp), ALLOCATABLE :: Vplume_2D_tot(:,:,:) ! total volume of 2-D plume in EU grid. cm3
  REAL(fp), ALLOCATABLE :: Vplume_1D_tot(:,:,:) ! total volume of 1-D plume in EU grid. cm3

  CHARACTER(LEN=16), ALLOCATABLE     :: spc_names_p_use(:) ! species used in plume-TOMAS model

  ! Diagnostic file definition, SO2 and Sulfate, bulk
  INTEGER               :: File_Smass_IU_2D       ! Diagnostic file nums, initialized in lagr_init
  INTEGER               :: File_spc_init_IU_2D    ! Diagnostic file nums, initialized in lagr_init
  INTEGER               :: File_Smass_IU_1D       ! Diagnostic file nums, initialized in lagr_init
  INTEGER               :: File_spc_init_IU_1D    ! Diagnostic file nums, initialized in lagr_init

  INTEGER               :: File_Plume_life_IU_2D  ! Diagnostic file nums, initialized in lagr_init
  INTEGER               :: File_Plume_life_IU_1D  ! Diagnostic file nums, initialized in lagr_init
  INTEGER               :: File_Plume_number_IU   ! Diagnostic file nums, initialized in lagr_init
  INTEGER               :: File_Plume_location_IU_2D   ! Diagnostic file location, initialized in lagr_init
  INTEGER               :: File_Plume_location_IU_1D   ! Diagnostic file location, initialized in lagr_init
  ! Variable to track sulfate mass (in molec S)
  ! _inj: injected
  ! _r1: release in entrainment/detrainment (check release in chem?)
  ! _r2: change due to chemistry
  ! _r3: mass release due to plume dissolve
  ! _r4: release due to 2D to 1D
  ! _1: mass in plume after injection
  ! _2: before physics
  ! _3: mass in Plume after Physical processes
  ! _4: before chemistry
  ! _5: mass in Plume after Chemical processes
  ! _6: before structure change
  ! _7: mass in plume after structure change
  ! _8: mass in plume core in 2-D plume box
  REAL(fp) 		          :: mass_S_SO2_inj_2D,   mass_S_SO2_r1_2D,   mass_S_SO2_r2_2D,   mass_S_SO2_r3_2D
  REAL(fp) 		          :: mass_S_SO2_r4
  REAL(fp) 		          :: mass_S_SO2_1_2D,     mass_S_SO2_2_2D,    mass_S_SO2_3_2D,    mass_S_SO2_4_2D
  REAL(fp) 		          :: mass_S_SO2_5_2D,     mass_S_SO2_6_2D,    mass_S_SO2_7_2D,    mass_S_SO2_8_2D

  REAL(fp) 		          :: mass_S_SO4_inj_2D,   mass_S_SO4_r1_2D,   mass_S_SO4_r2_2D,   mass_S_SO4_r3_2D
  REAL(fp) 		          :: mass_S_SO4_r4 
  REAL(fp) 		          :: mass_S_SO4_1_2D,     mass_S_SO4_2_2D,    mass_S_SO4_3_2D,    mass_S_SO4_4_2D
  REAL(fp) 		          :: mass_S_SO4_5_2D,     mass_S_SO4_6_2D,    mass_S_SO4_7_2D,    mass_S_SO4_8_2D

  REAL(fp) 		          :: mass_S_SO2_inj_1D,   mass_S_SO2_r1_1D,   mass_S_SO2_r2_1D,   mass_S_SO2_r3_1D
  REAL(fp) 		          :: mass_S_SO2_1_1D,     mass_S_SO2_2_1D,    mass_S_SO2_3_1D,    mass_S_SO2_4_1D
  REAL(fp) 		          :: mass_S_SO2_5_1D,     mass_S_SO2_6_1D,    mass_S_SO2_7_1D

  REAL(fp) 		          :: mass_S_SO4_inj_1D,   mass_S_SO4_r1_1D,   mass_S_SO4_r2_1D,   mass_S_SO4_r3_1D
  REAL(fp) 		          :: mass_S_SO4_1_1D,     mass_S_SO4_2_1D,    mass_S_SO4_3_1D,    mass_S_SO4_4_1D
  REAL(fp) 		          :: mass_S_SO4_5_1D,     mass_S_SO4_6_1D,    mass_S_SO4_7_1D

  REAL(fp) 		          :: mass_S_H2SO4_2D,     mass_S_H2SO4_1D 
  

  ! Variables to track plume values (used to convert mass into concentration)
  ! 1: befor phyhsics 2: after physics  3: after structure change;  4: Volume of plume core in 2-D plume box
  REAL(fp) 		          :: Vgrid_2D_tot_1, Vgrid_2D_tot_2, Vgrid_2D_tot_3, Vgrid_2D_tot_4
  REAL(fp) 		          :: Vgrid_1D_tot_1, Vgrid_1D_tot_2, Vgrid_1D_tot_3

  ! Variables to track initial mass read in plume (in molec), all tracers
  REAL(fp), ALLOCATABLE :: mass_spc_init_2D(:)
  REAL(fp), ALLOCATABLE :: mass_spc_init_1D(:)


  ! Diagnostic file definition, SO2 and Sulfate, size resolved tracer for TOMAS
  INTEGER               :: File_SF_bin_IU_2D ! Diagnostic file nums, initialized in lagr_init_tomas
  INTEGER               :: File_MK_bin_IU_2D ! Diagnostic file nums, initialized in lagr_init_tomas
  INTEGER               :: File_NK_bin_IU_2D ! Diagnostic file nums, initialized in lagr_init_tomas
  INTEGER               :: File_SF_bin_IU_1D ! Diagnostic file nums, initialized in lagr_init_tomas
  INTEGER               :: File_MK_bin_IU_1D ! Diagnostic file nums, initialized in lagr_init_tomas
  INTEGER               :: File_NK_bin_IU_1D ! Diagnostic file nums, initialized in lagr_init_tomas
  INTEGER               :: File_NK_bin_IU_2D_1D ! Diagnostic file nums, initialized in lagr_init_tomas
  INTEGER               :: File_SF_bin_IU_2D_1D ! Diagnostic file nums, initialized in lagr_init_tomas
  ! Variables to track tomas SF tracer (in molec), SF only
  REAL(fp), ALLOCATABLE :: mass_SF_bin_2D(:)
  REAL(fp), ALLOCATABLE :: mass_SF_bin_2D_1D(:)
  REAL(fp), ALLOCATABLE :: mass_SF_bin_1D(:)
  ! Variables to track tomas NK tracer (in molec)
  REAL(fp), ALLOCATABLE :: mass_NK_bin_2D(:)
  REAL(fp), ALLOCATABLE :: mass_NK_bin_2D_1D(:)
  REAL(fp), ALLOCATABLE :: mass_NK_bin_1D(:)
  ! Variables to track tomas total tracer (in molec), all tracer, for this case, mainly SF, NH4 and AW
  REAL(fp), ALLOCATABLE :: mass_MK_bin_2D(:)
  REAL(fp), ALLOCATABLE :: mass_MK_bin_1D(:)

  ! Diagnostic file names, initialized in lagr_init
  CHARACTER(LEN=255)    :: file_Smass_2D
  CHARACTER(LEN=255)    :: file_SF_bin_2D ! file to track tomas sf tracer
  CHARACTER(LEN=255)    :: file_NK_bin_2D ! file to track tomas nk tracer
  CHARACTER(LEN=255)    :: file_MK_bin_2D ! file to track tomas Mk tracer (all tracer include water)
  CHARACTER(LEN=255)    :: file_spc_init_2D
  CHARACTER(LEN=255)    :: file_Smass_1D
  CHARACTER(LEN=255)    :: file_SF_bin_1D ! file to track tomas sf tracer
  CHARACTER(LEN=255)    :: file_NK_bin_1D ! file to track tomas nk tracer
  CHARACTER(LEN=255)    :: file_MK_bin_1D ! file to track tomas Mk tracer (all tracer include water)
  CHARACTER(LEN=255)    :: file_spc_init_1D
  CHARACTER(LEN=255)    :: file_SF_bin_2D_1D ! file to track tomas sf tracer released from 2-D to 1-D
  CHARACTER(LEN=255)    :: file_NK_bin_2D_1D ! file to track tomas NK tracer released from 2-D to 1-D
  CHARACTER(LEN=255)    :: file_Plume_life_2D ! file to track 2-D plume lifetime
  CHARACTER(LEN=255)    :: file_Plume_life_1D ! file to track 1-D plume lifetime
  CHARACTER(LEN=255)    :: file_Plume_number  ! file to track number of plumes
  CHARACTER(LEN=255)    :: file_Plume_location_2D  ! file to track location of 2-D plumes
  CHARACTER(LEN=255)    :: file_Plume_location_1D  ! file to track location of 1-D plumes

  ! some parameter for sensitive test
  INTEGER, PARAMETER        :: N1_split           = 5            ! Cross-section splitting
  INTEGER, PARAMETER        :: N2_split           = 5            ! length splitting
  INTEGER, PARAMETER        :: Split_length       = 1            ! how many times of DX_GC
  REAL(fp), PARAMETER       :: Dissolve_criteria  = 10*0.01
  REAL(fp), PARAMETER       :: Volume_percent     = 30*0.01
  REAL(fp), PARAMETER       :: frac_mass = 0.95   ! fraction of tracer accumulated along the horizontal and verticle length scale
  ! REAL(fp), PARAMETER       :: Critical_day       = 28.0         ! [day] ! Move this into geoschem_config.yml
  REAL(fp), PARAMETER       :: Radical_exc_factor = 0.95_fp          ! BZ: Add for Radical exchange after chemistry, real number between 0 and 1
  REAL(fp), PARAMETER       :: Rb_min_shear       = 10.0_fp

  TYPE(Plume2d_list), POINTER :: Plume2d_tail, Plume2d_head
  TYPE(Plume1d_list), POINTER :: Plume1d_tail, Plume1d_head

  ! Diagnostic file, temporary
  INTEGER                :: file_2Dconc_SO4_ID, file_2Dconc_SO2_ID, file_2Dconc_OH_ID
  INTEGER                :: file_2Dconc_SO4_ID_1, file_2Dconc_SO2_ID_1, file_2Dconc_OH_ID_1
  INTEGER                :: file_2Dconc_SO4_ID_2, file_2Dconc_SO2_ID_2, file_2Dconc_OH_ID_2
  INTEGER                :: file_2Dconc_SO4_ID_3, file_2Dconc_SO2_ID_3, file_2Dconc_OH_ID_3

  INTEGER                :: file_2D_1D_conc_SO4_ID, file_2D_1D_conc_SO2_ID, file_2D_1D_conc_OH_ID
  INTEGER                :: file_1Dconc_SO4_ID, file_1Dconc_SO2_ID, file_1Dconc_OH_ID
  INTEGER                :: file_2D_conc_NKbin_ID, file_1D_conc_NKbin_ID
  INTEGER                :: file_2D_conc_SFbin_ID, file_1D_conc_SFbin_ID

  CHARACTER(LEN=255)     :: file_2Dconc_SO4, file_2Dconc_SO2, file_2Dconc_OH
  CHARACTER(LEN=255)     :: file_2Dconc_SO4_1, file_2Dconc_SO2_1, file_2Dconc_OH_1
  CHARACTER(LEN=255)     :: file_2Dconc_SO4_2, file_2Dconc_SO2_2, file_2Dconc_OH_2
  CHARACTER(LEN=255)     :: file_2Dconc_SO4_3, file_2Dconc_SO2_3, file_2Dconc_OH_3

  CHARACTER(LEN=255)     :: file_2D_1D_conc_SO4, file_2D_1D_conc_SO2, file_2D_1D_conc_OH
  CHARACTER(LEN=255)     :: file_1Dconc_SO4, file_1Dconc_SO2, file_1Dconc_OH
  CHARACTER(LEN=255)     :: file_2D_conc_NKbin, file_1D_conc_NKbin
  CHARACTER(LEN=255)     :: file_2D_conc_SFbin, file_1D_conc_SFbin

  ! Index of tested grid in 2-D, temporary, BZ
  INTEGER, PARAMETER     :: x_test = 3
  INTEGER, PARAMETER     :: y_test = 1
  INTEGER, PARAMETER     :: s_test = 3
  !-------------------------------------------------------------
  ! Some parameters to be retired 
  !integer               :: id_PASV_LA3, id_PASV_LA2, id_PASV_LA 
  !integer               :: id_PASV_EU2, id_PASV_EU
  !integer, parameter    :: i_tracer  = 1
  !integer, parameter    :: i_product = 2
  !real(fp), parameter       :: Kchem = 1.0e-20_fp ! chemical reaction rate
  ! use Kchem = 1.0e-20_fp for >1 year simulation
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Global variable for TOMAS-related processes
  INTEGER, PARAMETER   :: SRTSO4  = 1
  INTEGER, PARAMETER   :: SRTNACL = -1
  INTEGER, PARAMETER   :: SRTECOB = -1
  INTEGER, PARAMETER   :: SRTECIL = -1
  INTEGER, PARAMETER   :: SRTOCOB = -1
  INTEGER, PARAMETER   :: SRTOCIL = -1
  INTEGER, PARAMETER   :: SRTDUST = -1
  INTEGER, PARAMETER   :: SRTNH4  = 2
  INTEGER, PARAMETER   :: SRTH2O  = 3

  INTEGER, PARAMETER   :: ICOMPHARD =  3 ! number of used variable in TOMAS, not include Nk
  ! !PUBLIC DATA MEMBERS:
  INTEGER  :: bin_nuc = 1, tern_nuc = 1  ! Switches for nucleation type.
  INTEGER  :: act_nuc = 0 ! in BL
  INTEGER  :: ion_nuc = 0 ! 1 for modgil, 2 for Yu
  INTEGER  :: lowRH = 1    !This is to match AW more with AERONET (JKODROS 6/15)
  ! Arrays
  REAL(fp), SAVE,   ALLOCATABLE, TARGET :: Xk(:)
  ! REAL*4,   SAVE,   ALLOCATABLE :: MOLWT(:)
  REAL(fp), ALLOCATABLE         :: AVGMASS(:)       ! Average mass per particle mid-range of size bin [kg/no.]

CONTAINS
  SUBROUTINE lagrange_init_box(am_I_root, Input_Opt, State_Chm, State_Grid, State_Met, RC)

    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Met_Mod,   ONLY : MetState
    USE State_Chm_Mod,   ONLY : ChmState, Ind_
    USE State_Grid_Mod,  ONLY : GrdState
    USE Species_Mod,     ONLY : SpcConc
    USE TIME_MOD,        ONLY : GET_TS_DYN
    USE TIME_MOD,        ONLY : GET_YEAR, GET_MONTH, GET_DAY, GET_HOUR, GET_MINUTE, GET_SECOND
    USE TIME_MOD,        ONLY : ITS_TIME_FOR_CHEM, ITS_TIME_FOR_DYN, ITS_TIME_FOR_EXIT
    USE InquireMod,      ONLY : findFreeLun
    USE GcKpp_Parameters,           ONLY : NREACT, NSPEC

    LOGICAL,        INTENT(IN)            :: am_I_Root   ! Are we on the root CPU
    TYPE(MetState), INTENT(in)            :: State_Met
    TYPE(ChmState), INTENT(inout)         :: State_Chm
    TYPE(OptInput), INTENT(in)            :: Input_Opt
    TYPE(GrdState), INTENT(IN)            :: State_Grid  ! Grid State objectgg
    INTEGER,        INTENT(OUT)           :: RC         ! Success or failure
    ! Pointers
    TYPE(SpcConc), POINTER                :: Spc(:)
    
    INTEGER                       :: i_box, i_slab, ibin
    INTEGER                       :: ii, jj, kk
    INTEGER                       :: id_tracer, N
    INTEGER                       :: i_lon, i_lat, i_lev            !1:NX_GC
    INTEGER                       :: previous_units, previous_units_temp
    INTEGER                       :: IOS
    !INTEGER                       :: this_year, this_month, this_day, this_hour, this_minute, this_second

    !LOGICAL                       :: exe_dyn, exe_chem, exe_exit

    CHARACTER(LEN=255)            :: spc_name
    CHARACTER(LEN=255)            :: FILENAME, FileEntropy, File996
    CHARACTER(LEN=255)            :: FILENAME2, FILENAME3
    CHARACTER(LEN=255)            :: ErrMsg, ThisLoc
    
    REAL(fp)                      :: lon1, lat1, lon2, lat2
    REAL(fp)                      :: box_lon_edge, box_lat_edge
    REAL(fp)                      :: curr_lon, curr_lat, curr_lev
    REAL(fp)                      :: Dt_dyn
    REAL(fp)                      :: plume_len_deg, plume_rem_deg, add_deg 
    REAL(fp), PARAMETER           :: eps = 1.0e-10  


    !exe_dyn            =   ITS_TIME_FOR_DYN()
    !exe_chem           =   ITS_TIME_FOR_CHEM()
    !exe_exit           =   ITS_TIME_FOR_EXIT()

    !IF (exe_dyn)  WRITE(6, *)  "Debug (BZ): Initialization: Is time for dynamic"
    !IF (exe_chem) WRITE(6, *)  "Debug (BZ): Initialization: Is time for chem"
    !IF (exe_exit) WRITE(6, *)  "Debug (BZ): Initialization: Is time for exit"

    RC                 =   GC_SUCCESS
    ErrMsg             =   ''
    ThisLoc            =   ' -> at lagrange_init_box (in module GeosCore/lagrange_singlebox_mod.F90)'
    Spc                =>   State_Chm%Species

    Num_inject         =    0
    Num_Plume2d        =    0  ! Store the current number of 2D plume, change due to transfer or dissolve
    Num_Plume1d        =    0  ! Store the current number of 1D plume, change due to transfer or dissolve
    Num_Plume1d_new    =    0
    Num_dissolve_2D    =    0
    Num_dissolve_1D    =    0
    Num_transfer_2D    =    0
    Num_Plume2d_acc    =    0  ! Store the number of 2D plume created, no decrease due to transfer or dissolve
    Num_Plume1d_acc    =    0  ! Store the number of 1D plume created, no decrease due to transfer or dissolve
    time_elapsed       =    0.0_fp

    nspc_GC            =              State_Chm%nSpecies
    nspc_p             =              nspc_p_bulk
    NX_GC              =              State_Grid%NX
    NY_GC              =              State_Grid%NY
    NZ_GC              =              State_Grid%NZ
    DX_GC              =              State_Grid%DX
    DY_GC              =              State_Grid%DY
    X_edge             =>             State_Grid%XEdge(:,1) 
    Y_edge             =>             State_Grid%YEdge(1,:)
    P_edge             =>             State_Met%PEDGE(1,1,:) 
    X_edge2            =              X_edge(2)
    Y_edge2            =              Y_edge(2)
    X_mid              =>             State_Grid%XMid(:,1) ! Grid box longitude [degrees] ! XMID(:,1,1)   ! NX_GC ! new
    Y_mid              =>             State_Grid%YMid(1,:) ! Grid box latitude center [degree] ! YMID(1,:,1)
    P_mid              =>             State_Met%PMID(1,1,:)  ! Pressure at level centers (hPa)
   
    !Write (6, *) "Debug: (BZ) Num of reactions from KPP is: ", NREACT
    Write (6, *) "Debug: (BZ) Num of species in plume is : ", nspc_p
    ! create species name database
    ALLOCATE( spc_names_p_use(nspc_p), STAT=RC )
    IF ( RC /= 0 ) THEN
      errMsg = 'Error allocating spc_names_p_use'
      CALL ERROR_STOP( errMsg, thisLoc)
      RETURN
    ENDIF
    spc_names_p_use(1:nspc_p)=spc_names_p(1:nspc_p)

    ALLOCATE( RXNRATE_CONST_KPP(NX_GC, NY_GC, NZ_GC, NREACT ), STAT=RC )
    IF (RC /= 0) THEN
        errMsg = 'Error allocating RXNRATE_CONST_KPP'
        CALL ERROR_STOP( errMsg, thisLoc)
        RETURN
    ENDIF

   !  ALLOCATE( SpcConc_BEFORE_KPP(NX_GC, NY_GC, NZ_GC, NSPEC ), STAT=RC )
   !  IF (RC /= 0) THEN
   !      errMsg = 'Error allocating SpcConc_BEFORE_KPP'
   !      CALL ERROR_STOP( errMsg, thisLoc)
   !      RETURN
   !  ENDIF

   !  ALLOCATE( SpcConc_AFTER_KPP(NX_GC, NY_GC, NZ_GC, NSPEC ), STAT=RC )
   !  IF (RC /= 0) THEN
   !      errMsg = 'Error allocating SpcConc_AFTER_KPP'
   !      CALL ERROR_STOP( errMsg, thisLoc)
   !      RETURN
   !  ENDIF
    ALLOCATE( Vplume_2D_tot(NX_GC, NY_GC, NZ_GC ), STAT=RC )
    IF (RC /= 0) THEN
        errMsg = 'Error allocating Vplume_2D_tot'
        CALL ERROR_STOP( errMsg, thisLoc)
        RETURN
    ENDIF
    ALLOCATE( Vplume_1D_tot(NX_GC, NY_GC, NZ_GC ), STAT=RC )
    IF (RC /= 0) THEN
        errMsg = 'Error allocating Vplume_1D_tot'
        CALL ERROR_STOP( errMsg, thisLoc)
        RETURN
    ENDIF
   !  ALLOCATE( mass_spc_init_2D(10 + nspc_p_tomas_tracer * nBins), STAT=RC )
   !  IF (RC /= 0) THEN
   !      errMsg = 'Error allocating mass_spc_init_2D'
   !      CALL ERROR_STOP( errMsg, thisLoc)
   !      RETURN
   !  ENDIF
   !  ALLOCATE( mass_spc_init_1D(10 + nspc_p_tomas_tracer * nBins), STAT=RC )
   !  IF (RC /= 0) THEN
   !      errMsg = 'Error allocating mass_spc_init_1D'
   !      CALL ERROR_STOP( errMsg, thisLoc)
   !      RETURN
   !  ENDIF
    ! Copy all neccessary variable from geoschem.yml here
    use_lagrange        =             Input_Opt%LagrangianModel_Activate
    plume_inject_on     =             Input_Opt%PlumeInjection_Activate
    plume_diag          =             Input_Opt%PlumeInjection_Diag
    TROPP_sink          =             Input_Opt%TropSink_Activate
    Num_of_sources      =             Input_Opt%Plume_sources_num
    Length_init         =             Input_Opt%Initial_length*1000.0   ! km to m
    Aircraft_speed      =             Input_Opt%Aircraft_speed  !m/s
    Critical_day_2D     =             Input_Opt%Critical_day_2D    ! Day
    Critical_day_1D     =             Input_Opt%Critical_day_1D    ! Day
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
    n_slab_max                 =      (INT(n_y_max/4)+1)*4 ! close to n_y_max, number of slabs in 1D
    n_slab_max2                =      n_slab_max+2

   
   !  mass_S_SO2_1_2D               =      0.0_fp
   !  mass_S_SO2_2_2D               =      0.0_fp
   !  mass_S_SO2_3_2D               =      0.0_fp
   !  mass_S_SO2_inj_2D             =      0.0_fp
   !  mass_S_SO2_r1_2D              =      0.0_fp
   !  mass_S_SO2_r2_2D              =      0.0_fp
    
   !  mass_S_SO4_1_2D               =      0.0_fp
   !  mass_S_SO4_2_2D               =      0.0_fp
   !  mass_S_SO4_3_2D               =      0.0_fp
   !  mass_S_SO4_inj_2D             =      0.0_fp
   !  mass_S_SO4_r1_2D              =      0.0_fp
   !  mass_S_SO4_r2_2D              =      0.0_fp

   !  mass_S_SO2_1_1D               =      0.0_fp
   !  mass_S_SO2_2_1D               =      0.0_fp
   !  mass_S_SO2_3_1D               =      0.0_fp
   !  mass_S_SO2_inj_1D             =      0.0_fp
   !  mass_S_SO2_r1_1D              =      0.0_fp
   !  mass_S_SO2_r2_1D              =      0.0_fp
    
   !  mass_S_SO4_1_1D               =      0.0_fp
   !  mass_S_SO4_2_1D               =      0.0_fp
   !  mass_S_SO4_3_1D               =      0.0_fp
   !  mass_S_SO4_inj_1D             =      0.0_fp
   !  mass_S_SO4_r1_1D              =      0.0_fp
   !  mass_S_SO4_r2_1D              =      0.0_fp

   !  mass_S_SO2_r3                 =      0.0_fp
   !  mass_S_SO4_r3                 =      0.0_fp

    !this_year                     =      GET_YEAR()
    !this_month                    =      GET_MONTH()
    !this_day                      =      GET_DAY()
    !this_hour                     =      GET_HOUR()
    !this_minute                   =      GET_MINUTE()
    !this_second                   =      GET_SECOND()

    !WRITE(6,*) 'Debug (BZ): Initialization: Current simulation time is : ', this_year, this_month, &
    !this_day, this_hour, this_minute, this_second

    Dt_dyn = GET_TS_DYN()
    N_parcel = NINT(Aircraft_speed * Dt_dyn / Length_init)
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

    ! Creat file for sulfur mass diag
    File_Smass_IU_2D = findFreeLun()
    file_Smass_2D   = 'Plume_Sulfur_mass_2D.txt'
    OPEN( File_Smass_IU_2D, FILE=TRIM( file_Smass_2D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=IOS )
    WRITE(File_Smass_IU_2D,'(*(A,1X))' ) &
    'time_elapsed', 'mass_S_SO2_inj',  'mass_S_SO2_r1', 'mass_S_SO2_r2', 'mass_S_SO2_r3', &
         'mass_S_SO2_7', 'mass_S_SO2_8', &
         'mass_S_SO4_inj',  'mass_S_SO4_r1',  'mass_S_SO4_r2', 'mass_S_SO4_r3', & 
         'mass_S_SO4_7', 'mass_S_SO4_8', &
         'Vgrid_tot_3', 'Vgrid_tot_4', &
         'mass_S_SO2_r4', 'mass_S_SO4_r4'
    ! Return if there was an error opening the file
    IF ( IOS /= 0 ) THEN
        ! Define error message
        ErrMsg  = 'create Plume_Sulfur_mass_2D.txt (in lagrange_init_box)'
        CALL ERROR_STOP( ERRMSG, ThisLoc )
        RETURN
    ENDIF
    !CLOSE(File_Smass_IU_2D)

    ! Creat file for sulfur mass diag
    File_Smass_IU_1D = findFreeLun()
    file_Smass_1D   = 'Plume_Sulfur_mass_1D.txt'
    OPEN( File_Smass_IU_1D, FILE=TRIM( file_Smass_1D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=IOS )
    WRITE(File_Smass_IU_1D,'(*(A,1X))' ) &
    'time_elapsed', 'mass_S_SO2_inj',  'mass_S_SO2_r1', 'mass_S_SO2_r2', 'mass_S_SO2_r3', &
         'mass_S_SO2_7', &
         'mass_S_SO4_inj',  'mass_S_SO4_r1',  'mass_S_SO4_r2', 'mass_S_SO4_r3', & 
         'mass_S_SO4_7', 'Vgrid_tot_3'
    ! Return if there was an error opening the file
    IF ( IOS /= 0 ) THEN
        ! Define error message
        ErrMsg  = 'create Plume_Sulfur_mass_1D.txt (in lagrange_init_box)'
        CALL ERROR_STOP( ERRMSG, ThisLoc )
        RETURN
    ENDIF
    !CLOSE(File_Smass_IU_1D)

    File_spc_init_IU_2D = findFreeLun()
    file_spc_init_2D = 'Plume_species_initial_2D.txt'
    OPEN( File_spc_init_IU_2D, FILE=TRIM( file_spc_init_2D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=IOS )
    ! Return if there was an error opening the file
    IF ( IOS /= 0 ) THEN
        ! Define error message
        ErrMsg  = 'create Plume_species_initial_2D.txt (in lagrange_init_box)'
        CALL ERROR_STOP( ErrMsg, ThisLoc )
        RETURN
    ENDIF
    !CLOSE(File_spc_init_IU_2D)

    File_spc_init_IU_1D = findFreeLun()
    file_spc_init_1D = 'Plume_species_initial_1D.txt'
    OPEN( File_spc_init_IU_1D, FILE=TRIM( file_spc_init_1D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=IOS )
    ! Return if there was an error opening the file
    IF ( IOS /= 0 ) THEN
        ! Define error message
        ErrMsg  = 'create Plume_species_initial_1D.txt (in lagrange_init_box)'
        CALL ERROR_STOP( ErrMsg, ThisLoc )
        RETURN
    ENDIF
    !CLOSE(File_spc_init_IU_1D)

    File_Plume_life_IU_2D = findFreeLun()
    file_Plume_life_2D = 'Plume_lifetime_2D.txt'
    OPEN( File_Plume_life_IU_2D, FILE=TRIM( file_Plume_life_2D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=IOS )
    WRITE(File_Plume_life_IU_2D,'(A)') & 
    'PlumeLabel TransferLifetime DissolveLifetime'
    ! Return if there was an error opening the file
    IF ( IOS /= 0 ) THEN
        ! Define error message
        ErrMsg  = 'create Plume_lifetime_2D.txt (in lagrange_init_box)'
        CALL ERROR_STOP( ErrMsg, ThisLoc )
        RETURN
    ENDIF
    !CLOSE(File_Plume_life_IU_2D)

    File_Plume_life_IU_1D = findFreeLun()
    file_Plume_life_1D = 'Plume_lifetime_1D.txt'
    OPEN( File_Plume_life_IU_1D, FILE=TRIM( file_Plume_life_1D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=IOS )
    WRITE(File_Plume_life_IU_1D,'(A)') & 
    'PlumeLabel DissolveLifetime'
    ! Return if there was an error opening the file
    IF ( IOS /= 0 ) THEN
        ! Define error message
        ErrMsg  = 'create Plume_lifetime_1D.txt (in lagrange_init_box)'
        CALL ERROR_STOP( ErrMsg, ThisLoc )
        RETURN
    ENDIF
    !CLOSE(File_Plume_life_IU_1D)

    File_Plume_number_IU = findFreeLun()
    file_Plume_number = 'Plume_number.txt'
    OPEN( File_Plume_number_IU, FILE=TRIM( file_Plume_number ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=IOS )
    WRITE(File_Plume_number_IU,'(A)') &
    'time_elapsed Num_inject Num_Plume2d Num_Plume1d Num_dissolve_2D ' // &
    'Num_transfer_2D Num_dissolve_1D'
    ! Return if there was an error opening the file
    IF ( IOS /= 0 ) THEN
        ! Define error message
        ErrMsg  = 'create Plume_number.txt (in lagrange_init_box)'
        CALL ERROR_STOP( ErrMsg, ThisLoc )
        RETURN
    ENDIF
    !CLOSE(File_Plume_number_IU)

    File_Plume_location_IU_2D = findFreeLun()
    file_Plume_location_2D = 'Plume_location_2D.txt'
    OPEN( File_Plume_location_IU_2D, FILE=TRIM( file_Plume_location_2D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=IOS )
    WRITE(File_Plume_location_IU_2D,'(A)') & 
    'ElapsedTime PlumeLabel iLon iLat iLev Lon Lat Lev dx dy Length'
    ! Return if there was an error opening the file
    IF ( IOS /= 0 ) THEN
        ! Define error message
        ErrMsg  = 'create Plume_location_2D.txt (in lagrange_init_box)'
        CALL ERROR_STOP( ErrMsg, ThisLoc )
        RETURN
    ENDIF

    File_Plume_location_IU_1D = findFreeLun()
    file_Plume_location_1D = 'Plume_location_1D.txt'
    OPEN( File_Plume_location_IU_1D, FILE=TRIM( file_Plume_location_1D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=IOS )
    WRITE(File_Plume_location_IU_1D,'(A)') & 
    'ElapsedTime PlumeLabel iLon iLat iLev Lon Lat Lev Ra Rb Length'
    ! Return if there was an error opening the file
    IF ( IOS /= 0 ) THEN
        ! Define error message
        ErrMsg  = 'create Plume_location_1D.txt (in lagrange_init_box)'
        CALL ERROR_STOP( ErrMsg, ThisLoc )
        RETURN
    ENDIF
    
!90  FORMAT( /, A                 )
!95  FORMAT( A                    )
!100 FORMAT( A, L5                )
!105 FORMAT( A, I0                )
!110 FORMAT( A, A                 )
!120 FORMAT(A8,1X,"|", A8,1X,"|",A8,1X,"|",1X,A8,1X,"|",1X,A8,1X,"|",1X,A10,1X,"|",1X,A)
!121 FORMAT(I8,1X,"|",F8.3,1X,"|",F8.3,1X,"|",1X,F8.3,1X,"|",1X,F8.3,1X,"|",1X,F10.3,1X,"|",1X,A)

  END SUBROUTINE lagrange_init_box

  SUBROUTINE lagrange_init_box_tomas(am_I_root, Input_Opt, State_Chm, State_Grid, State_Met, RC)

    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Met_Mod,   ONLY : MetState
    USE State_Chm_Mod,   ONLY : ChmState, Ind_
    USE State_Grid_Mod,  ONLY : GrdState
    USE Species_Mod,     ONLY : SpcConc
    USE TIME_MOD,        ONLY : GET_TS_DYN
    USE TIME_MOD,        ONLY : GET_YEAR, GET_MONTH, GET_DAY, GET_HOUR, GET_MINUTE, GET_SECOND
    USE TIME_MOD,        ONLY : ITS_TIME_FOR_CHEM, ITS_TIME_FOR_DYN, ITS_TIME_FOR_EXIT
    USE InquireMod,      ONLY : findFreeLun
    !USE GcKpp_Parameters,           ONLY : NREACT, NSPEC

    LOGICAL,        INTENT(IN)            :: am_I_Root   ! Are we on the root CPU
    TYPE(MetState), INTENT(in)            :: State_Met
    TYPE(ChmState), INTENT(inout)         :: State_Chm
    TYPE(OptInput), INTENT(in)            :: Input_Opt
    TYPE(GrdState), INTENT(IN)            :: State_Grid  ! Grid State objectgg
    INTEGER,        INTENT(OUT)           :: RC         ! Success or failure
    
    INTEGER                               :: iBin, i_species
    INTEGER                               :: this_year, this_month, this_day, this_hour, this_minute, this_second
    
    LOGICAL                               :: exe_dyn, exe_chem, exe_exit

    REAL(fp)                              :: Mo
    
    CHARACTER(LEN = 16)                   :: first_tracer, this_tracer

    this_year                     =      GET_YEAR()
    this_month                    =      GET_MONTH()
    this_day                      =      GET_DAY()
    this_hour                     =      GET_HOUR()
    this_minute                   =      GET_MINUTE()
    this_second                   =      GET_SECOND()

    !WRITE(6,*) 'Debug (BZ): Initialization - TOMAS: Current simulation time is : ', this_year, this_month, &
    !this_day, this_hour, this_minute, this_second

    exe_dyn            =   ITS_TIME_FOR_DYN()
    exe_chem           =   ITS_TIME_FOR_CHEM()
    exe_exit           =   ITS_TIME_FOR_EXIT()

    !IF (exe_dyn)  WRITE(6, *)  "Debug (BZ): Initialization - TOMAS: Is time for dynamic"
    !IF (exe_chem) WRITE(6, *)  "Debug (BZ): Initialization - TOMAS: Is time for chem"
    !IF (exe_exit) WRITE(6, *)  "Debug (BZ): Initialization - TOMAS: Is time for exit"

    nBins              =              State_Chm%nTomasBins
    nspc_p             =              10 + nspc_p_tomas_tracer * nBins
    Write (6, *) "Debug: (BZ) Num of Bins in plume is : ", nBins
    Write (6, *) "Debug: (BZ) Num of species in plume changed to : ", nspc_p
    !-------------------------------------------------
    ! Initialization process for TOMAS
    ALLOCATE( Xk( nBins+1 ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'Xk [TOMAS] in plume' )
    Xk(:) = 0e+0_fp
    ALLOCATE( AVGMASS( nBins ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'AVGMASS [TOMAS] in plume' )
    AVGMASS(:) = 0e+0_fp

#if defined(TOMAS40)
    Mo = 1.0e-21_fp*2.e+0_fp**(-10)
#elif defined(TOMAS15)
    Mo = 1.0e-21_fp*4.e+0_fp**(-3)
#else
    Mo = 1.0e-21_fp
#endif
    ! Write (6, *) "Debug: (BZ) Mo =  : ", Mo
#if defined(TOMAS12) || defined(TOMAS15)
    DO ibin = 1, nBins + 1
       if(ibin.lt.nBins)then
          xk(ibin)=Mo * 4.e+0_fp**(ibin-1) !mass quadrupling
       else
          xk(ibin)=xk(ibin-1) * 32.e+0_fp
       endif
    ENDDO
#else
    DO ibin = 1, nBins + 1
       Xk( ibin ) = Mo * 2.e+0_fp ** ( ibin-1 )
    ENDDO
#endif

#ifdef TOMAS
    DO ibin = 1, nBins
       AVGMASS( ibin ) = sqrt(Xk(ibin)*Xk(ibin+1))
    ENDDO
#endif
    ! create species name database for TOMAS
    IF (ALLOCATED (spc_names_p_use)) THEN
      DEALLOCATE (spc_names_p_use)
      Write (6, *) "Debug: (BZ) Reallocate plume model database : spc_names_p_use"
    ENDIF
    ALLOCATE( spc_names_p_use(nspc_p), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'spc_names_p_use in plume' )
    spc_names_p_use(1:10)=spc_names_p(1:10)
    
     DO i_species = 1, nspc_p_tomas_tracer
         first_tracer = TRIM(spc_names_p(10+i_species))   
         DO ibin = 1, nBins
         write(this_tracer, '(A, I2.2)') TRIM(first_tracer(1:LEN_TRIM(first_tracer)-2)), ibin
         spc_names_p_use(10+i_species+(ibin-1)*nspc_p_tomas_tracer) = &
                  TRIM(this_tracer)
      ENDDO
    ENDDO
    ! Dubug printout
    Write (6,*) "Debug (BZ): Plume-TOMAS species database:"
    DO i_species = 1, nspc_p
       Write (6,*) "Species in Plume: Number", i_species, "; Name: ", TRIM(spc_names_p_use(i_species))
    ENDDO 

    ! Creat file for TOMAS tracer diag
    File_SF_bin_IU_2D = findFreeLun()
    file_SF_bin_2D   = 'Plume_SF_bin_2D.txt'
    OPEN( File_SF_bin_IU_2D, FILE=TRIM( file_SF_bin_2D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=RC )
    !CLOSE(File_SF_bin_IU_2D)

    File_SF_bin_IU_1D = findFreeLun()
    file_SF_bin_1D   = 'Plume_SF_bin_1D.txt'
    OPEN( File_SF_bin_IU_1D, FILE=TRIM( file_SF_bin_1D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=RC )
    !CLOSE(File_SF_bin_IU_1D)

    File_SF_bin_IU_2D_1D = findFreeLun()
    file_SF_bin_2D_1D   = 'Plume_SF_bin_2D_1D.txt'
    OPEN( File_SF_bin_IU_2D_1D, FILE=TRIM( file_SF_bin_2D_1D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=RC )
    !CLOSE(File_SF_bin_IU_2D_1D)

    
    File_NK_bin_IU_2D = findFreeLun()
    file_NK_bin_2D   = 'Plume_NK_bin_2D.txt'
    OPEN( File_NK_bin_IU_2D, FILE=TRIM( file_NK_bin_2D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=RC )
    !CLOSE(File_NK_bin_IU_2D)

    File_NK_bin_IU_1D = findFreeLun()
    file_NK_bin_1D   = 'Plume_NK_bin_1D.txt'
    OPEN( File_NK_bin_IU_1D, FILE=TRIM( file_NK_bin_1D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=RC )
    !CLOSE(File_NK_bin_IU_1D)

    File_NK_bin_IU_2D_1D = findFreeLun()
    file_NK_bin_2D_1D   = 'Plume_NK_bin_2D_1D.txt'
    OPEN( File_NK_bin_IU_2D_1D, FILE=TRIM( file_NK_bin_2D_1D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=RC )
    !CLOSE(File_NK_bin_IU_2D_1D)

    File_MK_bin_IU_2D = findFreeLun()
    file_MK_bin_2D   = 'Plume_MK_bin_2D.txt'
    OPEN( File_MK_bin_IU_2D, FILE=TRIM( file_MK_bin_2D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=RC )
    !CLOSE(File_MK_bin_IU_2D)

    File_MK_bin_IU_1D = findFreeLun()
    file_MK_bin_1D   = 'Plume_MK_bin_1D.txt'
    OPEN( File_MK_bin_IU_1D, FILE=TRIM( file_MK_bin_1D ), STATUS='REPLACE', &
        FORM='FORMATTED',  ACCESS='SEQUENTIAL',     IOSTAT=RC )
    !CLOSE(File_MK_bin_IU_1D)

    ! Creat arrays for TOMAS tracer diag
    ALLOCATE( mass_SF_bin_2D( nBins ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'mass_SF_bin_2D [TOMAS] in plume' )
    mass_SF_bin_2D(:) = 0e+0_fp

    ALLOCATE( mass_SF_bin_1D( nBins ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'mass_SF_bin_1D [TOMAS] in plume' )
    mass_SF_bin_1D(:) = 0e+0_fp

    ALLOCATE( mass_SF_bin_2D_1D( nBins ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'mass_SF_bin_2D_1D [TOMAS] in plume' )
    mass_SF_bin_2D_1D(:) = 0e+0_fp

    ALLOCATE( mass_NK_bin_2D( nBins ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'mass_NK_bin_2D [TOMAS] in plume' )
    mass_NK_bin_2D(:) = 0e+0_fp

    ALLOCATE( mass_NK_bin_1D( nBins ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'mass_NK_bin_1D [TOMAS] in plume' )
    mass_NK_bin_1D(:) = 0e+0_fp

    ALLOCATE( mass_NK_bin_2D_1D( nBins ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'mass_NK_bin_2D_1D [TOMAS] in plume' )
    mass_NK_bin_2D_1D(:) = 0e+0_fp

    ALLOCATE( mass_MK_bin_2D( nBins ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'mass_MK_bin_2D [TOMAS] in plume' )
    mass_MK_bin_2D(:) = 0e+0_fp

    ALLOCATE( mass_MK_bin_1D( nBins ), STAT=RC )
    IF ( RC /= 0 ) CALL ALLOC_ERR( 'mass_MK_bin_1D [TOMAS] in plume' )
    mass_MK_bin_1D(:) = 0e+0_fp

  END SUBROUTINE lagrange_init_box_tomas
  SUBROUTINE plume_inject_box(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)

    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Met_Mod,   ONLY : MetState
    USE State_Chm_Mod,   ONLY : ChmState, Ind_
    USE State_Grid_Mod,  ONLY : GrdState
    USE Species_Mod,     ONLY : SpcConc
    USE TIME_MOD,        ONLY : GET_TS_DYN, GET_TS_CHEM
    USE TIME_MOD,        ONLY : GET_YEAR, GET_MONTH, GET_DAY, GET_HOUR, GET_MINUTE, GET_SECOND
    USE TIME_MOD,        ONLY : ITS_TIME_FOR_CHEM, ITS_TIME_FOR_DYN, ITS_TIME_FOR_EXIT
    ! USE UnitConv_Mod !,    ONLY : Convert_Spc_Units, MOLECULES_SPECIES_PER_CM3
   

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


    INTEGER                :: i_box, i_lon, i_lat, i_lev, i_species, ibin
    INTEGER                :: nAdv, i_advect1
    INTEGER                :: id_tracer, id_tracer1
    INTEGER                :: previous_units, previous_units_temp
    INTEGER                :: Num_Stop
    INTEGER                :: this_year, this_month, this_day, this_hour, this_minute, this_second

    LOGICAL                :: exe_dyn, exe_chem, exe_exit

    REAL(fp)               :: box_lon, box_lat, box_lev
    REAL(fp), POINTER      :: PASV_EU
    REAL(fp)               :: MW_g
    REAL(fp)               :: Dt_dyn, Dt_chem
    REAL(fp)               :: Vgrid_2D, Vgrid_1D, Vgrid_EU
    REAL(fp)               :: Entropy2_Concnt, Entropy2_V, Entropy2
    REAL(fp)               :: tracer2_mol, air2_mol, mix2_ratio
    REAL(fp)               :: tracer0_mol, air0_mol, mix0_ratio
    !REAL(fp)               :: plume_length_deg, plane_route_deg, plane_loc_last

    !CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: spc_name
    CHARACTER(LEN=255)     :: ErrMsg
    CHARACTER(LEN=255)     :: FileEntropy
    CHARACTER(LEN=255)     :: ThisLoc

    ErrMsg                 =    ''
    Dt_dyn                 =    GET_TS_DYN()
    Dt_chem                =    GET_TS_CHEM()
    
    Spc                    =>   State_Chm%Species
    ThisLoc                =   ' -> at plume_inject_box (in module GeosCore/lagrange_singlebox_mod.F90)'

    id_SO4                 =    Ind_('SO4')
    id_SO2                 =    Ind_('SO2')
    id_OH                  =    Ind_('OH')
    id_HO2                 =    Ind_('HO2')
    id_PH2SO4              =    Ind_('PH2SO4')
    id_PASVLA              =    Ind_('PASVLA')


    mass_S_SO2_inj_2D      =    0.0_fp
    mass_S_SO4_inj_2D      =    0.0_fp
    mass_S_SO2_1_2D        =    0.0_fp
    mass_S_SO4_1_2D        =    0.0_fp

    mass_S_SO2_inj_1D      =    0.0_fp
    mass_S_SO4_inj_1D      =    0.0_fp
    mass_S_SO2_1_1D        =    0.0_fp
    mass_S_SO4_1_1D        =    0.0_fp
    
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
    !WRITE(6,'(a)') 'debug : Before plume injection, the unit is ' // TRIM(UNIT_STR(previous_units_temp))
    !WRITE(6,'(a)') 'debug : During plume injection, the unit is ' // TRIM(UNIT_STR(State_Chm%Species(id_SO4)%Units))
    
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

                i_lon = Find_iLonLat(box_lon, DX_GC, X_edge2)
                if(i_lon>NX_GC) i_lon=i_lon-NX_GC
                if(i_lon<1) i_lon=i_lon+NX_GC
                !WRITE(6,*) 'debug: ilon=', i_lon
                i_lat = Find_iLonLat(box_lat, DY_GC, Y_edge2)
                if(i_lat>NY_GC) i_lat=NY_GC
                if(i_lat<1) i_lat=1
                !WRITE(6,*) 'debug: ilat=', i_lat
                i_lev = Find_iPLev(box_lev,P_edge)
                !WRITE(6,*) 'debug: ilev1=', i_lev
                if(i_lev>NZ_GC) i_lev=NZ_GC
                ! instantly add injected species into Eulerian grid 
                DO i_species = 1, num_of_sources
                    spc_name = Plume_sources(i_species)%species
                    id_tracer   = Ind_(TRIM(spc_name))
                    !write(6,*) 'debug (BZ): species: ', TRIM(spc_name), 'id: ', id_tracer, 'conc before injection: ', Spc(id_tracer)%Conc(i_lon, i_lat, i_lev)
                    Spc(id_tracer)%Conc(i_lon, i_lat, i_lev) = Spc(id_tracer)%Conc(i_lon, i_lat, i_lev)  &
                            + (Length_init * Plume_sources(i_species)%rate / State_Chm%SpcData(id_tracer)%Info%MW_g * Avo) &
                            /(State_Met%AIRVOL(i_lon, i_lat, i_lev) * 1.0e+6_fp) ! molec/cm3
                    !write(6,*) 'debug (BZ): species: ', TRIM(spc_name), 'id: ', id_tracer, 'conc after injection: ', Spc(id_tracer)%Conc(i_lon, i_lat, i_lev)

                ENDDO
            ENDDO
            Num_inject = Num_Stop 

        ELSE
            DO i_box = Num_inject+1, Num_Stop, 1
                ALLOCATE(Plume2d_new)
                Num_Plume2d = Num_Plume2d + 1
                Num_Plume2d_acc = Num_Plume2d_acc + 1
                Plume2d_new%IsDissolve      = .False.
                Plume2d_new%IsTransfer      = .False.
                Plume2d_new%label           = Num_Plume2d_acc

                Plume2d_new%LON             = Plume_sources(1)%lon
                Plume2d_new%LAT             = GetInjectionLat(Plume_sources(1)%lat1,Plume_sources(1)%lat2, Length_init, i_box - 1)
                write(6,*) 'debug (BZ): injection lat: ', Plume2d_new%LAT, 'Num_inject: ', i_box
                Plume2d_new%LEV             = Plume_sources(1)%lev
                 
                i_lon = Find_iLonLat(Plume2d_new%LON, DX_GC, X_edge2)
                if(i_lon>NX_GC) i_lon=i_lon-NX_GC
                if(i_lon<1) i_lon=i_lon+NX_GC
                i_lat = Find_iLonLat(Plume2d_new%LAT, DY_GC, Y_edge2)
                if(i_lat>NY_GC) i_lat=NY_GC
                if(i_lat<1) i_lat=1
                i_lev = Find_iPLev(Plume2d_new%LEV,P_edge)
                if(i_lev>NZ_GC) i_lev=NZ_GC
                Plume2d_new%lat_ind = i_lat
                Plume2d_new%lon_ind = i_lon
                Plume2d_new%lev_ind = i_lev

                Plume2d_new%LENGTH = Length_init  !unit: m 
                Plume2d_new%ALPHA  = 0.0e+0_fp
                Plume2d_new%LIFE   = 0.0e+0_fp
                Plume2d_new%PDX    = Dx_init
                Plume2d_new%PDY    = Dy_init

                Vgrid_EU           = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ! [cm3]
                Vgrid_2D           = (Plume2d_new%PDX * Plume2d_new%PDY * Plume2d_new%LENGTH ) *1.E6_fp ! cm3
                ! print *, 'n_x_max=', n_x_max, 'n_y_max=', n_y_max, 'nspc_p=', nspc_p
                ! print *, 'allocated before?', allocated(Plume2d_new%CONCNT2d)
                ALLOCATE(Plume2d_new%CONCNT2d(n_x_max, n_y_max, nspc_p), STAT=RC)
                IF (RC /= 0) THEN
                     errMsg = 'Error allocating Plume2d_new%CONCNT2d'
                     !CALL GC_Error( errMsg, RC, thisLoc )
                     CALL ERROR_STOP( errMsg, thisLoc)
                     !RETURN
                 END IF
                ALLOCATE(Plume2d_new%MassRef2d(nspc_p))
                Plume2d_new%CONCNT2d = 0.0e+0_fp
                ! pass background concentration
                ! DO i_species = 1, n_species
                !    Plume2d_new%CONCNT2d(:,:,i_species) = Spc(i_species)%Conc(i_lon, i_lat, i_lev)  ! molec/cm3
                !    Plume2d_new%MassRef2d(i_species) = Spc(i_species)%Conc(i_lon, i_lat, i_lev) &
                !                                             * Vgrid_2D  *  n_x_max * n_y_max ! molec 
                !ENDDO
                !mass_S_SO2_b             =      mass_S_SO2_b + Plume2d_new%MassRef2d(id_SO2)
                !mass_S_SO4_b             =      mass_S_SO4_b + Plume2d_new%MassRef2d(id_SO4)

                ! pass background concentration
#ifdef TOMAS
                ! Read all bulk tracer and TOMAS tracer 01
                DO i_species = 1, 10 + nspc_p_tomas_tracer
                  !write(6,*) 'debug (BZ): i_species', i_species
                  spc_name = TRIM(spc_names_p(i_species))
                  ! Skip initializing these species:   
                  IF ((spc_name == 'PH2SO4').or. (spc_name == 'NH3') .or. &
                                 (spc_name == 'NH4') .or. (spc_name =='AW01') .or. &
                                  (spc_name =='H2O') .or. (spc_name == 'OH') .or. (spc_name == 'HO2')) CYCLE
                  !write(6,*) 'debug (BZ): spc_name in plume', TRIM(spc_name)
                  id_tracer   = Ind_(TRIM(spc_name))
                  !write(6,*) 'debug (BZ): species id in GEOS-Chem ', id_tracer
                  Plume2d_new%CONCNT2d(:,:,i_species) = Spc(id_tracer)%Conc(i_lon, i_lat, i_lev)
                  Spc(id_tracer)%Conc(i_lon, i_lat, i_lev) = Spc(id_tracer)%Conc(i_lon, i_lat, i_lev) * &
                                                              (1 - (Vgrid_2D*n_x_max*n_y_max)/Vgrid_EU)
                  Plume2d_new%MassRef2d(i_species) = Spc(id_tracer)%Conc(i_lon, i_lat, i_lev)*(n_x_max*n_y_max*Vgrid_2D)
                ENDDO
                ! Read TOMAS tracer 02, 03, ... nbins
                ! In GEOS-Chem the order is tracer1_01, tracer2_02,tracer3_03
                ! in plume model the order is designed to be tracer1_01, tracer2_01, ... , tracer1_02, tracer2_02...
                ! Skip initializing AW
                DO i_species = 1, nspc_p_tomas_tracer -1 
                  ! The first tracer of the same species
                  spc_name = spc_names_p(10 + i_species)
                  DO ibin = 2, nBins
                    !write(6,*) 'debug (BZ): species', TRIM(spc_name), ' ibin= ', ibin
                    ! Tracer in GEOS-Chem
                    id_tracer1   = Ind_(TRIM(spc_name))+ ibin-1
                    !write(6,*) 'debug (BZ): species id in GEOS-Chem ', id_tracer1
                    ! Tracer in plume
                    id_tracer  = 10 + i_species + nspc_p_tomas_tracer * (ibin-1)
                    !write(6,*) 'debug (BZ): species id in plume ', id_tracer
                    Plume2d_new%CONCNT2d(:,:,id_tracer) = Spc(id_tracer1)%Conc(i_lon, i_lat, i_lev)
                    Spc(id_tracer1)%Conc(i_lon, i_lat, i_lev) = Spc(id_tracer1)%Conc(i_lon, i_lat, i_lev) * &
                                                                (1 - (Vgrid_2D*n_x_max*n_y_max)/Vgrid_EU)
                    Plume2d_new%MassRef2d(id_tracer) = Spc(id_tracer1)%Conc(i_lon, i_lat, i_lev)*(n_x_max*n_y_max*Vgrid_2D)
                  ENDDO
                ENDDO
                  
#else
                DO i_species = 1, nspc_p
                  spc_name = spc_names_p(i_species)
                  IF ((spc_name == 'PH2SO4') .or. (spc_name == 'OH') .or. (spc_name == 'HO2'))  CYCLE
                  id_tracer   = Ind_(TRIM(spc_name))
                  Plume2d_new%CONCNT2d(:,:,i_species) = Spc(id_tracer)%Conc(i_lon, i_lat, i_lev)
                  Spc(id_tracer)%Conc(i_lon, i_lat, i_lev) = Spc(id_tracer)%Conc(i_lon, i_lat, i_lev) * &
                                                              (1 - (Vgrid_2D*n_x_max*n_y_max)/Vgrid_EU)
                  Plume2d_new%MassRef2d(i_species) = Spc(id_tracer)%Conc(i_lon, i_lat, i_lev)*(n_x_max*n_y_max*Vgrid_2D)
                ENDDO
#endif
                ! add injected mass
                DO i_species = 1, num_of_sources
                    spc_name = Plume_sources(i_species)%species
                    id_tracer = GET_PLUME_SPC_ID(spc_name)
                    id_tracer1 = Ind_(TRIM(spc_name))
                    IF(id_tracer<0) THEN
                      ErrMsg = 'Species: ' // TRIM(spc_name) // ' not exist in Plume model'
                      CALL ERROR_STOP (ErrMsg, ThisLoc)
                    ELSE
                    Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,id_tracer) = Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,id_tracer)   &
                                + (Plume2d_new%LENGTH * Plume_sources(i_species)%rate / State_Chm%SpcData(id_tracer1)%Info%MW_g * Avo) &
                                /(Plume2d_new%PDX * Plume2d_new%PDY *Plume2d_new%LENGTH*1.E6_fp ) ! molec/cm3
                    ENDIF                                
                ENDDO
                mass_S_SO2_inj_2D = mass_S_SO2_inj_2D + &
                        (Plume2d_new%LENGTH * Plume_sources(1)%rate) / 64.0_fp * Avo 
                mass_S_SO2_1_2D = mass_S_SO2_1_2D + SUM(Plume2d_new%CONCNT2d(:,:,id_SO2_p)) *Vgrid_2D
                mass_S_SO4_1_2D = mass_S_SO4_1_2D + SUM(Plume2d_new%CONCNT2d(:,:,id_SO4_p)) *Vgrid_2D

                OPEN( File_spc_init_IU_2D,      FILE=TRIM( file_spc_init_2D   ), STATUS='OLD',  &
                     POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
                WRITE(File_spc_init_IU_2D, '(*(g0,:,","))')  Plume2d_new%label, Plume2d_new%MassRef2d(1:nspc_p)

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
    !WRITE(6,'(a)') 'debug: after plume injection, the unit is ' // TRIM(UNIT_STR(State_Chm%Species(id_SO4)%Units))

  END SUBROUTINE plume_inject_box

  SUBROUTINE plume_model_box(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
  
    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Chm_Mod,   ONLY : ChmState, Ind_
    USE State_Met_Mod,   ONLY : MetState
    USE Species_Mod,     ONLY : SpcConc
    USE TIME_MOD,        ONLY : GET_TS_DYN, GET_TS_CHEM
    USE TIME_MOD,        ONLY : GET_YEAR, GET_MONTH, GET_DAY, GET_HOUR, GET_MINUTE, GET_SECOND
    USE TIME_MOD,        ONLY : ITS_TIME_FOR_CHEM, ITS_TIME_FOR_DYN, ITS_TIME_FOR_EXIT
    USE State_Grid_Mod,  ONLY : GrdState
    !USE State_Diag_Mod,           ONLY : DgnState
    !USE State_Diag_Mod,           ONLY : DgnMap
    ! USE UnitConv_Mod

    LOGICAL, INTENT(IN)           :: am_I_Root
    TYPE(MetState), INTENT(IN)    :: State_Met
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State objectgg
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    !TYPE(DgnState), INTENT(INOUT) :: State_Diag ! Diagnostics State object
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    INTEGER                :: previous_units, previous_units_temp
    INTEGER                :: ibin
    INTEGER                :: this_year, this_month, this_day, this_hour, this_minute, this_second

    REAL(fp)               :: Dt_dyn, Dt_chem
    !REAL(fp)               :: this_tau, this_taub
    LOGICAL                :: exe_dyn, exe_chem, exe_exit
    !CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: Datestr
    CHARACTER(LEN=255)     :: spc_name
    CHARACTER(LEN=255)     :: ErrMsg
    CHARACTER(LEN=255)     :: ThisLoc
    
    TYPE(SpcConc), POINTER :: Spc(:)
    !TYPE(Plume2d_list), POINTER :: Plume2d_tail=> NULL(), Plume2d_head=> NULL()

    this_year                     =      GET_YEAR()
    this_month                    =      GET_MONTH()
    this_day                      =      GET_DAY()
    this_hour                     =      GET_HOUR()
    this_minute                   =      GET_MINUTE()
    this_second                   =      GET_SECOND()
    WRITE(Datestr,'(I4.4,I2.2,I2.2,I2.2,I2.2,I2.2)') &
                  this_year, this_month, this_day, &
                  this_hour, this_minute, this_second
    !WRITE(6,*) 'Debug (BZ): Plume model: Current simulation time is : ', this_year, this_month, &
    !this_day, this_hour, this_minute, this_second

    exe_dyn            =   ITS_TIME_FOR_DYN()
    exe_chem           =   ITS_TIME_FOR_CHEM()
    exe_exit           =   ITS_TIME_FOR_EXIT()

    !IF (exe_dyn)  WRITE(6, *)  "Debug (BZ): Plume model: Is time for dynamic"
    !IF (exe_chem) WRITE(6, *)  "Debug (BZ): Plume model: Is time for chem"
    !IF (exe_exit) WRITE(6, *)  "Debug (BZ): Plume model: Is time for exit"

    ErrMsg                 =    ''
    Dt_dyn                 =    GET_TS_DYN()
    Dt_chem                 =   GET_TS_CHEM()
    Spc                    =>   State_Chm%Species
    ThisLoc                =   ' -> at plume_model_box (in module GeosCore/lagrange_singlebox_mod.F90)'

    id_SO4     =    Ind_('SO4')
    id_SO2     =    Ind_('SO2')
    id_OH      =    Ind_('OH')
    id_HO2     =    Ind_('HO2')
    id_PH2SO4  =    Ind_('PH2SO4')
    id_NH3     =    Ind_('NH3')
    id_NH4     =    Ind_('NH4')
    id_O3      =    Ind_('O3')
    id_H2O     =    Ind_('H2O')
    id_H2SO4   =    Ind_('H2SO4')
    id_NK01    =    Ind_('NK01')
    id_SF01    =    Ind_('SF01')
    id_AW01    =    Ind_('AW01')
    id_SO2pl   =    Ind_('SO2pl') 
    id_SO4pl   =    Ind_('SO4pl')
    !this_tau  = GET_TAU()
    !this_taub = GET_TAUb()
    !this_year = GET_YEAR()
    !this_month= GET_MONTH()

    
    !WRITE(6,'(a)') 'debug (BZ): TAU = ',  this_tau, 'TAUb = ', this_taub
    !WRITE(6,'(a)') 'debug (BZ): this_year = ',  this_year, 'this_month = ', this_month

   IF (use_lagrange .AND. plume_inject_on) THEN
      IF (.NOT. ASSOCIATED(Plume2d_head) .AND. &
               .NOT. ASSOCIATED(Plume1d_head)) RETURN
      IF (exe_dyn) THEN
         !WRITE(6, *)  "Debug (BZ): Plume model: Is time for dynamic"
         time_elapsed = time_elapsed + Dt_dyn
         ! In theory, before plume injection, the unit should be kg/kg
         CALL Convert_Spc_Units(                                            &
               Input_Opt      = Input_Opt,                                   &
               State_Chm      = State_Chm,                                   &
               State_Grid     = State_Grid,                                  &
               State_Met      = State_Met,                                   &
               new_units      = MOLECULES_SPECIES_PER_CM3,                   &
               previous_units = previous_units_temp,                              &
               RC             = RC                                          )
         ! WRITE(6,'(a)') 'debug : Before plume box model, the unit is ' // TRIM(UNIT_STR(previous_units_temp))
         ! WRITE(6,'(a)') 'debug : During plume box model, the unit is ' // TRIM(UNIT_STR(State_Chm%Species(id_SO4)%Units))

         ! Check that species units are in  molec/cm3
         IF ( Spc(id_SO4)%Units /= MOLECULES_SPECIES_PER_CM3 ) THEN
         ErrMsg = 'Incorrect species units: ' // TRIM(UNIT_STR(Spc(id_SO4)%Units))
         !CALL GC_Error( ErrMsg, RC, ThisLoc )
         CALL ERROR_STOP (ErrMsg, ThisLoc)
         ENDIF
         mass_S_SO2_r1_2D       =      0.0_fp
         mass_S_SO4_r1_2D       =      0.0_fp
         mass_S_SO2_r1_1D       =      0.0_fp
         mass_S_SO4_r1_1D       =      0.0_fp
         mass_S_SO2_r2_2D       =      0.0_fp
         mass_S_SO4_r2_2D       =      0.0_fp
         mass_S_SO2_r2_1D       =      0.0_fp
         mass_S_SO4_r2_1D       =      0.0_fp
         mass_S_SO2_r3_2D       =      0.0_fp
         mass_S_SO4_r3_2D       =      0.0_fp
         mass_S_SO2_r3_1D       =      0.0_fp
         mass_S_SO4_r3_1D       =      0.0_fp
         mass_S_SO2_r4          =      0.0_fp
         mass_S_SO4_r4          =      0.0_fp

         mass_S_SO2_7_2D        =      0.0_fp
         mass_S_SO4_7_2D        =      0.0_fp
         mass_S_SO2_8_2D        =      0.0_fp
         mass_S_SO4_8_2D        =      0.0_fp
         mass_S_SO2_7_1D        =      0.0_fp
         mass_S_SO4_7_1D        =      0.0_fp
         Vgrid_2D_tot_3         =      0.0_fp
         Vgrid_2D_tot_4         =      0.0_fp
         Vgrid_1D_tot_3         =      0.0_fp
#ifdef TOMAS
         mass_S_H2SO4_2D       = 0.0_fp
         mass_S_H2SO4_1D       = 0.0_fp
         mass_SF_bin_2D(:)     = 0.0_fp
         mass_NK_bin_2D(:)     = 0.0_fp
         mass_MK_bin_2D(:)     = 0.0_fp
         mass_SF_bin_1D(:)     = 0.0_fp
         mass_NK_bin_1D(:)     = 0.0_fp
         mass_MK_bin_1D(:)     = 0.0_fp
         ! mass_SF_bin_2D_1D(:)  = 0.0_fp
         ! mass_NK_bin_2D_1D(:)  = 0.0_fp
#endif

         Vplume_2D_tot(:,:,:)  = 0.0_fp
         Vplume_1D_tot(:,:,:)  = 0.0_fp

         CALL plume_physics(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
         ! Plume physical evolution: 
         ! Update plume lifetime 
         ! movement (advection); plume stretching; adiabatic volume change due to thermodynamics (p-T-V); 
         ! Update in-plume concentration due to volume change 
         ! In-plume species concentration: 
         ! advection, diffusion; entrainment


         IF (exe_chem) THEN
            !WRITE(6, *)  "Debug (BZ): Plume model: Is time for chem"
            ! BZ: See Chemistry_mod.F90
            ! Before Chemistry, do I need to set CO2 to 421ppm?This is set for Eulerian grid before chemistry
            ! This is necessary to reduce the error norm in KPP.
            ! See https://github.com/geoschem/geos-chem/issues/1529.
            CALL plume_chem_microphysics(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
            ! New module update chemistry and microphysics
         ENDIF

         CALL plume_structure_change(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
         ! Plume structure evolution:
         ! 2D to 1D
         ! 1D seg splitting
         ! Dissolve when meet the criteria 

         ! write diagnostic output
         CALL lagrange_write_std( time_elapsed )

         ! convert unit back
         CALL Convert_Spc_Units(                                            &
                  Input_Opt      = Input_Opt,                                   &
                  State_Chm      = State_Chm,                                   &
                  State_Grid     = State_Grid,                                  &
                  State_Met      = State_Met,                                   &
                  new_units      = previous_units_temp,                   &
                  RC             = RC                                          )
      ENDIF

   ELSE
      RETURN
   ENDIF
    
  END SUBROUTINE plume_model_box
  
  SUBROUTINE plume_mod_cleanup_box( RC)

     ! USE ErrCode_Mod

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

   !  IF ( ALLOCATED(SpcConc_BEFORE_KPP))  THEN
   !        DEALLOCATE(SpcConc_BEFORE_KPP, STAT=RC )
   !        CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:SpcConc_BEFORE_KPP', 2, RC )
   !        IF ( RC /= GC_SUCCESS ) RETURN
   !  ENDIF

   !  IF ( ALLOCATED(SpcConc_AFTER_KPP))  THEN
   !        DEALLOCATE(SpcConc_AFTER_KPP, STAT=RC )
   !        CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:SpcConc_AFTER_KPP', 2, RC )
   !        IF ( RC /= GC_SUCCESS ) RETURN
   !  ENDIF
    IF ( ALLOCATED(Vplume_2D_tot))  THEN
          DEALLOCATE(Vplume_2D_tot, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:Vplume_2D_tot', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    IF ( ALLOCATED(Vplume_1D_tot))  THEN
          DEALLOCATE(Vplume_1D_tot, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:Vplume_1D_tot', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    IF ( ALLOCATED(Plume_sources))  THEN
          DEALLOCATE(Plume_sources, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:Plume_sources', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    IF ( ALLOCATED(Xk))  THEN
          DEALLOCATE(Xk, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:Xk', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    IF ( ALLOCATED(AVGMASS))  THEN
          DEALLOCATE(AVGMASS, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:AVGMASS', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    IF ( ALLOCATED(spc_names_p_use))  THEN
          DEALLOCATE(spc_names_p_use, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:spc_names_p_use', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    IF ( ALLOCATED(mass_spc_init_2D))  THEN
          DEALLOCATE(mass_spc_init_2D, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:mass_spc_init_2D', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    IF ( ALLOCATED(mass_spc_init_1D))  THEN
          DEALLOCATE(mass_spc_init_1D, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:mass_spc_init_1D', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    IF ( ALLOCATED(mass_SF_bin_2D))  THEN
          DEALLOCATE(mass_SF_bin_2D, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:mass_SF_bin_2D', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    IF ( ALLOCATED(mass_SF_bin_1D))  THEN
          DEALLOCATE(mass_SF_bin_1D, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:mass_SF_bin_1D', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    IF ( ALLOCATED(mass_SF_bin_2D_1D))  THEN
          DEALLOCATE(mass_SF_bin_2D_1D, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:mass_SF_bin_2D_1D', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    IF ( ALLOCATED(mass_NK_bin_2D))  THEN
          DEALLOCATE(mass_NK_bin_2D, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:mass_NK_bin_2D', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF 
    IF ( ALLOCATED(mass_NK_bin_1D))  THEN
          DEALLOCATE(mass_NK_bin_1D, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:mass_NK_bin_1D', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    IF ( ALLOCATED(mass_NK_bin_2D_1D))  THEN
          DEALLOCATE(mass_NK_bin_2D_1D, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:mass_NK_bin_2D_1D', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF 
    IF ( ALLOCATED(mass_MK_bin_2D))  THEN
          DEALLOCATE(mass_MK_bin_2D, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:mass_MK_bin_2D', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    IF ( ALLOCATED(mass_MK_bin_1D))  THEN
          DEALLOCATE(mass_MK_bin_1D, STAT=RC )
          CALL GC_CheckVar( 'lagrange_singlebox_mod.F90:mass_MK_bin_1D', 2, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    
    ! Close all the files
    CLOSE(File_Smass_IU_2D)
    CLOSE(File_Smass_IU_1D)
    CLOSE(File_spc_init_IU_2D)
    CLOSE(File_spc_init_IU_1D)
    CLOSE(File_Plume_life_IU_2D)
    CLOSE(File_Plume_life_IU_1D)
    CLOSE(File_Plume_number_IU)
    CLOSE(File_Plume_location_IU_2D)
    CLOSE(File_Plume_location_IU_1D)
    CLOSE(File_SF_bin_IU_2D)
    CLOSE(File_SF_bin_IU_1D)
    CLOSE(File_SF_bin_IU_2D_1D)
    CLOSE(File_NK_bin_IU_2D)
    CLOSE(File_NK_bin_IU_1D)
    CLOSE(File_NK_bin_IU_2D_1D)
    CLOSE(File_MK_bin_IU_2D)
    CLOSE(File_MK_bin_IU_1D)
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
  ! All the use have been defined in host: Plume_box_model
  USE Input_Opt_Mod,   ONLY : OptInput, PlumeSource_t
  USE State_Chm_Mod,   ONLY : ChmState, Ind_
  USE State_Met_Mod,   ONLY : MetState
  USE Species_Mod,     ONLY : SpcConc, Species
  USE TIME_MOD,        ONLY : GET_TS_DYN
  USE State_Grid_Mod,  ONLY : GrdState
  USE InquireMod,      ONLY : findFreeLun
  !USE UnitConv_Mod

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
  INTEGER                :: i_species, i_tracer, ind_spc_GC, ind_spc_GC_bin1, ind_spc_p
  INTEGER                :: i_x, i_y
  INTEGER                :: i_slab
  INTEGER                :: Nt
  INTEGER                :: t1s
  INTEGER                :: OrigUnit
  INTEGER                :: box_label
  INTEGER                :: this_bin, this_tracer
  INTEGER                :: ibin

  REAL(fp)               :: MW_g
  REAL(fp)               :: Dt, Dt2
  REAL(fp)               :: ratio
  REAL(fp)               :: V_prev, V_new
  REAL(fp)               :: Vgrid_EU,  Vgrid_2D, Vgrid_1D, Vgrid_1D_new
  REAL(fp)               :: Ly ! Lyapunov exponent [s-1]
  REAL(fp)               :: length0
  REAL(fp)               :: Pdx, Pdy, Pdt
  REAL(fp)               :: box_lon, box_lat, box_lev                          
  REAL(fp)               :: box_length, box_alpha, theta_previous
  REAL(fp)               :: box_Ra, box_Rb, box_theta
  REAL(fp)               :: box_RA_tmp, box_Rb_tmp, box_theta_tmp
  REAL(fp)               :: box_RA_test, box_Rb_test, box_theta_test
  REAL(fp)               :: box_Ra_candidate, box_Rb_candidate, box_theta_candidate
  REAL(fp)               :: box_RA_final, box_Rb_final, box_theta_final
  !REAL(fp)               :: box_Ra(nspc_p), box_Rb(nspc_p), box_theta(nspc_p)
  REAL(fp)               :: box_life
  REAL(fp)               :: box_x_PS, box_y_PS
  REAL(fp)               :: box_u, box_v, box_omeg
  REAL(fp)               :: dbox_lon, dbox_lat, dbox_lev
  REAL(fp)               :: dbox_x_PS, dbox_y_PS
  REAL(fp)               :: curr_lon, curr_lat, curr_pressure
  REAL(fp)               :: curr_u, curr_v, curr_omeg, curr_ptemp
  REAL(fp)               :: curr_u_PS, curr_v_PS
  REAL(fp)               :: curr_T1, next_T2 
  REAL(fp)               :: RK_Dx_PS, RK_Dy_PS
  REAL(fp)               :: RK_x_PS, RK_y_PS
  REAL(fp)               :: RK_lon, RK_lat
  REAL(fp)               :: wind_s_shear, Ptemp_shear
  REAL(fp)               :: Cv, Ch, Omega_N, N_BV
  REAL(fp)               :: eddy_v, eddy_h, eddy_D
  REAL(fp)               :: CFL, CFL_adv, CFL_dif_h, CFL_dif_v, max_u, CFL_dif_D, adv_term 
  REAL(fp)               :: Pc_middle, Pc_bottom, Pc_top, Pc_left, Pc_right, Pc_update
  REAL(fp)               :: TOMAS_Scale
  REAL(fp)               :: mass_before_clip, mass_after_clip
  REAL(fp)               :: mass_before_clip_1D(nspc_p), mass_after_clip_1D(nspc_p)
  REAL(fp)               :: background_conc, background_conc_NK, background_conc_SF
  REAL(fp)               :: mass_plume, mass_plume_new, D_mass_plume, mass_plume_scale
  REAL(fp)               :: background_mass, background_mass_new
  REAL(fp)               :: excess_mass
  ! real(fp)               :: Massref_species
  REAL(fp)               :: RK_Dt(5)
  REAL(fp)               :: RK_u(4), RK_v(4), RK_omeg(4)
  REAL(fp)               :: RK_Dlon(4), RK_Dlat(4), RK_Dlev(4)
  REAL(fp)               :: Pu(n_x_max2,n_y_max2) ! for 2D advection
  !real(fp)               :: Pc(n_x_max,n_y_max), Pc2(n_x_max,n_y_max) !, Ec(n_x_max,n_y_max)
  !real(fp)               :: Pc_bdy(n_x_max2,n_y_max2)
  REAL(fp)               :: C2d_prev(n_x_max,n_y_max) !, C2d_prev_extra(n_x_max,n_y_max)
  REAL(fp)               :: C2d_new(n_x_max2,n_y_max2)
  REAL(fp)               :: C2d_bg(n_x_max2,n_y_max2)
  REAL(fp)               :: C1d_bg(n_slab_max2), C1d_new(n_slab_max2), C1d_prev(n_slab_max)

#ifdef TOMAS
  REAL(fp)               :: C2d_prev_NK(n_x_max,n_y_max), C2d_prev_SF(n_x_max,n_y_max) 
  REAL(fp)               :: C2d_new_NK(n_x_max,n_y_max), C2d_new_SF(n_x_max,n_y_max)
  REAL(fp)               :: C1d_prev_NK(n_slab_max), C1d_prev_SF(n_slab_max) 
  REAL(fp)               :: C1d_new_NK(n_slab_max), C1d_new_SF(n_slab_max)
#endif
  REAL(fp)               :: Xscale, Yscale
  

  REAL(fp), dimension(:,:,:), allocatable :: box_concnt_2D
  REAL(fp), dimension(:,:), allocatable :: box_concnt_1D, box_concnt_1D_prev
  REAL(fp), POINTER      :: u(:,:,:)
  REAL(fp), POINTER      :: v(:,:,:)
  REAL(fp), POINTER      :: omeg(:,:,:)
  REAL(fp), POINTER      :: Ptemp(:,:,:)
  REAL(fp), POINTER      :: T1(:,:,:)
  REAL(fp), POINTER      :: T2(:,:,:)
  REAL(fp), POINTER      :: P_BXHEIGHT(:,:,:)
  REAL(fp), POINTER      :: TROPP(:,:)

  LOGICAL                :: CFL_ok

  TYPE(SpcConc), POINTER        :: Spc(:)
  TYPE(Species), POINTER        :: SpcInfo
  TYPE(Plume2d_list), POINTER   :: Plume2d_new, Plume2d_curr, Plume2d_prev
  TYPE(Plume1d_list), POINTER   :: Plume1d_new, Plume1d_curr, Plume1d_prev


  
  CHARACTER(LEN=255)     :: ErrMsg
  CHARACTER(LEN=255)     :: ThisLoc
  CHARACTER(LEN=255)     :: spc_name
  

  ! Variables for diagnostic file
  INTEGER                :: file_2Dconc_NK01_ID_1,    file_2Dconc_NK01_ID_2
  CHARACTER(LEN=255)     :: file_2Dconc_NK01_1,       file_2Dconc_NK01_2

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
  SpcInfo                =>   NULL()
  Dt = GET_TS_DYN()
  ThisLoc                =   ' -> at plume_physics (in module GeosCore/lagrange_singlebox_mod.F90)'
  RC     =  GC_SUCCESS
  ErrMsg = ''
  NULLIFY(Plume2d_new, Plume2d_curr, Plume2d_prev)
  NULLIFY(Plume1d_new, Plume1d_curr, Plume1d_prev)
  !IF(Stop_inject==1) GOTO 400 ! deallocate and nullify -> exit
  
  ALLOCATE(box_concnt_2D(n_x_max, n_y_max, nspc_p )) !
  ALLOCATE(box_concnt_1D(n_slab_max, nspc_p ))
  ALLOCATE(box_concnt_1D_prev(n_slab_max, nspc_p ))

  RK_Dt(1) = 0.0
  RK_Dt(2) = 0.5*Dt
  RK_Dt(3) = 0.5*Dt
  RK_Dt(4) = Dt
  RK_Dt(5) = 0.0
  
  u                  =>     State_Met%U   ! m/s
  v                  =>     State_Met%V   ! m/s
  omeg               =>     State_Met%OMEGA     ! Updraft velocity [Pa/s]
  Ptemp              =>     State_Met%THETA     ! Potential temperature [K]
  T1                 =>     State_Met%TMPU1     ! Temperature at start of timestep [K]
  T2                 =>     State_Met%TMPU2     ! Temperature at end of timestep [K]
  P_BXHEIGHT         =>     State_Met%BXHEIGHT  ![NX_GC,NY_GC,KKPAR]
  TROPP              =>     State_Met%TROPP     ! tropopause pressure, hpa


  !=======================================================================
  ! for 2D plume: Run Lagrangian trajectory-track HERE
  !=======================================================================

  IF(.NOT.ASSOCIATED(Plume2d_head)) GOTO 401

  Plume2d_curr => Plume2d_head
  DO WHILE(ASSOCIATED(Plume2d_curr))
   
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
    
    Vgrid_2D       = Pdx*Pdy*box_length*1.0e+6_fp ! [cm3]
   !  Vgrid_2D_tot_1 = Vgrid_2D_tot_1 + Vgrid_2D * n_x_max * n_y_max
    !write(6,*) 'debug (BZ): solve plume physics in 2-D plume box: ', box_label
   !  mass_S_SO2_2_2D = mass_S_SO2_2_2D + & 
   !                      SUM(Plume2d_curr%CONCNT2d(:, :, id_SO2_p)) * Vgrid_2D 
               
   !  mass_S_SO4_2_2D = mass_S_SO4_2_2D +  &
   !       SUM(Plume2d_curr%CONCNT2d(:,:, id_SO4_p)) * Vgrid_2D
   !  Write (6, *) 'Debug: BZ:  (2-D phys before - 1), SO2 mass = ', SUM(box_concnt_2D(:,:, id_SO2_p)) * Vgrid_2D, &
   !             'SO4 mass = ',  SUM(box_concnt_2D(:,:, id_SO4_p)) * Vgrid_2D
   !  Write (6, *) 'Debug: BZ:  (2-D phys before - 2), SO2 mass = ', SUM(Plume2d_curr%CONCNT2d(:,:, id_SO2_p)) * Vgrid_2D, &
   !             'SO4 mass = ',  SUM(Plume2d_curr%CONCNT2d(:,:, id_SO4_p)) * Vgrid_2D 
    !write(6,*) 'debug (BZ): SO2mass before physics: ', SUM(box_concnt_2D(:,:,id_SO2_p)) *  Vgrid_2D
    !write(6,*) 'debug (BZ): SO2conc before physics: ', SUM(box_concnt_2D(:,:,id_SO2_p)) /  (n_x_max*n_y_max)
    !write(6,*) 'debug (BZ): SO2conc Euleria Grid: ', Spc(id_SO2)%Conc(i_lon,i_lat,i_lev)
! #ifdef TOMAS
!    Write (6, *) 'Debug: BZ:  (2-D test grid) before physics, SO4 (molec)= ', box_concnt_2D(x_test,y_test, id_SO4_p)
!    Write (6, *) 'Debug: BZ: (2-D test grid) before physics, Nk bin 15 (molec)= ', box_concnt_2D(x_test,y_test, 53)
!    Write (6, *) 'Debug: BZ: (2-D test grid) before physics, Nk bin 15 (molec)= ', Plume2d_curr%CONCNT2d(x_test,y_test, 53)
!    Write (6, *) 'Debug: BZ: (2-D test grid) before physics SF bin 15 (molec)= ', box_concnt_2D(x_test,y_test, 54)
! #endif
    ! mass_S_SO2_1_1  = mass_S_SO2_1_1 + SUM(box_concnt_2D(:,:,id_SO2)) *  Vgrid_2D
    ! mass_S_SO4_1_1  = mass_S_SO4_1_1 + SUM(box_concnt_2D(:,:,id_SO4)) *  Vgrid_2D
    
    box_life = box_life + Dt

    curr_lon      = box_lon
    curr_lat      = box_lat
    curr_pressure = box_lev
    DO Ki = 1,4,1
      !------------------------------------------------------------------
      ! For vertical wind speed:
      ! pay attention for the polar region * * *
      !------------------------------------------------------------------
      if(abs(curr_lat)>Y_mid(NY_GC))then
        curr_omeg = Interplt_wind_RLL_polar(omeg, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
      else
        curr_omeg = Interplt_wind_RLL(omeg, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
      endif

      RK_omeg(Ki) = curr_omeg
      RK_Dlev(Ki)   = Dt * curr_omeg / 100.0     ! Pa => hPa

      curr_pressure = box_lev + RK_Dt(Ki+1) * curr_omeg / 100.0

      if(curr_pressure<P_mid(NZ_GC)) &
            curr_pressure = P_mid(NZ_GC) !+ ( P_mid(NZ_GC) - curr_pressure )
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

        if(abs(curr_lat)>Y_mid(NY_GC))then 
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
    do while (box_lat > Y_edge(NY_GC+1))
        box_lat = Y_edge(NY_GC+1) - ( box_lat-Y_edge(NY_GC+1) )
    end do
    do while (box_lat < Y_edge(1))
        box_lat = Y_edge(1) - ( box_lat-Y_edge(1) )
    end do

    do while (box_lon > X_edge(NX_GC+1))
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
    if(abs(curr_lat)>Y_mid(NY_GC))then
        curr_T1 = Interplt_wind_RLL_polar(T1, i_lon, i_lat, &
                                i_lev, curr_lon, curr_lat, curr_pressure)
    else
        curr_T1 = Interplt_wind_RLL(T1, i_lon, i_lat, i_lev, &
                                      curr_lon, curr_lat, curr_pressure)
    endif

    ! update index based on new location
    i_lon = Find_iLonLat(box_lon, DX_GC, X_edge2)
    if(i_lon>NX_GC) i_lon=i_lon-NX_GC
    if(i_lon<1) i_lon=i_lon+NX_GC

    i_lat = Find_iLonLat(box_lat, DY_GC, Y_edge2)
    if(i_lat>NY_GC) i_lat=NY_GC
    if(i_lat<1) i_lat=1

    i_lev = Find_iPLev(box_lev,P_edge)
    if(i_lev>NZ_GC) i_lev=NZ_GC


    if(abs(box_lat)>Y_mid(NY_GC))then
        next_T2 = Interplt_wind_RLL_polar(T2, i_lon, i_lat, i_lev, box_lon, box_lat, box_lev)
    else
        next_T2 = Interplt_wind_RLL(T2, i_lon, i_lat, i_lev, box_lon, box_lat, box_lev)
    endif

    ! PV=nRT, V2 = T2/P2 : T1/P1 * V1
    ratio = SQRT( (next_T2/box_lev)/(curr_T1/curr_pressure) )

    ! assume the volume change mainly apply to the cross-section, 
    ! the box_length would not change 

    V_prev = Pdx*Pdy*box_length*1.0e+6_fp
    V_new  = V_prev*ratio*ratio

    Pdx = Pdx *ratio
    Pdy = Pdy *ratio

    DO i_species = 1, nspc_p, 1
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
    IF((box_u**2.0_fp+box_v**2.0_fp)==0)then
        box_alpha = 0.0
    ELSE
      IF(box_v>=0)THEN
        box_alpha = ACOS( box_u/SQRT(box_u**2.0_fp+box_v**2.0_fp) ) 
      ELSE
        box_alpha = 2.0_fp*PI - ACOS( box_u/SQRT(box_u**2.0_fp+box_v**2.0_fp) )
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
    Plume2d_curr%LON         = box_lon
    Plume2d_curr%LAT         = box_lat
    Plume2d_curr%LEV         = box_lev
    Plume2d_curr%lat_ind     = i_lat
    Plume2d_curr%lon_ind     = i_lon
    Plume2d_curr%lev_ind     = i_lev
    Plume2d_curr%LENGTH      = box_length
    Plume2d_curr%ALPHA       = box_alpha
    Plume2d_curr%label       = box_label
    Plume2d_curr%LIFE        = box_life
    Plume2d_curr%PDX         = Pdx
    Plume2d_curr%PDY         = Pdy
    Plume2d_curr%CONCNT2d    = box_concnt_2D

    Vplume_2D_tot(i_lon, i_lat, i_lev) = &
               Vplume_2D_tot(i_lon, i_lat, i_lev)  +  &
               (Pdx*Pdy*box_length*n_x_max*n_y_max) * 1.0e+6_fp ! [cm3]
      ! file_2Dconc_SO2_ID_1 = findFreeLun()
      ! WRITE(file_2Dconc_SO2_1,'("Plume-2D_SO2_conc_",I0,"_1.txt")') NINT(time_elapsed)
      ! CALL PLUME_CONC_DIAG_FILES_2D(file_2Dconc_SO2_ID_1, file_2Dconc_SO2_1, Plume2d_curr%CONCNT2d(:,:, id_SO2_p), RC)

      ! file_2Dconc_SO4_ID_1 = findFreeLun()
      ! WRITE(file_2Dconc_SO4_1,'("Plume-2D_SO4_conc_",I0,"_1.txt")') NINT(time_elapsed)
      ! CALL PLUME_CONC_DIAG_FILES_2D(file_2Dconc_SO4_ID_1, file_2Dconc_SO4_1, Plume2d_curr%CONCNT2d(:,:, id_SO4_p), RC)


      ! file_2Dconc_OH_ID_1 = findFreeLun()
      ! WRITE(file_2Dconc_OH_1,'("Plume-2D_OH_conc_",I0,"_1.txt")') NINT(time_elapsed)
      ! OPEN(file_2Dconc_OH_ID_1, FILE=TRIM(file_2Dconc_OH_1), STATUS='REPLACE', &
      !    FORM='FORMATTED', ACCESS='SEQUENTIAL', IOSTAT=RC)
      ! DO i_x = 1, n_x_max
      !    WRITE(file_2Dconc_OH_ID_1,'(*(ES12.4,1X))') &
      !       (Plume2d_curr%CONCNT2d(i_x,i_y,id_OH_p), i_y = 1, n_y_max)
      ! ENDDO
      ! CLOSE(file_2Dconc_OH_ID_1)

    ! write(6,*) 'debug (BZ): SO2mass middle physics: ', SUM(box_concnt_2D(:,:,id_SO2)) * V_grid_2D
    ! write(6,*) 'debug (BZ): SO4mass middle physics: ', SUM(box_concnt_2D(:,:,id_SO4)) * V_grid_2D
    !write(6,*) 'debug (BZ): SO2mass after plume shape change: ', SUM(box_concnt_2D(:,:,id_SO2_p)) *  Vgrid_2D
    !write(6,*) 'debug (BZ): SO4mass after plume shape change: ', SUM(box_concnt_2D(:,:,id_SO4_p)) *  Vgrid_2D
    !write(6,*) 'debug (BZ): SO2conc after plume shape change: ', SUM(box_concnt_2D(:,:,id_SO2_p)) /  (n_x_max*n_y_max)
    !write(6,*) 'debug (BZ): SO4conc after plume shape change: ', SUM(box_concnt_2D(:,:,id_SO4_p)) /  (n_x_max*n_y_max)

! #ifdef TOMAS
!    Write (6, *) 'Debug: BZ:  (2-D test grid) after plume shape change, SO4 (molec)= ', box_concnt_2D(x_test,y_test, id_SO4_p)
!    Write (6, *) 'Debug: BZ: (2-D test grid)after plume shape change, Nk bin 15 (molec)= ', box_concnt_2D(x_test,y_test, 53)
!    Write (6, *) 'Debug: BZ: (2-D test grid)after plume shape change, Nk bin 15 (molec)= ', Plume2d_curr%CONCNT2d(x_test,y_test, 53)
!    Write (6, *) 'Debug: BZ: (2-D test grid) after plume shape change, SF bin 15 (molec)= ', box_concnt_2D(x_test,y_test, 54)
! #endif
    !-------- solve in-plume concentration here --------
    ! advection and diffusion

    curr_lon         = box_lon
    curr_lat         = box_lat
    curr_pressure    = box_lev      ! hPa
    
    Vgrid_EU         = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ! [cm3]
    Vgrid_2D         = Pdx*Pdy*box_length*1.0e+6_fp ! [cm3]
    ! Vgrid_2D_tot_2   = Vgrid_2D_tot_2 + Vgrid_2D * n_x_max * n_y_max
    !====================================================================
    ! calculate the wind shear along plume corss-section
    ! clock-wise 90 degree from the plume length direction (box_alpha)
    ! calculate the diffusivity in horizontal and vertical direction
    !====================================================================

    ! calculate the wind_s shear along pressure direction
    wind_s_shear = Wind_shear_s(u, v, P_BXHEIGHT, box_alpha, i_lon, &
                      i_lat, i_lev,curr_lon, curr_lat, curr_pressure)
    !Write(6, *) "(Debug: BZ) Wind shear horizontal is:  ", wind_s_shear
    ! Calculate vertical eddy diffusivity (U.Schumann, 2012) :
    Cv = 0.2
    Omega_N = 0.1
    Ptemp_shear = Vertical_shear(Ptemp, P_BXHEIGHT, i_lon, i_lat, &
                          i_lev,curr_lon, curr_lat, curr_pressure)
    !Write(6, *) "(Debug: BZ) Wind shear verticle is:  ", Ptemp_shear
    !--------------------------------------------------------------------
    ! interpolate potential temperature for Plume module:
    !--------------------------------------------------------------------
    IF(abs(curr_lat)>Y_mid(NY_GC))then
      curr_Ptemp = Interplt_wind_RLL_polar(Ptemp, i_lon, i_lat, &
                              i_lev, curr_lon, curr_lat, curr_pressure)
    ELSE
      curr_Ptemp = Interplt_wind_RLL(Ptemp, i_lon, i_lat, i_lev, &
                                      curr_lon, curr_lat, curr_pressure)
    ENDIF

    
    IF(N_BV<=0.001) N_BV = 0.001
    N_BV = SQRT(Ptemp_shear*g0/curr_Ptemp)
    ! diffusivity unit: [m2/s]
    eddy_v = Cv * Omega_N**2 / N_BV
    eddy_h = 10.0
    !Write(6, *) "(Debug: BZ) eddy_v is:  ", eddy_v

    ! Define the wind field based on wind shear
    DO i_y = 1, n_y_max2
      Pu(:,i_y) = (i_y-n_y_mid2)*Pdy * wind_s_shear ! [m s-1]
    ENDDO
   !  Write(6, *) "(Debug: BZ) wind field based on shear is: Pu(x_test, y_test) = ", Pu(x_test+1, y_test+1), &
   !    "Pu(x_test-1, y_test) =", Pu(x_test, y_test+1), &
   !    "Pu(x_test+1, y_test) =", Pu(x_test+2, y_test+1), &
   !    "Pu(x_test, y_test-1) =", Pu(x_test+1, y_test), &
   !    "Pu(x_test, y_test+1) =", Pu(x_test+1, y_test+2)

    !--------------------------------------------------------------
    ! if Pdx become smaller than half Dx_init,
    ! combine 9 grids into 1 grids. 
    ! Update: box_concnt_2D(), Pdx(), Pdy(),  Extra_mass_2D()
    ! (BZ): The logic here need to be clarified, temporiraly disable
    ! resize and throw an error when not meet CFL condition
    !--------------------------------------------------------------
   !  IF(Pdx<=0.5*Dx_init) THEN
   !    errMsg = 'Size are too small to meet the CFL condition! '
   !    Write(6, *) "(Debug: BZ) Pdx = :  ", Pdx, "Pdy = :  ", Pdy, "length = ", box_length
   !        !CALL GC_Error( errMsg, RC, thisLoc )
   !    CALL ERROR_STOP( errMsg, thisLoc)
   !  ENDIF
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
    Nt = CEILING(Dt/120)
    Pdt = Dt/REAL(Nt, fp) ! Dt=600 ! FLOOR(Pdx/Pu(1,1)/10)*10

     ! Find the best Pdt to meet CFL condition:
!700     CONTINUE
      !CFL = Pdt*Pu(1,1)/Pdx
      !IF(MAX( ABS(CFL), ABS(2*eddy_h*Pdt/(Pdx**2)), &
      !                  ABS(2*eddy_v*Pdt/(Pdy**2)) ) > 0.8)THEN

      !Nt = Nt+1
      !Pdt = Dt/Nt
      !GOTO 700

      !ENDIF
      DO 
         max_u = MAXVAL( ABS(Pu(2:n_x_max2-1,2:n_y_max2-1)) )
         CFL_adv = max_u * Pdt / Pdx
         CFL_dif_h = eddy_h * Pdt / Pdx**2
         CFL_dif_v = eddy_v * Pdt / Pdy**2
         !IF (MAX(CFL_adv, 2.0_fp*(CFL_dif_h+CFL_dif_v))<=0.8_fp) EXIT
         IF (CFL_adv + 2.0_fp*(CFL_dif_h+CFL_dif_v)<=0.8_fp) EXIT
         
         Nt = Nt+1
         IF (Nt .gt. 300) THEN
            Nt = 300
            !errMsg = '1-D plume diffusion: cannot satisfy CFL condition by substepping, exceed maximum 500'
            WRITE(6,*) "Debug (BZ): 2-D plume box cannot satisfy CFL condition, dissolve the box! "
            Write(6, *) "(Debug: BZ) Pdx = :  ", Pdx, "Pdy = :  ", Pdy, "length = ", box_length, &
            "Pdt= ", Pdt, "eddy_h= ", eddy_h, "eddy_v= ", eddy_v, "max_u=", max_u, &
            "CFL_adv= ", CFL_adv, "CFL_dif_h= ", CFL_dif_h, "CFL_dif_v= ", CFL_dif_v 
            Plume2d_curr%IsDissolve = .True.
            !CALL GC_Warning( ErrMsg, RC, ThisLoc )
            !CALL ERROR_STOP(errMsg, thisLoc)
            EXIT
            
         ENDIF
         Pdt = Dt/Nt
      ENDDO
      ! WRITE(6,*) 'Debug (BZ): (2-D plume) Numbers of substep for advection & diffusion, Nt= ', Nt, &
      ! 'time interval Pdt(s) = ', Pdt, "eddy_h= ", eddy_h, "eddy_v= ", eddy_v, "max_u=", max_u, &
      ! "CFL_adv= ", CFL_adv, "CFL_dif_h= ", CFL_dif_h, "CFL_dif_v= ", CFL_dif_v, &
      ! 'box = ', box_label, "Pdx = :  ", Pdx, "Pdy = :  ", Pdy, "length = ", box_length

      ! file_2Dconc_NK01_ID_1 = findFreeLun()
      ! WRITE(file_2Dconc_NK01_1,  &
      !       '("Plume-2D_NK01_mass_ID_",I0,"_time_", I0,"_beforephys.txt")') &
      !       Plume2d_curr%LABEL, NINT(time_elapsed)
      ! CALL PLUME_CONC_DIAG_FILES_2D(   &
      !       file_2Dconc_NK01_ID_1, file_2Dconc_NK01_1,   &
      !       Plume2d_curr%CONCNT2d(:,:, id_NK01_p)*Vgrid_2D, RC)
    !-------------------------------------------------------------------
    ! Calculate the advection-diffusion in 2D grids
    !   - Consider flux-limited / positive conservative scheme to enforce 
    !     donor-cell mass constraints and prevent negative concentrations.
    !     require multiple additional loops, substantially higher computational cost
    !
    !   - Cuurently update SF and NK tracer based on bulk sulfate dilution
    !   - This assume advection/diffusion changes total sulfate but does not directly
    !     modify the sulfate size distribution; size redistribution is handled later
    !     by TOMAS microphysics. (method 1)

    !   - This implementation treats TOMAS tracer separately, but cap the values
    !     at each substep to avoid 
    !     physically impossible negative values (method 2)
    !-------------------------------------------------------------------
      IF (Plume2d_curr%IsDissolve .OR. Plume2d_curr%IsTransfer) THEN
         GOTO 1110
      ENDIF
      ! Method 2
      DO i_species= 1, nspc_p, 1
         !Skip physics of these species
         IF ((i_species == id_NH3_p) .or.  (i_species == id_NH4_p) .or. &
                (i_species ==id_H2O_p).or.(i_species ==id_PH2SO4_p) .or. &
                 (i_species ==id_AW01_p).or. (i_species ==id_HO2_p) .or. &
                 (i_species ==id_OH_p)) CYCLE
         !(i_species == id_OH_p) .or. (i_species == id_HO2_p) .or.  &
         
         IF (i_species < 11) THEN
            spc_name = TRIM(spc_names_p(i_species))
            ind_spc_GC = Ind_(TRIM(spc_name))
            background_conc = Spc(ind_spc_GC)%Conc(i_lon, i_lat, i_lev)
         ELSEIF (i_species .lt. nspc_p) THEN
            ! TOMAS tracer
            i_tracer = MOD((i_species - 10), nspc_p_tomas_tracer)
            IF (i_tracer == 0) CYCLE ! Skip TOMAS AW tracer ! i_tracer = nspc_p_tomas_tracer
            spc_name = TRIM(spc_names_p(10+i_tracer))
            ind_spc_GC_bin1 = Ind_(TRIM(spc_name))
            ibin = (i_species - 11) / nspc_p_tomas_tracer + 1
            ind_spc_GC =  ind_spc_GC_bin1 +ibin - 1
            background_conc = Spc(ind_spc_GC)%Conc(i_lon, i_lat, i_lev)
         ENDIF
         SpcInfo => State_Chm%SpcData(ind_spc_GC)%Info
         C2d_prev(1:n_x_max,1:n_y_max) = &
                                 box_concnt_2D(1:n_x_max,1:n_y_max,i_species)
         C2d_bg(:,:)  = background_conc
         C2d_bg(2:n_x_max2-1,2:n_y_max2-1) =C2d_prev(1:n_x_max,1:n_y_max)
         C2d_new (:,:) = C2d_bg (:,:)
         DO t1s = 1, Nt 
            C2d_bg(:,:)  = background_conc
            C2d_bg(2:n_x_max2-1,2:n_y_max2-1) =C2d_new(2:n_x_max2-1,2:n_y_max2-1) 
            ! Only calculate the vertical half 2D domain
            !$OMP PARALLEL DO           &
            !$OMP DEFAULT( SHARED     ) &
            !$OMP PRIVATE(i_y,i_x,CFL,Pc_middle,Pc_top,Pc_bottom,Pc_right,Pc_left, Pc_update)
            DO i_y = 2, n_y_mid2, 1
            DO i_x = 2, n_x_max2-1, 1
              Pc_middle = C2d_bg( i_x,   i_y  )
              Pc_top    = C2d_bg( i_x,   i_y+1)
              Pc_bottom = C2d_bg( i_x,   i_y-1)
              Pc_right  = C2d_bg( i_x+1, i_y  )
              Pc_left   = C2d_bg( i_x-1, i_y  )
             
              CFL       = Pdt*Pu(i_x,i_y)/Pdx
              IF (Pu(i_x,i_y) >= 0.0_fp) THEN
                  adv_term = -CFL * (Pc_middle - Pc_left)
              ELSE
                  adv_term = -CFL * (Pc_right - Pc_middle)
              ENDIF
              Pc_update = Pc_middle + adv_term          &
                        + Pdt*( eddy_h*( Pc_right -2*Pc_middle +Pc_left   ) /(Pdx**2) &
                        + eddy_v*( Pc_top   -2*Pc_middle +Pc_bottom ) /(Pdy**2) )
            !   Pc_update = Pc_middle           &
            !     - 0.5 * CFL    * ( Pc_right - Pc_left )   &
            !     + 0.5 * CFL**2 * ( Pc_right - 2*Pc_middle + Pc_left )         &
            !     + Pdt*( eddy_h*( Pc_right -2*Pc_middle +Pc_left   ) /(Pdx**2) &
            !            +eddy_v*( Pc_top   -2*Pc_middle +Pc_bottom ) /(Pdy**2) )
              !C2d_new(i_x, i_y) = MAX(Pc_update, 0.0_fp)
              ! (BZ) If did not cap, might still have negative values due to rounding issue?
              C2d_new(i_x, i_y) = Pc_update
              ! update the other half based on vertical symmetry 
              C2d_new(n_x_max2+1-i_x, n_y_max2+1-i_y) = C2d_new(i_x, i_y)
            ENDDO
            ENDDO
            !$OMP END PARALLEL DO
         ENDDO ! DO t1s = 1, NINT(Dt/Pdt)
         ! IF (i_species == id_NK01_p)THEN
         !    file_2Dconc_NK01_ID_2 = findFreeLun()
         !    WRITE(file_2Dconc_NK01_2,  &
         !          '("Plume-2D_NK01_mass_ID_",I0,"_time_", I0,"_afterphys.txt")') &
         !          Plume2d_curr%LABEL, NINT(time_elapsed)
         !    CALL PLUME_CONC_DIAG_FILES_2D(   &
         !          file_2Dconc_NK01_ID_2, file_2Dconc_NK01_2,   &
         !          C2d_new(2:n_x_max2-1,2:n_y_max2-1)*Vgrid_2D, RC)
         ! ENDIF
         !================================================================
         ! Calculate the mass exchange of plume to background cell
         ! Update the concentration in the background and plume
         ! accordingly
         !================================================================
         ! the boundary always represents the background concentration
         mass_plume      = Vgrid_2D * SUM(C2d_prev(:,:))! molec
         mass_plume_new  = Vgrid_2D * SUM(C2d_new(2:n_x_max2-1,2:n_y_max2-1))

         D_mass_plume    = mass_plume_new - mass_plume
         background_mass     = background_conc * Vgrid_EU
         excess_mass = D_mass_plume - background_mass

         IF ((background_mass <= 0.0_fp).AND. (D_mass_plume.GT. 0.0_fp)) THEN
            WRITE (6, *) "(Debug: BZ) Mass enter the plume but background is 0, Revise to the prevous conc"
            WRITE (6, *) "(2-D) Plume num: ", Plume2d_curr%label, '; Species: ', TRIM(SpcInfo%Name), &
                        'D_mass_plume: ', D_mass_plume, 'background_mass: ', background_mass
            C2d_new(2:n_x_max2-1,2:n_y_max2-1) = C2d_prev
            mass_plume_new = mass_plume
            D_mass_plume = 0.0_fp
         ELSEIF (D_mass_plume .GT. background_mass) THEN
            WRITE (6, *) "(Debug: BZ) Mass enter the plume larger than mass in background, scaled mass entered"
            WRITE (6, *)  "(2-D) Plume num: ", Plume2d_curr%label, '; Species: ', TRIM(SpcInfo%Name), &
                           'D_mass_plume: ', D_mass_plume, 'background_mass: ', background_mass
            
            C2d_new(2:n_x_max2-1,2:n_y_max2-1) = C2d_new(2:n_x_max2-1,2:n_y_max2-1) * &
                           (mass_plume + background_mass) / mass_plume_new
            mass_plume_new = Vgrid_2D * SUM( C2d_new(2:n_x_max2-1,2:n_y_max2-1))
            D_mass_plume   = mass_plume_new - mass_plume
         ENDIF
         background_mass_new = MAX( 0.0_fp, background_mass - D_mass_plume )
         box_concnt_2D(:,:,i_species) = C2d_new(2:n_x_max2-1,2:n_y_max2-1)
         Spc(ind_spc_GC)%Conc(i_lon,i_lat,i_lev) = background_mass_new / Vgrid_EU

         IF (i_species .eq. id_SO2_p) THEN
            mass_S_SO2_r1_2D = mass_S_SO2_r1_2D + D_mass_plume
            ! WRITE (6, *) "(Debug: BZ) (2-D) SO2 Mass enter plume num: ", &
            !           Plume2d_curr%label, 'is D_mass_plume= ', D_mass_plume
         ENDIF
         IF (i_species .eq. id_SO4_p) THEN
            ! WRITE (6, *) "(Debug: BZ) (2-D) SO4 Mass enter plume num: ", &
            !           Plume2d_curr%label, 'is D_mass_plume= ', D_mass_plume          
            mass_S_SO4_r1_2D = mass_S_SO4_r1_2D + D_mass_plume
         ENDIF
      ENDDO
      
      Plume2d_curr%CONCNT2d    = box_concnt_2D
      
! #ifdef TOMAS
!    Write (6, *) 'Debug: BZ:  (2-D test grid) after plume physics, SO4 (molec)= ', box_concnt_2D(x_test,y_test, id_SO4_p)
!    Write (6, *) 'Debug: BZ: (2-D test grid)after plume physics, Nk bin 15 (molec)= ', box_concnt_2D(x_test,y_test, 53)
!    Write (6, *) 'Debug: BZ: (2-D test grid)after plume physics, Nk bin 15 (molec)= ', Plume2d_curr%CONCNT2d(x_test,y_test, 53)
!    Write (6, *) 'Debug: BZ: (2-D test grid) after plume physics, SF bin 15 (molec)= ', box_concnt_2D(x_test,y_test, 54)
! #endif
      ! mass_S_SO2_3_2D = mass_S_SO2_3_2D + SUM(box_concnt_2D(:,:,id_SO2_p)) * Vgrid_2D
      ! mass_S_SO4_3_2D = mass_S_SO4_3_2D + SUM(box_concnt_2D(:,:,id_SO4_p)) * Vgrid_2D
      ! Write (6, *) 'Debug: BZ:  (2-D phys after - 1), SO2 mass = ', SUM(box_concnt_2D(:,:, id_SO2_p)) * Vgrid_2D, &
      !          'SO4 mass = ',  SUM(box_concnt_2D(:,:, id_SO4_p)) * Vgrid_2D
      ! Write (6, *) 'Debug: BZ:  (2-D phys after - 2), SO2 mass = ', SUM(Plume2d_curr%CONCNT2d(:,:, id_SO2_p)) * Vgrid_2D, &
      !          'SO4 mass = ',  SUM(Plume2d_curr%CONCNT2d(:,:, id_SO4_p)) * Vgrid_2D
      ! file_2Dconc_SO2_ID_2 = findFreeLun()
      ! WRITE(file_2Dconc_SO2_2,'("Plume-2D_SO2_conc_",I0,"_2.txt")') NINT(time_elapsed)
      ! CALL PLUME_CONC_DIAG_FILES_2D(file_2Dconc_SO2_ID_2, file_2Dconc_SO2_2, Plume2d_curr%CONCNT2d(:,:, id_SO4_p), RC)

      ! file_2Dconc_SO4_ID_2 = findFreeLun()
      ! WRITE(file_2Dconc_SO4_2,'("Plume-2D_SO4_conc_",I0,"_2.txt")') NINT(time_elapsed)
      ! CALL PLUME_CONC_DIAG_FILES_2D(file_2Dconc_SO4_ID_2, file_2Dconc_SO4_2, Plume2d_curr%CONCNT2d(:,:, id_SO4_p), RC)
      


      ! --------------------------------------------------------------
      ! Decide if this plume will be dissolve or transfer:
      ! --------------------------------------------------------------
      ! Fourth Judge:
      ! If the plume touch the tropopause, dissolve the plume
      IF (box_lev>TROPP(i_lon,i_lat).AND.(TROPP_sink)) THEN
         Plume2d_curr%IsDissolve = .True.
      ENDIF
      ! IF (Plume2d_curr%LIFE>Critical_day_2D * 3600.0_fp * 24.0_fp) THEN
      !    Plume2d_curr%IsDissolve = .True.
      ! ENDIF
      ! Change from 2D to 1D, 
      ! once the tilting degree is bigger than 88 deg (88/180*3.14)
      IF(Plume2d_curr%LIFE>8.0*3600.0)THEN
         Xscale = Get_XYscale(Plume2d_curr%CONCNT2d(:,:,id_SO2_p), Pdx, Pdy, frac_mass, 2)
         Yscale = Get_XYscale(Plume2d_curr%CONCNT2d(:,:,id_SO2_p), Pdx, Pdy, frac_mass, 1)
         box_theta = ATAN( Xscale/Yscale )
         Write(6, * ) "Debug (BZ): box_theta for SO2: box_label= ", box_label, &
                        "Xscale= ", Xscale, " Yscale= ", Yscale,  "box_theta= ",box_theta, &
                        "box_alpha= ",box_alpha, &
                        'Pdx = ', Plume2d_curr%Pdx, 'Pdy = ', Plume2d_curr%Pdy, 'Length = ',  Plume2d_curr%length
         ! box_theta = ATAN( Xscale/Yscale )
         IF ((Xscale/Yscale .gt. 25.0_fp) .OR.(Plume2d_curr%LIFE > Critical_day_2D * 3600.0_fp * 24.0_fp)) THEN
            Plume2d_curr%IsTransfer = .True.
         ENDIF 
      ENDIF
1110  CONTINUE
      Plume2d_curr => Plume2d_curr%next
   ENDDO  ! DO WHILE(ASSOCIATED(Plume2d))

      ! file_2Dconc_OH_ID_2 = findFreeLun()
      ! WRITE(file_2Dconc_OH_2,'("Plume-2D_OH_conc_",I0,"_2.txt")') NINT(time_elapsed)
      ! OPEN(file_2Dconc_OH_ID_2, FILE=TRIM(file_2Dconc_OH_2), STATUS='REPLACE', &
      !    FORM='FORMATTED', ACCESS='SEQUENTIAL', IOSTAT=RC)
      ! DO i_x = 1, n_x_max
      !    WRITE(file_2Dconc_OH_ID_2,'(*(ES12.4,1X))') &
      !       (Plume2d_curr%CONCNT2d(i_x,i_y,id_OH_p), i_y = 1, n_y_max)
      ! ENDDO
      ! CLOSE(file_2Dconc_OH_ID_2)
        !write(6,*) 'debug (BZ): SO2mass after inplume conc: ', SUM(box_concnt_2D(:,:,id_SO2_p)) *  Vgrid_2D
        !write(6,*) 'debug (BZ): SO4mass after inplume conc: ', SUM(box_concnt_2D(:,:,id_SO4_p)) *  Vgrid_2D
        !write(6,*) 'debug (BZ): SO2conc after inplume conc: ', SUM(box_concnt_2D(:,:,id_SO2_p)) /  (n_x_max*n_y_max)
        !write(6,*) 'debug (BZ): SO4conc after inplume conc: ', SUM(box_concnt_2D(:,:,id_SO4_p)) /  (n_x_max*n_y_max)
!#ifdef TOMAS
   !Write (6, *) 'Debug: BZ: after inplume conc, SO4 (molec)= ', box_concnt_2D(1,2, id_SO4_p)
   !Write (6, *) 'Debug: BZ: after inplume conc, Nk bin 14 (molec)= ', box_concnt_2D(1,2, 10+1+(14-1)*nspc_p_tomas_tracer)
   ! Write (6, *) 'Debug: BZ: after inplume conc, Nk bin sum (molec)= ', SUM(box_concnt_2D(:,:, 10+1+(14-1)*nspc_p_tomas_tracer))
   !Write (6, *) 'Debug: BZ: after inplume conc Mk(SO4) bin 1 (molec)= ', box_concnt_2D(1,2, 12)
!#endif
    
    !write(6,*) 'debug (BZ): SO2mass after physics: ', SUM(box_concnt_2D(:,:,id_SO2)) * V_grid_2D
    !write(6,*) 'debug (BZ): SO4mass after physics: ', SUM(box_concnt_2D(:,:,id_SO4)) * V_grid_2D

401 CONTINUE
  !=======================================================================
  ! For 1d plume: Run Lagrangian trajectory-track HERE
  !=======================================================================
  IF(.NOT.ASSOCIATED(Plume1d_head)) GOTO 400
  Plume1d_curr => Plume1d_head
  DO WHILE(ASSOCIATED(Plume1d_curr))
    box_lon       = Plume1d_curr%LON
    box_lat       = Plume1d_curr%LAT
    box_lev       = Plume1d_curr%LEV
    i_lat         = Plume1d_curr%lat_ind
    i_lon         = Plume1d_curr%lon_ind
    i_lev         = Plume1d_curr%lev_ind

    box_length    = Plume1d_curr%LENGTH
    box_alpha     = Plume1d_curr%ALPHA
    box_label     = Plume1d_curr%label
    box_life      = Plume1d_curr%LIFE
    box_Ra        = Plume1d_curr%RA
    box_Rb        = Plume1d_curr%RB
    box_theta     = Plume1d_curr%THETA

    box_concnt_1D = Plume1d_curr%CONCNT1d

    Vgrid_1D      = Plume1d_curr%RA * Plume1d_curr%RB * Plume1d_curr%LENGTH *1.0e+6_fp ! [cm3]
    ! Vgrid_1D_tot_1= Vgrid_1D_tot_1+Vgrid_1D * n_slab_max
    !write(6,*) 'debug (BZ): solve plume physics in 1-D plume box: ', box_label
   !  mass_S_SO2_2_1D = mass_S_SO2_2_1D + SUM(box_concnt_1D(:,id_SO2_p)) * Vgrid_1D
   !  mass_S_SO4_2_1D = mass_S_SO4_2_1D + SUM(box_concnt_1D(:,id_SO4_p)) * Vgrid_1D
   !  Write (6, *) 'Debug: BZ:  (1-D phys before - 1), SO2 mass = ', SUM(box_concnt_1D(:,id_SO2_p)) * Vgrid_1D, &
   !             'SO4 mass = ',  SUM(box_concnt_1D(:, id_SO4_p)) * Vgrid_1D
   !  Write (6, *) 'Debug: BZ:  (1-D phys before - 2), SO2 mass = ', SUM(Plume1d_curr%CONCNT1d(:, id_SO2_p)) * Vgrid_1D, &
   !             'SO4 mass = ',  SUM(Plume1d_curr%CONCNT1d(:, id_SO4_p)) * Vgrid_1D
    box_life = box_life + Dt

    curr_lon      = box_lon
    curr_lat      = box_lat
    curr_pressure = box_lev
    
    DO Ki = 1,4,1
      !------------------------------------------------------------------
      ! For vertical wind speed:
      ! pay attention for the polar region * * *
      !------------------------------------------------------------------
      if(abs(curr_lat)>Y_mid(NY_GC))then
         curr_omeg = Interplt_wind_RLL_polar(omeg, i_lon, i_lat, i_lev, &
                                       curr_lon, curr_lat, curr_pressure)
      else
         curr_omeg = Interplt_wind_RLL(omeg, i_lon, i_lat, i_lev, &
                                       curr_lon, curr_lat, curr_pressure)
      endif

      RK_omeg(Ki) = curr_omeg
      RK_Dlev(Ki)   = Dt * curr_omeg / 100.0     ! Pa => hPa
      curr_pressure = box_lev + RK_Dt(Ki+1) * curr_omeg / 100.0

      if(curr_pressure<P_mid(NZ_GC)) &
            curr_pressure = P_mid(NZ_GC) !+ ( P_mid(LLPAR) - curr_pressure )
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

         if(abs(curr_lat)>Y_mid(NY_GC))then
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
      ! make sure the location is not out of range
      do while (box_lat > Y_edge(NY_GC+1))
         box_lat = Y_edge(NY_GC+1) - ( box_lat-Y_edge(NY_GC+1) )
      end do
      do while (box_lat < Y_edge(1))
         box_lat = Y_edge(1) - ( box_lat-Y_edge(1) )
      end do
      do while (box_lon > X_edge(NX_GC+1))
         box_lon = box_lon - 360.0
      end do
      do while (box_lon < X_edge(1))
         box_lon = box_lon + 360.0
      end do
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
      if(abs(curr_lat)>Y_mid(NY_GC))then
         curr_T1 = Interplt_wind_RLL_polar(T1, i_lon, i_lat, &
                                 i_lev, curr_lon, curr_lat, curr_pressure)
      else
         curr_T1 = Interplt_wind_RLL(T1, i_lon, i_lat, i_lev, &
                                        curr_lon, curr_lat, curr_pressure)
      endif

      ! update index based on new location
      i_lon = Find_iLonLat(box_lon, DX_GC, X_edge2)
      if(i_lon>NX_GC) i_lon=i_lon-NX_GC
      if(i_lon<1) i_lon=i_lon+NX_GC

      i_lat = Find_iLonLat(box_lat, DY_GC, Y_edge2)
      if(i_lat>NY_GC) i_lat=NY_GC
      if(i_lat<1) i_lat=1

      i_lev = Find_iPLev(box_lev,P_edge)
      if(i_lev>NZ_GC) i_lev=NZ_GC

      if(abs(box_lat)>Y_mid(NY_GC))then
         next_T2 = Interplt_wind_RLL_polar(T2, i_lon, i_lat, i_lev, &
                           box_lon, box_lat, box_lev)
      else
         next_T2 = Interplt_wind_RLL(T2, i_lon, i_lat, i_lev, &
                           box_lon, box_lat, box_lev)
      endif

      ! PV=nRT, V2 = T2/P2 : T1/P1 * V1
      ratio = SQRT( (next_T2/box_lev)/(curr_T1/curr_pressure) )

      ! assume the volume change mainly apply to the cross-section, 
      ! the box_length would not change 

      !V_prev = box_Ra*box_Rb*box_length*1.0e+6_fp
      !V_new = V_prev * ratio * ratio
      box_Ra = box_Ra *ratio
      box_Rb = box_Rb *ratio

      box_concnt_1D(:,:) = box_concnt_1D(:,:)/(ratio **2.0_fp)

      !------------------------------------------------------------------
      ! calcualte the box_alpha [0,2*PI)
      !------------------------------------------------------------------
      IF((box_u**2.0_fp+box_v**2.0_fp)==0)then
          box_alpha = 0.0
      ELSEIF(box_v>=0.0_fp)THEN
          box_alpha = ACOS( box_u/SQRT(box_u**2.0_fp+box_v**2.0_fp) )
      ELSE
          box_alpha = 2.0_fp*PI - ACOS( box_u/SQRT(box_u**2.0_fp+box_v**2.0_fp) )
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
      box_Ra = box_Ra*SQRT(length0/box_length)
      box_Rb = box_Rb*SQRT(length0/box_length)
      ! ----------------------------------------------------------------
      ! update     
      ! ----------------------------------------------------------------
      Plume1d_curr%LON    = box_lon
      Plume1d_curr%LAT    = box_lat
      Plume1d_curr%LEV    = box_lev
      Plume1d_curr%lon_ind = i_lon
      Plume1d_curr%lat_ind = i_lat
      Plume1d_curr%lev_ind = i_lev

      Plume1d_curr%LENGTH = box_length
      Plume1d_curr%ALPHA  = box_alpha
      Plume1d_curr%label  = box_label
      Plume1d_curr%LIFE   = box_life
      Plume1d_curr%RA     = box_Ra
      Plume1d_curr%RB     = box_Rb

      Plume1d_curr%CONCNT1d = box_concnt_1D
      Vplume_1D_tot(i_lon, i_lat, i_lev) = &
               Vplume_1D_tot(i_lon, i_lat, i_lev)  +  &
               (box_Ra*box_Rb*box_length*n_slab_max) * 1.0e+6_fp ! [cm3]
      !-------- solve in-plume concentration here --------
      ! advection and diffusion
      curr_lon         = box_lon
      curr_lat         = box_lat
      curr_pressure    = box_lev      ! hPa

      Vgrid_EU         = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ! [cm3]
      Vgrid_1D         = Plume1d_curr%RA * Plume1d_curr%RB * Plume1d_curr%LENGTH *1.0e+6_fp ! [cm3]
      ! Vgrid_1D_tot_2   = Vgrid_1D_tot_2 + Vgrid_1D * n_slab_max
      !====================================================================
      ! calculate the wind shear along plume corss-section
      ! calculate the diffusivity in horizontal and vertical direction
      !====================================================================
      ! calculate the wind_s shear along pressure direction
      wind_s_shear = Wind_shear_s(u, v, P_BXHEIGHT, box_alpha, i_lon, &
                         i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
      
      ! Calculate vertical eddy diffusivity (U.Schumann, 2012) :
      Cv = 0.2_fp
      Omega_N = 0.1_fp
      Ptemp_shear = Vertical_shear(Ptemp, P_BXHEIGHT, i_lon, i_lat, &
                           i_lev, curr_lon, curr_lat, curr_pressure)
      !--------------------------------------------------------------------
      ! interpolate potential temperature for Plume module:
      !--------------------------------------------------------------------
      if(abs(curr_lat)>Y_mid(NY_GC))then
         curr_Ptemp = Interplt_wind_RLL_polar(Ptemp, i_lon, i_lat, &
                                 i_lev, curr_lon, curr_lat, curr_pressure)
      else
         curr_Ptemp = Interplt_wind_RLL(Ptemp, i_lon, i_lat, i_lev, &
                                       curr_lon, curr_lat, curr_pressure)
      endif

      
      IF(N_BV<=0.001) N_BV = 0.001
      N_BV = SQRT(Ptemp_shear*g0/curr_Ptemp)
      ! diffusivity unit: [m2/s]
      eddy_v = Cv * Omega_N**2 / N_BV
      eddy_h = 10.0
      !Write(6, *) "(Debug: BZ) 1-D: wind_s_shear = ", wind_s_shear, "; Ptemp_shear=", Ptemp_shear, &
      !   "eddy_v=", eddy_v
      ! ============================================
      ! Decide the time step based on CFL condition
      ! 2*k*Dt/Dr<1 (or Dr-2*k*Dt>0) for diffusion
      !=============================================
      ! ignore the diffusivity in long radius direction
      
      !Nt = CEILING(Dt / 120.0_fp)
      !Nt = MAX(1, Nt)
      Nt = 15
      
      DO 
         Dt2 = Dt / REAL(Nt, fp)
         box_theta_test = Plume1d_curr%THETA
         box_Ra_test = Plume1d_curr%RA
         box_Rb_test = Plume1d_curr%RB
         CFL_ok = .TRUE.
         DO t1s = 1, Nt
            theta_previous   = box_theta_test
            box_theta_candidate = ATAN( TAN(theta_previous) + wind_s_shear*Dt2 )
            box_Ra_candidate = box_Ra_test  * SQRT(TAN(box_theta_candidate)**2+1.0_fp) &
                              / SQRT(TAN(theta_previous)**2+1.0_fp)
            box_Rb_candidate = Vgrid_1D/(box_Ra_candidate*box_length*1.0e+6_fp)
            !IF (box_Rb_candidate >= Rb_min_shear) THEN
               box_theta_test = box_theta_candidate
               box_Ra_test= box_Ra_candidate
               box_Rb_test = box_Rb_candidate
            !ENDIF

            eddy_D = eddy_v*SIN(box_theta_test)**2 + &
                     eddy_h*COS(box_theta_test)**2 
            CFL_dif_D = eddy_D * Dt2 / box_Rb_test**2
            !WRITE(6,*) "Debug (BZ): t1s = ", t1s, "; eddy_D = ", eddy_D, "; Dt2 = ", Dt2, &
            !            "; CFL_dif_D = ", CFL_dif_D, "; Nt = ", Nt, &
            !            "; box_Rb = ", box_Rb_test, "; box_Ra = ", box_Ra_test, "; box_theta = ", box_theta_test
                   
            IF (CFL_dif_D > 0.4_fp) THEN
               CFL_ok = .FALSE.
               EXIT
            ENDIF
         ENDDO

         IF (CFL_ok) EXIT
         
         Nt = Nt + 5

         IF (Nt .gt. 300) THEN
            Nt = 300
            !errMsg = '1-D plume diffusion: cannot satisfy CFL condition by substepping, exceed maximum 500'
            WRITE(6,*) "Debug (BZ): cannot satisfy CFL condition, dissolve the box! "
            WRITE(6,*) "Debug (BZ): box_Rb = ", Plume1d_curr%Rb, "; box_Ra = ", Plume1d_curr%Ra, "; box_theta = ", Plume1d_curr%THETA
            WRITE(6,*) "Debug (BZ): box_Rb_test = ", box_Rb_test, "; box_Ra_test = ", box_Ra_test, "; box_theta_test = ", box_theta_test
            Plume1d_curr%IsDissolve = .True.
            !CALL GC_Warning( ErrMsg, RC, ThisLoc )
            !CALL ERROR_STOP(errMsg, thisLoc)
            EXIT
            
         ENDIF
      ENDDO
      IF (Plume1d_curr%IsDissolve) THEN
         GOTO 1111
      ENDIF
      Dt2 = Dt / REAL(Nt, fp)
      box_Ra_final = box_Ra_test
      box_Rb_final = box_Rb_test
      box_theta_final = box_theta_test
      !WRITE(6,*) 'Debug (BZ): (1-D plume) Numbers of substep for advection & diffusion, Nt= ', Nt, &
      !'time interval Dt2(s) = ', Dt2, 'box = ', box_label

      box_concnt_1D_prev = box_concnt_1D
      box_RA_tmp = Plume1d_curr%RA     
      box_Rb_tmp = Plume1d_curr%RB
      box_theta_tmp = Plume1d_curr%THETA
      DO i_species = 1, nspc_p
         !Skip physics of these species
         IF ((i_species == id_OH_p) .or. (i_species == id_HO2_p) .or.  &
            (i_species == id_NH3_p) .or.  (i_species == id_NH4_p) .or. &
            (i_species ==id_H2O_p).or.(i_species ==id_PH2SO4_p) .or. &
            (i_species ==id_AW01_p)) CYCLE

         !Vgrid_1D = box_Ra(i_species) * box_Rb(i_species) *box_length*1.0e+6_fp
         box_RA = box_RA_tmp
         box_RB = box_Rb_tmp
         box_theta = box_theta_tmp

         IF (i_species < 11) THEN
            spc_name = TRIM(spc_names_p(i_species))
            ind_spc_GC = Ind_(TRIM(spc_name))
            background_conc = Spc(ind_spc_GC)%Conc(i_lon, i_lat, i_lev)
         ELSEIF (i_species .lt. nspc_p) THEN
            ! TOMAS tracer
            i_tracer = MOD((i_species - 10), nspc_p_tomas_tracer)
            IF (i_tracer == 0)  CYCLE ! Skip TOMAS AW tracer ! i_tracer = nspc_p_tomas_tracer
            spc_name = TRIM(spc_names_p(10+i_tracer))
            ind_spc_GC_bin1 = Ind_(TRIM(spc_name))
            ibin = (i_species - 11) / nspc_p_tomas_tracer + 1
            ind_spc_GC =  ind_spc_GC_bin1 +ibin - 1
            background_conc = Spc(ind_spc_GC)%Conc(i_lon, i_lat, i_lev)
         ENDIF
         SpcInfo => State_Chm%SpcData(ind_spc_GC)%Info

         C1d_prev(1:n_slab_max) = box_concnt_1D(1:n_slab_max,i_species)
         C1d_bg(:)  = MAX(0.0_fp, background_conc)
         C1d_bg(2:n_slab_max2-1) =C1d_prev(1:n_slab_max)
         C1d_new(:) = C1d_bg(:)

         DO t1s = 1, Nt

            C1d_bg(:)  = MAX(0.0_fp, background_conc)
            C1d_bg(2:n_slab_max2-1) =C1d_new(2:n_slab_max2-1)

            theta_previous   = box_theta
            box_theta_candidate = ATAN( TAN(theta_previous) + wind_s_shear*Dt2 )
            box_Ra_candidate = box_Ra  * SQRT(TAN(box_theta_candidate)**2+1.0_fp) &
                              / sqrt(TAN(theta_previous)**2+1.0_fp)
            box_Rb_candidate = Vgrid_1D/(box_Ra_candidate*box_length*1.0e+6_fp)
            ! Apply shear distortion only if Rb remains above threshold
            !IF (box_Rb_candidate >= Rb_min_shear) THEN
               box_theta = box_theta_candidate
               box_Ra    = box_Ra_candidate
               box_Rb    = box_Rb_candidate
            !ENDIF
            eddy_D = eddy_v*SIN(box_theta)**2 + eddy_h*COS(box_theta)**2 ! b
            CFL_dif_D = eddy_D * Dt2 / box_Rb**2
            ! WRITE(6,*) "Debug (BZ): t1s = ", t1s, "; eddy_D = ", eddy_D, "; Dt2 = ", Dt2, "; box_Rb = ", box_Rb, &
            !          "; CFL_dif_D = ", CFL_dif_D, "species = ", TRIM(spc_names_p_use(i_species))

            C1d_new(2:n_slab_max2-1) = C1d_bg(2:n_slab_max2-1)     &
               + CFL_dif_D*(   C1d_bg(3:n_slab_max2)   &
                        -2.0_fp*C1d_bg(2:n_slab_max2-1) &
                        +  C1d_bg(1:n_slab_max2-2) ) 
            !C1d_new(2:n_slab_max2-1) = MAX(C1d_new(2:n_slab_max2-1), 0.0_fp)
         ENDDO

         !Vgrid_1D_new    =  box_Ra(i_species) * box_Rb(i_species) *box_length*1.0e+6_fp
         ! Vgrid_1D_new    =  box_Ra * box_Rb *box_length*1.0e+6_fp

         mass_plume      = Vgrid_1D * SUM(C1d_prev(:))! molec
         mass_plume_new  = Vgrid_1D * SUM(C1d_new(2:n_slab_max2-1))

         D_mass_plume    = mass_plume_new - mass_plume
         background_mass     = background_conc * Vgrid_EU
         excess_mass = D_mass_plume - background_mass

         IF ((background_mass <= 0.0_fp).AND. (D_mass_plume.GT. 0.0_fp)) THEN
            WRITE (6, *) "(Debug: BZ) (1-D): Mass enter the plume but background is 0, Revise to the prevous conc"
            WRITE (6, *) "Plume num: ", Plume1d_curr%label, '; Species: ', TRIM(SpcInfo%Name), &
                        '; D_mass_plume: ', D_mass_plume, '; background_mass: ', background_mass
            C1d_new(2:n_slab_max2-1) = C1d_prev
            D_mass_plume    = 0.0_fp
            mass_plume_new = mass_plume
         ELSEIF (D_mass_plume .GT. background_mass) THEN
            WRITE (6, *) "(Debug: BZ) (1-D) Mass enter the plume larger than mass in background, scaled mass entered"
            WRITE (6, *)  "Plume num: ", Plume1d_curr%label, '; Species: ', TRIM(SpcInfo%Name), &
                           '; D_mass_plume: ', D_mass_plume, '; background_mass: ', background_mass
            
            C1d_new(2:n_slab_max2-1) = C1d_new(2:n_slab_max2-1) * &
                           (mass_plume + background_mass) / mass_plume_new
            mass_plume_new = Vgrid_1D * SUM(C1d_new(2:n_slab_max2-1))
            D_mass_plume   = mass_plume_new - mass_plume
         ENDIF

         background_mass_new = MAX( 0.0_fp, background_mass - D_mass_plume )
         Spc(ind_spc_GC)%Conc(i_lon,i_lat,i_lev) = background_mass_new / Vgrid_EU
         box_concnt_1D(:,i_species) = C1d_new(2:n_slab_max2-1)

         IF (i_species .eq. id_SO2_p) THEN
            mass_S_SO2_r1_1D = mass_S_SO2_r1_1D + D_mass_plume
            !WRITE (6, *) "(Debug: BZ) (1-D) SO2 Mass enter plume num: ", &
            !         Plume1d_curr%label, 'is D_mass_plume= ', D_mass_plume
            !mass_S_SO2_2_1D = mass_S_SO2_2_1D + SUM(box_concnt_1D(:,i_species)) * Vgrid_1D
         ENDIF
         IF (i_species .eq. id_SO4_p) THEN
            !WRITE (6, *) "(Debug: BZ) (1-D) SO4 Mass enter plume num: ", &
            !         Plume1d_curr%label, 'is D_mass_plume= ', D_mass_plume        
            mass_S_SO4_r1_1D = mass_S_SO4_r1_1D + D_mass_plume
            !mass_S_SO4_2_1D = mass_S_SO4_2_1D + SUM(box_concnt_1D(:,i_species)) * Vgrid_1D
         ENDIF

      ENDDO
      Plume1d_curr%CONCNT1D = box_concnt_1D
      Plume1d_curr%RA = box_Ra_final
      Plume1d_curr%RB = box_Rb_final
      Plume1d_curr%THETA = box_theta_final
     
      ! --------------------------------------------------------------
      ! Decide if this plume will be dissolve
      ! Location + lifetime
      ! --------------------------------------------------------------
      IF ((box_lev>TROPP(i_lon,i_lat).AND.(TROPP_sink)).OR.(box_life>Critical_day_1D * 3600.0_fp*24.0_fp )) THEN
         Plume1d_curr%IsDissolve = .True.
      ENDIF
      ! --------------------------------------------------------------
      ! Decide if we need to split the plume
      ! --------------------------------------------------------------

! 1111 mass_S_SO2_3_1D = mass_S_SO2_3_1D + SUM(Plume1d_curr%CONCNT1D(:,id_SO2_p)) * Vgrid_1D
!      mass_S_SO4_3_1D = mass_S_SO4_3_1D + SUM(Plume1d_curr%CONCNT1D(:,id_SO4_p)) * Vgrid_1D
!      Write (6, *) 'Debug: BZ:  (1-D phys aftere - 1), SO2 mass = ', SUM(box_concnt_1D(:,id_SO2_p)) * Vgrid_1D, &
!                'SO4 mass = ',  SUM(box_concnt_1D(:, id_SO4_p)) * Vgrid_1D
!      Write (6, *) 'Debug: BZ:  (1-D phys after - 2), SO2 mass = ', SUM(Plume1d_curr%CONCNT1d(:, id_SO2_p)) * Vgrid_1D, &
!                'SO4 mass = ',  SUM(Plume1d_curr%CONCNT1d(:, id_SO4_p)) * Vgrid_1D
1111  CONTINUE
      Plume1d_curr => Plume1d_curr%next
  ENDDO

400 CONTINUE
  !------------------------------------------------------------------
  ! Everything is done, clean up pointers
  !------------------------------------------------------------------
  ! deallocate unused space
  IF(allocated(box_concnt_2D)) deallocate(box_concnt_2D)
  IF(allocated(box_concnt_1D)) deallocate(box_concnt_1D)
  IF(allocated(box_concnt_1D_prev)) deallocate(box_concnt_1D_prev)
  ! Nullify pointers
  IF(ASSOCIATED(u)) nullify(u)
  IF(ASSOCIATED(v)) nullify(v)
  IF(ASSOCIATED(omeg)) nullify(omeg)
  IF(ASSOCIATED(Ptemp)) nullify(Ptemp)
  IF(ASSOCIATED(T1)) nullify(T1)
  IF(ASSOCIATED(T2)) nullify(T2)
  IF(ASSOCIATED(P_BXHEIGHT)) nullify(P_BXHEIGHT)
  IF(ASSOCIATED(TROPP)) nullify(TROPP)
  IF(ASSOCIATED(Spc)) nullify(Spc)
  !IF(ASSOCIATED(X_edge)) nullify(X_edge)
  !IF(ASSOCIATED(Y_edge)) nullify(Y_edge)

  !IF(ASSOCIATED(PASV_EU)) nullify(PASV_EU)

  IF(ASSOCIATED(Plume2d_new)) nullify(Plume2d_new)
  IF(ASSOCIATED(Plume2d_curr)) nullify(Plume2d_curr)
  IF(ASSOCIATED(Plume2d_prev)) nullify(Plume2d_prev)
  !IF(ASSOCIATED(Plume1d_new)) nullify(Plume1d_new)
  IF(ASSOCIATED(Plume1d_curr)) nullify(Plume1d_curr)
  !IF(ASSOCIATED(Plume1d_prev)) nullify(Plume1d_prev)

  END SUBROUTINE plume_physics

  SUBROUTINE plume_chem_microphysics(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
    ! All the use have been defined in host: Plume_box_model
    USE Input_Opt_Mod,            ONLY : OptInput, PlumeSource_t
    USE State_Chm_Mod,            ONLY : ChmState, Ind_
    USE State_Met_Mod,            ONLY : MetState
    USE Species_Mod,              ONLY : SpcConc, Species
    USE TIME_MOD,                 ONLY : GET_TS_DYN, GET_TS_CHEM
    USE State_Grid_Mod,           ONLY : GrdState
    USE UnitConv_Mod
    USE ERROR_MOD
    USE InquireMod,      ONLY : findFreeLun
    !USE State_Diag_Mod,           ONLY : DgnState
    !USE State_Diag_Mod,           ONLY : DgnMap
    ! KPP-related module
!#ifdef KPP_INTEGRATOR_AUTOREDUCE
!    USE fullchem_AutoReduceFuncs, ONLY : fullchem_AR_KeepHalogensActive
!    USE fullchem_AutoReduceFuncs, ONLY : fullchem_AR_SetKeepActive
!    USE fullchem_AutoReduceFuncs, ONLY : fullchem_AR_UpdateKppDiags
!    USE fullchem_AutoReduceFuncs, ONLY : fullchem_AR_SetIntegratorOptions
!#endif
!    USE GcKpp_Global
    ! USE GcKpp_Parameters,           ONLY : NREACT, NSPEC
!    USE Gckpp_Monitor,            ONLY : SPC_NAMES, Eqn_Names, Fam_Names
!    USE GcKpp_Rates,              ONLY : UPDATE_RCONST, RCONST
!    USE GcKpp_Integrator,         ONLY : Integrate
!    USE GcKpp_Function
    ! GEOS-Chem related module
!    USE UCX_MOD,                  ONLY : SO4_PHOTFRAC

    LOGICAL, INTENT(IN)           :: am_I_Root
    TYPE(MetState), INTENT(IN)    :: State_Met
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State objectgg
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    !TYPE(DgnState), INTENT(INOUT) :: State_Diag ! Diagnostics State object
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    
    INTEGER                       :: N
    INTEGER                       :: ind_spc_p
    INTEGER                       :: i_box, i_lon, i_lat, i_lev
    INTEGER                       :: i_species, i_species_1,i_phot, i_kpp, i_rxn
    INTEGER                       :: i_x, i_y, i_slab
    !INTEGER                       :: n_species 
    INTEGER                       :: Thread, IERR,  P, F, errorCount
    INTEGER                       :: SpcID, KppID
    INTEGER                       :: RXN_O3_1, RXN_O3_2
    !INTEGER                       :: ind_SO2, ind_SO4, ind_OH
    
    INTEGER                       :: ISTATUS(20)
    INTEGER                       :: ICNTRL (20)
    INTEGER                       :: chem_status
    INTEGER                       :: n_low_SO2, n_low_OH, n_low_rate
    INTEGER                       :: n_debug_chem, n_debug_chem_max
    INTEGER, ALLOCATABLE          :: debug_ix(:), debug_iy(:), debug_ibox(:), debug_islab(:)
    INTEGER, ALLOCATABLE          :: debug_status(:)
    
    REAL(fp)                      :: Vgrid_2D, Vgrid_1D, Vgrid_EU, Vgrid_1D_SO4, Vgrid_1D_SO2, Vgrid_1D_PH2SO4
    REAL(fp)                      :: mass_OH, mass_HO2
    REAL(fp)                      :: Dt
    REAL(fp)                      :: SO4_FRAC,   SR,        LWC
    !REAL(dp)                      :: RCNTRL (20)
    !REAL(dp)                      :: RSTATE (20)
    !REAL(dp)                      :: C_before_integrate(NSPEC)
    !REAL(dp)                      :: local_RCONST(NREACT)
    REAL(fp)                      :: H2SO4_RATE_2d(n_x_max,n_y_max) ! H2SO4 prod rate [kg s-1]
    REAL(fp)                      :: H2SO4_RATE_1d(n_slab_max) ! H2SO4 prod rate [kg s-1]
    REAL(fp)                      :: PSO4AQ_RATE_2d(n_x_max,n_y_max) ! Cld chem sulfate prod rate [kg s-1]
    REAL(fp)                      :: K_SO2_OH
    REAL(fp)                      :: full_exchange_conc
    REAL(fp)                      :: chem_debug_value
    REAL(fp)                      :: C_before_Chem(nspc_p),C_after_Chem(nspc_p)
    REAL(fp)                      :: Core_Mean ! Mean concentration of plume core in 2-D concentration array
    REAL(fp)                      :: Conc_background_Mean ! mean conc in Eulerian background 
    REAL(fp)                      :: Vplume_frac
    
    REAL(fp), ALLOCATABLE         :: box_concnt_2D(:,:,:), box_concnt_2D_prev(:,:,:)
    REAL(fp), ALLOCATABLE         :: box_concnt_1D(:,:), box_concnt_1D_prev(:,:)
    REAL(fp), ALLOCATABLE         :: debug_value(:)
    REAL(fp), ALLOCATABLE         :: mass_OH_consum_plume(:,:,:), mass_HO2_consum_plume(:,:,:) ! molec, record the total amount of OH/HO2 consume during plume chemistry
    REAL(fp), ALLOCATABLE         :: mass_NH3_consum_plume(:,:,:), mass_NH4_consum_plume(:,:,:) ! molec, record the total amount of NH3/NH4 consume during plume microphysics



    LOGICAL                       :: Failed2x,  Size_Res, doSuppress

    CHARACTER(LEN=255)            :: ErrMsg
    CHARACTER(LEN=255)            :: ThisLoc

    TYPE(SpcConc), POINTER        :: Spc(:)
    TYPE(Species), POINTER        :: SpcInfo

    TYPE(Plume2d_list), POINTER :: Plume2d_next, Plume2d_curr, Plume2d_prev
    TYPE(Plume1d_list), POINTER :: Plume1d_next, Plume1d_curr, Plume1d_prev
    
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
    ! TOMAS in plume notes: Bingqing Zhang
    ! - Temporary define TOMAS related variable here so that it can be moved
    ! to a separate module if necessary
    ! - For TOMAS, species need to be in unit: kg, convert before feeding to TOMAS array
    ! - TOMAS species order (in Mk, Gc) is different from that used in Plume and GC
    ! hard coded below: 1 SO4, 2, NH4, 3, H2O
    ! - Temporarily set unused variable = -1
    !========================================================================
    
    !REAL(fp)             :: Mo
    

    LOGICAL              :: COND, COAG, NUCL !<step5.1> switch for each process (win 4/8/06)
    LOGICAL              :: PRINTNEG  !<step4.0-temp> (win, 3/24/05)
    LOGICAL              :: ERRORSWITCH  !<step4.2> To see where mnfix found negative value (win, 9/12/05)
    LOGICAL              :: ERRSPOT   !<step4.4> To see where so4cond found errors (win, 9/21/05)
    LOGICAL              :: PRINTDEBUG !<step4.3> Print out for debugging (win, 9/16/05)

    INTEGER              :: ibin, i_L
    INTEGER              :: num_iter
    INTEGER              :: id_tracer, id_tracer_1
    
    REAL(fp)             :: BOXVOL,  BOXMASS, TEMPTMS,   PRES,   RHTOMAS
    REAL(fp)             :: surf_area     ! aerosol surface area [micon^2 cm^-3]
    REAL(fp)             :: h2so4rate_o ! H2SO4rate for the specific grid cell
    REAL(fp)             :: fn  ! nucleation rate of clusters cm-3 s-1
    REAL(fp)             :: fn1 ! formation rate of particles to first size bin cm-3 s-1
    REAL(fp)             :: nucrate(State_Grid%NY,State_Grid%NZ)
    REAL(fp)             :: nucrate1(State_Grid%NY,State_Grid%NZ)
    REAL(fp)             :: molwt_spc
    REAL(fp)             :: NH4bulk, NH3_to_NH4, CEPS
    REAL(fp)             :: ionrate  ! ion pair formation rate [ion pairs cm^-3 s^-1]
    REAL(fp)             :: tot_s_1
    REAL(fp)             :: tot_n_1
    REAL(fp)             :: TOT_MK, TOT_NK
    REAL(fp)             :: Nk(nBins), Nkd(nBins), Nkout(nBins), Nknuc(nBins), Nkcond(nBins)
    REAL(fp)             :: Mk(nBins, ICOMPHARD), Mkd(nBins, ICOMPHARD), Mkout(nBins,ICOMPHARD), Mknuc(nBins,ICOMPHARD), Mkcond(nBins,ICOMPHARD)
    REAL(fp)             :: Gc(ICOMPHARD), Gcd(ICOMPHARD), Gcout(ICOMPHARD)
    !REAL(fp)             :: 
    REAL(fp)             :: TRANSFER(nbins)
    REAL(fp)             :: mass_NH3, mass_NH4
    parameter ( CEPS=1.e-17_fp )
    ! Arguments for CHECK_VALUE; avoids array temporaries (bmy, 1/28/14)
    CHARACTER(LEN=255) :: ERR_VAR
    CHARACTER(LEN=255) :: ERR_MSG
    INTEGER            :: ERR_IND(4)

    CHARACTER(LEN=255)     :: spc_name
    ! This variable is intend to use as a reference, when exchange species from TOMAS to plume
    ! convert to original unit mole/cm3 and exchange with box_concnt_2D
    ! REAL(fp), ALLOCATABLE         :: box_concnt_2D_kg(:,:,:) 
    ! Initialization
    ! Initialize switches for each microphysical process
    COND                   = .TRUE.
    COAG                   = .TRUE.
    NUCL                   = .TRUE.

    ! Initialize debugging and error-signal switches
    PRINTNEG               = .FALSE.
    ERRORSWITCH            = .FALSE.
    PRINTDEBUG             = .FALSE.
    ERRSPOT                = .FALSE.
    !========================================================================
    ! plume_chem_microphysics begins here!
    ! New version solve simplified chemistry independently
    ! Old version info below: 
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
    NULLIFY(Plume1d_next, Plume1d_curr, Plume1d_prev)

    Spc                    =>   State_Chm%Species
    SpcInfo                =>   NULL()
    ThisLoc                =    ' -> at plume_chem_microphysics (in module GeosCore/lagrange_singlebox_mod.F90)'
    RC                     =    GC_SUCCESS
    ErrMsg                 =    ''
    i_box                  =    0
    ! n_species              =    State_Chm%nSpecies
    Thread                 =    1
    errorCount             =    0
    Failed2x               =   .FALSE.
    doSuppress             =   .FALSE.
    !id_SO2                 =    Ind_('SO2')
    !id_SO4                 =    Ind_('SO4')
    !id_OH                  =    Ind_('OH')
    H2SO4_RATE_2d          =    0.0_fp
    H2SO4_RATE_1d          =    0.0_fp
    PSO4AQ_RATE_2d         =    0.0_fp
    ! RXN_O3_1 specifies: O3 + hv -> O2 + O
    ! RXN_O3_2 specifies: O3 + hv -> O2 + O(1D)
    ! (BZ): For debug purpose
    RXN_O3_1              = State_Chm%Phot%RXN_O3_1
    RXN_O3_2              = State_Chm%Phot%RXN_O3_2
    ! Noted that in GEOS-Chem, the default chemistry timestep is 20min, dynamic time step is 10min
    ! Here we implement chemsitry timestep using 10min
    Dt                     =    GET_TS_CHEM()
    n_debug_chem_max       =    n_x_max * n_y_max
    

   !  mass_S_SO2_4_2D    = 0.0_fp
   !  mass_S_SO4_4_2D    = 0.0_fp
   !  mass_S_SO2_5_2D    = 0.0_fp
   !  mass_S_SO4_5_2D    = 0.0_fp
   !  mass_S_SO2_4_1D    = 0.0_fp
   !  mass_S_SO4_4_1D    = 0.0_fp
   !  mass_S_SO2_5_1D    = 0.0_fp
   !  mass_S_SO4_5_1D    = 0.0_fp

!     mass_S_SO2_r2_2D   = 0.0_fp
!     mass_S_SO4_r2_2D   = 0.0_fp

!     mass_S_SO2_r2_1D   = 0.0_fp
!     mass_S_SO4_r2_1D   = 0.0_fp
! #ifdef TOMAS
!     mass_S_H2SO4_2D       = 0.0_fp
!     mass_S_H2SO4_1D       = 0.0_fp
! #endif
    ALLOCATE(box_concnt_2D(n_x_max, n_y_max, nspc_p ))
    ALLOCATE(box_concnt_2D_prev(n_x_max, n_y_max, nspc_p ))
    ALLOCATE(box_concnt_1D(n_slab_max, nspc_p ))
    ALLOCATE(box_concnt_1D_prev(n_slab_max, nspc_p ))
    ALLOCATE(debug_ix(n_debug_chem_max))
    ALLOCATE(debug_iy(n_debug_chem_max))
    ALLOCATE(debug_islab(n_debug_chem_max))
    ALLOCATE(debug_ibox(n_debug_chem_max))
    ALLOCATE(debug_status(n_debug_chem_max))
    ALLOCATE(debug_value(n_debug_chem_max))

    ALLOCATE(mass_OH_consum_plume(NX_GC, NY_GC, NZ_GC ))
    ALLOCATE(mass_HO2_consum_plume(NX_GC, NY_GC, NZ_GC ))
    mass_OH_consum_plume   =    0.0_fp
    mass_HO2_consum_plume  =    0.0_fp

#ifdef TOMAS
    !ALLOCATE(box_concnt_2D_kg(n_x_max, n_y_max, nspc_p ))
    ALLOCATE(mass_NH3_consum_plume(NX_GC, NY_GC, NZ_GC ))
    ALLOCATE(mass_NH4_consum_plume(NX_GC, NY_GC, NZ_GC ))
    mass_NH3_consum_plume   =    0.0_fp
    mass_NH4_consum_plume  =    0.0_fp
#endif
    ! Set up integration convergence conditions and timesteps
    ! This is defined in gckpp_Global and set to be public
    ! ATOL = State_Chm%KPP_AbsTol   ! Absolute tolerance
    ! RTOL = State_Chm%KPP_RelTol   ! Relative tolerance

    ! IF (State_Diag%Archive_RxnConst        ) Write(6, *) "Debug: (BZ) ; Archive_RxnRate", State_Diag%Archive_RxnRate
    ! For debug process, print rate constant 202: SO2 + OH {+M} = SO4 + HO2 + PH2SO4 :
    !Write (6, *) "Debug: (BZ): In Plume  (Before Plume Chem): rate constant for RXN SO2_OH_RXN_ID = ", &
    !  RXNRATE_CONST_KPP(23, 40, 39 ,SO2_OH_RXN_ID)
   !  write(6,*) 'debug (BZ): Before Plume Chem, species ID --- '
   !  write(6,*) 'id_OH = ', id_OH, 'id_OH_p = ', id_OH_p, &
   !    'id_HO2 = ', id_HO2, 'id_HO2_p = ', id_HO2_p, &
   !    'id_SO2 = ', id_SO2, 'id_SO2_p = ', id_SO2_p, &
   !    'id_SO4 = ', id_SO4, 'id_SO4_p = ', id_SO4_p

    IF(.NOT.ASSOCIATED(Plume2d_head)) GOTO 401
    Plume2d_curr => Plume2d_head
    DO WHILE(ASSOCIATED(Plume2d_curr))
      
      i_box         = Plume2d_curr%label
      i_lon         = Plume2d_curr%lon_ind
      i_lat         = Plume2d_curr%lat_ind
      i_lev         = Plume2d_curr%lev_ind
      !write(6,*) 'debug (BZ): solve plume Chemistry in plume 2-D box: ', i_box

      Vgrid_2D        =  Plume2d_curr%Pdx * Plume2d_curr%Pdy * Plume2d_curr%length * 1.0e+6_fp ! [cm3]
      Vgrid_EU        =  State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp
      K_SO2_OH        =  RXNRATE_CONST_KPP(i_lon,i_lat,i_lev , SO2_OH_RXN_ID)
      box_concnt_2D   =  Plume2d_curr%CONCNT2d
      
      !Test if we need to do the chemistry for box (I,J,L), otherwise move onto the next box.
      ! MaxChemLev = MaxStratLev = 59 
      ! Hard coded in GeosUtil/gc_grid_mod.F90
      ! BZ: Maybe if plume reach out of chem grid, directly release species to Eulerian grid and delete the plume box?
      IF ( .not. State_Met%InChemGrid(i_lon,i_lat,i_lev) ) THEN
        WRITE(6,*) 'Debug (BZ): 2-D plume outside chem grid: (i_box, i_lon, i_lat, i_lev): ',    &
        Plume2d_curr%label, i_lon, i_lat, i_lev
        Plume2d_curr%IsDissolve = .True.
      ENDIF
      IF (Plume2d_curr%IsDissolve .OR. Plume2d_curr%IsTransfer) THEN
         GOTO 1112
      ENDIF


      ! Exchange OH/HO2 with background before chemistry 
      ! write(6,*) 'debug (BZ): (2-D) Before exchange, background OH conc =  ', Spc(id_OH)%Conc(i_lon,i_lat,i_lev), &
      ! '; background HO2 conc = ', Spc(id_HO2)%Conc(i_lon,i_lat,i_lev), 'Plume number: ', Plume2d_curr%label, &
      ! 'Plume location (X, Y, L) = ', i_lon, i_lat, i_lev
      ! mass_OH     = Spc(id_OH)%Conc(i_lon,i_lat,i_lev)*Vgrid_EU +      &
      !                     SUM(box_concnt_2D(:,:,id_OH_p))*Vgrid_2D
      ! mass_HO2    = Spc(id_HO2)%Conc(i_lon,i_lat,i_lev)*Vgrid_EU +     &
      !                     SUM(box_concnt_2D(:,:,id_HO2_p))*Vgrid_2D
      !full_exchange_conc = mass_OH  / (Vgrid_EU + Vgrid_2D*n_x_max*n_y_max)

      mass_OH     = Spc(id_OH)%Conc(i_lon,i_lat,i_lev)*Vgrid_EU
      mass_HO2    = Spc(id_HO2)%Conc(i_lon,i_lat,i_lev)*Vgrid_EU
      
      ! full_exchange_conc = mass_OH  / (Vgrid_EU + Vplume_2D_tot(i_lon, i_lat, i_lev) + Vplume_1D_tot(i_lon, i_lat, i_lev))
      full_exchange_conc =Spc(id_OH)%Conc(i_lon,i_lat,i_lev)
      ! Vplume_frac = (Vplume_2D_tot(i_lon, i_lat, i_lev) + Vplume_1D_tot(i_lon, i_lat, i_lev)) /  &
      !       (Vgrid_EU + Vplume_2D_tot(i_lon, i_lat, i_lev) + Vplume_1D_tot(i_lon, i_lat, i_lev))
      Vplume_frac =  (Vplume_2D_tot(i_lon, i_lat, i_lev) + Vplume_1D_tot(i_lon, i_lat, i_lev)) /  Vgrid_EU
      Write(6, *) "debug (BZ): plume volume fraction = ", Vplume_frac, " at location (X, Y, L) = ", i_lon, i_lat, i_lev

      box_concnt_2D(:,:,id_OH_p)   = full_exchange_conc * Radical_exc_factor
      mass_OH_consum_plume(i_lon, i_lat, i_lev)  =  mass_OH_consum_plume(i_lon, i_lat, i_lev) + SUM(box_concnt_2D(:,:,id_OH_p))*Vgrid_2D
      ! Spc(id_OH)%Conc(i_lon,i_lat,i_lev) = (mass_OH -SUM(box_concnt_2D(:,:,id_OH_p))*Vgrid_2D ) / &
      !                                         Vgrid_EU
      
      ! full_exchange_conc = mass_HO2  / (Vgrid_EU + Vplume_2D_tot(i_lon, i_lat, i_lev) + Vplume_1D_tot(i_lon, i_lat, i_lev))
      full_exchange_conc = Spc(id_HO2)%Conc(i_lon,i_lat,i_lev)
      box_concnt_2D(:,:,id_HO2_p)   = full_exchange_conc * Radical_exc_factor
      ! Spc(id_HO2)%Conc(i_lon,i_lat,i_lev) = (mass_HO2 -SUM(box_concnt_2D(:,:,id_HO2_p))*Vgrid_2D ) / &
      !                                         Vgrid_EU
      mass_HO2_consum_plume(i_lon, i_lat, i_lev)  =  mass_HO2_consum_plume(i_lon, i_lat, i_lev) + SUM(box_concnt_2D(:,:,id_HO2_p))*Vgrid_2D
      ! write(6,*) 'debug (BZ): (2-D) After exchange, before chem, background OH conc =  ', Spc(id_OH)%Conc(i_lon,i_lat,i_lev), &
      ! '; background HO2 conc = ', Spc(id_HO2)%Conc(i_lon,i_lat,i_lev), 'Plume number: ', Plume2d_curr%label, &
      ! 'Plume location (X, Y, L) = ', i_lon, i_lat, i_lev

      ! mass_S_SO2_4_2D = mass_S_SO2_4_2D + SUM(Plume2d_curr%CONCNT2D(:,:,id_SO2_p)) * Vgrid_2D
      ! mass_S_SO4_4_2D = mass_S_SO4_4_2D + SUM(Plume2d_curr%CONCNT2D(:,:,id_SO4_p)) * Vgrid_2D
      
!#ifdef TOMAS
!      Write (6, *) 'Debug: BZ: before chemistry, Nk bin 14 (molec)= ', box_concnt_2D(1,2, 10+1+(14-1)*nspc_p_tomas_tracer)
!      Write (6, *) 'Debug: BZ: before chemistry Mk(SO4) bin 1 (molec)= ', box_concnt_2D(1,2, 12)
!#endif
      !WRITE(6,*) 'Debug (BZ): Euleria grid: Conc of SO2 ', Spc(id_SO2)%Conc(i_lon,i_lat,i_lev)
      !WRITE(6,*) 'Debug (BZ): Euleria grid: Conc of SO4 ', Spc(id_SO4)%Conc(i_lon,i_lat,i_lev)
      !WRITE(6,*) 'Debug (BZ): Euleria grid: Conc of OH ', Spc(id_OH)%Conc(i_lon,i_lat,i_lev)

      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc before Chem: Conc of SO2 ', SUM(Plume2d_curr%CONCNT2d(:,:,id_SO2_p))/(n_x_max*n_y_max)
      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc before Chem: Conc of SO4 ', SUM(Plume2d_curr%CONCNT2d(:,:,id_SO4_p))/(n_x_max*n_y_max)
      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc before Chem: Conc of OH ', SUM(Plume2d_curr%CONCNT2d(:,:,id_OH_p))/(n_x_max*n_y_max)

      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc before Chem: Conc of SO2 ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_SO2_p)
      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc before Chem: Conc of SO4 ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_SO4_p)
      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc before Chem: Conc of OH ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_OH_p)
      ! Write (6, *) 'Debug: BZ:  (2-D chem before), SO2 mass = ', SUM(box_concnt_2D(:,:, id_SO2_p)) * Vgrid_2D, &
      !                'SO4 mass = ',  SUM(box_concnt_2D(:,:, id_SO4_p)) * Vgrid_2D
      ! Write (6, *) 'Debug: BZ:  (2-D chem before - 2), SO2 mass = ', SUM(Plume2d_curr%CONCNT2d(:,:, id_SO2_p)) * Vgrid_2D, &
      !                'SO4 mass = ',  SUM(Plume2d_curr%CONCNT2d(:,:, id_SO4_p)) * Vgrid_2D 
      n_debug_chem     =  0
      debug_ix         =  0
      debug_iy         =  0
      debug_ibox       =  0
      debug_status     =  0
      debug_value      =  0.0_fp
      box_concnt_2D_prev = box_concnt_2D
      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE(i_y,i_x)&
      !$OMP PRIVATE(chem_status, chem_debug_value)&
      !$OMP PRIVATE(C_before_Chem, C_after_Chem)&
      !$OMP COLLAPSE( 2                    )
      !!$OMP SCHEDULE( DYNAMIC, 24          )&
      !!$OMP REDUCTION( +:n_debug_chem        )

      DO i_y = 1, n_y_max, 1
        DO i_x = 1, n_x_max, 1
          
          chem_status     =  0
          chem_debug_value = 0.0_fp

          C_before_Chem(id_SO2_p)     =  box_concnt_2D(i_x,i_y,id_SO2_p)
          C_before_Chem(id_SO4_p)     =  box_concnt_2D(i_x,i_y,id_SO4_p)
          C_before_Chem(id_OH_p)      =  box_concnt_2D(i_x,i_y,id_OH_p)
          C_before_Chem(id_HO2_p)     =  box_concnt_2D(i_x,i_y,id_HO2_p)
          C_before_Chem(id_PH2SO4_p)  =  0.0_fp

          C_after_Chem              =  C_before_Chem
          ! CHEM_SO2_OH_PLUME (dt, K, SO2, OH, SO4, HO2, PH2SO4, i_x, i_y, i_box) 
          CALL CHEM_SO2_OH_PLUME( Dt, K_SO2_OH, C_after_Chem(id_SO2_p),       &
                                     C_after_Chem(id_OH_p), C_after_Chem(id_SO4_p), &
                                     C_after_Chem(id_HO2_p), C_after_Chem(id_PH2SO4_p), &
                                     chem_status,chem_debug_value )

          IF ( chem_status /= 0 ) THEN
            !$OMP CRITICAL(debug_store)
            IF ( n_debug_chem < n_debug_chem_max ) THEN
                n_debug_chem = n_debug_chem + 1
                debug_ix(n_debug_chem)     = i_x
                debug_iy(n_debug_chem)     = i_y
                debug_ibox(n_debug_chem)   = i_box
                debug_status(n_debug_chem) = chem_status
                debug_value(n_debug_chem)  = chem_debug_value
            ENDIF
            !$OMP END CRITICAL(debug_store)
          ENDIF

          box_concnt_2D(i_x,i_y,id_SO2_p)    = REAL( C_after_Chem(id_SO2_p), kind=fp )
          box_concnt_2D(i_x,i_y,id_SO4_p)    = REAL( C_after_Chem(id_SO4_p), kind=fp )
          box_concnt_2D(i_x,i_y,id_HO2_p)    = REAL( C_after_Chem(id_HO2_p), kind=fp )
          box_concnt_2D(i_x,i_y,id_OH_p)     = REAL( C_after_Chem(id_OH_p), kind=fp )
          box_concnt_2D(i_x,i_y,id_PH2SO4_p) = REAL( C_after_Chem(id_PH2SO4_p), kind=fp )
          H2SO4_RATE_2D(i_x,i_y) = C_after_Chem(id_PH2SO4_p) / AVO * 98.e-3_fp * &
                           (Vgrid_2D) / Dt  ! kg s-1 box-1

        ENDDO ! DO i_x = 1, n_x_max, 1
      ENDDO ! i_y = 1, n_y_max, 1
      !$OMP END PARALLEL DO
      
      ! Chemistry debug output
      DO N = 1, n_debug_chem
        SELECT CASE ( debug_status(N) )
        CASE (1)
          WRITE(6,*) "Skip chem in 2-D Plume: low SO2 at 2-D segment [x,y,box]= ", &
                    debug_ix(N), debug_iy(N), debug_ibox(N), &
               " value = ", debug_value(n)
        CASE (2)
          WRITE(6,*) "Skip chem in 2-D Plume: low OH at 2-D segment [x,y,box]= ", &
                    debug_ix(N), debug_iy(N), debug_ibox(N), &
               " value = ", debug_value(n)
        CASE (3)
          WRITE(6,*) "Skip chem in 2-D Plume: low K[SO2][OH]dt at 2-D segment [x,y,box] = ", &
                    debug_ix(N), debug_iy(N), debug_ibox(N), &
                    " value = ", debug_value(N)
        END SELECT
      ENDDO
!#ifdef TOMAS
   !Write (6, *) 'Debug: BZ:  (2-D test grid) after plume chemistry, SO4 (molec)= ', box_concnt_2D(x_test,y_test, id_SO4_p)
   !Write (6, *) 'Debug: BZ: (2-D test grid) after plume chemistry, Nk bin 15 (molec)= ', box_concnt_2D(x_test,y_test, 53)
   !Write (6, *) 'Debug: BZ: (2-D test grid) after plume chemistry, SF bin 15 (molec)= ', box_concnt_2D(x_test,y_test, 54)
!#endif

      
      
      mass_S_SO2_r2_2D = mass_S_SO2_r2_2D + SUM(box_concnt_2D(:,:, id_SO2_p) - box_concnt_2D_prev(:,:, id_SO2_p)) * Vgrid_2D
      mass_S_SO4_r2_2D = mass_S_SO4_r2_2D + SUM(box_concnt_2D(:,:, id_SO4_p) - box_concnt_2D_prev(:,:, id_SO4_p)) * Vgrid_2D

      ! Subtract remaining radical
      mass_OH_consum_plume(i_lon, i_lat, i_lev)  = mass_OH_consum_plume(i_lon, i_lat, i_lev)   - SUM(box_concnt_2D(:,:,id_OH_p))*Vgrid_2D
      mass_HO2_consum_plume(i_lon, i_lat, i_lev) = mass_HO2_consum_plume(i_lon, i_lat, i_lev)  - SUM(box_concnt_2D(:,:,id_HO2_p))*Vgrid_2D
      ! ! Release remaining radical to the background
      ! mass_OH = SUM(box_concnt_2D(:,:, id_OH_p)) * Vgrid_2D
      ! Spc(id_OH)%Conc(i_lon,i_lat,i_lev)  =  Spc(id_OH)%Conc(i_lon,i_lat,i_lev) +  &
      !          mass_OH/Vgrid_EU
      ! box_concnt_2D(:,:, id_OH_p)  = 0.0_fp
      
      ! mass_HO2 = SUM(box_concnt_2D(:,:, id_HO2_p)) * Vgrid_2D
      ! Spc(id_OH)%Conc(i_lon,i_lat,i_lev)  =  Spc(id_HO2)%Conc(i_lon,i_lat,i_lev) +  &
      !          mass_HO2/Vgrid_EU
      ! box_concnt_2D(:,:, id_HO2_p)  = 0.0_fp
      ! write(6,*) 'debug (BZ): (2-D) After Chemistry, after exchange, background OH conc =  ', Spc(id_OH)%Conc(i_lon,i_lat,i_lev), &
      ! '; background HO2 conc = ', Spc(id_HO2)%Conc(i_lon,i_lat,i_lev), 'Plume number: ', Plume2d_curr%label, &
      ! 'Plume location (X, Y, L) = ', i_lon, i_lat, i_lev
      ! Write (6, *) 'Debug: BZ:  (2-D chem after), SO2 mass = ', SUM(box_concnt_2D(:,:, id_SO2_p)) * Vgrid_2D, &
      !                'SO4 mass = ',  SUM(box_concnt_2D(:,:, id_SO4_p)) * Vgrid_2D
      
      ! file_2Dconc_SO2_ID_3 = findFreeLun()
      ! WRITE(file_2Dconc_SO2_3,'("Plume-2D_SO2_conc_",I0,"_3.txt")') NINT(time_elapsed)
      ! CALL PLUME_CONC_DIAG_FILES_2D(file_2Dconc_SO2_ID_3, file_2Dconc_SO2_3, box_concnt_2D(:,:, id_SO2_p), RC)

      ! file_2Dconc_SO4_ID_3 = findFreeLun()
      ! WRITE(file_2Dconc_SO4_3,'("Plume-2D_SO4_conc_",I0,"_3.txt")') NINT(time_elapsed)
      ! CALL PLUME_CONC_DIAG_FILES_2D(file_2Dconc_SO4_ID_3, file_2Dconc_SO4_3, box_concnt_2D(:,:, id_SO4_p), RC)

      ! file_2Dconc_OH_ID_3 = findFreeLun()
      ! WRITE(file_2Dconc_OH_3,'("Plume-2D_OH_conc_",I0,"_3.txt")') NINT(time_elapsed)
      ! OPEN(file_2Dconc_OH_ID_3, FILE=TRIM(file_2Dconc_OH_3), STATUS='REPLACE', &
      !    FORM='FORMATTED', ACCESS='SEQUENTIAL', IOSTAT=RC)
      ! DO i_x = 1, n_x_max
      !    WRITE(file_2Dconc_OH_ID_3,'(*(ES12.4,1X))') &
      !       (box_concnt_2D(i_x,i_y,id_OH_p), i_y = 1, n_y_max)
      ! ENDDO
      ! CLOSE(file_2Dconc_OH_ID_3)
#ifdef TOMAS

      ! Before doing chemistry, read background NH3/NH4
      mass_NH3 = Spc(id_NH3)%Conc(i_lon,i_lat,i_lev)* Vgrid_EU 
      ! full_exchange_conc = mass_NH3  / (Vgrid_EU + Vplume_2D_tot(i_lon, i_lat, i_lev) + Vplume_1D_tot(i_lon, i_lat, i_lev))
      full_exchange_conc =Spc(id_NH3)%Conc(i_lon,i_lat,i_lev)
      box_concnt_2D(:,:,id_NH3_p) = full_exchange_conc
      mass_NH3_consum_plume(i_lon, i_lat, i_lev) =    &
               mass_NH3_consum_plume(i_lon, i_lat, i_lev) + SUM(box_concnt_2D(:,:,id_NH3_p)) *Vgrid_2D

      ! Spc(id_NH3)%Conc(i_lon,i_lat,i_lev) = mass_NH3 / &
      !          (Vgrid_2D*n_x_max*n_y_max + Vgrid_EU)

      mass_NH4 = Spc(id_NH4)%Conc(i_lon,i_lat,i_lev)* Vgrid_EU
      ! full_exchange_conc = mass_NH4  / (Vgrid_EU + Vplume_2D_tot(i_lon, i_lat, i_lev) + Vplume_1D_tot(i_lon, i_lat, i_lev))
      full_exchange_conc = Spc(id_NH4)%Conc(i_lon,i_lat,i_lev)
      box_concnt_2D(:,:,id_NH4_p) = full_exchange_conc
      mass_NH4_consum_plume(i_lon, i_lat, i_lev) =    &
               mass_NH4_consum_plume(i_lon, i_lat, i_lev) + SUM(box_concnt_2D(:,:,id_NH4_p)) *Vgrid_2D
      ! Spc(id_NH4)%Conc(i_lon,i_lat,i_lev) = mass_NH4 / &
      !       (Vgrid_2D*n_x_max*n_y_max + Vgrid_EU)
      !write(6,*) 'debug (BZ): solve plume microphysics in plume 2-D box: ', i_box
      PRES    = State_Met%PMID(i_lon,i_lat,i_lev)*100.0 ! in Pa
      TEMPTMS = State_Met%T(i_lon,i_lat,i_lev)
      BOXMASS = State_Met%AD(i_lon,i_lat,i_lev)*Vgrid_2D/Vgrid_EU ! Dry air mass, kg
      RHTOMAS = State_Met%RH(i_lon,i_lat,i_lev)/ 1.e2
      IF ( RHTOMAS > 0.99 ) RHTOMAS = 0.99
      BOXVOL  = Vgrid_2D  !cm3
      !Write (6, *) 'Debug (BZ): PRES(Pa)=', PRES, ' TEMPTMS(K)=', TEMPTMS, ' BOXMASS=', BOXMASS, ' RHTOMAS(0-1)=', RHTOMAS
      !Write (6, *) 'Debug: BZ: after chemistry, Nk bin 14 (molec)= ', box_concnt_2D(1,2, 10+1+(14-1)*nspc_p_tomas_tracer)
      !Write (6, *) 'Debug: BZ: after chemistry Mk(SO4) bin 1 (molec)= ', box_concnt_2D(1,2, 12)
      !write(6,*) 'debug (BZ): solve plume microphysics in plume 1-D box: ', i_box

      ! TOMAS species unit should be in kg
      ! GC, MK, NK
      !$OMP PARALLEL DO         &
      !$OMP DEFAULT( SHARED )   &
      !$OMP PRIVATE( i_y, i_x, i_species, i_species_1, ibin, id_tracer, id_tracer_1 )  &
      !$OMP PRIVATE(spc_name, molwt_spc, tot_n_1, tot_s_1, TOT_NK, TOT_MK)    &
      !$OMP PRIVATE(Gc, Gcd, GCout, Nk, Nkd, Nkout, Nknuc, Nkcond, MK, Mkd,Mkout, Mknuc, Mkcond)        &
      !$OMP PRIVATE( TRANSFER, H2SO4rate_o, NH4bulk, NH3_to_NH4, fn, fn1,num_iter, ionrate, surf_area)       & 
      !$OMP PRIVATE(ERRORSWITCH, PRINTNEG, PRINTDEBUG, ERRSPOT)  &
      !$OMP PRIVATE( ERR_VAR, ERR_MSG, ERR_IND )        &    
      !$OMP COLLAPSE( 2                    )  &         
      !$OMP SCHEDULE( DYNAMIC )
      DO i_y = 1, n_y_max, 1
         DO i_x = 1, n_x_max, 1
            
            ERRORSWITCH = .FALSE.
            ! Initialize all condensible gas values to zero
            ! Gc(srtso4) will remain zero until within cond_nuc where the
            ! pseudo steady state H2SO4 concentration will be put in this place.
            Gc(:) = 0.0_fp
            ! Swap Spc into Nk, Mk, Gc arrays, unit kg
            Nk(:) = 0.0_fp
            Mk(:,:) = 0.0_fp
            DO ibin = 1, nBins
               ! Mol weight for Nk is 1
               i_species_1 = id_NK01_p+(ibin-1)*nspc_p_tomas_tracer
               NK(ibin) = box_concnt_2D(i_x,i_y,i_species_1)*BOXVOL/Avo/1.0E+3_fp
               DO i_species = 1, ICOMPHARD-2 ! skip NH4 and H2O
                  i_species_1 = id_NK01_p+i_species+(ibin-1)*nspc_p_tomas_tracer
                  spc_name = spc_names_p(id_NK01_p+i_species)
                  id_tracer   = Ind_(TRIM(spc_name))
                  molwt_spc       = State_Chm%SpcData(id_tracer)%Info%MW_g
                  MK(ibin,i_species) = box_concnt_2D(i_x,i_y,i_species_1)*BOXVOL/Avo*molwt_spc/1.0E+3_fp
                  !IF( IT_IS_NAN( MK(ibin,i_species) ) ) THEN
                  !    PRINT *,'+++++++ Found NaN in AEROPHYS ++++++++'
                  !    PRINT *,'Location (i_plume, i_x, i_y):',i_box,i_x,i_y,'Bin',ibin,'comp',spc_name
                  !ENDIF
               ENDDO
               id_tracer   = Ind_('AW01')
               molwt_spc       = State_Chm%SpcData(id_tracer)%Info%MW_g
               MK(ibin,SRTH2O) = box_concnt_2D(i_x,i_y,id_AW01_p+(ibin-1)*nspc_p_tomas_tracer)*BOXVOL/Avo*molwt_spc/1.0E+3_fp
            ENDDO
            !IF ((i_x.eq.x_test) .AND. (i_y.eq.y_test)) THEN
            !   WRITE(6, *) "Debug (BZ): [TOMAS] NK 15 in = ", NK(15)
            !ENDIF
            ! Get NH4 mass from the bulk mass and scale to bin with sulfate
            IF ( SRTNH4 > 0 ) THEN

               NH4bulk = box_concnt_2D(i_x,i_y,id_NH4_p)*BOXVOL/Avo*18.0_fp/1.0E+3_fp
               CALL NH4BULKTOBIN( MK(:,SRTSO4), NH4bulk, TRANSFER )
               MK(1:nbins,SRTNH4) = TRANSFER(1:nbins)
               Gc(SRTNH4) = box_concnt_2D(i_x,i_y,id_NH3_p)*BOXVOL/Avo*17.0_fp/1.0E+3_fp

            ENDIF
            ! Give it the pseudo-steady state value instead later (win,9/30/08)
            !GC(SRTSO4) = Spc(id_H2SO4)%Conc(I,J,L)
          
            H2SO4rate_o = H2SO4_RATE_2D(i_x, i_y)  ! [kg s-1]
            IF ( H2SO4rate_o .lt. 0.e0 ) THEN
                Print*, 'Debug TOMAS: (2-D) H2SO4RATE = ', H2SO4rate_o, 'i_box = ', i_box, &
                    'i_x = ', i_x, 'i_y = ', i_y
                H2SO4rate_o = 0.e+0_fp
            ENDIF
            
            ! nitrogen and sulfur mass checks
            ! get the total mass of N
            tot_n_1 = Gc(srtnh4)*14.e+0_fp/17.e+0_fp
            do ibin=1,nbins
               tot_n_1 = tot_n_1 + Mk(ibin,srtnh4)*14.e+0_fp/18.e+0_fp
            enddo

            ! get the total mass of S
            tot_s_1 = H2SO4rate_o*Dt*32.e+0_fp/98.e+0_fp
            do ibin=1,nbins
               tot_s_1 = tot_s_1 + Mk(ibin,srtso4)*32.e+0_fp/96.e+0_fp
            enddo
            !!$OMP CRITICAL(EZWATEREQM_TEST)
            !Do water eqm at appropriate times
            CALL EZWATEREQM( MK, RHTOMAS )
            !!$OMP END CRITICAL(EZWATEREQM_TEST)
            !IF ((i_x .eq. 1) .AND. (i_y.eq.2)) THEN
            !      Write (6, *) 'Debug: BZ: before microphysics, Nk bin 14 (kg)= ', Nk(14)
            !      Write (6, *) 'Debug: BZ: before microphysics Mk(SO4) bin 1 (kg) = ', Mk(1,1)
            !ENDIF
            !Fix any inconsistencies in M/N distribution (because of advection)
            !IF ((i_x .eq. x_test) .AND. (i_y.eq.y_test)) THEN
            !   WRITE (6, *) "Debug (BZ), AEROPHYS-MNFIX (1), Nk is: "
            !   WRITE(*,'(15(1X,ES12.4))') (Nk(ibin), ibin=1,15)
            !   WRITE (6, *) "Debug (BZ), AEROPHYS-MNFIX (1), SF is: "
            !   WRITE(*,'(15(1X,ES12.4))') (Mk(ibin,SRTSO4), ibin=1,15)
            !ENDIF
            CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
            CALL MNFIX( NK, MK, ERRORSWITCH )
            IF ( ERRORSWITCH ) THEN
               PRINT *,'Aerophys 1: (2-D) MNFIX found error at (i_box, i_x, i_y) = ',i_box, i_x, i_y
               WRITE (6, *) "Debug (BZ), AEROPHYS-MNFIX (1), Nk is: "
               WRITE(*,'(15(1X,ES12.4))') (Nk(ibin), ibin=1,15)
               WRITE (6, *) "Debug (BZ), AEROPHYS-MNFIX (1), SF is: "
               WRITE(*,'(15(1X,ES12.4))') (Mk(ibin,SRTSO4), ibin=1,15)
               WRITE (6, *) "Debug (BZ), AEROPHYS-MNFIX (1), H2O is: "
               WRITE(*,'(15(1X,ES12.4))') (Mk(ibin,SRTH2O), ibin=1,15)
               CALL ERROR_STOP('AEROPHYS-MNFIX (1)','Enter microphys')
            ENDIF
            !IF ((i_x.eq.x_test) .AND. (i_y.eq.y_test)) THEN
            !   WRITE(6, *) "Debug (BZ): [TOMAS] NK 15 [Aerophys 1] = ", NK(15), 'H2SO4rate = ', H2SO4rate_o, &
            !      'Nkdtot = ', SUM(Nkd(:)), 'NKtot = ', SUM(Nk(:)), &
            !      'SFdtot = ', SUM(Mkd(:,srtso4)), 'SFtot = ', SUM(Mk(:,srtso4))
            !ENDIF
            !IF ((i_x .eq. 1) .AND. (i_y.eq.2)) THEN
            !      Write (6, *) 'Debug: BZ: before microphysics 1, Nk bin 14 (kg)= ', Nk(14)
            !      Write (6, *) 'Debug: BZ: before microphysics Mk(SO4) bin 1 (kg) = ', Mk(1,1)
            !ENDIF

            ! Before doing any cond/nucl/coag, check if there's any aerosol in
            ! the current box
            TOT_NK = SUM(NK)
            IF(TOT_NK .lt. 1.e-5_fp) THEN
               IF( .NOT. SPINUP(5.0)) THEN
                  print *,'No aerosol in box (i_box, i_x, i_y) ',i_box, i_x, i_y,'-->SKIP'
               ENDIF
               CYCLE
            ENDIF

            !---------------------------------------
            ! Condensation and nucleation (coupled)
            !---------------------------------------
            IF ( COND .AND. NUCL .AND. H2SO4rate_o > 0.e0_fp) THEN

               !if(printdebug .and. i==iob.and.j==job.and.l==lob) ERRORSWITCH =.TRUE.

               CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
               !IF ((i_x .eq. 1) .AND. (i_y.eq.2)) THEN
               !   Write (6, *) 'Debug: BZ: before cond_nuc, Nk bin 14 (kg)= ', Nk(14)
               !   Write (6, *) 'Debug: BZ: before cond_nuc, Mk(SO4) bin 1 (kg)= ', Mk(1,1)
               !ENDIF
               CALL COND_NUC(Nk,Mk,Gc,Nkout,Mkout,Gcout,fn,fn1, &
                              H2SO4rate_o,Dt,num_iter,Nknuc,Mknuc,Nkcond,Mkcond, &
                              ionrate, surf_area, BOXVOL, BOXMASS, TEMPTMS, PRES, &
                              RHTOMAS, ERRORSWITCH, i_lev)
               
                  
               IF ( ERRORSWITCH ) THEN
                  !PRINT *,'Aerophys: found error at (i_box, i_x, i_y)',i_box, i_x, i_y
                  WRITE(errMsg, '(A,I0,A,I0,A,I0)') &
                     'Aerophys, after cond_nuc (2-D): found error at i_box=', i_box, &
                     ' i_x=', i_x, &
                     ' i_y=', i_y
                  CALL ERROR_STOP( errMsg, thisLoc)
                  !CALL ERROR_STOP('AEROPHYS','After cond_nuc')
               ENDIF
              
               ERR_VAR = 'Gcout'
               ERR_MSG = 'After COND_NUC'
               ! check for NaN and Inf (win, 10/4/08)
               DO i_species = 1, ICOMPHARD
                  ERR_IND(1) = i_box
                  ERR_IND(2) = i_x
                  ERR_IND(3) = i_y
                  ERR_IND(4) = 0
      !             IF (SPINUP(14.0) .and. Gcout(jc) /= Gcout(jc) ) THEN
                  IF( SPINUP(14.0) .AND. IT_IS_NAN( Gcout(i_species) ) ) THEN
                     Gcout(i_species) = 0.0e+0_fp ! reset Nan to zero during spinup, bc 18/12/23
                     print*,'(2-D): Reset Gcout NaN to zero at ',i_box, i_x, i_y
                  ELSEIF ( SPINUP(14.0) .AND. .not. IT_IS_FINITE( Gcout(i_species) ) ) THEN
                     Gcout(i_species) = 0.0e+0_fp ! reset Inf to zero during spinup, bc 18/12/23
                     print*,'(2-D): Reset Gcout Inf to zero at ',i_box, i_x, i_y
                  ELSE
                  ! call check_value( Gcout(N), ERR_IND, ERR_VAR, ERR_MSG )
                  ENDIF
                  !if( IT_IS_FINITE(Gcout(jc))) then
                  !   print *,'xxxxxxxxx Found Inf in Gcout xxxxxxxxxxxxxx'
                  !   print *,'Location ',I,J,L, 'comp',jc
                  !   call debugprint( Nkout, Mkout, i,j,l,'After COND_NUC')
                  !   stop
                  !endif
               ENDDO

               !get nucleation diagnostic
               !DO ibin = 1, nbins
               !  NK(ibin) = NKnuc(ibin)
               !  DO N = 1, ICOMPHARD
               !      MK(ibin,N) = MKnuc(ibin,N)
               !  ENDDO
               !ENDDO

               !get condensation diagnostic
               !DO ibin = 1, nBins
               !  NK(ibin) = NKcond(ibin)
               !  DO N = 1, ICOMPHARD
               !      MK(ibin,N) = MKcond(ibin,N)
               !  ENDDO
               !ENDDO

               ! Update GC, NK, Mk
               Gc(srtnh4)=Gcout(srtnh4)
               Gc(srtso4)=Gcout(srtso4)
               !nucrate(j,l)=nucrate(j,l)+fn
               !nucrate1(j,l)=nucrate1(j,l)+fn1

               DO ibin = 1, nBins
                  NK(ibin) = NKout(ibin)
                  DO i_species = 1, ICOMPHARD
                     MK(ibin,i_species) = MKout(ibin,i_species)
                  ENDDO
               ENDDO
            ENDIF ! end of cond and nuc !
            ! IF ((i_x.eq.x_test) .AND. (i_y.eq.y_test)) THEN
            !       WRITE(6, *) "Debug (BZ): [TOMAS] NK 15 [After COND_NUC] = ", NK(15)
            ! ENDIF
            !IF ((i_x .eq. 1) .AND. (i_y.eq.2)) THEN
            !   WRITE (6, *) 'Debug: BZ: after cond_nuc, Nk bin 14 (kg)= ', Nk(14)
            !   WRITE (6, *) 'Debug: BZ: after cond_nuc, Mk(SO4) bin 1 (kg)= ', Mk(1,1)
            !ENDIF
            ! nitrogen and sulfur mass checks and fix
            !tot_n_1a = Gc(srtnh4)*14.e+0_fp/17.e+0_fp
            !do k=1,ibins
            !    tot_n_1a = tot_n_1a + Mk(k,srtnh4)*14.e+0_fp/18.e+0_fp
            !enddo
            !tot_s_1a = 0.e+0_fp
            !do k=1,ibins
            !    tot_s_1a = tot_s_1a + Mk(k,srtso4)*32.e+0_fp/96.e+0_fp
            !enddo

            CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
            !print *, 'mnfix in tomas_mod:677'
            CALL MNFIX( Nk, Mk, ERRORSWITCH )
            IF ( ERRORSWITCH ) THEN
               PRINT *,'Aerophys 2: (2-D): MNFIX found error at (i_box, i_x, i_y)',i_box, i_x, i_y
               IF( .not. SPINUP(14.0) ) THEN
                  CALL ERROR_STOP('AEROPHYS-MNFIX (2)','After cond/nucl')
               ELSE
                  PRINT *,'Let error go during spin up'
               ENDIF
            ENDIF
            ! IF ((i_x.eq.x_test) .AND. (i_y.eq.y_test)) THEN
            !       WRITE(6, *) "Debug (BZ): [TOMAS] NK 15 [Aerophys 2] = ", NK(15), &
            !       'Nkdtot = ', SUM(Nkd(:)), 'NKtot = ', SUM(Nk(:)), &
            !       'SFdtot = ', SUM(Mkd(:,srtso4)), 'SFtot = ', SUM(Mk(:,srtso4))
            ! ENDIF
            !-----------------------------
            ! Coagulation
            !-----------------------------

            !if(printdebug .and. i==iob.and.j==job.and.l==lob) ERRORSWITCH =.TRUE.
            IF( COAG )  THEN
               CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
               !IF ((i_x .eq. 1) .AND. (i_y.eq.2)) THEN
               !   Write (6, *) 'Debug: BZ: before coagulation, Nk bin 14 (kg)= ', Nk(14)
               !   Write (6, *) 'Debug: BZ: before coagulation, Mk(SO4) bin 1 (kg)= ', Mk(1,1)
               !ENDIF
               CALL MULTICOAG( Dt, Nk, Mk, BOXVOL, PRES, TEMPTMS, ERRORSWITCH )
               IF ( ERRORSWITCH ) THEN
                  print*,'(2-D) error after coagulation at (i_box, i_x, i_y) = ',i_box, i_x, i_y
                  !CALL DEBUGPRINT( Nk, Mk, I, J, L,'After coagulation' )
               ENDIF
               
               !Fix any inconsistency after coagulation (win, 4/18/06)
               CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
               !if(printdebug .and. i==iob.and.j==job.and.l==lob) &
            !     ERRORSWITCH=.true. !4/18/06 win
               !print *, 'mnfix in tomas_mod:719'
               CALL MNFIX( NK, MK, ERRORSWITCH )

               IF ( ERRORSWITCH ) THEN
                  PRINT *,'Aerophys 3: (2-D): MNFIX found error at (i_box, i_x, i_y) = ',i_box, i_x, i_y
                  IF( .not. SPINUP(14.0) ) THEN
                     CALL ERROR_STOP('AEROPHYS-MNFIX (3)', 'After COAGULATION'  )
                  ELSE
                     PRINT *,'Let error go during spin up'
                  ENDIF
               ENDIF

            ENDIF ! end of coagulation
            ! IF ((i_x.eq.x_test) .AND. (i_y.eq.y_test)) THEN
            !       WRITE(6, *) "Debug (BZ): [TOMAS] NK 15 [after coagulation] = ", NK(15), &
            !       'Nkdtot = ', SUM(Nkd(:)), 'NKtot = ', SUM(Nk(:)), &
            !       'SFdtot = ', SUM(Mkd(:,srtso4)), 'SFtot = ', SUM(Mk(:,srtso4))
            ! ENDIF
            !IF ((i_x .eq. 1) .AND. (i_y.eq.2)) THEN
            !      Write (6, *) 'Debug: BZ: after coagulation, Nk bin 14 (kg)= ', Nk(14)
            !      Write (6, *) 'Debug: BZ: after coagulation, Mk(SO4) bin 1 (kg)= ', Mk(1,1)
            !ENDIF
            ! Do water eqm at appropriate times
            CALL EZNH3EQM( Gc, Mk )
            CALL EZWATEREQM ( MK, RHTOMAS )
            !****************************
            ! End of aerosol dynamics
            !****************************
            !Fix any inconsistencies in M/N distribution (because of advection)
            CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)

            ! Make sure anything that leaves AEROPHYS is free of any error
            ! This MNFIX call could be temporary (?) or just leave it here and
            ! monitor if the error fixed is significantly large meaning some
            ! serious problem needs to be investigated
            !if(printdebug .and. i==iob.and.j==job.and.l==lob) ERRORSWITCH =.true.
            !print *, 'mnfix in tomas_mod:758'
            CALL MNFIX(NK,MK,ERRORSWITCH)
            IF ( ERRORSWITCH ) THEN
               PRINT *,'End of Aerophys 4 (2=D): MNFIX found error at (i_box, i_x, i_y)',i_box, i_x, i_y
               IF( .not. SPINUP(14.0) ) THEN
                  CALL ERROR_STOP('AEROPHYS-MNFIX (4)', 'End of microphysics')
               ELSE
                  PRINT *,'Let error go during spin up'
               ENDIF
            ENDIF
            ! IF ((i_x.eq.x_test) .AND. (i_y.eq.y_test)) THEN
            !       WRITE(6, *) "Debug (BZ): [TOMAS] NK 15 [Aerophys 4] = ", NK(15), &
            !       'Nkdtot = ', SUM(Nkd(:)), 'NKtot = ', SUM(Nk(:)), &
            !       'SFdtot = ', SUM(Mkd(:,srtso4)), 'SFtot = ', SUM(Mk(:,srtso4))
            ! ENDIF
            ! Swap Nk, Mk, and Gc arrays back to conc array
            ! convert unit from kg to molec/cm3
            DO ibin = 1, nBins
               id_tracer = id_NK01_p + (ibin-1) * nspc_p_tomas_tracer
               box_concnt_2D(i_x,i_y,id_tracer) = Nk(ibin)*Avo*1.0E+3_fp/BOXVOL
               DO i_species = 1, ICOMPHARD -2 
                  id_tracer = id_NK01_p+i_species+(ibin-1)*nspc_p_tomas_tracer
                  spc_name = spc_names_p(id_NK01_p+i_species)
                  id_tracer_1   = Ind_(TRIM(spc_name))
                  molwt_spc       = State_Chm%SpcData(id_tracer_1)%Info%MW_g
                  box_concnt_2D(i_x,i_y,id_tracer) = MK(ibin,i_species)*1.0E+3_fp/molwt_spc * Avo/BOXVOL
               ENDDO
               id_tracer = id_AW01_p + (ibin-1) * nspc_p_tomas_tracer
               molwt_spc       =  State_Chm%SpcData(id_AW01)%Info%MW_g
               box_concnt_2D(i_x,i_y,id_tracer) = MK(ibin,SRTH2O)*1.0E+3_fp/molwt_spc * Avo/BOXVOL
            ENDDO
            id_tracer = GET_PLUME_SPC_ID('H2SO4')
            molwt_spc       =  State_Chm%SpcData(Ind_('H2SO4'))%Info%MW_g
            box_concnt_2D(i_x,i_y,id_tracer) = GC(SRTSO4)*1.0E+3_fp/molwt_spc * Avo/BOXVOL

            ! Calculate NH3 gas lost to aerosol phase as NH4
            molwt_spc       =  State_Chm%SpcData(Ind_('NH3'))%Info%MW_g
            NH3_to_NH4 = box_concnt_2D(i_x,i_y,id_NH3_p)-GC(SRTNH4)*1.0E+3_fp/molwt_spc * Avo/BOXVOL
            ! Update the bulk NH4 aerosol species
            IF ( NH3_to_NH4 > 0e+0_fp ) THEN
               !Spc(id_NH4)%Conc(I,J,L) = Spc(id_NH4)%Conc(I,J,L) + &
               !                    NH3_to_NH4/17.e+0_fp*18.e+0_fp
               box_concnt_2D(i_x,i_y,id_NH4_p) =  box_concnt_2D(i_x,i_y,id_NH4_p) + &
                                 NH3_to_NH4
            ENDIF
            ! Update NH3 gas species (win, 10/6/08)
            ! plus tiny amount CEPS in case zero causes some problem
            molwt_spc       =  State_Chm%SpcData(Ind_('NH3'))%Info%MW_g
            box_concnt_2D(i_x,i_y,id_NH3_p) = GC(SRTNH4)*1.0E+3_fp/molwt_spc * Avo/BOXVOL  + CEPS !MUST CHECK THIS!! (win,9/26/08)
            !IF ((i_x .eq. 1) .AND. (i_y.eq.2)) THEN
            !      Write (6, *) 'Debug: BZ: after microphysics, Nk bin 14 (molec)= ', box_concnt_2D(i_x,i_y, 10+1+(14-1)*nspc_p_tomas_tracer)
            !      Write (6, *) 'Debug: BZ: after microphysics Mk(SO4) bin 1 (molec)= ', box_concnt_2D(i_x,i_y, 12)
            !ENDIF
         ENDDO
      ENDDO
      !$OMP END PARALLEL DO
      mass_S_H2SO4_2D = mass_S_H2SO4_2D + &
                SUM(box_concnt_2D(:,:, id_H2SO4_p)) * Vgrid_2D
      ! After TOMAS, exchange background NH3/NH4
      ! mass_NH3 = Spc(id_NH3)%Conc(i_lon,i_lat,i_lev)* Vgrid_EU + &
      !          SUM(box_concnt_2D(:,:,id_NH3_p))*Vgrid_2D
      ! box_concnt_2D(:,:,id_NH3_p) = mass_NH3  / &
      !                      (Vgrid_2D*n_x_max*n_y_max + Vgrid_EU)
      ! Spc(id_NH3)%Conc(i_lon,i_lat,i_lev) = mass_NH3 / &
      !                      (Vgrid_2D*n_x_max*n_y_max + Vgrid_EU)

      ! mass_NH4 = Spc(id_NH4)%Conc(i_lon,i_lat,i_lev)* Vgrid_EU + &
      !            SUM(box_concnt_2D(:,:,id_NH4_p))*Vgrid_2D
      ! box_concnt_2D(:,:,id_NH4_p) = mass_NH4 / &
      !                      (Vgrid_2D*n_x_max*n_y_max + Vgrid_EU)
      ! Spc(id_NH4)%Conc(i_lon,i_lat,i_lev) = mass_NH4 / &
      !                (Vgrid_2D*n_x_max*n_y_max + Vgrid_EU)
      mass_NH3_consum_plume(i_lon, i_lat, i_lev)  =   &
      mass_NH3_consum_plume(i_lon, i_lat, i_lev)  -   SUM(box_concnt_2D(:,:,id_NH3_p))*Vgrid_2D

      mass_NH4_consum_plume(i_lon, i_lat, i_lev)  =   &
      mass_NH4_consum_plume(i_lon, i_lat, i_lev)  -   SUM(box_concnt_2D(:,:,id_NH4_p))*Vgrid_2D
#endif
      Plume2d_curr%CONCNT2d = box_concnt_2D
      ! Concentration criteria: if SO2 concentration no larger than background, then dissolve the plume
      Core_Mean = Compute_plume_core_conc (Plume2d_curr%CONCNT2d(:,:,id_SO2_p))
      Conc_background_Mean = Spc(id_SO2)%Conc(i_lon,i_lat,i_lev)
      IF (Core_Mean .lt. Conc_background_Mean ) THEN
         Plume2d_curr%IsDissolve = .True.
      ENDIF
#ifdef TOMAS
   !Write (6, *) 'Debug: BZ:  (2-D test grid) after plume microphysics, SO4 (molec)= ', box_concnt_2D(x_test,y_test, id_SO4_p)
   !Write (6, *) 'Debug: BZ: (2-D test grid) after plume microphysics, Nk bin 15 (molec)= ', box_concnt_2D(x_test,y_test, 53)
   !Write (6, *) 'Debug: BZ: (2-D test grid) after plume microphysics, Nk bin 15 (molec)= ', Plume2d_curr%CONCNT2d(x_test,y_test, 53)
   !Write (6, *) 'Debug: BZ: (2-D test grid) after plume microphysics, SF bin 15 (molec)= ', box_concnt_2D(x_test,y_test, 54)
#endif 
1112  CONTINUE
      ! mass_S_SO2_5_2D = mass_S_SO2_5_2D + SUM(Plume2d_curr%CONCNT2D(:,:,id_SO2_p)) * Vgrid_2D
      ! mass_S_SO4_5_2D = mass_S_SO4_5_2D + SUM(Plume2d_curr%CONCNT2D(:,:,id_SO4_p)) * Vgrid_2D
      ! Write (6, *) 'Debug: BZ:  (2-D chem after), SO2 mass = ', SUM(box_concnt_2D(:,:, id_SO2_p)) * Vgrid_2D, &
      !                'SO4 mass = ',  SUM(box_concnt_2D(:,:, id_SO4_p)) * Vgrid_2D
      ! Write (6, *) 'Debug: BZ:  (2-D chem after - 2), SO2 mass = ', SUM(Plume2d_curr%CONCNT2d(:,:, id_SO2_p)) * Vgrid_2D, &
      !                'SO4 mass = ',  SUM(Plume2d_curr%CONCNT2d(:,:, id_SO4_p)) * Vgrid_2D
      Plume2d_curr => Plume2d_curr%next
   ENDDO

401 CONTINUE

   IF(.NOT.ASSOCIATED(Plume1d_head)) GOTO 400
   Plume1d_curr => Plume1d_head
   DO WHILE(ASSOCIATED(Plume1d_curr))
      
      i_box         = Plume1d_curr%label
      i_lon         = Plume1d_curr%lon_ind
      i_lat         = Plume1d_curr%lat_ind
      i_lev         = Plume1d_curr%lev_ind
      !WRITE(6,*) 'debug (BZ): solve plume Chemistry in plume 1-D box: ', i_box
      !Test if we need to do the chemistry for box (I,J,L), otherwise move onto the next box.
      ! MaxChemLev = MaxStratLev = 59 
      ! Hard coded in GeosUtil/gc_grid_mod.F90
      ! BZ: Maybe if plume reach out of chem grid, directly release species to Eulerian grid and delete the plume box?
      

      ! Vgrid_1D_SO4        =  Plume1d_curr%Ra(id_SO4_p) * Plume1d_curr%Rb(id_SO4_p) * Plume1d_curr%length * 1.0e+6_fp ! [cm3]
      ! Vgrid_1D_SO2        =  Plume1d_curr%Ra(id_SO2_p) * Plume1d_curr%Rb(id_SO2_p) * Plume1d_curr%length * 1.0e+6_fp ! [cm3]
      ! Vgrid_1D_PH2SO4        =  Plume1d_curr%Ra(id_PH2SO4_p) * Plume1d_curr%Rb(id_PH2SO4_p) * Plume1d_curr%length * 1.0e+6_fp ! [cm3]
      Vgrid_1D        =  Plume1d_curr%Ra * Plume1d_curr%Rb * Plume1d_curr%length * 1.0e+6_fp ! [cm3]
      Vgrid_EU        =  State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp
      K_SO2_OH        =  RXNRATE_CONST_KPP(i_lon,i_lat,i_lev , SO2_OH_RXN_ID)
      box_concnt_1D   =  Plume1d_curr%CONCNT1d
      ! mass_S_SO2_4_1D = mass_S_SO2_4_1D + SUM(Plume1d_curr%CONCNT1D(:,id_SO2_p)) * Vgrid_1D
      ! mass_S_SO4_4_1D = mass_S_SO4_4_1D + SUM(Plume1d_curr%CONCNT1D(:,id_SO4_p)) * Vgrid_1D
      ! Write (6, *) 'Debug: BZ:  (1-D chem before - 1), SO2 mass = ', SUM(box_concnt_1D(:,id_SO2_p)) * Vgrid_1D, &
      !          'SO4 mass = ',  SUM(box_concnt_1D(:, id_SO4_p)) * Vgrid_1D
      ! Write (6, *) 'Debug: BZ:  (1-D chem before - 2), SO2 mass = ', SUM(Plume1d_curr%CONCNT1d(:, id_SO2_p)) * Vgrid_1D, &
      !          'SO4 mass = ',  SUM(Plume1d_curr%CONCNT1d(:, id_SO4_p)) * Vgrid_1D

      IF ( .not. State_Met%InChemGrid(i_lon,i_lat,i_lev) ) THEN
         WRITE(6,*) 'Debug (BZ): 1-D plume segment outside chem grid: (i_box, i_lon, i_lat, i_lev): ',    &
         Plume1d_curr%label, i_lon, i_lat, i_lev
         Plume1d_curr%IsDissolve = .True.
      ENDIF
      IF (Plume1d_curr%IsDissolve) THEN
         GOTO 1113
      ENDIF

      ! write(6,*) 'debug (BZ): (1-D) Before exchange, background OH conc =  ', Spc(id_OH)%Conc(i_lon,i_lat,i_lev), &
      ! '; background HO2 conc = ', Spc(id_HO2)%Conc(i_lon,i_lat,i_lev), 'Plume number: ', Plume1d_curr%label, &
      ! 'Plume location (X, Y, L) = ', i_lon, i_lat, i_lev
      
      ! Exchange OH/HO2 with background before chemistry
      ! mass_OH     = Spc(id_OH)%Conc(i_lon,i_lat,i_lev)*Vgrid_EU +      &
      !                   SUM(box_concnt_1D(:,id_OH_p))*Vgrid_1D
      ! mass_HO2    = Spc(id_HO2)%Conc(i_lon,i_lat,i_lev)*Vgrid_EU +     &
      !                   SUM(box_concnt_1D(:,id_HO2_p))*Vgrid_1D

      ! full_exchange_conc = mass_OH  / (Vgrid_EU + Vgrid_1D*n_slab_max)

      ! box_concnt_1D(:,id_OH_p)   = full_exchange_conc * Radical_exc_factor
      ! Spc(id_OH)%Conc(i_lon,i_lat,i_lev) = (mass_OH -SUM(box_concnt_1D(:,id_OH_p))*Vgrid_1D ) / &
      !                                        Vgrid_EU
      
      ! full_exchange_conc = mass_HO2  / (Vgrid_EU + Vgrid_1D*n_slab_max)

      ! box_concnt_1D(:,id_HO2_p)   = full_exchange_conc * Radical_exc_factor
      ! Spc(id_HO2)%Conc(i_lon,i_lat,i_lev) = (mass_HO2 -SUM(box_concnt_1D(:,id_HO2_p))*Vgrid_1D ) / &
      !                                        Vgrid_EU
      mass_OH     = Spc(id_OH)%Conc(i_lon,i_lat,i_lev)*Vgrid_EU
      mass_HO2    = Spc(id_HO2)%Conc(i_lon,i_lat,i_lev)*Vgrid_EU

      ! full_exchange_conc = mass_OH  / (Vgrid_EU + Vplume_2D_tot(i_lon, i_lat, i_lev) + Vplume_1D_tot(i_lon, i_lat, i_lev))
      full_exchange_conc =  Spc(id_OH)%Conc(i_lon,i_lat,i_lev)
      box_concnt_1D(:,id_OH_p)   = full_exchange_conc * Radical_exc_factor
      mass_OH_consum_plume(i_lon, i_lat, i_lev)  =  mass_OH_consum_plume(i_lon, i_lat, i_lev) + SUM(box_concnt_1D(:,id_OH_p))*Vgrid_1D

      ! full_exchange_conc = mass_HO2  / (Vgrid_EU + Vplume_2D_tot(i_lon, i_lat, i_lev) + Vplume_1D_tot(i_lon, i_lat, i_lev))
      full_exchange_conc = Spc(id_HO2)%Conc(i_lon,i_lat,i_lev)
      box_concnt_1D(:,id_HO2_p)   = full_exchange_conc * Radical_exc_factor
      mass_HO2_consum_plume(i_lon, i_lat, i_lev)  =  mass_HO2_consum_plume(i_lon, i_lat, i_lev) + SUM(box_concnt_1D(:,id_HO2_p))*Vgrid_1D
      ! write(6,*) 'debug (BZ): (1-D) After exchange, background OH conc =  ', Spc(id_OH)%Conc(i_lon,i_lat,i_lev), &
      ! '; background HO2 conc = ', Spc(id_HO2)%Conc(i_lon,i_lat,i_lev), 'Plume number: ', Plume1d_curr%label, &
      ! 'Plume location (X, Y, L) = ', i_lon, i_lat, i_lev
!#ifdef TOMAS
      !Write (6, *) 'Debug: BZ: before chemistry, Nk bin 14 (molec)= ', box_concnt_2D(1,2, 10+1+(14-1)*nspc_p_tomas_tracer)
      !Write (6, *) 'Debug: BZ: before chemistry Mk(SO4) bin 1 (molec)= ', box_concnt_2D(1,2, 12)
!#endif
      !WRITE(6,*) 'Debug (BZ): Euleria grid: Conc of SO2 ', Spc(id_SO2)%Conc(i_lon,i_lat,i_lev)
      !WRITE(6,*) 'Debug (BZ): Euleria grid: Conc of SO4 ', Spc(id_SO4)%Conc(i_lon,i_lat,i_lev)
      !WRITE(6,*) 'Debug (BZ): Euleria grid: Conc of OH ', Spc(id_OH)%Conc(i_lon,i_lat,i_lev)

      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc before Chem: Conc of SO2 ', SUM(Plume2d_curr%CONCNT2d(:,:,id_SO2_p))/(n_x_max*n_y_max)
      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc before Chem: Conc of SO4 ', SUM(Plume2d_curr%CONCNT2d(:,:,id_SO4_p))/(n_x_max*n_y_max)
      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc before Chem: Conc of OH ', SUM(Plume2d_curr%CONCNT2d(:,:,id_OH_p))/(n_x_max*n_y_max)

      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc before Chem: Conc of SO2 ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_SO2_p)
      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc before Chem: Conc of SO4 ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_SO4_p)
      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc before Chem: Conc of OH ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_OH_p)
      box_concnt_1D_prev = box_concnt_1D
      n_debug_chem     =  0
      debug_islab      =  0
      debug_ibox       =  0
      debug_status     =  0
      debug_value      =  0.0_fp
      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE(i_slab) &
      !$OMP PRIVATE(chem_status, chem_debug_value)&
      !$OMP PRIVATE(C_before_Chem, C_after_Chem)
      !!$OMP COLLAPSE( 2                    )
      !!$OMP SCHEDULE( DYNAMIC, 24          )&
      !!$OMP REDUCTION( +:n_debug_chem        )

      DO i_slab = 1, n_slab_max, 1
         chem_status     =  0
         chem_debug_value = 0.0_fp

         C_before_Chem(id_SO2_p)     =  box_concnt_1D(i_slab,id_SO2_p)
         C_before_Chem(id_SO4_p)     =  box_concnt_1D(i_slab,id_SO4_p)
         C_before_Chem(id_OH_p)      =  box_concnt_1D(i_slab,id_OH_p)
         C_before_Chem(id_HO2_p)     =  box_concnt_1D(i_slab,id_HO2_p)
         C_before_Chem(id_PH2SO4_p)  =  0.0_fp

         C_after_Chem              =  C_before_Chem
         ! CHEM_SO2_OH_PLUME (dt, K, SO2, OH, SO4, HO2, PH2SO4, i_x, i_y, i_box) 
         CALL CHEM_SO2_OH_PLUME( Dt, K_SO2_OH, C_after_Chem(id_SO2_p),       &
                                    C_after_Chem(id_OH_p), C_after_Chem(id_SO4_p), &
                                    C_after_Chem(id_HO2_p), C_after_Chem(id_PH2SO4_p), &
                                    chem_status,chem_debug_value )

         IF ( chem_status /= 0 ) THEN
            !$OMP CRITICAL(debug_store)
            IF ( n_debug_chem < n_debug_chem_max ) THEN
               n_debug_chem = n_debug_chem + 1
               debug_islab(n_debug_chem)     = i_slab
               debug_ibox(n_debug_chem)   = i_box
               debug_status(n_debug_chem) = chem_status
               debug_value(n_debug_chem)  = chem_debug_value
            ENDIF
            !$OMP END CRITICAL(debug_store)
         ENDIF

         box_concnt_1D(i_slab,id_SO2_p)    = REAL( C_after_Chem(id_SO2_p), kind=fp )
         box_concnt_1D(i_slab,id_SO4_p)    = REAL( C_after_Chem(id_SO4_p), kind=fp )
         box_concnt_1D(i_slab,id_HO2_p)    = REAL( C_after_Chem(id_HO2_p), kind=fp )
         box_concnt_1D(i_slab,id_OH_p)     = REAL( C_after_Chem(id_OH_p), kind=fp )
         box_concnt_1D(i_slab,id_PH2SO4_p) = REAL( C_after_Chem(id_PH2SO4_p), kind=fp )
         H2SO4_RATE_1D(i_slab) = C_after_Chem(id_PH2SO4_p) / AVO * 98.e-3_fp * &
                           (Vgrid_1D) / Dt  ! kg s-1 box-1
      ENDDO
      !$OMP END PARALLEL DO
      
      ! Chemistry debug output
      DO N = 1, n_debug_chem
         SELECT CASE ( debug_status(N) )
         CASE (1)
            WRITE(6,*) "Skip chem: low SO2 at 1-D segment [slab,box]= ", &
                     debug_islab(N), debug_ibox(N), &
                  " value = ", debug_value(n)
         CASE (2)
            WRITE(6,*) "Skip chem: low OH at 1-D segment [slab,box]= ", &
                     debug_islab(N), debug_ibox(N), &
                  " value = ", debug_value(n)
         CASE (3)
            WRITE(6,*) "Skip chem: low K[SO2][OH]dt at 1-D segment [slab,box] = ", &
                     debug_islab(N), debug_ibox(N), &
                     " value = ", debug_value(N)
         END SELECT
      ENDDO

      
      !box_concnt_1D(:, id_PH2SO4_p) = 0.0_fp
      mass_S_SO2_r2_1D = mass_S_SO2_r2_1D + SUM(box_concnt_1D(:,id_SO2_p)-box_concnt_1D_prev(:,id_SO2_p))* Vgrid_1D
      mass_S_SO4_r2_1D = mass_S_SO4_r2_1D + SUM(box_concnt_1D(:,id_SO4_p)-box_concnt_1D_prev(:,id_SO4_p))* Vgrid_1D

      ! ! Exchange OH/HO2 with background after chemistry
      ! mass_OH     = Spc(id_OH)%Conc(i_lon,i_lat,i_lev)*Vgrid_EU +      &
      !                   SUM(box_concnt_1D(:,id_OH_p))*Vgrid_1D
      ! mass_HO2    = Spc(id_HO2)%Conc(i_lon,i_lat,i_lev)*Vgrid_EU +     &
      !                   SUM(box_concnt_1D(:,id_HO2_p))*Vgrid_1D

      ! full_exchange_conc = mass_OH  / (Vgrid_EU + Vgrid_1D*n_slab_max)

      ! box_concnt_1D(:,id_OH_p)   = full_exchange_conc * Radical_exc_factor
      ! Spc(id_OH)%Conc(i_lon,i_lat,i_lev) = (mass_OH -SUM(box_concnt_1D(:,id_OH_p))*Vgrid_1D ) / &
      !                                        Vgrid_EU
      
      ! full_exchange_conc = mass_HO2  / (Vgrid_EU + Vgrid_1D*n_slab_max)

      ! box_concnt_1D(:,id_HO2_p)   = full_exchange_conc * Radical_exc_factor
      ! Spc(id_HO2)%Conc(i_lon,i_lat,i_lev) = (mass_HO2 -SUM(box_concnt_1D(:,id_HO2_p))*Vgrid_1D ) / &
      !                                        Vgrid_EU
      ! write(6,*) 'debug (BZ): (1-D) After plume chem, after exchange, background OH conc =  ', Spc(id_OH)%Conc(i_lon,i_lat,i_lev), &
      ! '; background HO2 conc = ', Spc(id_HO2)%Conc(i_lon,i_lat,i_lev), 'Plume number: ', Plume1d_curr%label, &
      ! 'Plume location (X, Y, L) = ', i_lon, i_lat, i_lev
      
      ! After Chemistry
      ! Subtract remaining radical
      mass_OH_consum_plume(i_lon, i_lat, i_lev)  = mass_OH_consum_plume(i_lon, i_lat, i_lev)   - SUM(box_concnt_1D(:,id_OH_p))*Vgrid_1D
      mass_HO2_consum_plume(i_lon, i_lat, i_lev) = mass_HO2_consum_plume(i_lon, i_lat, i_lev)  - SUM(box_concnt_1D(:,id_HO2_p))*Vgrid_1D

#ifdef TOMAS
      !Write (6, *) 'Debug: BZ: after chemistry, Nk bin 14 (molec)= ', box_concnt_2D(1,2, 10+1+(14-1)*nspc_p_tomas_tracer)
      !Write (6, *) 'Debug: BZ: after chemistry Mk(SO4) bin 1 (molec)= ', box_concnt_2D(1,2, 12)
      !write(6,*) 'debug (BZ): solve plume microphysics in plume 1-D box: ', i_box
      ! Before TOMAS, read background NH3/NH4, assuming same size with SO4 grid
      mass_NH3 = Spc(id_NH3)%Conc(i_lon,i_lat,i_lev)* Vgrid_EU 
      ! full_exchange_conc = mass_NH3 / (Vgrid_EU + Vplume_2D_tot(i_lon, i_lat, i_lev) + Vplume_1D_tot(i_lon, i_lat, i_lev))
      full_exchange_conc =  Spc(id_NH3)%Conc(i_lon,i_lat,i_lev)
      box_concnt_1D(:,id_NH3_p) = full_exchange_conc
      mass_NH3_consum_plume(i_lon,i_lat,i_lev) = &
      mass_NH3_consum_plume(i_lon,i_lat,i_lev) + SUM(box_concnt_1D(:,id_NH3_p))*Vgrid_1D
      ! Spc(id_NH3)%Conc(i_lon,i_lat,i_lev) = mass_NH3 / &
      !          (Vgrid_1D*n_slab_max + Vgrid_EU)
      mass_NH4 = Spc(id_NH4)%Conc(i_lon,i_lat,i_lev)* Vgrid_EU
      ! full_exchange_conc = mass_NH4 / (Vgrid_EU + Vplume_2D_tot(i_lon, i_lat, i_lev) + Vplume_1D_tot(i_lon, i_lat, i_lev)) 
      full_exchange_conc = Spc(id_NH4)%Conc(i_lon,i_lat,i_lev)
      box_concnt_1D(:,id_NH4_p) = full_exchange_conc
      mass_NH4_consum_plume(i_lon,i_lat,i_lev) = &
      mass_NH4_consum_plume(i_lon,i_lat,i_lev) + SUM(box_concnt_1D(:,id_NH4_p))*Vgrid_1D
      ! Spc(id_NH4)%Conc(i_lon,i_lat,i_lev) =  mass_NH4 / &
      !          (Vgrid_1D*n_slab_max + Vgrid_EU)
      PRES    = State_Met%PMID(i_lon,i_lat,i_lev)*100.0 ! in Pa
      TEMPTMS = State_Met%T(i_lon,i_lat,i_lev)
      BOXMASS = State_Met%AD(i_lon,i_lat,i_lev)*Vgrid_1D/Vgrid_EU ! Dry air mass, kg
      RHTOMAS = State_Met%RH(i_lon,i_lat,i_lev)/ 1.e2
      IF ( RHTOMAS > 0.99 ) RHTOMAS = 0.99
      BOXVOL  = Vgrid_1D  !cm3
      ! Temporary use NK grid volume 
      !Vgrid_1D = Plume1d_curr%Ra(id_NK01_p) * Plume1d_curr%Rb(id_NK01_p) * Plume1d_curr%length * 1.0e+6_fp ! [cm3]
      ! TOMAS species unit should be in kg
      ! GC, MK, NK
      !$OMP PARALLEL DO         &
      !$OMP DEFAULT( SHARED )   &
      !$OMP PRIVATE( i_slab, i_species, i_species_1, ibin, id_tracer, id_tracer_1 )  &
      !$OMP PRIVATE(spc_name, molwt_spc, tot_n_1, tot_s_1, TOT_NK, TOT_MK)    &
      !$OMP PRIVATE(Gc, Gcd, GCout, Nk, Nkd, Nkout, Nknuc, Nkcond, MK, Mkd,Mkout, Mknuc, Mkcond)        &
      !$OMP PRIVATE( TRANSFER, H2SO4rate_o, NH4bulk, NH3_to_NH4, fn, fn1,num_iter, ionrate, surf_area)       & 
      !$OMP PRIVATE(ERRORSWITCH, PRINTNEG, PRINTDEBUG, ERRSPOT)  &
      !$OMP PRIVATE( ERR_VAR, ERR_MSG, ERR_IND )  &                   
      !$OMP SCHEDULE( DYNAMIC )
      DO i_slab = 1, n_slab_max, 1
         ERRORSWITCH = .FALSE.
         ! Initialize all condensible gas values to zero
         ! Gc(srtso4) will remain zero until within cond_nuc where the
         ! pseudo steady state H2SO4 concentration will be put in this place.
         Gc(:)       = 0.0_fp
          
         ! Swap Spc into Nk, Mk, Gc arrays, unit kg
         Nk(:)       = 0.0_fp
         Mk(:,:)     = 0.0_fp
         DO ibin = 1, nBins
            ! Mol weight for Nk is 1
            i_species_1 = id_NK01_p+(ibin-1)*nspc_p_tomas_tracer
            NK(ibin) = box_concnt_1D(i_slab,i_species_1)*BOXVOL/Avo/1.0E+3_fp
            DO i_species = 1, ICOMPHARD-2 ! skip NH4 and H2O
               i_species_1 = id_NK01_p+i_species+(ibin-1)*nspc_p_tomas_tracer
               spc_name = spc_names_p(id_NK01_p+i_species)
               id_tracer   = Ind_(TRIM(spc_name))
               molwt_spc       = State_Chm%SpcData(id_tracer)%Info%MW_g
               MK(ibin,i_species) = box_concnt_1D(i_slab,i_species_1)*BOXVOL/Avo*molwt_spc/1.0E+3_fp

               !IF( IT_IS_NAN( MK(ibin,i_species) ) ) THEN
               !    PRINT *,'+++++++ Found NaN in AEROPHYS ++++++++'
               !    PRINT *,'Location (i_plume, i_x, i_y):',i_box,i_x,i_y,'Bin',ibin,'comp',spc_name
               !ENDIF
            ENDDO
            id_tracer   = Ind_('AW01')
            molwt_spc       = State_Chm%SpcData(id_tracer)%Info%MW_g
            MK(ibin,SRTH2O) = box_concnt_1D(i_slab,id_AW01_p+(ibin-1)*nspc_p_tomas_tracer)*BOXVOL/Avo*molwt_spc/1.0E+3_fp
         ENDDO

         ! Get NH4 mass from the bulk mass and scale to bin with sulfate
         IF ( SRTNH4 > 0 ) THEN
            NH4bulk = box_concnt_1D(i_slab,id_NH4_p)*BOXVOL/Avo*18.0_fp/1.0E+3_fp
            CALL NH4BULKTOBIN( MK(:,SRTSO4), NH4bulk, TRANSFER )
            MK(1:nbins,SRTNH4) = TRANSFER(1:nbins)
            Gc(SRTNH4) = box_concnt_1D(i_slab,id_NH3_p)*BOXVOL/Avo*17.0_fp/1.0E+3_fp
         ENDIF
          ! Give it the pseudo-steady state value instead later (win,9/30/08)
          !GC(SRTSO4) = Spc(id_H2SO4)%Conc(I,J,L)
          
         H2SO4rate_o = H2SO4_RATE_1D(i_slab)  ! [kg s-1]

         ! nitrogen and sulfur mass checks
         ! get the total mass of N
         tot_n_1 = Gc(srtnh4)*14.e+0_fp/17.e+0_fp
         DO ibin=1,nbins
            tot_n_1 = tot_n_1 + Mk(ibin,srtnh4)*14.e+0_fp/18.e+0_fp
         ENDDO

         ! get the total mass of S
         tot_s_1 = H2SO4rate_o*Dt*32.e+0_fp/98.e+0_fp
         DO ibin=1,nbins
            tot_s_1 = tot_s_1 + Mk(ibin,srtso4)*32.e+0_fp/96.e+0_fp
         ENDDO

         !Do water eqm at appropriate times
         CALL EZWATEREQM( MK, RHTOMAS )

          !Fix any inconsistencies in M/N distribution (because of advection)
         CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
         CALL MNFIX( NK, MK, ERRORSWITCH )
         IF ( ERRORSWITCH ) THEN
            PRINT *,'Aerophys 1 (1-D): MNFIX found error at 1-D plume segment (i_box, i_slab) = ',i_box, i_slab
            CALL ERROR_STOP('AEROPHYS-MNFIX (1)','Enter microphys')
         ENDIF

         ! Before doing any cond/nucl/coag, check if there's any aerosol in
         ! the current box
         TOT_NK = SUM(NK)
         IF(TOT_NK .lt. 1.e-5_fp) THEN
            IF( .NOT. SPINUP(5.0)) THEN
                  print *,'No aerosol in 1-D plume segment (i_box, i_slab) ',i_box, i_slab, '-->SKIP'
            ENDIF
            CYCLE
         ENDIF

         !---------------------------------------
         ! Condensation and nucleation (coupled)
         !---------------------------------------
         IF ( COND .AND. NUCL .AND. H2SO4rate_o > 0.e0_fp) THEN

            !if(printdebug .and. i==iob.and.j==job.and.l==lob) ERRORSWITCH =.TRUE.
            CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
            CALL COND_NUC(Nk,Mk,Gc,Nkout,Mkout,Gcout,fn,fn1, &
                           H2SO4rate_o,Dt,num_iter,Nknuc,Mknuc,Nkcond,Mkcond, &
                           ionrate, surf_area, BOXVOL, BOXMASS, TEMPTMS, PRES, &
                           RHTOMAS, ERRORSWITCH, i_lev)
         
            IF ( ERRORSWITCH ) THEN
               !PRINT *,'Aerophys: found error at (i_box, i_x, i_y)',i_box, i_x, i_y
               WRITE(errMsg, '(A,I0,A,I0,A,I0)') &
                  'Aerophys, 1-D plume segment, after cond_nuc: found error at i_box=', i_box, &
                  ' i_s;ab=', i_slab
               CALL ERROR_STOP( errMsg, thisLoc)
               !CALL ERROR_STOP('AEROPHYS','After cond_nuc')
            ENDIF

            ERR_VAR = 'Gcout'
            ERR_MSG = 'After COND_NUC'
            ! check for NaN and Inf (win, 10/4/08)
            DO i_species = 1, ICOMPHARD
               !ERR_IND(1) = i_box
               !ERR_IND(2) = i_slab
               !ERR_IND(3) = 0
               !ERR_IND(4) = 0
   !             IF (SPINUP(14.0) .and. Gcout(jc) /= Gcout(jc) ) THEN
               IF( SPINUP(14.0) .AND. IT_IS_NAN( Gcout(i_species) ) ) THEN
                  Gcout(i_species) = 0.0e+0_fp ! reset Nan to zero during spinup, bc 18/12/23
                  print*,'(1-D) Reset Gcout NaN to zero at ',i_box, i_slab
               ELSEIF ( SPINUP(14.0) .AND. .not. IT_IS_FINITE( Gcout(i_species) ) ) THEN
                  Gcout(i_species) = 0.0e+0_fp ! reset Inf to zero during spinup, bc 18/12/23
                  print*,'(1-D) Reset Gcout Inf to zero at ',i_box, i_slab
               ELSE
               ! call check_value( Gcout(N), ERR_IND, ERR_VAR, ERR_MSG )
               ENDIF
               !if( IT_IS_FINITE(Gcout(jc))) then
               !   print *,'xxxxxxxxx Found Inf in Gcout xxxxxxxxxxxxxx'
               !   print *,'Location ',I,J,L, 'comp',jc
               !   call debugprint( Nkout, Mkout, i,j,l,'After COND_NUC')
               !   stop
               !endif
            ENDDO

            !get nucleation diagnostic
            !DO ibin = 1, nbins
            !  NK(ibin) = NKnuc(ibin)
            !  DO N = 1, ICOMPHARD
            !      MK(ibin,N) = MKnuc(ibin,N)
            !  ENDDO
            !ENDDO

            !get condensation diagnostic
            !DO ibin = 1, nBins
            !  NK(ibin) = NKcond(ibin)
            !  DO N = 1, ICOMPHARD
            !      MK(ibin,N) = MKcond(ibin,N)
            !  ENDDO
            !ENDDO

            ! Update GC, NK, Mk
            Gc(srtnh4)=Gcout(srtnh4)
            Gc(srtso4)=Gcout(srtso4)
            !nucrate(j,l)=nucrate(j,l)+fn
            !nucrate1(j,l)=nucrate1(j,l)+fn1

            DO ibin = 1, nBins
               NK(ibin) = NKout(ibin)
               DO i_species = 1, ICOMPHARD
                  MK(ibin,i_species) = MKout(ibin,i_species)
               ENDDO
            ENDDO

         ENDIF ! end of cond and nuc !

         ! nitrogen and sulfur mass checks and fix
         !tot_n_1a = Gc(srtnh4)*14.e+0_fp/17.e+0_fp
         !do k=1,ibins
         !    tot_n_1a = tot_n_1a + Mk(k,srtnh4)*14.e+0_fp/18.e+0_fp
         !enddo
         !tot_s_1a = 0.e+0_fp
         !do k=1,ibins
         !    tot_s_1a = tot_s_1a + Mk(k,srtso4)*32.e+0_fp/96.e+0_fp
         !enddo

         CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
         !print *, 'mnfix in tomas_mod:677'
         CALL MNFIX( Nk, Mk, ERRORSWITCH )
         IF ( ERRORSWITCH ) THEN
            PRINT *,'Aerophys 2 (1-D): MNFIX found error at (i_box, i_slab)',i_box, i_slab
            IF( .not. SPINUP(14.0) ) THEN
               CALL ERROR_STOP('AEROPHYS-MNFIX (2), 1-D','After cond/nucl')
            ELSE
               PRINT *,'Let error go during spin up'
            ENDIF
         ENDIF
         !-----------------------------
         ! Coagulation
         !-----------------------------

         !if(printdebug .and. i==iob.and.j==job.and.l==lob) ERRORSWITCH =.TRUE.
         IF( COAG )  THEN
            CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
            CALL MULTICOAG( Dt, Nk, Mk, BOXVOL, PRES, TEMPTMS, ERRORSWITCH )
            IF ( ERRORSWITCH ) THEN
               print*,'(1-D): error after coagulation at (i_box, i_slab) = ',i_box, i_slab
               !CALL DEBUGPRINT( Nk, Mk, I, J, L,'After coagulation' )
            ENDIF
         
            !Fix any inconsistency after coagulation (win, 4/18/06)
            CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)
            !if(printdebug .and. i==iob.and.j==job.and.l==lob) &
            !     ERRORSWITCH=.true. !4/18/06 win
            !print *, 'mnfix in tomas_mod:719'
            CALL MNFIX( NK, MK, ERRORSWITCH )

            IF ( ERRORSWITCH ) THEN
               PRINT *,'Aerophys 3 (1-D): MNFIX found error at (i_box, i_slab) = ',i_box, i_slab
               IF( .not. SPINUP(14.0) ) THEN
                  CALL ERROR_STOP('AEROPHYS-MNFIX (3)', 'After COAGULATION'  )
               ELSE
                  PRINT *,'Let error go during spin up'
               ENDIF
            ENDIF

         ENDIF ! end of coagulation

         ! Do water eqm at appropriate times
         CALL EZNH3EQM( Gc, Mk )
         CALL EZWATEREQM ( MK, RHTOMAS )
         !****************************
         ! End of aerosol dynamics
         !****************************
         !Fix any inconsistencies in M/N distribution (because of advection)
         CALL STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd)

         ! Make sure anything that leaves AEROPHYS is free of any error
         ! This MNFIX call could be temporary (?) or just leave it here and
         ! monitor if the error fixed is significantly large meaning some
         ! serious problem needs to be investigated
         !if(printdebug .and. i==iob.and.j==job.and.l==lob) ERRORSWITCH =.true.
         !print *, 'mnfix in tomas_mod:758'
         CALL MNFIX(NK,MK,ERRORSWITCH)
         IF ( ERRORSWITCH ) THEN
            PRINT *,'End of Aerophys (1-D): MNFIX found error at (i_box, i_slab)',i_box, i_slab
            IF( .not. SPINUP(14.0) ) THEN
               CALL ERROR_STOP('AEROPHYS-MNFIX (4)', 'End of microphysics')
            ELSE
               PRINT *,'Let error go during spin up'
            ENDIF
         ENDIF

         ! Swap Nk, Mk, and Gc arrays back to conc array
         ! convert unit from kg to molec/cm3
         DO ibin = 1, nBins
            id_tracer = id_NK01_p + (ibin-1) * nspc_p_tomas_tracer
            box_concnt_1D(i_slab,id_tracer) = Nk(ibin)*Avo*1.0E+3_fp/BOXVOL
            DO i_species = 1, ICOMPHARD -2 
               id_tracer = id_NK01_p+i_species+(ibin-1)*nspc_p_tomas_tracer
               spc_name = spc_names_p(id_NK01_p+i_species)
               id_tracer_1   = Ind_(TRIM(spc_name))
               molwt_spc       = State_Chm%SpcData(id_tracer_1)%Info%MW_g
               box_concnt_1D(i_slab, id_tracer) = MK(ibin,i_species)*1.0E+3_fp/molwt_spc * Avo/BOXVOL
            ENDDO
            id_tracer = id_AW01_p + (ibin-1) * nspc_p_tomas_tracer
            molwt_spc       =  State_Chm%SpcData(id_AW01)%Info%MW_g
            box_concnt_1D(i_slab,id_tracer) = MK(ibin,SRTH2O)*1.0E+3_fp/molwt_spc * Avo/BOXVOL
         ENDDO

         id_tracer = GET_PLUME_SPC_ID('H2SO4')
         molwt_spc       =  State_Chm%SpcData(Ind_('H2SO4'))%Info%MW_g
         box_concnt_1D(i_slab,id_tracer) = GC(SRTSO4)*1.0E+3_fp/molwt_spc * Avo/BOXVOL

         ! Calculate NH3 gas lost to aerosol phase as NH4
         molwt_spc       =  State_Chm%SpcData(Ind_('NH3'))%Info%MW_g
         NH3_to_NH4 = box_concnt_1D(i_slab,id_NH3_p)-GC(SRTNH4)*1.0E+3_fp/molwt_spc * Avo/BOXVOL
         ! Update the bulk NH4 aerosol species
         IF ( NH3_to_NH4 > 0e+0_fp ) THEN
         !Spc(id_NH4)%Conc(I,J,L) = Spc(id_NH4)%Conc(I,J,L) + &
         !                    NH3_to_NH4/17.e+0_fp*18.e+0_fp
         box_concnt_1D(i_slab, id_NH4_p) =  box_concnt_1D(i_slab, id_NH4_p) + &
                              NH3_to_NH4
         ENDIF
         ! Update NH3 gas species (win, 10/6/08)
         ! plus tiny amount CEPS in case zero causes some problem
         molwt_spc       =  State_Chm%SpcData(Ind_('NH3'))%Info%MW_g
         box_concnt_1D(i_slab,id_NH3_p) = GC(SRTNH4)*1.0E+3_fp/molwt_spc * Avo/BOXVOL  + CEPS !MUST CHECK THIS!! (win,9/26/08)
      ENDDO
      !$OMP END PARALLEL DO
      
      mass_S_H2SO4_1D = mass_S_H2SO4_1D + &
                SUM(box_concnt_1D(:, id_H2SO4_p)) * Vgrid_1D
      ! After TOMAS, read background NH3/NH4, assuming same size with SO4 grid
      ! mass_NH3 = Spc(id_NH3)%Conc(i_lon,i_lat,i_lev)* Vgrid_EU + &
      !          SUM(box_concnt_1D(:,id_NH3_p))*Vgrid_1D
      ! box_concnt_1D(:,id_NH3_p) = mass_NH3 / &
      !          (Vgrid_1D*n_slab_max + Vgrid_EU)
      ! Spc(id_NH3)%Conc(i_lon,i_lat,i_lev) = mass_NH3 / &
      !          (Vgrid_1D*n_slab_max + Vgrid_EU)
      ! mass_NH4 = Spc(id_NH4)%Conc(i_lon,i_lat,i_lev)* Vgrid_EU + &
      !          SUM(box_concnt_1D(:,id_NH4_p))*Vgrid_1D         
      ! box_concnt_1D(:,id_NH4_p) = mass_NH4 / &
      !          (Vgrid_1D*n_slab_max + Vgrid_EU)
      ! Spc(id_NH4)%Conc(i_lon,i_lat,i_lev) =  mass_NH4 / &
      !          (Vgrid_1D*n_slab_max + Vgrid_EU)
      mass_NH3_consum_plume(i_lon, i_lat, i_lev)  =   &
      mass_NH3_consum_plume(i_lon, i_lat, i_lev)  -   SUM(box_concnt_1D(:,id_NH3_p))*Vgrid_1D

      mass_NH4_consum_plume(i_lon, i_lat, i_lev)  =   &
      mass_NH4_consum_plume(i_lon, i_lat, i_lev)  -   SUM(box_concnt_1D(:,id_NH4_p))*Vgrid_1D
#endif
      Plume1d_curr%CONCNT1d = box_concnt_1D
      Core_Mean = SUM(box_concnt_1D(:,id_SO2_p))/REAL(n_slab_max, fp)
      Conc_background_Mean = Spc(id_SO2)%Conc(i_lon,i_lat,i_lev)
      IF (Core_Mean .lt. Conc_background_Mean ) THEN
         Plume1d_curr%IsDissolve = .True.
      ENDIF
1113  CONTINUE
      ! mass_S_SO2_5_1D = mass_S_SO2_5_1D + SUM(Plume1d_curr%CONCNT1D(:,id_SO2_p)) * Vgrid_1D
      ! mass_S_SO4_5_1D = mass_S_SO4_5_1D + SUM(Plume1d_curr%CONCNT1D(:,id_SO4_p)) * Vgrid_1D
      ! Write (6, *) 'Debug: BZ:  (1-D chem after - 1), SO2 mass = ', SUM(box_concnt_1D(:,id_SO2_p)) * Vgrid_1D, &
      !          'SO4 mass = ',  SUM(box_concnt_1D(:, id_SO4_p)) * Vgrid_1D
      ! Write (6, *) 'Debug: BZ:  (1-D chem after - 2), SO2 mass = ', SUM(Plume1d_curr%CONCNT1d(:, id_SO2_p)) * Vgrid_1D, &
      !          'SO4 mass = ',  SUM(Plume1d_curr%CONCNT1d(:, id_SO4_p)) * Vgrid_1D
      Plume1d_curr => Plume1d_curr%next
   ENDDO ! DO WHILE(ASSOCIATED(Plume1d_curr))         
      !========================================================================
      ! Temporary implement TOMAS microphysics here so it can be moved to 
      ! a separate module if necessary
      ! TOMAS module: AEROPHYS                            
      !========================================================================
      !                                           
!       DO i_y = 1, n_y_max, 1
!         DO i_x = 1, n_x_max, 1
!       ! Some species (counter or diagnostic species) needed to be zero out before solving KPP chemistry
!       ! See details in Do_Chemistry in Fullchem_mod.F90
!       DO i_species = 1, n_species
!         ! Get info about this species from the species database
!         SpcInfo      => State_Chm%SpcData(i_species)%Info

!         ! isoprene oxidation counter species
!         IF ( TRIM( SpcInfo%Name ) == 'LISOPOH' .or. &
!               TRIM( SpcInfo%Name ) == 'LISOPNO3' ) THEN
!             Plume2d_curr%CONCNT2d(:,:,i_species) = 0.0_fp
!         ENDIF

!        ! aromatic oxidation counter species
!         IF ( Input_Opt%LSOA .or. Input_Opt%LSVPOA ) THEN
!             SELECT CASE ( TRIM( SpcInfo%Name ) )
!               CASE ( 'LBRO2H', 'LBRO2N', 'LTRO2H', 'LTRO2N', &
!                       'LXRO2H', 'LXRO2N', 'LNRO2H', 'LNRO2N' )
!                   Plume2d_curr%CONCNT2d(:,:,i_species) = 0.0_fp
!             END SELECT
!         ENDIF

!         ! Sulfate gas/cloud prod diagnostic species
!         IF ( TRIM( SpcInfo%Name ) == 'PH2SO4' .or. &
!               TRIM( SpcInfo%Name ) == 'PSO4AQ' ) THEN
!             Plume2d_curr%CONCNT2d(:,:,i_species) = 0.0_fp
!         ENDIF

!         ! Free pointer
!         SpcInfo => NULL()
!       ENDDO

!       !!! (BZ) Test output: Will Photolysis rate array be initalized every dynamic timestep?
!       ! State_Chm%Phot%ZPJ(i_lev,n_photoRxn,i_lat,i_lon)
!       ! WRITE(6,*) 'Debug (BZ): Do_Chemistry  (In plume): PhotoRxn: O3 + hv -> O2 + O; rate: ', State_Chm%Phot%ZPJ(39,RXN_O3_1,23,40)
!       ! WRITE(6,*) 'Debug (BZ): Do_Chemistry  (In plume): PhotoRxn:  O3 + hv -> O2 + O(1D); rate: ', State_Chm%Phot%ZPJ(39,RXN_O3_2,23,40)

!       !========================================================================
!       ! MAIN LOOP: Compute reaction rates and call chemical solver
!       !========================================================================
      
!       !$OMP PARALLEL DO           &
!       !$OMP DEFAULT( SHARED     ) &
!       !$OMP PRIVATE(i_y,i_x, Thread)&
!       !$OMP PRIVATE( SO4_FRAC, IERR,     RCNTRL,  ISTATUS,   RSTATE     )&
!       !$OMP PRIVATE( ICNTRL,   C_before_integrate )&
!       !$OMP PRIVATE( SpcID,    KppID,    F,       P            )&
!       !$OMP PRIVATE( SR               )&
!       !$OMP PRIVATE( SIZE_RES, LWC                                            )&
!       !$OMP COLLAPSE( 2                    )&
!       !$OMP SCHEDULE( DYNAMIC, 24          )&
!       !$OMP REDUCTION( +:errorCount        )
!       DO i_y = 1, n_y_max, 1
!         DO i_x = 1, n_x_max, 1
!           ! Skip to the end of the loop if we have failed integration twice
!           IF ( Failed2x ) CYCLE
!            ! Initialize private loop variables for each (i_x, i_y)
!           IERR      = 0                        ! KPP success or failure flag
!           ISTATUS   = 0.0_dp                   ! Rosenbrock output
!           ICNTRL    = 0                        ! Rosenbrock input (integer)
!           RCNTRL    = 0.0_fp                   ! Rosenbrock input (real)
!           RSTATE    = 0.0_dp                   ! Rosenbrock output
!           SO4_FRAC  = 0.0_fp                   ! Frac of SO4 avail for photolysis
!           P         = 0                        ! GEOS-Chem photolyis species ID
!           !LCH4      = 0.0_fp                   ! P/L diag: Methane loss rate
!           !PCO_TOT   = 0.0_fp                   ! P/L diag: Total P(CO)
!           !PCO_CH4   = 0.0_fp                   ! P/L diag: P(CO) from CH4
!           !PCO_NMVOC = 0.0_fp                   ! P/L diag: P(CO) from NMVOC
!           SR        = 0.0_fp                   ! Enhancement to O2 catalysis rate
!           LWC       = 0.0_fp                   ! Liquid water content
!           SIZE_RES  = .FALSE.                  ! Size resolved calculation?
!           C         = 0.0_dp                   ! KPP species conc's
!           RCONST    = 0.0_dp                   ! KPP rate constants
!           PHOTOL    = 0.0_dp                   ! Photolysis array for KPP
!           K_CLD     = 0.0_dp                   ! Sulfur in-cloud rxn het rates
!           K_MT      = 0.0_dp                   ! Sulfur sea salt rxn het rates
!           CFACTOR   = 1.0_dp                   ! KPP conversion factor
!           SRO3      = 0.0_dp                   ! Enhanced sulfate production of
!           SRHOBr    = 0.0_dp                   !  O3, HOBr, HCl in size-resolved
!           SRHOCl    = 0.0_dp                   !  cloud droplets
! #ifdef MODEL_CLASSIC
! #ifndef NO_OMP
!           Thread    = OMP_GET_THREAD_NUM() + 1 ! OpenMP thread number
! #endif
! #endif
! #ifdef KPP_INTEGRATOR_AUTOREDUCE
!        ! Per discussions for Lin et al., force keepActive throughout the
!        ! atmosphere if keepActive option is enabled. (hplin, 2/9/22)
!        CALL fullchem_AR_SetKeepActive( option=.TRUE. )
! #endif
!           ! Get photolysis rates (daytime only)
!           ! Update SUNCOSmid threshold from 0 to cos(98 degrees)
!           ! Loop over the FAST-JX photolysis species
!           ! IF ( State_Met%SUNCOSmid(i_lon,i_lat) > -0.1391731e+0_fp ) THEN
!           !  ! Only proceed if doing photolysis
!           !  IF ( Input_Opt%Do_Photolysis ) THEN
!           !    DO i_phot = 1, State_Chm%Phot%nMaxPhotRxns

!                 ! Copy photolysis rate from FAST_JX into KPP PHOTOL array
!           !      PHOTOL(i_phot) = State_Chm%Phot%ZPJ(i_lev,i_phot,i_lon,i_lat)
!                 ! In GEOS-Chem, maybe useful to archieve instantaneous photolysis rate [s -1] and noon time photolysis rate [s -1]
!                 ! The mapping between the GEOS-Chem photolysis species and
!                 ! the FAST-JX photolysis species is contained in the lookup
!                 ! table in input file FJX_j2j.dat.
!                 ! Some GEOS-Chem photolysis species may have multiple
!                 ! branches for photolysis reactions.  These will be
!                 ! represented by multiple entries in the FJX_j2j.dat
!                 ! lookup table.
!                 !    NOTE: For convenience, we have stored the GEOS-Chem
!                 !    photolysis species index (range: 1..State_Chm%nPhotol)
!                 !    for each of the FAST-JX photolysis species (range;
!                 !    1..State_Chm%Phot%nMaxPhotRxns) in the GC_PHOTO_ID array

!                 ! GC photolysis species index
!             !    P = State_Chm%Phot%GC_Photo_Id(i_phot)
!                 ! Below is diagnostic steps and might be ignored for now
!                 ! If this FAST_JX photolysis species maps to a valid
!                 ! GEOS-Chem photolysis species (for this simulation)...
!                 !IF ( P > 0 .and. P <= State_Chm%nPhotol ) THEN
!                   !
!                 !ELSE IF ( P == State_Chm%nPhotol+1 ) THEN
!                   ! J(O3_O1D).  This used to be stored as the nPhotol+1st
!                   ! diagnostic in Jval, but needed to be broken off
!                   ! to facilitate cleaner diagnostic indexing (bmy, 6/3/20)
!                 !ELSE IF ( P == State_Chm%nPhotol+2 ) THEN
!                   ! J(O3_O3P).  This used to be stored as the nPhotol+2nd
!                   ! diagnostic in Jval, but needed to be broken off
!                   ! to facilitate cleaner diagnostic indexing (bmy, 6/3/20)
!                 !ENDIF

!              ! ENDDO
!             !ENDIF
!           !ENDIF
          
!           ! Initialize the KPP "C" vector of species concentrations [molec/cm3]
!           DO i_species = 1, NSPEC
!             SpcID = State_Chm%Map_KppSpc(i_species)
!             C(i_species)  = 0.0_dp
!             IF ( SpcId > 0 ) C(i_species) =  Plume2d_curr%CONCNT2d(i_x,i_y,i_species)
!           ENDDO

!           !=====================================================================
!           ! CHEMISTRY MECHANISM INITIALIZATION (#1)
!           !
!           ! Populate KPP global variables and arrays in gckpp_global.F90
!           !
!           ! NOTE: This has to be done before Set_Sulfur_Chem_Rates, so that
!           ! the NUMDEN and SR_TEMP KPP variables will be populated first.
!           ! Otherwise this can lead to differences in output that are evident
!           ! when running with different numbers of OpenMP cores.
!           ! See https://github.com/geoschem/geos-chem/issues/1157
!           !    -- Bob Yantosca (08 Mar 2022)
!           !=====================================================================
!           ! Copy values into the various KPP global variables
!           CALL Set_inPlume_2d_Kpp_GridBox_Values( I_EU          = i_lon,                          &
!                                         J_EU          = i_lat,                          &
!                                         L_EU         = i_lev,                          &
!                                         I_LA         =   i_x,                     &
!                                         J_LA        = i_y,                       &
!                                         Input_Opt  = Input_Opt,                  &
!                                         State_Chm  = State_Chm,                  &
!                                         State_Grid = State_Grid,                 &
!                                         State_Met  = State_Met,                  &
!                                         RC         = RC                         )
!           !=====================================================================
!           ! CHEMISTRY MECHANISM INITIALIZATION (#2)
!           !
!           ! Update reaction rates [1/s] for sulfur chemistry in cloud and on
!           ! seasalt.  These will be passed to the KPP chemical solver.
!           !
!           ! NOTE: This has to be done before fullchem_SetStateHet so that
!           ! State_Chm%HSO3_aq and State_Chm%SO3_aq will be populated first.
!           ! These are copied into State_Het%HSO3_aq and State_Het%SO3_aq.
!           ! See https://github.com/geoschem/geos-chem/issues/1157
!           !    -- Bob Yantosca (08 Mar 2022)
!           !=====================================================================
!           ! Compute sulfur chemistry reaction rates [1/s]
!           ! If size_res = T, we'll call fullchem_HetDropChem below.
!           ! Could remove the use of State_Diag
!           ! defines the variables State_Chm%HSO3_aq and State_Chm%SO3aq
!           ! Therefore, we must call Set_Sulfur_Chem_Rates after
!           !  Set_KPP_GridBox_Values, but before fullchem_SetStateHet.  Otherwise we
!           !  will not be able to copy State_Chm%HSO3_aq to State_Het%HSO3_aq and
!           !  State_Chm%SO3_aq to State_Het%SO3_aq properly.
!           !CALL Set_Sulfur_Chem_Rates( I          = i_lon,                           &
!           !                            J          = i_lat,                           &
!           !                            L          = i_lev,                           &
!           !                            Input_Opt  = Input_Opt,                   &
!           !                            State_Chm  = State_Chm,                   &
!           !                            State_Diag = State_Diag,                  &
!           !                            State_Grid = State_Grid,                  &
!           !                            State_Met  = State_Met,                   &
!           !                            size_res   = size_res,                    &
!           !                            RC         = RC                          )

!           !=====================================================================
!           ! CHEMISTRY MECHANISM INITIALIZATION (#3)
!           !
!           ! Populate the various fields of the State_Het object.
!           !
!           ! NOTE: This has to be done after fullchem_SetStateHet so that
!           ! State_Chm%HSO3_aq and State_Chm%SO3_aq will be populated first.
!           ! These are copied into State_Het%HSO3_aq and State_Het%SO3_aq.
!           ! See https://github.com/geoschem/geos-chem/issues/1157
!           !    -- Bob Yantosca (08 Mar 2022)
!           !=====================================================================

!           ! Populate fields of the State_Het object
!           ! These values are used in the heterogenous chemistry reaction rate computations.
!           ! Ignore heteogeneous reaction for now and complete later
!           !CALL fullchem_SetStateHet( I         = I,                             &
!           !                            J         = J,                             &
!           !                            L         = L,                             &
!           !                            id_SALA   = id_SALA,                       &
!           !                            id_SALAAL = id_SALAAL,                     &
!           !                            id_SALC   = id_SALC,                       &
!           !                            id_SALCAL = id_SALCAL,                     &
!           !                            Input_Opt = Input_Opt,                     &
!           !                            State_Chm = State_Chm,                     &
!           !                            State_Met = State_Met,                     &
!           !                            H         = State_Het,                     &
!           !                            RC        = RC                            )

!           !=====================================================================
!           ! (ignore for now) CHEMISTRY MECHANISM INITIALIZATION (#5)
!           !
!           ! Call Het_Drop_Chem (formerly located in sulfate_mod.F90) to
!           ! estimate the in-cloud sulfate production rate in heterogeneous
!           ! cloud droplets based on the Yuen et al., 1996 parameterization.
!           ! Code by Becky Alexander (2011) with updates by Mike Long and Bob
!           ! Yantosca (2021).
!           !
!           ! We will only call Het_Drop_Chem if:
!           ! (1) It is at least 0.01% cloudy
!           ! (2) We are doing a size-resolved computation
!           ! (3) The grid box is over water
!           ! (4) The temperature is above -5C
!           ! (5) Liquid water content is nonzero
!           !=====================================================================

!           !=====================================================================
!           ! Prepare arrays
!           !=====================================================================

!           ! Zero out dummy species index in KPP
!           !! Since we did not initialize PL_Kpp_Id, need to define dummy species from another way

!           DO i_kpp = 1, NFAM
!               KppID = Ind_( TRIM ( Fam_Names(i_kpp) ), 'K' )
!               ! Exit if an invalid ID is encountered
!               !IF ( KppId <= 0 ) THEN
!               !    ErrMsg = 'Invalid KPP ID for prod/loss species: '            // &
!               !         TRIM( Fam_Names(i_kpp) )
!               !    CALL GC_Error( ErrMsg, RC, ThisLoc )
!               !    RETURN
!               !ENDIF
!               IF ( KppID > 0 ) C(KppID) = 0.0_dp
!           ENDDO

!           !=====================================================================
!           ! Update reaction rates
!           !=====================================================================

!           ! Update the array of rate constants
!           ! Mannually Set K_CLoud and K_MT = 0, 
!           ! Set REACTION RELATED TO SEA SALT AND IODONE =0
!           ! Read from last timestep
!           ! CALL Update_RCONST()
!             RCONST = RXNRATE_CONST_KPP (i_lon, i_lat, i_lev, :)
!           !=====================================================================
!           ! HISTORY (aka netCDF diagnostics)
!           !
!           ! Archive KPP reaction rates [molec cm-3 s-1]
!           ! See gckpp_Monitor.F90 for a list of chemical reactions
!           !=====================================================================

! #ifdef KPP_INTEGRATOR_AUTOREDUCE
!           !=====================================================================
!           ! This is currently an empty function, might be discarded?
!           ! Set options for the KPP integrator in vectors ICNTRL and RCNTRL
!           ! This now needs to be done within the parallel loop
!           !=====================================================================
!           CALL fullchem_AR_SetIntegratorOptions( Input_Opt, State_Chm,          &
!                                                  State_Met, FirstChem,          &
!                                                  i_lon,    i_lat,  i_lev,       &
!                                                  ICNTRL,    RCNTRL             )
!           ! BZ, this needs to be modified from KPP/fullchem/fullchem_AutoReduceFuncs.F90
!           ! Initialize Hstart (the starting value of the integration step
!           ! size with the value of Hnew (the last predicted but not yet 
!           ! taken timestep)  saved to the the restart file.
!           RCNTRL(3) = State_Chm%KPPHvalue(i_lon,i_lat,i_lev)
!           !---------------------------------------------------------------------
!           ! Auto-reduce threshold, Method 1: Pressure-dependent
!           !                                            
!           !   Actual_Threshold =
!           !                                           Mid-Pressure at Level
!           !     AUTOREDUCE_THRESHOLD (at surface) * --------------------------
!           !                                          "Mid-Pressure" at Sfc.
!           !
!           !---------------------------------------------------------------------
!           IF ( .not. Input_Opt%AUTOREDUCE_IS_KEY_THRESHOLD ) THEN
!             IF ( Input_Opt%AUTOREDUCE_IS_PRS_THRESHOLD ) THEN
!                 RCNTRL(12) = Input_Opt%AUTOREDUCE_THRESHOLD                        & 
!                           * State_Met%PMID(i_lon,i_lat,i_lev)                                 & 
!                           / State_Met%PMID(i_lon,i_lat,1)
!             ENDIF
            
!             IF ( .not. Input_Opt%AUTOREDUCE_IS_PRS_THRESHOLD ) THEN
!                 RCNTRL(12) = Input_Opt%AUTOREDUCE_THRESHOLD
!             ENDIF
!           ENDIF

!           !---------------------------------------------------------------------
!           ! Auto-reduce threshold, Method 2: Determine threshold 
!           ! dynamically by scaling rates of key species.
!           !---------------------------------------------------------------------
!           IF ( Input_Opt%AUTOREDUCE_IS_KEY_THRESHOLD ) THEN

!             !--------------------------------
!             ! Daytime target species (OH)
!             !--------------------------------
!             ICNTRL(14) = ind_OH
!             RCNTRL(14) = Input_Opt%AUTOREDUCE_TUNING_OH
            
!             !--------------------------------
!             ! Nighttime target species (NO2)
!             !--------------------------------
!             ! COMMENTS BY HAIPENG LIN:
!             ! 1e6 daytime conc...testing shows 5e-5 as an offset here works best.
!             ! Use JNO2 as night determination.
!             ! RXN_NO2: NO2 + hv --> NO  + O
!             ! JNO2 ranges from 0 to 0.02 and is order ~ 1e-4 at the terminator. 
!             ! We set this threshold to be slightly relaxed so it captures the 
!             ! terminator, but this needs some tweaking.
!             !
!             ! For some reason, RXN_NO2 as a proxy fails to propagate the sunset 
!             ! terminator even though all diagnostics seem fine, and after a while 
!             ! only the OH scheme applies.  Use SUNCOSmid as a proxy to fix this. 
!             ! (hplin, 4/20/22)
!             ! IF(ZPJ(L,RXN_NO2,I,J) .eq. 0.0_fp) THEN
!             !
!             IF( State_Met%SUNCOSmid(i_lon, i_lat) .le. -0.1391731e+0_dp ) THEN
!                 ICNTRL(14) = ind_NO2
!                 RCNTRL(14) = Input_Opt%AUTOREDUCE_TUNING_NO2
!             ENDIF
!           ENDIF
! #endif
!           ! ICNTRL(16) option
!           ! 0 -> do nothing.
!           ! 1 -> set negative values to zero
!           ! 2 -> return with error code
!           ! 3 -> stop at negative
!           ICNTRL(16) = 1
!           ICNTRL(15) =  -1 ! ICNTRL(15) = -1 ! Do not call Update_* functions within the integrator
!           !=====================================================================
!           ! Integrate the box forwards
!           !=====================================================================
!           ! BZ manually set H2, O2, N2 equals to background
!           ! C(345) = 
!           C_before_integrate = C
!           CALL Integrate( 0.0_dp, DT, ICNTRL, RCNTRL, ISTATUS, RSTATE, IERR )
!           ! Print grid box indices to screen if integrate failed
!           IF ( IERR < 0 ) THEN

!               ! Turn off error output after a certain limit is reached
!               IF ( .not. doSuppress ) THEN
!                 WRITE( 6, * ) '### INTEGRATE RETURNED ERROR AT: ', i_x, i_y, Plume2d_curr%label
!                 errorCount = errorCount + 1
!                 IF ( errorCount > INTEGRATE_FAIL_TOGGLE ) THEN
!                     WRITE( 6, '(a)' ) &
!                       '### Further error output has been switched off'
!                     doSuppress = .TRUE.
!                 ENDIF
!               ENDIF
!             ENDIF
!           !=====================================================================
!           ! HISTORY: Archive KPP solver diagnostics
!           !=====================================================================
          
!           !=====================================================================
!           ! Try another time if it failed
!           !=====================================================================
!           IF ( IERR < 0 ) THEN
!             ! Zero the first time step (Hstart).  Also reset C with
!             ! concentrations prior to the 1st call to "Integrate".
!             RCNTRL(3) = 0.0_dp
!             C         = C_before_integrate

!             ! Disable auto-reduce solver for the second iteration for safety
!             IF ( Input_Opt%Use_AutoReduce ) THEN
!               RCNTRL(12) = -1.0_dp ! without using ICNTRL
!             ENDIF

!             ! Update rates again
!             ! CALL Update_RCONST()

!             ! Call the integrator
!             ! NOTE: Some integrators (like LSODE) will overwrite the TIN value
!             ! upon exit.  To prevent this, pass 0.0, DT as the 1st 2 arguments.
!             CALL Integrate( 0.0_dp, DT, ICNTRL, RCNTRL, ISTATUS, RSTATE, IERR )
            
!             !==================================================================
!             ! Exit upon the second failure
!             !==================================================================
!             IF ( IERR < 0 ) THEN
!               ! Print error message
!               WRITE(6,     '(a   )' ) '## INTEGRATE FAILED TWICE !!! '
!               WRITE(ERRMSG,'(a,i3)' ) 'Integrator error code :', IERR
             
!              !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!              ! Make sure only one thread at a time executes this block
!              !$OMP CRITICAL
!              !
!              ! Set a flag to break out of loop gracefully
!              ! NOTE: You can set a GDB breakpoint here to examine the error
!              ! BZ: Temperoray let the code continue  
!              ! Failed2x = .TRUE.

!              ! Print concentrations at failure grid box
!              PRINT*, REPEAT( '#', 79 )
!              PRINT*, '### KPP DEBUG OUTPUT!'
!              PRINT*, '### Species concentrations at problem box ',i_x, i_y, Plume2d_curr%label
!              PRINT*, REPEAT( '#', 79 )
!              DO i_species = 1, NSPEC
!                 IF (C(i_species) .lt. 0.0_dp) THEN
!                   C(i_species) = C_before_integrate(i_species)
!                   Write (6, *) "Revert Species num: ", i_species
!                 ENDIF
!              ENDDO
!              ! DO i_species = 1, n_species
!              !   PRINT*, C(i_species), TRIM( ADJUSTL( SPC_NAMES(i_species) ) )
!              !ENDDO

!              ! Print rate constants at failure grid box
!              PRINT*, REPEAT( '#', 79 )
!              PRINT*, '### KPP DEBUG OUTPUT!'
!              PRINT*, '### Reaction rates at problem box ', i_x, i_y, Plume2d_curr%label
!              PRINT*, REPEAT( '#', 79 )
!              ! DO i_rxn = 1, NREACT
!               !  PRINT*, RCONST(i_rxn), TRIM( ADJUSTL( EQN_NAMES(i_rxn) ) )
!              ! ENDDO
!              !
!              !$OMP END CRITICAL
!              !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!              ! Start skipping to end of loop upon 2 failures in a row
!              CYCLE
!             ENDIF
!           ENDIF

!           !=====================================================================
!           ! Check we have no negative values and copy the concentrations
!           ! calculated from the C array back into Plume conc array
!           !=====================================================================

!           ! Loop over KPP species
!           DO i_species = 1, NSPEC

!               ! GEOS-Chem species ID
!               SpcID = State_Chm%Map_KppSpc(i_species)

!               ! Skip if this is not a GEOS-Chem species
!               IF ( SpcID <= 0 ) CYCLE


!               ! Set negative concentrations to zero
!               C(i_species) = MAX( C(i_species), 0.0_dp )

!               ! Copy concentrations back into State_Chm%Species
!                Plume2d_curr%CONCNT2d(i_x,i_y,i_species) = REAL( C(i_species), kind=fp )

!           ENDDO
! #ifdef TOMAS
!               !-----------------------------------------------------------------
!               ! FOR TOMAS MICROPHYSICS:
!               !
!               ! Obtain P/L with a unit [kg S] for tracing
!               ! gas-phase sulfur species production (SO2, SO4, MSA)
!               ! (win, 8/4/09)
!               !
!               ! TODO: Abstract this to a subroutine, to simplify DO_FULLCHEM
!               !-----------------------------------------------------------------
!               ! Calculate H2SO4 production rate [kg s-1] in each
!               ! time step (win, 8/4/09)
!               H2SO4_RATE_2d(i_x, i_y)= C(ind_PH2SO4) / AVO * 98.e-3_fp * &
!                            State_Met%AIRVOL(i_lon,i_lat,i_lev)    * &
!                            1.0e+6_fp / DT  ! kg s-1 box-1
        
!               IF ( H2SO4_RATE_2d(i_x, i_y) < 0.0d0) THEN
!                 !ErrMsg = "H2SO4_RATE_2D negative in (Plumeid, x, y):", &
!                 !    Plume2d_curr%label, i_x, i_y, "was:", H2SO4_RATE_2d(i_x, i_y), "  setting to 0.0d0"
!                 !CALL GC_Warning( ErrMsg, RC, ThisLoc )
!                 WRITE(6, *) "H2SO4_RATE_2D negative in (Plumeid, x, y):", &
!                     Plume2d_curr%label, i_x, i_y, "was:", H2SO4_RATE_2d(i_x, i_y), "  setting to 0.0d0"
!                 H2SO4_RATE_2d(i_x, i_y) = 0.0d0
!               ENDIF

!               PSO4AQ_RATE_2d(i_x, i_y) = C(ind_PSO4AQ) / AVO * 98.e-3_fp * &
!                             State_Met%AIRVOL(i_lon,i_lat,i_lev)    * &
!                             1.0e+6_fp ! kg per timestep box-1

!               IF ( PSO4AQ_RATE_2d(i_x, i_y) < 0.0d0) THEN
!                 !ErrMsg = "PSO4AQ_RATE_2D negative in (Plumeid, x, y):", &
!                 !    Plume2d_curr%label, i_x, i_y, "was:", PSO4AQ_RATE_2d(i_x, i_y), "  setting to 0.0d0"
                
!                 !CALL GC_Warning( ErrMsg, RC, ThisLoc )
!                 WRITE(6, *) "PSO4AQ_RATE_2D negative in (Plumeid, x, y):", &
!                       Plume2d_curr%label, i_x, i_y, "was:", PSO4AQ_RATE_2d(i_x, i_y), "  setting to 0.0d0"
!                 PSO4AQ_RATE_2d(i_x, i_y) = 0.0d0
!               ENDIF
! #endif
!               !====================================================================
!               ! Archieve KPP diagnostic output and write diagnostics files
!               ! e.g., Chemical production and loss
!               ! OH reactivity: inverse of its life-time
!               ! GcKpp_Util -> Get_OHreactivity
!               !====================================================================
!         ENDDO
!       ENDDO
!       !$OMP END PARALLEL DO
      
      
      !=======================================================================
      ! Return gracefully if integration failed 2x anywhere
      ! (as we cannot break out of a parallel DO loop!)
      !=======================================================================
      !IF ( Failed2x ) THEN
      !  ErrMsg = 'KPP failed to converge after 2 iterations!'
      !  CALL GC_Error( ErrMsg, RC, ThisLoc )
      !  RETURN
      !ENDIF

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

      !Plume2d_curr%CONCNT2d = box_concnt_2D
      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc after Chem: Conc of SO2 ', SUM(Plume2d_curr%CONCNT2d(:,:,id_SO2_p))/(n_x_max*n_y_max)
      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc after Chem: Conc of SO4 ', SUM(Plume2d_curr%CONCNT2d(:,:,id_SO4_p))/(n_x_max*n_y_max)
      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' ave conc after Chem: Conc of OH ', SUM(Plume2d_curr%CONCNT2d(:,:,id_OH_p))/(n_x_max*n_y_max)

      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc after Chem: Conc of SO2 ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_SO2_p)
      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc after Chem: Conc of SO4 ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_SO4_p)
      !WRITE(6,*) 'Debug (BZ): Plume box 2d : ',Plume2d_curr%label,' center conc after Chem: Conc of OH ', Plume2d_curr%CONCNT2d(n_x_mid,n_y_mid,id_OH_p)
      

    !Write (6, *) "Debug: (BZ): In Plume  (After Plume Chem): rate constant for SO2_OH_RXN_ID = ", &
    !    State_Diag%RxnConst(23, 40, 39 ,SO2_OH_RXN_ID)
400 CONTINUE
    ! Exchange OH and HO2 with background when 1) in chemical grid 2) if turn on tropp sink, only do exchange in stratosphere
   ! Write (6, *) "Debug: (BZ): Plume chem before exchange: SUM(Spc(id_OH)%Conc) =  ", SUM(Spc(id_OH)%Conc), &
   !             "SUM(Spc(id_HO2)%Conc) = ", SUM(Spc(id_HO2)%Conc), &
   !             "SUM(Spc(id_NH3)%Conc) = ", SUM(Spc(id_NH3)%Conc), "SUM(Spc(id_NH4)%Conc) = ", SUM(Spc(id_NH4)%Conc)
    !    State_Diag%RxnConst(23, 40, 39 ,SO2_OH_RXN_ID)
    DO i_lev = 1, NZ_GC
      DO i_lat = 1, NY_GC
         DO i_lon = 1, NX_GC
            IF (State_Met%InChemGrid(i_lon,i_lat,i_lev).AND.  &
                        ( (.NOT. TROPP_sink) .OR. (State_Met%InStratosphere(i_lon,i_lat,i_lev)) ) ) THEN

               Vgrid_EU        =  State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ! cm3
               Spc(id_OH)%Conc(i_lon,i_lat,i_lev) = Spc(id_OH)%Conc(i_lon,i_lat,i_lev) - &
                     (mass_OH_consum_plume(i_lon, i_lat, i_lev) / Vgrid_EU)
               Spc(id_HO2)%Conc(i_lon,i_lat,i_lev) = Spc(id_HO2)%Conc(i_lon,i_lat,i_lev) - &
                     (mass_HO2_consum_plume(i_lon, i_lat, i_lev) / Vgrid_EU)
#ifdef TOMAS
               Spc(id_NH3)%Conc(i_lon,i_lat,i_lev) = Spc(id_NH3)%Conc(i_lon,i_lat,i_lev) - &
                     (mass_NH3_consum_plume(i_lon, i_lat, i_lev) / Vgrid_EU)
               Spc(id_NH4)%Conc(i_lon,i_lat,i_lev) = Spc(id_NH4)%Conc(i_lon,i_lat,i_lev) - &
                     (mass_NH4_consum_plume(i_lon, i_lat, i_lev) / Vgrid_EU)
#endif
            ENDIF
         ENDDO
      ENDDO
   ENDDO



   ! Write (6, *) "Debug: (BZ): Plume chem after exchange: SUM(Spc(id_OH)%Conc) =  ", SUM(Spc(id_OH)%Conc), &
   !             "SUM(Spc(id_HO2)%Conc) = ", SUM(Spc(id_HO2)%Conc), &
   !             "SUM(Spc(id_NH3)%Conc) = ", SUM(Spc(id_NH3)%Conc), "SUM(Spc(id_NH4)%Conc) = ", SUM(Spc(id_NH4)%Conc)

    ! deallocate unused space
    IF(allocated(box_concnt_2D)) deallocate(box_concnt_2D)
    IF(allocated(box_concnt_2D_prev)) deallocate(box_concnt_2D_prev)
    !IF(allocated(box_concnt_2D_kg)) deallocate(box_concnt_2D_kg)
    IF(allocated(debug_ix)) deallocate(debug_ix)
    IF(allocated(debug_iy)) deallocate(debug_iy)
    IF(allocated(debug_ibox)) deallocate(debug_ibox)
    IF(allocated(debug_islab)) deallocate(debug_islab)
    IF(allocated(debug_status)) deallocate(debug_status)
    IF(allocated(debug_value)) deallocate(debug_value)
    IF(ALLOCATED(mass_OH_consum_plume)) DEALLOCATE(mass_OH_consum_plume)
    IF(ALLOCATED(mass_HO2_consum_plume)) DEALLOCATE(mass_HO2_consum_plume)
    IF(ALLOCATED(mass_NH3_consum_plume)) DEALLOCATE(mass_NH3_consum_plume)
    IF(ALLOCATED(mass_NH4_consum_plume)) DEALLOCATE(mass_NH4_consum_plume)
    
    IF(ASSOCIATED(Spc)) nullify(Spc)

    IF(allocated(box_concnt_1D)) deallocate(box_concnt_1D)
    IF(allocated(box_concnt_1D_prev)) deallocate(box_concnt_1D_prev)

    IF(ASSOCIATED(Plume2d_next)) nullify(Plume2d_next)
    IF(ASSOCIATED(Plume2d_curr)) nullify(Plume2d_curr)
    IF(ASSOCIATED(Plume2d_prev)) nullify(Plume2d_prev)
    IF(ASSOCIATED(Plume1d_next)) nullify(Plume1d_next)
    IF(ASSOCIATED(Plume1d_curr)) nullify(Plume1d_curr)
    IF(ASSOCIATED(Plume1d_prev)) nullify(Plume1d_prev)
  END SUBROUTINE plume_chem_microphysics

  SUBROUTINE plume_structure_change(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
    ! All the use have been defined in host: Plume_box_model
    USE Input_Opt_Mod,   ONLY : OptInput, PlumeSource_t
    USE State_Chm_Mod,   ONLY : ChmState, Ind_
    USE State_Met_Mod,   ONLY : MetState
    USE Species_Mod,     ONLY : SpcConc
    USE TIME_MOD,        ONLY : GET_TS_DYN
    USE TIME_MOD,        ONLY : ITS_TIME_FOR_EXIT
    USE State_Grid_Mod,  ONLY : GrdState
    ! For creating diagnostic files
    USE InquireMod,      ONLY : findFreeLun
    !USE UnitConv_Mod

  !    USE GC_GRID_MOD,   ONLY : XEDGE, YEDGE
  !    USE CMN_SIZE_Mod,  ONLY : DLAT, DLON !new

    LOGICAL,        INTENT(IN)    :: am_I_Root
    TYPE(MetState), INTENT(IN)    :: State_Met
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State objectgg
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    TYPE(SpcConc), POINTER        :: Spc(:)

    INTEGER                       :: i_box, i_lon, i_lat, i_lev, ibin
    INTEGER                       :: i_species, n_species, i_tracer
    INTEGER                       :: ind_spc_GC, ind_spc_p, ind_spc_GC_bin1
    INTEGER                       :: Stop_loop
    INTEGER                       :: N_core

    LOGICAL                       :: exe_exit

    REAL(fp)                      :: box_lon, box_lat, box_lev
    REAL(fp)                      :: box_length, box_alpha
    REAL(fp)                      :: box_extra, box_life, box_label
   ! REAL(fp)                      :: box_RA(nspc_p), box_Rb(nspc_p), box_theta(nspc_p)
    REAL(fp)                      :: box_Ra, box_Rb, box_theta
    REAL(fp)                      :: Pdx, Pdy
    REAL(fp)                      :: Vgrid_EU, Vgrid_2D, Vgrid_1D
    REAL(fp)                      :: Vgrid_1D_temp
    REAL(fp)                      :: conc_background, mass_release, background_mass
    REAL(fp)                      :: Dt
    REAL(fp)                      :: ConcSlab(n_slab_max)
    REAL(fp)                      :: mass_plume_1D, mass_plume_2D, mass_plume_diff, mass_plume_bg
    REAL(fp)                      :: Xscale, Yscale
    REAL(fp)                      :: Xslab(128, n_slab_max),  Yslab(128, n_slab_max) ! intermediate product of interpolating from 2-D to 1-D, size 128 is defined in subroutine: Define_Slab_Grid
    REAL(fp)                      :: C99, Core_Threshold
    
    CHARACTER(LEN=255)            :: spc_name
    CHARACTER(LEN=255)            :: ErrMsg
    CHARACTER(LEN=255)            :: ThisLoc

    REAL(fp), DIMENSION(:,:,:), ALLOCATABLE :: box_concnt_2D
    REAL(fp), DIMENSION(:,:), ALLOCATABLE   :: box_concnt_1D

    TYPE(Plume2d_list), POINTER :: Plume2d_next, Plume2d_curr, Plume2d_prev
    TYPE(Plume1d_list), POINTER :: Plume1d_next, Plume1d_curr, Plume1d_prev, Plume1d_new

    ! Diagnostic file variable
    INTEGER                :: i_x, i_y, i_slab
    
    INTEGER                :: file_2Dmass_NK01_ID, file_1Dmass_NK01_ID
    CHARACTER(LEN=255)     :: file_2Dmass_NK01,  file_1Dmass_NK01



    Spc                    =>  State_Chm%Species
    n_species              =   State_Chm%nSpecies
    Dt                     =   GET_TS_DYN()
    exe_exit               =   ITS_TIME_FOR_EXIT()
    ThisLoc                =   ' -> at plume_structure_change (in module GeosCore/lagrange_singlebox_mod.F90)'
    RC                     =   GC_SUCCESS
    ErrMsg                 =   ''
    
    Spc(id_SO2pl)%Conc(:,:,:)    =   0.0_fp
    Spc(id_SO4pl)%Conc(:,:,:)    =   0.0_fp
   !  mass_S_SO2_r3_2D       =   0.0_fp
   !  mass_S_SO4_r3_2D       =   0.0_fp
   !  mass_S_SO2_r3_1D       =   0.0_fp
   !  mass_S_SO4_r3_1D       =   0.0_fp
   !  mass_S_SO2_r4          =   0.0_fp
   !  mass_S_SO4_r4          =   0.0_fp
   !  mass_S_SO2_6_2D        =   0.0_fp
   !  mass_S_SO4_6_2D        =   0.0_fp
   !  mass_S_SO2_6_1D        =   0.0_fp
   !  mass_S_SO4_6_1D        =   0.0_fp
   !  mass_S_SO2_7_2D        =   0.0_fp
   !  mass_S_SO4_7_2D        =   0.0_fp
   !  mass_S_SO2_7_1D        =   0.0_fp
   !  mass_S_SO4_7_1D        =   0.0_fp
   !  Vgrid_2D_tot_3         =   0.0_fp
   !  Vgrid_1D_tot_3         =   0.0_fp
    
! #ifdef TOMAS
!     mass_SF_bin_2D(:)     = 0.0_fp
!     mass_NK_bin_2D(:)     = 0.0_fp
!     mass_MK_bin_2D(:)     = 0.0_fp
!     mass_SF_bin_1D(:)     = 0.0_fp
!     mass_NK_bin_1D(:)     = 0.0_fp
!     mass_MK_bin_1D(:)     = 0.0_fp
!     mass_SF_bin_2D_1D(:)  = 0.0_fp
!     mass_NK_bin_2D_1D(:)  = 0.0_fp
! #endif

    NULLIFY(Plume2d_next, Plume2d_curr, Plume2d_prev)
    NULLIFY(Plume1d_next, Plume1d_curr, Plume1d_prev, Plume1d_new)
    ! Debug for dissolving plume at exit time
    IF (ITS_TIME_FOR_EXIT()) THEN
     Write (6, *) "Debug (BZ): time for exit detected in plume model"
    ENDIF
    !ALLOCATE(box_concnt_2D(n_x_max, n_y_max, n_species))
    !---------------------------------------------------------------------
    ! convert 2D seg to 1D seg if:
    ! 1) Plume2d%IsTransfer=1, or;
    ! 2) Lifetime larger than a critical value

    IF(.NOT.ASSOCIATED(Plume2d_head)) GOTO 401
    Plume2d_curr => Plume2d_head

    DO WHILE(ASSOCIATED(Plume2d_curr))

      box_lon       = Plume2d_curr%LON
      box_lat       = Plume2d_curr%LAT
      box_lev       = Plume2d_curr%LEV
      box_length    = Plume2d_curr%LENGTH
      box_alpha     = Plume2d_curr%ALPHA
      box_label     = Plume2d_curr%label
      box_life      = Plume2d_curr%LIFE
      Pdx           = Plume2d_curr%PDX
      Pdy           = Plume2d_curr%PDY

      i_lon         = Plume2d_curr%lon_ind
      i_lat         = Plume2d_curr%lat_ind
      i_lev         = Plume2d_curr%lev_ind

      Vgrid_EU     = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ! [cm3]
      Vgrid_2D       = Pdx*Pdy*box_length*1.0e+6_fp

      ! mass_S_SO2_6_2D = mass_S_SO2_6_2D + SUM(Plume2d_curr%CONCNT2D(:,:, id_SO2_p)) * Vgrid_2D
      ! mass_S_SO4_6_2D = mass_S_SO4_6_2D + SUM(Plume2d_curr%CONCNT2D(:,:, id_SO4_p)) * Vgrid_2D

      Plume2d_next => Plume2d_curr%next
      WRITE(File_Plume_location_IU_2D,'(5(I0,1X),6(F10.3, 1X))') &
            NINT(time_elapsed), Plume2d_curr%LABEL, &
            Plume2d_curr%lon_ind, Plume2d_curr%lat_ind, Plume2d_curr%lev_ind, &
            Plume2d_curr%LON, Plume2d_curr%LAT, Plume2d_curr%LEV, &
            Plume2d_curr%PDX, Plume2d_curr%PDY,  Plume2d_curr%LENGTH
      ! Below print 2-D conc matrix in each plume, turn off to avoid massive output files
      !file_2Dconc_SO4_ID = findFreeLun()
      !WRITE(file_2Dconc_SO4,'("Plume-2D_SO4_conc_",I0,".txt")') NINT(time_elapsed)
      !CALL PLUME_CONC_DIAG_FILES_2D(file_2Dconc_SO4_ID, file_2Dconc_SO4, Plume2d_curr%CONCNT2d(:,:,id_SO4_p), RC)
      !WRITE(6, *) "(Debug: BZ), 2-D SO4 at test grid =", Plume2d_curr%CONCNT2d(x_test,y_test, id_SO4_p)

      !file_2Dconc_SO2_ID = findFreeLun()
      !WRITE(file_2Dconc_SO2,'("Plume-2D_SO2_conc_",I0,".txt")') NINT(time_elapsed)
      !CALL PLUME_CONC_DIAG_FILES_2D(file_2Dconc_SO2_ID, file_2Dconc_SO2, Plume2d_curr%CONCNT2d(:,:,id_SO2_p), RC)

      !file_2Dconc_OH_ID = findFreeLun()
      !WRITE(file_2Dconc_OH,'("Plume-2D_OH_conc_",I0,".txt")') NINT(time_elapsed)
      !CALL PLUME_CONC_DIAG_FILES_2D(file_2Dconc_OH_ID, file_2Dconc_OH, Plume2d_curr%CONCNT2d(:,:,id_OH_p), RC)
      ! DO ibin = 1, nbins
      !   ind_spc_p = 10+1+(ibin-1)*nspc_p_tomas_tracer
      !   file_2D_conc_NKbin_ID = findFreeLun()
      !   WRITE(file_2D_conc_NKbin,'("Plume-2D_NKbin_",I0,"_", I0,".txt")') ibin, NINT(time_elapsed)
      !   CALL PLUME_CONC_DIAG_FILES_2D(file_2D_conc_NKbin_ID, file_2D_conc_NKbin, Plume2d_curr%CONCNT2d(:,:,ind_spc_p), RC)
      !   IF (ibin .eq. 15) THEN
      !    WRITE(6, *) "(Debug: BZ), 2-D NK bin 15 at test grid =", Plume2d_curr%CONCNT2d(x_test,y_test, ind_spc_p)
      !   ENDIF

      !   ind_spc_p = 10+2+(ibin-1)*nspc_p_tomas_tracer
      !   file_2D_conc_SFbin_ID = findFreeLun()
      !   WRITE(file_2D_conc_SFbin,'("Plume-2D_SFbin_",I0,"_", I0,".txt")') ibin, NINT(time_elapsed)
      !   CALL PLUME_CONC_DIAG_FILES_2D(file_2D_conc_SFbin_ID, file_2D_conc_SFbin, Plume2d_curr%CONCNT2d(:,:,ind_spc_p), RC)
      
      ! ENDDO
      ! --------------------------------------------------------------------
      ! If plume2d_curr%IsTransfer is true, create 1D segment
      ! If plume2d_curr%IsTransfer is true or If plume2d_curr%IsDissolve is true
      ! Dissolve the 2D seg
      ! --------------------------------------------------------------------
      IF(Plume2d_curr%IsTransfer) THEN
         Num_transfer_2D = Num_transfer_2D + 1
         WRITE(File_Plume_life_IU_2D,'(I0,2(1X,F10.1))') & 
            Plume2d_curr%LABEL, Plume2d_curr%LIFE, 0.0_fp
         ! file_2Dmass_NK01_ID = findFreeLun()
         ! WRITE(file_2Dmass_NK01,  &
         !    '("Plume-2D_NK01_mass_ID_",I0,"_time_", I0,".txt")') &
         !    Plume2d_curr%LABEL, NINT(time_elapsed)
         ! CALL PLUME_CONC_DIAG_FILES_2D(   &
         !    file_2Dmass_NK01_ID, file_2Dmass_NK01,   &
         !    Plume2d_curr%CONCNT2d(:,:, id_NK01_p)*Vgrid_2D, RC)
         ! Creating 1D list
         Num_Plume1d = Num_Plume1d + 1
         Num_Plume1d_acc = Num_Plume1d_acc +1
         ALLOCATE(Plume1d_new)
         Plume1d_new%LABEL = Num_Plume1d_acc
         Plume1d_new%LON = Plume2d_curr%LON
         Plume1d_new%LAT = Plume2d_curr%LAT
         Plume1d_new%LEV = Plume2d_curr%LEV
         Plume1d_new%lon_ind = Plume2d_curr%lon_ind
         Plume1d_new%lat_ind = Plume2d_curr%lat_ind
         Plume1d_new%lev_ind = Plume2d_curr%lev_ind
         Plume1d_new%IsDissolve = .False.
         Plume1d_new%LENGTH = Plume2d_curr%LENGTH
         Plume1d_new%ALPHA  = Plume2d_curr%ALPHA
         Plume1d_new%LIFE   = Plume2d_curr%LIFE
         ALLOCATE(Plume1d_new%CONCNT1d(n_slab_max,nspc_p))
         ! ALLOCATE(Plume1d_new%RA(nspc_p))
         ! ALLOCATE(Plume1d_new%RB(nspc_p))
         ! ALLOCATE(Plume1d_new%THETA(nspc_p))
         NULLIFY(Plume1d_new%next)

         ! Update plume size and shape based on SO2
         ! Xscale = Get_XYscale(Plume2d_curr%CONCNT2d(:,:,id_SO2_p), Pdx, Pdy, frac_mass, 2)
         ! Yscale = Get_XYscale(Plume2d_curr%CONCNT2d(:,:,id_SO2_p), Pdx, Pdy, frac_mass, 1)
         ! box_theta = ATAN( Xscale/Yscale )
         ! CALL Slab_init_bilinear(Pdx, Pdy, box_theta, Plume2d_curr%CONCNT2d(:,:,id_SO2_p), &
         !             Yscale, box_Ra, box_Rb, ConcSlab)

         ! Plume1d_new%RA = box_Ra
         ! Plume1d_new%RB = box_Rb
         ! Plume1d_new%THETA = box_theta
         ! Vgrid_1D = Plume1d_new%RA*Plume1d_new%RB*Plume1d_new%LENGTH*1.0e+6_fp

         ! mass_plume_2D = Vgrid_2D * SUM(Plume2d_curr%CONCNT2d(:,:,id_SO2_p))
         ! mass_plume_1D = Vgrid_1D * SUM(ConcSlab)
         ! mass_plume_diff = mass_plume_1D - mass_plume_2D
         ! Plume1d_new%CONCNT1d(:,id_SO2_p) = ConcSlab/mass_plume_1D*mass_plume_2D

         ! IF(abs(mass_plume_diff)/mass_plume_2D>0.01)THEN
         !       ! BZ need to skip PH2SO4, OH/HO2, and think abuout how to treat TOMAS tracers
         !       Write (6, *) "Debug (BZ): More than 1% mass change in plume from 2-D to 1-D at i_box = ", Plume2d_curr%label, &
         !             "Species: SO2; mass_plume_2D= ", mass_plume_2D, &
         !             'mass_plume_1D= ', mass_plume_1D
         ! ENDIF
         ! Update RA, RB, THETA, CONCNT1D
         ! Use size based on SO4 conc
         Xscale = Get_XYscale(Plume2d_curr%CONCNT2d(:,:,id_SO2_p), Pdx, Pdy, frac_mass, 2)
         Yscale = Get_XYscale(Plume2d_curr%CONCNT2d(:,:,id_SO2_p), Pdx, Pdy, frac_mass, 1)
         box_theta = ATAN( Xscale/Yscale )
         
         CALL Define_Slab_Grid(                              &
               Pdx, Pdy, box_theta, Yscale,                   &
               box_Ra, box_Rb, Xslab, Yslab                   &
         )
         Plume1d_new%THETA = box_theta
         Plume1d_new%RA = box_Ra
         Plume1d_new%RB = box_Rb
         Vgrid_1D = box_Ra*box_Rb*Plume1d_new%LENGTH*1.0e+6_fp
         !Write(6, * ) "Debug (BZ): Xscale= ", Xscale, " Yscale= ", Yscale, "Pdx= ", Pdx, "Pdy= ", Pdy, &
         !                  "Box_length= ",box_length, "Box_theta=", box_theta
         DO i_species = 1, nspc_p
            IF ((i_species == id_OH_p) .or. (i_species == id_HO2_p) .or.  &
               (i_species == id_NH3_p) .or. (i_species == id_NH4_p).or. (i_species == id_H2O_p) ) CYCLE
            ! Xscale = Get_XYscale(Plume2d_curr%CONCNT2d(:,:,i_species), Pdx, Pdy, frac_mass, 2)
            ! Yscale = Get_XYscale(Plume2d_curr%CONCNT2d(:,:,i_species), Pdx, Pdy, frac_mass, 1)
            ! box_theta(i_species) = ATAN( Xscale/Yscale )
            
            ! CALL Slab_init_bilinear(Pdx, Pdy, box_theta, Plume2d_curr%CONCNT2d(:,:,i_species), &
            !          Yscale, box_Ra, box_Rb, ConcSlab)
            ConcSlab(:) = 0.0_fp
            CALL Interpolate_To_Slab(Pdx, Pdy, Plume2d_curr%CONCNT2d(:,:,i_species), &
                                     box_Ra, box_Rb, Xslab, Yslab, ConcSlab )
            !Vgrid_1D_temp = box_Ra*box_Rb*Plume1d_new%LENGTH*1.0e+6_fp

            mass_plume_2D = Vgrid_2D * SUM(Plume2d_curr%CONCNT2d(:,:,i_species))
            mass_plume_1D = Vgrid_1D * SUM(ConcSlab)
            mass_plume_diff = mass_plume_1D - mass_plume_2D
            IF (i_species .eq. id_SO2_p) THEN
               mass_S_SO2_r4 = mass_S_SO2_r4 + mass_plume_diff
            ENDIF
            IF (i_species .eq. id_SO4_p) THEN
               mass_S_SO4_r4 = mass_S_SO4_r4 + mass_plume_diff
            ENDIF
            ! IF (MOD((i_species - 10), nspc_p_tomas_tracer) .eq. 1) THEN
            !    ! TOMAS tracer, NK
            !    ibin = FLOOR((i_species - 10.0_fp)/nspc_p_tomas_tracer) + 1
            !    mass_NK_bin_2D_1D(ibin) = mass_NK_bin_2D_1D(ibin) + mass_plume_diff
            ! ENDIF
            ! IF (MOD((i_species - 10), nspc_p_tomas_tracer) .eq. 2) THEN
            !    ! TOMAS tracer, SF
            !    ibin = FLOOR((i_species - 10.0_fp)/nspc_p_tomas_tracer) + 1
            !    mass_SF_bin_2D_1D(ibin) = mass_SF_bin_2D_1D(ibin) + mass_plume_diff
            ! ENDIF
            ! ConcSlab = ConcSlab/mass_plume_1D*mass_plume_2D
            IF(mass_plume_diff>0.0_fp)THEN
                  Write (6, *) "Debug (BZ): Mass enter the plume from 2-D to 1-D at i_box = ", Plume2d_curr%label, &
                        "Species: ", TRIM(spc_names_p_use(i_species)), " mass_plume_2D= ", mass_plume_2D, &
                        'mass_plume_1D= ', mass_plume_1D
                  errMsg = 'Mass enter the plume from 2-D to 1-D grid! '
                  CALL ERROR_STOP( errMsg, thisLoc)
            ELSE
               IF(abs(mass_plume_diff)/mass_plume_2D>0.01)THEN
                     Write (6, *) "Debug (BZ): More than 1% mass change in plume from 2-D to 1-D at i_box = ", Plume2d_curr%label, &
                           "Species: ", TRIM(spc_names_p_use(i_species)), " mass_plume_2D= ", mass_plume_2D, &
                           'mass_plume_1D= ', mass_plume_1D
               ENDIF
               ! Release the species to the background
               IF (i_species < 11) THEN
                  spc_name = TRIM(spc_names_p(i_species))
                  ind_spc_GC = Ind_(TRIM(spc_name))
                  background_mass = Spc(ind_spc_GC)%Conc(i_lon, i_lat, i_lev) * Vgrid_EU
               ELSEIF (i_species .lt. nspc_p) THEN
                  ! TOMAS tracer
                  i_tracer = MOD((i_species - 10), nspc_p_tomas_tracer)
                  IF (i_tracer == 0) CYCLE ! Skip TOMAS AW tracer ! i_tracer = nspc_p_tomas_tracer
                  spc_name = TRIM(spc_names_p(10+i_tracer))
                  ind_spc_GC_bin1 = Ind_(TRIM(spc_name))
                  ibin = FLOOR((i_species - 10.0_fp)/nspc_p_tomas_tracer) + 1
                  ind_spc_GC =  ind_spc_GC_bin1 +ibin - 1
                  background_mass = Spc(ind_spc_GC)%Conc(i_lon, i_lat, i_lev) * Vgrid_EU
               ENDIF
               Spc(ind_spc_GC)%Conc(i_lon, i_lat, i_lev) = (background_mass - mass_plume_diff) / Vgrid_EU
            ENDIF
            Plume1d_new%CONCNT1d(:,i_species) = ConcSlab
            ! Project from species Ra, Rb grid to SO2 Ra, Rb grid
            !IF (mass_plume_1D .GT. 1e-10_fp) THEN
            !   Plume1d_new%CONCNT1d(:,i_species) = ConcSlab*Vgrid_1D_temp/Vgrid_1D
            !ENDIF
         ENDDO
         ! file_1Dmass_NK01_ID = findFreeLun()
         ! WRITE(file_1Dmass_NK01,  &
         !       '("Plume-1D_NK01_mass_ID_",I0,"_time_", I0,".txt")') &
         !       Plume1d_new%LABEL, NINT(time_elapsed)
         ! CALL PLUME_CONC_DIAG_FILES_1D( &
         !       file_1Dmass_NK01_ID, file_1Dmass_NK01, &
         !       Plume1d_new%CONCNT1d(:,id_NK01_p) * &
         !       Plume1d_new%RA*Plume1d_new%RB*Plume1d_new%LENGTH*1.0e+6_fp, &
         !       RC)
         

         ! If no existing 1D segment, creating the first node
         IF(.NOT.ASSOCIATED(Plume1d_head))THEN
            Plume1d_head => Plume1d_new
            Plume1d_tail => Plume1d_new
            WRITE(6,*)'*** Ceate the first Plume1d node: label = ', Plume1d_new%label
         ELSE
            Plume1d_tail%next => Plume1d_new
            Plume1d_tail      => Plume1d_new
         ENDIF
         NULLIFY(Plume1d_new)

      ENDIF

      !IF((box_life .GT. Critical_day*24.0*60.0*60).OR. (ITS_TIME_FOR_EXIT())) THEN
      IF(Plume2d_curr%IsDissolve .OR. Plume2d_curr%IsTransfer) THEN
         Num_Plume2d = Num_Plume2d - 1
        !mass_S_SO2_r2_2D = mass_S_SO2_r2_2D + SUM(Plume2d_curr%CONCNT2d(:,:,id_SO2_p))  * Vgrid_2D 
                      !- Plume2d_curr%MassRef2d(id_SO2)
        !mass_S_SO4_r2_2D = mass_S_SO4_r2 + SUM(Plume2d_curr%CONCNT2d(:,:,id_SO4_p))  * Vgrid_2D 
                      ! Plume2d_curr%MassRef2d(id_SO4)
        !mass_S_SO2_4 = mass_S_SO2_3 - (SUM(Plume2d_curr%CONCNT2d(:,:,id_SO2))  * Vgrid_2D - &
                      ! Plume2d_curr%MassRef2d(id_SO2))
        !mass_S_SO4_4 = mass_S_SO4_3 - (SUM(Plume2d_curr%CONCNT2d(:,:,id_SO4))  * Vgrid_2D - &
                      ! Plume2d_curr%MassRef2d(id_SO4))

        ! Release mass in Eulerian grid only if IsDissolve is true
         IF (Plume2d_curr%IsDissolve) THEN
            Num_dissolve_2D = Num_dissolve_2D + 1
            DO i_species = 1, nspc_p
               IF ((i_species == id_OH_p) .or. (i_species == id_HO2_p) .or.  &
               (i_species == id_NH3_p) .or. &
                        (i_species == id_NH4_p) .or. (i_species == id_H2O_p) ) CYCLE
               spc_name     = TRIM(spc_names_p_use(i_species))
               ind_spc_GC   = Ind_(TRIM(spc_name))
               mass_release = SUM(Plume2d_curr%CONCNT2d(:,:,i_species))  * Vgrid_2D
               Spc(ind_spc_GC)%Conc(i_lon, i_lat, i_lev) =       &
                              MAX(0.0_fp, (Spc(ind_spc_GC)%Conc(i_lon, i_lat, i_lev)* Vgrid_EU + &
                              mass_release) / Vgrid_EU)
            ENDDO
            mass_S_SO2_r3_2D = mass_S_SO2_r3_2D - SUM(Plume2d_curr%CONCNT2d(:,:,id_SO2_p))  * Vgrid_2D 
            mass_S_SO4_r3_2D = mass_S_SO4_r3_2D - SUM(Plume2d_curr%CONCNT2d(:,:,id_SO4_p))  * Vgrid_2D
            WRITE(File_Plume_life_IU_2D,'(I0,2(1X,F10.1))') & 
            Plume2d_curr%LABEL, 0.0_fp, Plume2d_curr%LIFE
         ENDIF

        ! Delete the head node
        IF (.NOT. ASSOCIATED(Plume2d_prev)) THEN
            WRITE(6,*)'*** (2-D): Deleting head node: label = ', Plume2d_curr%label
            Plume2d_head => Plume2d_next
         ! Delete middle or tail node
        ELSE 
             Plume2d_prev%next => Plume2d_next
             ! If delete the tail node, update the tail pointer
             IF (.NOT.ASSOCIATED(Plume2d_next))THEN
                Plume2d_tail => Plume2d_prev
                WRITE(6,*)'*** (2-D): Deleting tail node: label = ', Plume2d_curr%label
             ELSE
                WRITE(6,*)'*** (2-D): Deleting middle node: label = ', Plume2d_curr%label
             ENDIF
         ENDIF
         DEALLOCATE(Plume2d_curr)
         Plume2d_curr => Plume2d_next

      ELSE
         mass_S_SO2_7_2D = mass_S_SO2_7_2D + SUM(Plume2d_curr%CONCNT2d(:,:,id_SO2_p))*Vgrid_2D
         mass_S_SO4_7_2D = mass_S_SO4_7_2D + SUM(Plume2d_curr%CONCNT2d(:,:,id_SO4_p))*Vgrid_2D
         Vgrid_2D_tot_3  = Vgrid_2D_tot_3 + Vgrid_2D*n_x_max*n_y_max

         C99                 =     Compute_C99_2D(Plume2d_curr%CONCNT2d(:,:,id_SO2_p))
         Core_Threshold      =     0.01_fp * C99
         N_Core              =     COUNT(Plume2d_curr%CONCNT2d(:,:,id_SO2_p) >= Core_Threshold)
         mass_S_SO2_8_2D     =     mass_S_SO2_8_2D + SUM(                                               &
                                   Plume2d_curr%CONCNT2d(:,:,id_SO2_p),                                 &
                                   MASK = Plume2d_curr%CONCNT2d(:,:,id_SO2_p) >= Core_Threshold         &
                                   ) *Vgrid_2D
         ! Define the plume core based on SO2 concentration
         mass_S_SO4_8_2D     =     mass_S_SO4_8_2D + SUM(                                               &
                                   Plume2d_curr%CONCNT2d(:,:,id_SO4_p),                                 &
                                   MASK = Plume2d_curr%CONCNT2d(:,:,id_SO2_p) >= Core_Threshold         &
                                   ) *Vgrid_2D
         Vgrid_2D_tot_4  = Vgrid_2D_tot_4 + Vgrid_2D*N_Core
         Spc(id_SO2pl)%Conc(i_lon, i_lat, i_lev) =   &
               (  Spc(id_SO2pl)%Conc(i_lon, i_lat, i_lev) * Vgrid_EU     +     &
                  SUM(Plume2d_curr%CONCNT2d(:,:,id_SO2_p))*Vgrid_2D     )  /   Vgrid_EU
         Spc(id_SO4pl)%Conc(i_lon, i_lat, i_lev) =   &
               (  Spc(id_SO4pl)%Conc(i_lon, i_lat, i_lev) * Vgrid_EU     +     &
                  SUM(Plume2d_curr%CONCNT2d(:,:,id_SO4_p))*Vgrid_2D     )  /   Vgrid_EU
#ifdef TOMAS
         DO ibin = 1, nbins
            mass_SF_bin_2D(ibin) = mass_SF_bin_2D(ibin) + &
                  SUM(Plume2d_curr%CONCNT2d(:,:, 10+2+(ibin-1)*nspc_p_tomas_tracer)) * Vgrid_2D
            mass_NK_bin_2D(ibin) = mass_NK_bin_2D(ibin) + &
                  SUM(Plume2d_curr%CONCNT2d(:,:, 10+1+(ibin-1)*nspc_p_tomas_tracer)) * Vgrid_2D
            mass_MK_bin_2D(ibin) = mass_MK_bin_2D(ibin) + &
                  SUM(Plume2d_curr%CONCNT2d(:,:, 10+2+(ibin-1)*nspc_p_tomas_tracer)) * Vgrid_2D * 1.375_fp + &
                  SUM(Plume2d_curr%CONCNT2d(:,:, 10+3+(ibin-1)*nspc_p_tomas_tracer)) * Vgrid_2D 
         ENDDO
#endif
         Plume2d_prev => Plume2d_curr
         Plume2d_curr => Plume2d_next 
      ENDIF
   ENDDO
   ! If all nodes were deleted
   IF ( .NOT. ASSOCIATED(Plume2d_head) ) THEN
      NULLIFY(Plume2d_tail)
   ENDIF

      

401 CONTINUE
    IF(.NOT.ASSOCIATED(Plume1d_head)) GOTO 400
    ! Start to dissolve Plume1d
    Plume1d_curr => Plume1d_head

    DO WHILE(ASSOCIATED(Plume1d_curr))

      ! file_1Dconc_SO4_ID = findFreeLun()
      ! WRITE(file_1Dconc_SO4,'("Plume-1D_SO4_conc_",I0,".txt")') NINT(time_elapsed)
      ! CALL PLUME_CONC_DIAG_FILES_1D(file_1Dconc_SO4_ID, file_1Dconc_SO4, Plume1d_curr%CONCNT1d(:,id_SO4_p), RC)

      ! file_1Dconc_SO2_ID = findFreeLun()
      ! WRITE(file_1Dconc_SO2,'("Plume-1D_SO2_conc_",I0,".txt")') NINT(time_elapsed)
      ! CALL PLUME_CONC_DIAG_FILES_1D(file_1Dconc_SO2_ID, file_1Dconc_SO2, Plume1d_curr%CONCNT1d(:,id_SO2_p), RC)

      ! ! file_1Dconc_OH_ID = findFreeLun()
      ! ! WRITE(file_1Dconc_OH,'("Plume-1D_OH_conc_",I0,".txt")') NINT(time_elapsed)
      ! ! CALL PLUME_CONC_DIAG_FILES_1D(file_1Dconc_OH_ID, file_1Dconc_OH, Plume1d_curr%CONCNT1d(:,id_OH_p), RC)
      
      ! DO ibin = 1, nbins
      !    ind_spc_p = 10+1+(ibin-1)*nspc_p_tomas_tracer
      !    file_1D_conc_NKbin_ID = findFreeLun()
      !    WRITE(file_1D_conc_NKbin,'("Plume-1D_NKbin_",I0,"_", I0,".txt")') ibin, NINT(time_elapsed)
      !    CALL PLUME_CONC_DIAG_FILES_1D(file_1D_conc_NKbin_ID, file_1D_conc_NKbin, Plume1d_curr%CONCNT1d(:,ind_spc_p), RC)
      
      !    ind_spc_p = 10+2+(ibin-1)*nspc_p_tomas_tracer
      !    file_1D_conc_SFbin_ID = findFreeLun()
      !    WRITE(file_1D_conc_SFbin,'("Plume-1D_SFbin_",I0,"_", I0,".txt")') ibin, NINT(time_elapsed)
      !    CALL PLUME_CONC_DIAG_FILES_1D(file_1D_conc_SFbin_ID, file_1D_conc_SFbin, Plume1d_curr%CONCNT1d(:,ind_spc_p), RC)

      ! ENDDO
      
      i_lon = Plume1d_curr%lon_ind
      i_lat = Plume1d_curr%lat_ind
      i_lev = Plume1d_curr%lev_ind
      box_RA =Plume1d_curr%RA
      box_RB =Plume1d_curr%RB
      box_length = Plume1d_curr%LENGTH
      Vgrid_1D     = box_RA*box_RB*box_length*1.0e+6_fp
      Vgrid_EU     = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ! [cm3]
      ! mass_S_SO2_6_1D = mass_S_SO2_6_1D + SUM(Plume1d_curr%CONCNT1D(:,id_SO2_p)) * Vgrid_1D
      ! mass_S_SO4_6_1D = mass_S_SO4_6_1D + SUM(Plume1d_curr%CONCNT1D(:,id_SO4_p)) * Vgrid_1D

      Plume1d_next => Plume1d_curr%next

      WRITE(File_Plume_location_IU_1D,'(5(I0,1X),3(F10.3, 1X))') &
            NINT(time_elapsed), Plume1d_curr%LABEL, &
            Plume1d_curr%lon_ind, Plume1d_curr%lat_ind, Plume1d_curr%lev_ind, &
            Plume1d_curr%LON, Plume1d_curr%LAT, Plume1d_curr%LEV, &
            Plume1d_curr%RA, Plume1d_curr%RB, Plume1d_curr%Length

      IF(Plume1d_curr%IsDissolve) THEN
         Num_Plume1d = Num_Plume1d -1 
         Num_dissolve_1D = Num_dissolve_1D + 1
         WRITE(File_Plume_life_IU_1D,'(I0,1X,F10.1)') &
            Plume1d_curr%LABEL, Plume1d_curr%LIFE
        ! Release the mass
         DO i_species = 1, nspc_p
            IF ((i_species == id_OH_p) .or. (i_species == id_HO2_p) .or.  &
               (i_species == id_NH3_p) .or. &
                        (i_species == id_NH4_p).or. &
               (i_species ==id_H2O_p) ) CYCLE
            !Vgrid_1D     = box_RA(i_species)*box_RB(i_species)*box_length*1.0e+6_fp
            spc_name     = TRIM(spc_names_p_use(i_species))
            ind_spc_GC   = Ind_(TRIM(spc_name))
            mass_release = SUM(Plume1d_curr%CONCNT1d(:,i_species))  * Vgrid_1D
            Spc(ind_spc_GC)%Conc(i_lon, i_lat, i_lev) =       &
                           MAX(0.0_fp, (Spc(ind_spc_GC)%Conc(i_lon, i_lat, i_lev)* Vgrid_EU + &
                           mass_release) / Vgrid_EU)
         ENDDO
         mass_S_SO2_r3_1D = mass_S_SO2_r3_1D - SUM(Plume1d_curr%CONCNT1d(:,id_SO2_p))  * Vgrid_1D 
         mass_S_SO4_r3_1D = mass_S_SO4_r3_1D - SUM(Plume1d_curr%CONCNT1d(:,id_SO4_p))  * Vgrid_1D 
        ! Delete the head node
        IF (.NOT. ASSOCIATED(Plume1d_prev)) THEN
            
            WRITE(6,*)'*** (1-D): Deleting head node: label = ', Plume1d_curr%label
            Plume1d_head => Plume1d_next
         ! Delete middle or tail node
        ELSE 
             Plume1d_prev%next => Plume1d_next
             ! If delete the tail node, update the tail pointer
             IF (.NOT.ASSOCIATED(Plume1d_next))THEN
                Plume1d_tail => Plume1d_prev
                WRITE(6,*)'*** (1-D): Deleting tail node: label = ', Plume1d_curr%label
             ELSE
                WRITE(6,*)'*** (1-D): Deleting middle node: label = ', Plume1d_curr%label
             ENDIF
         ENDIF
         DEALLOCATE(Plume1d_curr)
         Plume1d_curr => Plume1d_next
      !ELSEIF (splitting plume cross section)

      !ELSEIF (splitting plume length)
         
      ELSE
         mass_S_SO2_7_1D = mass_S_SO2_7_1D + SUM(Plume1d_curr%CONCNT1d(:,id_SO2_p))*Vgrid_1D
         mass_S_SO4_7_1D = mass_S_SO4_7_1D + SUM(Plume1d_curr%CONCNT1d(:,id_SO4_p))*Vgrid_1D
         Vgrid_1D_tot_3  = Vgrid_1D_tot_3 + Vgrid_1D*n_slab_max
         Spc(id_SO2pl)%Conc(i_lon, i_lat, i_lev) =   &
               (  Spc(id_SO2pl)%Conc(i_lon, i_lat, i_lev) * Vgrid_EU     +     &
                  SUM(Plume1d_curr%CONCNT1d(:,id_SO2_p))*Vgrid_1D     )  /   Vgrid_EU
         Spc(id_SO4pl)%Conc(i_lon, i_lat, i_lev) =   &
               (  Spc(id_SO4pl)%Conc(i_lon, i_lat, i_lev) * Vgrid_EU     +     &
                  SUM(Plume1d_curr%CONCNT1d(:,id_SO4_p))*Vgrid_1D     )  /   Vgrid_EU
#ifdef TOMAS
         DO ibin = 1, nbins 
            mass_SF_bin_1D(ibin) = mass_SF_bin_1D(ibin) + &
                  SUM(Plume1d_curr%CONCNT1d(:, 10+2+(ibin-1)*nspc_p_tomas_tracer)) * Vgrid_1D
            mass_NK_bin_1D(ibin) = mass_NK_bin_1D(ibin) + &
                  SUM(Plume1d_curr%CONCNT1d(:, 10+1+(ibin-1)*nspc_p_tomas_tracer)) * Vgrid_1D
            mass_MK_bin_1D(ibin) = mass_MK_bin_1D(ibin) + &
                  SUM(Plume1d_curr%CONCNT1d(:, 10+2+(ibin-1)*nspc_p_tomas_tracer)) * Vgrid_1D * 1.375_fp + &
                  SUM(Plume1d_curr%CONCNT1d(:, 10+3+(ibin-1)*nspc_p_tomas_tracer)) * Vgrid_1D 
         ENDDO
         ! WRITE (6, *) "Debug (BZ), 1-D-1, Nk is: "
         ! WRITE(*,'(15(1X,ES12.4))') (mass_NK_bin_1D(ibin), ibin=1,15)
         ! WRITE (6, *) "Debug (BZ), 1-D-1, SF is: "
         ! WRITE(*,'(15(1X,ES12.4))') (mass_SF_bin_1D(ibin), ibin=1,15)
#endif
         Plume1d_prev => Plume1d_curr
         Plume1d_curr => Plume1d_next 
      ENDIF
    ENDDO
    ! If all nodes were deleted
   IF ( .NOT. ASSOCIATED(Plume1d_head) ) THEN
      NULLIFY(Plume1d_tail)
   ENDIF
400 CONTINUE


    ! Cleanup pointer
    !IF (ALLOCATED(box_concnt_2D)) DEALLOCATE(box_concnt_2D)
   IF(ASSOCIATED(Plume2d_next)) nullify(Plume2d_next)
   IF(ASSOCIATED(Plume2d_curr)) nullify(Plume2d_curr)
   IF(ASSOCIATED(Plume2d_prev)) nullify(Plume2d_prev)
   IF(ASSOCIATED(Plume1d_new)) nullify(Plume1d_new)
   IF(ASSOCIATED(Plume1d_next)) nullify(Plume1d_next)
   IF(ASSOCIATED(Plume1d_curr)) nullify(Plume1d_curr)
   IF(ASSOCIATED(Plume1d_prev)) nullify(Plume1d_prev)
  END SUBROUTINE plume_structure_change

  SUBROUTINE lagrange_write_std(time_elapsed)

   !LOGICAL,          INTENT(IN   ) :: am_I_Root   ! root CPU?
   !INTEGER,          INTENT(INOUT) :: RC          ! Failure or success
   REAL(fp),          INTENT(INOUT) :: time_elapsed

   INTEGER                       :: ibin
   
   !OPEN( File_Smass_IU_2D,      FILE=TRIM( file_Smass_2D   ), STATUS='OLD',  &
   !      POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
   ! WRITE(File_Smass_IU_2D,'(*(ES12.4,1X))') time_elapsed, mass_S_SO2_inj_2D, &
   !       mass_S_SO2_r1_2D,  mass_S_SO2_r2_2D, mass_S_SO2_r3_2D, &
   !       mass_S_SO2_1_2D, mass_S_SO2_2_2D, mass_S_SO2_3_2D, mass_S_SO2_4_2D, &
   !        mass_S_SO2_5_2D, mass_S_SO2_6_2D, mass_S_SO2_7_2D, &
   !       mass_S_SO4_inj_2D,  mass_S_SO4_r1_2D,  mass_S_SO4_r2_2D, mass_S_SO4_r3_2D, & 
   !       mass_S_SO4_1_2D, mass_S_SO4_2_2D, mass_S_SO4_3_2D, mass_S_SO4_4_2D, &
   !       mass_S_SO4_5_2D, mass_S_SO4_6_2D, mass_S_SO4_7_2D, &
   !       Vgrid_2D_tot_1, Vgrid_2D_tot_2, Vgrid_2D_tot_3, &
   !       mass_S_SO2_r4, mass_S_SO4_r4
      WRITE(File_Smass_IU_2D,'(*(ES12.4,1X))') time_elapsed, mass_S_SO2_inj_2D, &
         mass_S_SO2_r1_2D,  mass_S_SO2_r2_2D, mass_S_SO2_r3_2D, &
         mass_S_SO2_7_2D, mass_S_SO2_8_2D, &
         mass_S_SO4_inj_2D,  mass_S_SO4_r1_2D,  mass_S_SO4_r2_2D, mass_S_SO4_r3_2D, & 
         mass_S_SO4_7_2D, mass_S_SO4_8_2D, &
         Vgrid_2D_tot_3, Vgrid_2D_tot_4, &
         mass_S_SO2_r4, mass_S_SO4_r4

   !OPEN( File_Smass_IU_1D,      FILE=TRIM( file_Smass_1D   ), STATUS='OLD',  &
   !      POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
   ! WRITE(File_Smass_IU_1D,'(*(ES12.4,1X))') time_elapsed, mass_S_SO2_inj_1D, &
   !       mass_S_SO2_r1_1D,  mass_S_SO2_r2_1D, mass_S_SO2_r3_1D, & 
   !       mass_S_SO2_1_1D, mass_S_SO2_2_1D, mass_S_SO2_3_1D, mass_S_SO2_4_1D, &
   !       mass_S_SO2_5_1D, mass_S_SO2_6_1D, mass_S_SO2_7_1D, &
   !       mass_S_SO4_inj_1D,  mass_S_SO4_r1_1D,  mass_S_SO4_r2_1D, mass_S_SO4_r3_1D, & 
   !       mass_S_SO4_1_1D, mass_S_SO4_2_1D, mass_S_SO4_3_1D, mass_S_SO4_4_1D, &
   !       mass_S_SO4_5_1D, mass_S_SO4_6_1D, mass_S_SO4_7_1D, &
   !       Vgrid_1D_tot_1, Vgrid_1D_tot_2, Vgrid_1D_tot_3
   WRITE(File_Smass_IU_1D,'(*(ES12.4,1X))') time_elapsed, mass_S_SO2_inj_1D, &
         mass_S_SO2_r1_1D,  mass_S_SO2_r2_1D, mass_S_SO2_r3_1D, & 
         mass_S_SO2_7_1D, &
         mass_S_SO4_inj_1D,  mass_S_SO4_r1_1D,  mass_S_SO4_r2_1D, mass_S_SO4_r3_1D, & 
         mass_S_SO4_7_1D, &
         Vgrid_1D_tot_3
   !OPEN( File_Plume_number_IU, FILE=TRIM(file_Plume_number), STATUS='OLD',  &
   !      POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
   WRITE(File_Plume_number_IU, '(ES12.4,6(1X,I0))') time_elapsed, Num_inject, &
   Num_Plume2d, Num_Plume1d, Num_dissolve_2D, Num_transfer_2D, Num_dissolve_1D
#ifdef TOMAS
      !OPEN( File_SF_bin_IU_2D,      FILE=TRIM( file_SF_bin_2D   ), STATUS='OLD',  &
      !    POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
      WRITE(File_SF_bin_IU_2D, '(ES12.4, *(1X,",",1X, ES12.4))') time_elapsed, (mass_SF_bin_2D(ibin), ibin= 1, nBins), &
                  mass_S_H2SO4_2D

      !OPEN( File_NK_bin_IU_2D,      FILE=TRIM( file_NK_bin_2D   ), STATUS='OLD',  &
      !    POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
      WRITE(File_NK_bin_IU_2D, '(ES12.4, *(1X,",",1X, ES12.4))') time_elapsed, (mass_NK_bin_2D(ibin), ibin= 1, nBins)

      !OPEN( File_MK_bin_IU_2D,      FILE=TRIM( file_MK_bin_2D   ), STATUS='OLD',  &
      !    POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
      WRITE(File_MK_bin_IU_2D, '(ES12.4, *(1X,",",1X, ES12.4))') time_elapsed, (mass_MK_bin_2D(ibin), ibin= 1, nBins)

      ! OPEN( File_NK_bin_IU_2D_1D,      FILE=TRIM( file_NK_bin_2D_1D   ), STATUS='OLD',  &
      !     POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
      ! WRITE(File_NK_bin_IU_2D_1D, '(ES12.4, *(1X,",",1X, ES12.4))') time_elapsed, (mass_NK_bin_2D_1D(ibin), ibin= 1, nBins)

      ! OPEN( File_SF_bin_IU_2D_1D,      FILE=TRIM( file_SF_bin_2D_1D   ), STATUS='OLD',  &
      !     POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
      ! WRITE(File_SF_bin_IU_2D_1D, '(ES12.4, *(1X,",",1X, ES12.4))') time_elapsed, (mass_SF_bin_2D_1D(ibin), ibin= 1, nBins)

      OPEN( File_SF_bin_IU_1D,      FILE=TRIM( file_SF_bin_1D   ), STATUS='OLD',  &
          POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
      WRITE(File_SF_bin_IU_1D, '(ES12.4, *(1X,",",1X, ES12.4))') time_elapsed, (mass_SF_bin_1D(ibin), ibin= 1, nBins), &
                  mass_S_H2SO4_1D

      OPEN( File_NK_bin_IU_1D,      FILE=TRIM( file_NK_bin_1D   ), STATUS='OLD',  &
          POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
      WRITE(File_NK_bin_IU_1D, '(ES12.4, *(1X,",",1X, ES12.4))') time_elapsed, (mass_NK_bin_1D(ibin), ibin= 1, nBins)

      OPEN( File_MK_bin_IU_1D,      FILE=TRIM( file_MK_bin_1D   ), STATUS='OLD',  &
          POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
      ! WRITE(File_MK_bin_IU_1D, '(ES12.4, *(1X,",",1X, ES12.4))') time_elapsed, (mass_MK_bin_1D(ibin), ibin= 1, nBins)
      ! WRITE (6, *) "Debug (BZ), 1-D-2, Nk is: "
      ! WRITE(*,'(15(1X,ES12.4))') (mass_NK_bin_1D(ibin), ibin=1,15)
      ! WRITE (6, *) "Debug (BZ), 1-D-2, SF is: "
      ! WRITE(*,'(15(1X,ES12.4))') (mass_SF_bin_1D(ibin), ibin=1,15)
#endif

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
          DX_m   = DX_GC/360.0 * 2*PI *Re*COS(box_lat/180*PI)
        else
          next_i_lon = i_lon - 1
          DX_m   = -1*DX_GC/360.0 * 2*PI *Re*COS(box_lat/180*PI)
        endif

        next_i_lat = i_lat

        if(next_i_lon>NX_GC) next_i_lon=next_i_lon-NX_GC
        if(next_i_lon<1)     next_i_lon=next_i_lon+NX_GC
        if(next_i_lat>NY_GC) next_i_lat=NY_GC
        if(next_i_lat<1)     next_i_lat=1

        D_wind = u(next_i_lon, next_i_lat, i_lev)-u(i_lon, i_lat, i_lev)
        Ly     = D_wind/DX_m

      ELSEIF(box_alpha>=1.25*PI)THEN
        if(box_lat>=Y_mid(i_lat))then
          next_i_lat = i_lat + 1
          DY_m   = DY_GC/360.0 * 2*PI*Re
        else
          next_i_lat = i_lat - 1
          DY_m   = -1*DY_GC/360.0 * 2*PI*Re
        endif

        next_i_lon = i_lon

        if(next_i_lon>NX_GC) next_i_lon=next_i_lon-NX_GC
        if(next_i_lon<1)     next_i_lon=next_i_lon+NX_GC
        if(next_i_lat>NY_GC) next_i_lat=NY_GC
        if(next_i_lat<1)     next_i_lat=1

        D_wind = v(next_i_lon, next_i_lat, i_lev)-v(i_lon, i_lat, i_lev)
        Ly     = D_wind/DY_m

      ELSEIF(box_alpha>=0.75*PI)THEN
        if(box_lon>=X_mid(i_lon))then
          next_i_lon = i_lon + 1
          DX_m   = DX_GC/360.0 * 2*PI *Re*COS(box_lat/180*PI)
        else
          next_i_lon = i_lon - 1
          DX_m   = -1*DX_GC/360.0 * 2*PI *Re*COS(box_lat/180*PI)
        endif

        next_i_lat = i_lat

        if(next_i_lon>NX_GC) next_i_lon=next_i_lon-NX_GC
        if(next_i_lon<1)     next_i_lon=next_i_lon+NX_GC
        if(next_i_lat>NY_GC) next_i_lat=NY_GC
        if(next_i_lat<1)     next_i_lat=1

        D_wind = u(next_i_lon, next_i_lat, i_lev)-u(i_lon, i_lat, i_lev)
        Ly     = D_wind/DX_m

      ELSEIF(box_alpha>=0.25*PI)THEN
        if(box_lat>=Y_mid(i_lat))then
          next_i_lat = i_lat + 1
          DY_m   = DY_GC/360.0 * 2*PI*Re
        else
          next_i_lat = i_lat - 1
          DY_m   = -1*DY_GC/360.0 * 2*PI*Re
        endif

        next_i_lon = i_lon

        if(next_i_lon>NX_GC) next_i_lon=next_i_lon-NX_GC
        if(next_i_lon<1)     next_i_lon=next_i_lon+NX_GC
        if(next_i_lat>NY_GC) next_i_lat=NY_GC
        if(next_i_lat<1)     next_i_lat=1

        D_wind = v(next_i_lon, next_i_lat, i_lev)-v(i_lon, i_lat, i_lev)
        Ly     = D_wind/DY_m

      ELSE
        if(box_lon>=X_mid(i_lon))then
          next_i_lon = i_lon + 1
          DX_m   = DX_GC/360.0 * 2*PI *Re*COS(box_lat/180*PI)
        else
          next_i_lon = i_lon - 1
          DX_m   = -1*DX_GC/360.0 * 2*PI *Re*COS(box_lat/180*PI)
        endif

        next_i_lat = i_lat

        if(next_i_lon>NX_GC) next_i_lon=next_i_lon-NX_GC
        if(next_i_lon<1)     next_i_lon=next_i_lon+NX_GC
        if(next_i_lat>NY_GC) next_i_lat=NY_GC
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
    if(init_lev==NZ_GC) init_lev = NZ_GC-1


    ! calculate the distance between particle and grid point
    do i = 1,2
    do j = 1,2
      ii = i + init_lon - 1
      jj = j + init_lat - 1

      ! For some special circumstance:
      if(ii==0)then
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii+NX_GC), Y_mid(jj))
      else if(ii==(NX_GC+1))then
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii-NX_GC), Y_mid(jj))
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
            wind_lonlat(k) =  Weight(1,1) * wind(NX_GC,init_lat,kk) &
                          + Weight(1,2) * wind(NX_GC,init_lat+1,kk) &
                          + Weight(2,1) * wind(init_lon+1,init_lat,kk) &
                          + Weight(2,2) * wind(init_lon+1,init_lat+1,kk)
        else if(init_lon==NX_GC)then
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
      if(init_lev==NZ_GC) init_lev = NZ_GC-1


      ! calculate the distance between particle and grid point
      j = 1
      do i = 1,2
        ii = i + init_lon - 1
        jj = j + init_lat - 1

        ! For some special circumstance:
        if(ii==0)then
        distance(i)= Distance_Circle(curr_lon, curr_lat, X_mid(ii+NX_GC), Y_mid(jj))
        else if(ii==(NX_GC+1))then
        distance(i)= Distance_Circle(curr_lon, curr_lat, X_mid(1), Y_mid(jj))
        else
        distance(i)= Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
        endif

      enddo


    if(ii==0)then
      distance(3)= Distance_Circle(curr_lon, curr_lat, X_mid(ii+NX_GC), 90.0e+0_fp)
    else if(ii==(NX_GC+1))then
      distance(3)= Distance_Circle(curr_lon, curr_lat, X_mid(1), 90.0e+0_fp)
    else
      distance(3)= Distance_Circle(curr_lon, curr_lat, X_mid(ii), 90.0e+0_fp)
    endif


    IF(distance(3)==0.0)THEN
        do k=1,2
          kk = k + init_lev - 1
          wind_lonlat(k) = SUM(wind(:,init_lat,kk))/NX_GC
        enddo
    ELSE
        ! Calculate the inverse distance weight
        do i=1,3
            Weight(i) = 1.0/distance(i) / sum( 1.0/distance(:) )
        enddo

        do k=1,2
          kk = k + init_lev - 1

          wind_polar = SUM(wind(:,init_lat,kk))/NX_GC      

          if(init_lon==0)then    
              wind_lonlat(k) =  Weight(1)*wind(NX_GC,init_lat,kk)   &
                              + Weight(2)*wind(init_lon+1,init_lat,kk)   &
                              + Weight(3)*wind_polar
          else if(init_lon==NX_GC)then
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
    real(fp)          :: uv_polars(NX_GC)
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
    if(init_lev==NZ_GC) init_lev = NZ_GC-1


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
       if(ii==NX_GC+1)then
          ii = 1
       endif
       if(ii==0)then
          ii = NX_GC
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
      do ii = 1,NX_GC
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

            uv_PS(3,k) = SUM(uv_polars)/NX_GC
          ENDIF

          IF(i_uv==0)THEN ! for v
          do ii = 1,NX_GC
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
             uv_PS(3,k) = SUM(uv_polars)/NX_GC
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
    if(init_lev==NZ_GC) init_lev = NZ_GC-1

    
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
      if(jj==NY_GC+1)then
        jj = jj-1
      endif

    
      ! For lon=180 deg:
      if(ii==NX_GC+1)then
        ii = 1
      endif
      if(ii==0)then
        ii = NX_GC
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
    ! for lon: Xedge_Sec - DX_GC = Xedge_first
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
    Find_iPLev = MAX(1, MIN(Find_iPLev, NZ_GC))
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
    init_lat = MAX(1, MIN(init_lat, NY_GC-1))
    ! For pressure level, P_mid(1) is about surface pressure
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif
    if(init_lev==0) init_lev = 1
    if(init_lev==NZ_GC) init_lev = NZ_GC-1

    do i = 1,2
    do j = 1,2

      ii = i + init_lon - 1
      jj = j + init_lat - 1

      ! For some special circumstance:
      if(ii==0)then
        distance(i,j) = &
           Distance_Circle(curr_lon, curr_lat, X_mid(ii+NX_GC), Y_mid(jj))
      else if(ii==(NX_GC+1))then
        distance(i,j) = &
           Distance_Circle(curr_lon, curr_lat, X_mid(ii-NX_GC), Y_mid(jj))
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
          u_lonlat(k) =  Weight(1,1) * u(NX_GC,init_lat,kk) &
                       + Weight(1,2) * u(NX_GC,init_lat+1,kk) &
                       + Weight(2,1) * u(init_lon+1,init_lat,kk) &
                       + Weight(2,2) * u(init_lon+1,init_lat+1,kk)
          v_lonlat(k) =  Weight(1,1) * v(NX_GC,init_lat,kk) &
                       + Weight(1,2) * v(NX_GC,init_lat+1,kk) &
                       + Weight(2,1) * v(init_lon+1,init_lat,kk) &
                       + Weight(2,2) * v(init_lon+1,init_lat+1,kk)
      else if(init_lon==NX_GC)then
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
    if(init_lon==0) init_lon=NX_GC
    if(init_lon==NX_GC+1) init_lon=1
    ! Index issue of P_BXHEIGHT, need to revisit later, BZ
    Delt_height = Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev),    &
                          P_edge(init_lev), P_edge(init_lev+1), 1 ) &   
                + Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev+1),   &
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
    init_lat = MAX(1, MIN(init_lat, NY_GC-1))
    ! For pressure level, P_mid(1) is about surface pressure
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif
    if(init_lev==0) init_lev = 1
    if(init_lev==NZ_GC) init_lev = NZ_GC-1

    do i = 1,2
    do j = 1,2

      ii = i + init_lon - 1
      jj = j + init_lat - 1

      ! For some special circumstance:
      if(ii==0)then
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii+NX_GC), Y_mid(jj))
      else if(ii==(NX_GC+1))then
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii-NX_GC), Y_mid(jj))
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
      ! Index issue of var, need to revisit later, BZ
      IF(init_lon==0)THEN
        var_lonlat(k) =  Weight(1,1) *var(NX_GC,init_lat,kk)   &
                       + Weight(1,2) *var(NX_GC,init_lat+1,kk)   &
                       + Weight(2,1) *var(1,init_lat,kk) &
                       + Weight(2,2) *var(1,init_lat+1,kk)
      ELSE IF(init_lon==NX_GC)THEN
        var_lonlat(k) =  Weight(1,1) *var(NX_GC,init_lat,kk)   &
                       + Weight(1,2) *var(NX_GC,init_lat+1,kk)   &
                       + Weight(2,1) *var(1,init_lat,kk) &
                       + Weight(2,2) *var(1,init_lat+1,kk)
      ELSE
        var_lonlat(k) =  Weight(1,1) *var(init_lon,init_lat,kk)   &
                       + Weight(1,2) *var(init_lon,init_lat+1,kk)   &
                       + Weight(2,1) *var(init_lon+1,init_lat,kk) &
                       + Weight(2,2) *var(init_lon+1,init_lat+1,kk)
      ENDIF
    enddo


    ! second vertical shear of wind
    
    if(init_lon==0)then
      Delt_height = Pa2meter( P_BXHEIGHT(NX_GC,init_lat,init_lev),    &
                            P_edge(init_lev), P_edge(init_lev+1), 1 ) &
                 + Pa2meter( P_BXHEIGHT(NX_GC,init_lat,init_lev+1),   &
                            P_edge(init_lev), P_edge(init_lev+1), 0 )
    else if(init_lon==NX_GC)then
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
!===================================================================
  FUNCTION Compute_plume_core_conc(Conc2D) Result(Core_Mean)
   ! This function compute average concentration of 2-D concentration array
   ! in plume core only, excluding surrounding cells with background conc
    REAL(fp), INTENT(IN)     ::    Conc2D(:,:)
    REAL(fp)                 ::    Core_Mean
    REAL(fp)                 ::    C99, Core_Threshold
    INTEGER                  ::    N_core

    C99                 =     Compute_C99_2D(Conc2D)
    Core_Threshold      =     0.01_fp * C99
    N_Core              =     COUNT(Conc2D >= Core_Threshold)
    IF (N_Core > 0) THEN

    Core_Mean = SUM(                                    &
               Conc2D,                                  &
               MASK = Conc2D >= Core_Threshold          &
               ) / REAL(N_Core, kind = fp)
    ELSE
    Core_Mean = 0.0_fp

    END IF
  END FUNCTION Compute_plume_core_conc
!===================================================================
  FUNCTION Compute_C99_2D(Conc2D) Result(C99)
   USE ERROR_MOD,    ONLY : ERROR_STOP, IT_IS_NAN
   ! This function calculates the 99th percentile of 2-D concentration
   ! array, which is then used to determine the average concentration 
   !inside the plume (consider grid cell with concentration
   ! larger than 0.01*C99 only to exclude cells with background 
   !concentration level )
   REAL(fp), INTENT(IN)     ::    Conc2D(:,:)
   REAL(fp)                 ::    C99
   REAL(fp), ALLOCATABLE    ::    Conc1D(:)

   IF (SIZE(Conc2D) == 0) THEN
        CALL ERROR_STOP('input array is empty','Compute_C99_2D')
   END IF

   ALLOCATE(Conc1D(SIZE(Conc2D)))
   Conc1D = RESHAPE(Conc2D, [SIZE(Conc2D)])
   C99 = Compute_C99_1D(Conc1D)
   DEALLOCATE(Conc1D)
  END FUNCTION Compute_C99_2D
!===================================================================

!===================================================================
  FUNCTION Compute_C99_1D( Conc1D ) RESULT( C99 )
    USE ERROR_MOD,    ONLY : ERROR_STOP, IT_IS_NAN
    REAL(fp), INTENT(IN)     ::    Conc1D(:)
    REAL(fp)                 ::    C99
    REAL(fp), ALLOCATABLE    ::    Sorted_Conc1D(:)
    INTEGER                  ::    N
    INTEGER                  ::    Idx
     
    
    N = SIZE(Conc1D)
    IF (N == 0) THEN
        CALL ERROR_STOP('input array is empty','Compute_C99_1D')
    END IF

    ALLOCATE(Sorted_Conc1D(N))
    Sorted_Conc1D = Conc1D
    CALL QuickSort_Real(Sorted_Conc1D, 1, N)
    Idx = CEILING(0.99_fp * REAL(N, kind=fp))
    Idx = MAX(1, MIN(Idx, N))
    C99 = Sorted_Conc1D(Idx)
    DEALLOCATE(Sorted_Conc1D)
  END FUNCTION Compute_C99_1D
!======================================================================
! Sort a real-valued array in ascending order
!======================================================================
RECURSIVE SUBROUTINE QuickSort_Real( Values, First, Last )

    REAL(fp),     INTENT(INOUT)     :: Values(:)
    INTEGER,      INTENT(IN)        :: First
    INTEGER,      INTENT(IN)        :: Last

    INTEGER                         :: I
    INTEGER                         :: J
    REAL(fp)                        :: Pivot
    REAL(fp)                        :: Temporary

    IF (First >= Last) RETURN

    I = First
    J = Last

    Pivot = Values((First + Last) / 2)

    DO

        DO WHILE (Values(I) < Pivot)
            I = I + 1
        END DO

        DO WHILE (Values(J) > Pivot)
            J = J - 1
        END DO

        IF (I <= J) THEN

            Temporary = Values(I)
            Values(I) = Values(J)
            Values(J) = Temporary

            I = I + 1
            J = J - 1

        END IF

        IF (I > J) EXIT

    END DO

    IF (First < J) CALL QuickSort_Real(Values, First, J)
    IF (I < Last)  CALL QuickSort_Real(Values, I, Last)

END SUBROUTINE QuickSort_Real

!======================================================================  
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


	! Check whether target grid cell area is equal to Dx_GC*Dy_gc
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
						box_Ra, box_Rb, Cslab)
 
    IMPLICIT NONE

    REAL(fp)    :: Pdx, Pdy, theta1, Height1
    REAL(fp), INTENT(INOUT)  :: box_Ra, box_Rb
    REAL(fp), INTENT(INOUT)  :: Cslab(n_slab_max)
    REAL(fp)    :: Pc_2D(n_x_max,n_y_max) !, Ec_2D(n_x_max,n_y_max)
    ! INTEGER, INTENT(IN)      :: n_slab_max

    !INTEGER, parameter     :: Nb = n_slab_max
    INTEGER, parameter     :: Na = 128
    !INTEGER                :: Nb  
    
    INTEGER     :: Nb_mid, Na_mid

    REAL(fp)    :: X2d(Na,n_slab_max), Y2d(Na,n_slab_max), C2d(Na,n_slab_max) !, Extra_C2d(Na,Nb)

    REAL(fp)    :: LenB, LenA
    REAL(fp)    :: Adx, Ady, Bdx, Bdy
    REAL(fp)    :: Prod, M, Lb, La

    REAL(fp)    :: Cslab_tmp(n_slab_max) !, Extra_slab(Nb)

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

        DO j=1,n_slab_max,1
          Cslab_tmp(j)     = SUM(C2d(:,j)) *(LenA*LenB)/ (La*Lb)
        ENDDO


      ! Return:
      box_Ra = La
      box_Rb = Lb
      Cslab(1:n_slab_max) = Cslab_tmp(1:n_slab_max)

  END SUBROUTINE Slab_init_bilinear

  ! Definea 1-D grid based on one specific species distribution, 
  ! and interpolate all species to the same grid
  SUBROUTINE Define_Slab_Grid(                                  &
               Pdx, Pdy, theta1, SO2_Yscale,                    &
               box_Ra, box_Rb, X2d, Y2d                         &
             )

    IMPLICIT NONE

    INTEGER, PARAMETER :: Na = 128

    REAL(fp), INTENT(IN)  :: Pdx
    REAL(fp), INTENT(IN)  :: Pdy
    REAL(fp), INTENT(IN)  :: theta1
    REAL(fp), INTENT(IN)  :: SO2_Yscale

    REAL(fp), INTENT(OUT) :: box_Ra
    REAL(fp), INTENT(OUT) :: box_Rb

    REAL(fp), INTENT(OUT) :: X2d(Na,n_slab_max)
    REAL(fp), INTENT(OUT) :: Y2d(Na,n_slab_max)

    INTEGER :: i, j
    INTEGER :: Na_mid, Nb_mid

    REAL(fp) :: LenA, LenB
    REAL(fp) :: Adx, Ady
    REAL(fp) :: Bdx, Bdy

    Na_mid = Na / 2
    Nb_mid = n_slab_max / 2

    ! Define the physical geometry using SO2 only.
    !
    ! Retain your original geometry relation here if that is
    ! the intended model definition.
    box_Ra = SO2_Yscale * TAN(theta1)

    ! Resolution of each slab in the cross-plume d direction
    box_Rb = Pdy

    ! Sampling interval along the uniform b direction
    LenA = box_Ra / REAL(Na,fp)

    ! Slab resolution along the resolved d direction
    LenB = box_Rb

    ! Direction b:
    ! theta = 90 degrees minus the angle between b and s.
    Adx = LenA * SIN(theta1)
    Ady = LenA * COS(theta1)

    ! Direction d, perpendicular to b
    Bdx = LenB * COS(theta1)
    Bdy = LenB * SIN(theta1)

    ! Construct the center slab.
    DO i = 1, Na
        X2d(i,Nb_mid) =                              &
            -LenA*REAL(Na_mid,fp)*SIN(theta1)         &
            + LenA*(REAL(i,fp)-0.5_fp)*SIN(theta1)

        Y2d(i,Nb_mid) =                              &
            -LenA*REAL(Na_mid,fp)*COS(theta1)         &
            + LenA*(REAL(i,fp)-0.5_fp)*COS(theta1)
    END DO

    ! Offset by half a slab so that the coordinates represent
    ! slab centers consistently.
    X2d(:,Nb_mid) = X2d(:,Nb_mid) + 0.5_fp*Bdx
    Y2d(:,Nb_mid) = Y2d(:,Nb_mid) - 0.5_fp*Bdy

    DO j = Nb_mid + 1, n_slab_max
        X2d(:,j) = X2d(:,j-1) - Bdx
        Y2d(:,j) = Y2d(:,j-1) + Bdy
    END DO

    DO j = Nb_mid - 1, 1, -1
        X2d(:,j) = X2d(:,j+1) + Bdx
        Y2d(:,j) = Y2d(:,j+1) - Bdy
    END DO

  END SUBROUTINE Define_Slab_Grid
  
  SUBROUTINE Interpolate_To_Slab(                               &
               Pdx, Pdy, Pc_2D, box_Ra, box_Rb,                 &
               X2d, Y2d, Cslab                                  &
             )

    IMPLICIT NONE

    INTEGER, PARAMETER :: Na = 128

    REAL(fp), INTENT(IN)  :: Pdx
    REAL(fp), INTENT(IN)  :: Pdy
    REAL(fp), INTENT(IN)  :: Pc_2D(n_x_max,n_y_max)

    REAL(fp), INTENT(IN)  :: box_Ra
    REAL(fp), INTENT(IN)  :: box_Rb

    REAL(fp), INTENT(IN)  :: X2d(Na,n_slab_max)
    REAL(fp), INTENT(IN)  :: Y2d(Na,n_slab_max)

    REAL(fp), INTENT(OUT) :: Cslab(n_slab_max)

    REAL(fp) :: C2d(Na,n_slab_max)

    INTEGER :: i, j

    !$OMP PARALLEL DO COLLAPSE(2)                     &
    !$OMP DEFAULT(SHARED) PRIVATE(i,j)
    DO j = 1, n_slab_max
        DO i = 1, Na
            C2d(i,j) = Interplt_2D(                    &
                Pdx, Pdy, X2d(i,j), Y2d(i,j), Pc_2D, 2 &
            )
        END DO
    END DO
    !$OMP END PARALLEL DO

    ! Concentration is assumed uniform along b.
    ! Therefore average the sampled concentrations along b.
    DO j = 1, n_slab_max
        Cslab(j) = SUM(C2d(:,j)) / REAL(Na,fp)
    END DO

  END SUBROUTINE Interpolate_To_Slab


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

SUBROUTINE CHEM_SO2_OH_PLUME( dt, K, SO2, OH, SO4, HO2, PH2SO4, chem_status, chem_debug_value)
  !-------------------------------------------------------------------
  ! Do one implicit chemistry step for:
  !
  !   SO2 + OH -> SO4 + HO2 + PH2SO4
  !
  ! Rate = k * [SO2] * [OH]
  !
  ! Backward Euler:
  !   x = dt * k * (SO2_new) * (OH_new)
  ! with
  !   SO2_new = SO2_old - x
  !   OH_new  = OH_old  - x
  !
  ! This gives a quadratic equation for x:
  !
  !   a*x^2 - b*x + c = 0
  !
  ! where
  !   a = dt*k
  !   b = 1 + dt*k*(SO2_old + OH_old)
  !   c = dt*k*SO2_old*OH_old
  !
  ! We choose the smaller root to keep 0 <= x <= min(SO2,OH).
  !-------------------------------------------------------------------
  ! Inputs
  REAL(fp),        INTENT(IN)  :: dt    ! Chemistry timestep
  REAL(fp),        INTENT(IN)  :: K    ! Reaction rate constant
  !INTEGER,         INTENT(IN)  :: i_x
  !INTEGER,         INTENT(IN)  :: i_y
  !INTEGER,         INTENT(IN)  :: i_box ! Print debug info

  ! Outputs
  INTEGER,         INTENT(OUT) :: chem_status
  REAL(fp),        INTENT(OUT) :: chem_debug_value

  ! In/Out concentrations, molec/cm3 
  REAL(fp),        INTENT(INOUT) :: SO2
  REAL(fp),        INTENT(INOUT) :: OH
  REAL(fp),        INTENT(INOUT) :: SO4
  REAL(fp),        INTENT(INOUT) :: HO2
  REAL(fp),        INTENT(INOUT) :: PH2SO4
  REAL(fp), PARAMETER            :: SO2_min = 1.0e-6_fp
  REAL(fp), PARAMETER            :: OH_min  = 1.0e-6_fp
  REAL(fp), PARAMETER            :: x_min   = 1.0e-12_fp

  ! Local variables
  REAL(fp)                       :: SO2_old, OH_old
  REAL(fp)                       :: x, a, b, c, disc
  !REAL(fp), PARAMETER            :: tiny = 1.0d-300
  REAL(fp), PARAMETER            :: tiny = 100.0_fp * epsilon(1.0_fp)

  chem_status  = 0
  ! Save old values
  SO2_old      = max(SO2, 0.0_fp)
  OH_old       = max(OH , 0.0_fp)

  IF (K <= 0.0_fp) RETURN

  ! (Double check) Typical background level SO2 is ~1E6, [OH] is ~1E1, K=1E-13
  IF (SO2_old <= SO2_min) THEN
    !Write (6, *) "Debug: (BZ): At plume grid [x, y, box]", i_x, i_y, i_box
    !Write (6, *) "Debug: (BZ) Skip plume chem due to low SO2 conc: ", SO2_old
    chem_debug_value = SO2_old
    chem_status = 1
    RETURN
  ENDIF
  IF (OH_old  <= OH_min ) THEN
    !Write (6, *) "Debug: (BZ): At plume grid [x, y, box]", i_x, i_y, i_box
    !Write (6, *) "Debug: (BZ) Skip plume chem due to low OH conc: ", OH_old
    chem_debug_value = OH_old
    chem_status = 2
    RETURN
  ENDIF
  IF (dt * K * SO2_old * OH_old <= x_min) THEN
    !Write (6, *) "Debug: (BZ): At plume grid [x, y, box]", i_x, i_y, i_box
    !Write (6, *) "Debug: (BZ) Skip plume chem due to low K[SO2][OH]dt: ", dt * K * SO2_old * OH_old
    chem_debug_value = dt * K * SO2_old * OH_old
    chem_status = 3
    RETURN
  ENDIF
  ! Backward Euler coefficients
  a = dt * K
  b = 1.0_fp + a * (SO2_old + OH_old)
  c = a * SO2_old * OH_old

  ! If a is extremely small, fall back to explicit Euler
  if (abs(a) < tiny) then
    x = dt * k * SO2_old * OH_old
  else
    disc = b*b - 4.0_fp*a*c
    disc = max(disc, 0.0_fp)

    ! Smaller root
    x = (2.0_fp*c) / (b + sqrt(disc))
  end if

  ! Enforce physical bounds
  x = max(x, 0.0_fp)
  x = min(x, SO2_old)
  x = min(x, OH_old)

  ! Update reactants
  SO2 = SO2_old - x
  OH  = OH_old  - x

  ! Update products using the stoichiometry exactly as written
  SO4    = max( SO4 + x,    0.0_fp )
  HO2    = max( HO2 + x,    0.0_fp )
  PH2SO4 = max( PH2SO4 + x, 0.0_fp )
ENDSUBROUTINE CHEM_SO2_OH_PLUME
FUNCTION GET_PLUME_SPC_ID(name) RESULT(id)
  USE CharPak_Mod,     ONLY : To_UpperCase
  
  CHARACTER(LEN=*), INTENT(IN) :: name
  CHARACTER(LEN=LEN(name))     :: name_std
  INTEGER                      :: id, i
  
  name_std =  To_UpperCase(ADJUSTL(TRIM(name)))        
  id       = -1

  DO i = 1, nspc_p
    IF ( name_std == To_UpperCase(ADJUSTL(TRIM(spc_names_p(i)))) ) THEN
      id = i
      RETURN
    END IF
  END DO
END FUNCTION GET_PLUME_SPC_ID

SUBROUTINE PLUME_CONC_DIAG_FILES_2D(file_id, file_name_str, conc_2D, RC)
   
   !USE InquireMod,      ONLY : findFreeLun

   ! Input:
   INTEGER,           INTENT(IN)    :: file_id
   REAL(fp),          INTENT(IN)    :: conc_2D(n_x_max, n_y_max)
   CHARACTER(LEN=255),  INTENT(IN)  :: file_name_str
   ! Output:
   INTEGER,           INTENT(INOUT) :: RC

   ! Other vars
   INTEGER                          :: i_x, i_y
   CHARACTER(LEN=255)               :: file_name_str_1

   !file_id = findFreeLun()
   file_name_str_1 = ADJUSTL(file_name_str)
   OPEN(file_id, FILE=TRIM(file_name_str_1), STATUS='REPLACE', &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', IOSTAT=RC)
   DO i_x = 1, n_x_max
      WRITE(file_id,'(*(ES12.4,1X))') &
         (conc_2D(i_x,i_y), i_y = 1, n_y_max)
   ENDDO
   CLOSE(file_id)

END SUBROUTINE

SUBROUTINE PLUME_CONC_DIAG_FILES_1D(file_id, file_name_str, conc_1D, RC)
   
   !USE InquireMod,      ONLY : findFreeLun

   ! Input:
   INTEGER,           INTENT(IN)    :: file_id
   REAL(fp),          INTENT(IN)    :: conc_1D(n_slab_max)
   CHARACTER(LEN=255),  INTENT(IN)    :: file_name_str
   ! Output:
   INTEGER,           INTENT(INOUT) :: RC

   ! Other vars
   INTEGER                          :: i_slab
   CHARACTER(LEN=255)               :: file_name_str_1

   !file_id = findFreeLun()
   file_name_str_1 = ADJUSTL(file_name_str)
   OPEN(file_id, FILE=TRIM(file_name_str_1), STATUS='REPLACE', &
         FORM='FORMATTED', ACCESS='SEQUENTIAL', IOSTAT=RC)

   WRITE(file_id,'(*(ES12.4,1X))') &
      (conc_1D(i_slab), i_slab = 1, n_slab_max)

   CLOSE(file_id)

END SUBROUTINE PLUME_CONC_DIAG_FILES_1D

! TOMAS related subroutine and functions below
! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
! Copied from GEOS-Chem TOMAS_mod.F90
! !IROUTINE: nh4bulktobin
!
! !DESCRIPTION: Subroutine NH4BULKTOBIN takes the bulk ammonium aerosol from
!  GEOS-Chem and fraction it to each bin according to sulfate mole fraction in
!  each bin
!  Written by Win Trivityanurak, Sep 26, 2008
!  .
!  Make sure that we work with mass or mass conc.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NH4BULKTOBIN( MSULF, NH4B, MAMMO )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN)   :: MSULF(nbins)  ! size-resolved sulfate [kg]
    REAL(fp),  INTENT(IN)   :: NH4B          ! Bulk NH4 mass [kg]
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),  INTENT(OUT)  :: MAMMO(nbins)  ! size-resolved NH4 [kg]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                 :: K
    REAL(fp)                :: TOTMASS, NH4TEMP

    !=================================================================
    ! NH4BULKTOBIN begins here
    !=================================================================

    MAMMO(:) = 0.e+0_fp

    ! Sum the total sulfate
    TOTMASS = 0.e+0_fp
    DO K = 1, nbins
       TOTMASS = TOTMASS + MSULF(K)
    ENDDO

    IF ( TOTMASS .eq. 0.e+0_fp ) RETURN

    ! Limit the amount of NH4 entering TOMAS calculation
    ! if it is very NH4-rich, just limit the amount to balance
    ! existing 30-bin-summed SO4 assuming (NH4)2SO4 in such case
    !  (NH4)2 mass = (SO4)mass / 96. * 2. * 18. = 0.375*(SO4)mass
    ! (win, 9/28/08)
    NH4TEMP = NH4B
    IF ( NH4B/TOTMASS > 0.375e+0_fp ) &  !make sure we use mass ratio
         NH4TEMP = 0.375e+0_fp * TOTMASS

    ! Calculate ammonium aerosol scale to each bin
    DO K = 1, nbins
       MAMMO(K) = MSULF(K) / TOTMASS * NH4TEMP
    ENDDO

    !write(777,*) NH4B/TOTMASS

    RETURN

  END SUBROUTINE NH4BULKTOBIN
!EOC

!  Copy from GEOS-Chem tomas_mod.F90
!BOP
!
! !IROUTINE: ezwatereqm
!
! !DESCRIPTION: WRITTEN BY Peter Adams, March 2000
!     .
!     This routine uses the current RH to calculate how much water is
!     in equilibrium with the aerosol.  Aerosol water concentrations
!     are assumed to be in equilibrium at all times and the array of
!     concentrations is updated accordingly.
!     .
!     Introduced to GEOS-CHEM by Win Trivitayanurak. May 8, 2006.
!     This file is replacing the old ezwatereqm that was not compatible
!     with multicomponent aerosols.  The new ezwatereqm use external
!     functions to do ISORROPIA-result curve fitting for each aerosol
!     component.
!     WARNING :
!      *** Watch out for the new aerosol species added in the future!
!     .
!     This version of the routine works for sulfate and sea salt
!     particles.  They are assumed to be externally mixed and their
!     associated water is added up to get total aerosol water.
!     wr is the ratio of wet mass to dry mass of a particle.  Instead
!     of calling a thermodynamic equilibrium code, this routine uses a
!     simple curve fits to estimate wr based on the current humidity.
!     The curve fit is based on ISORROPIA/HETP results for ammonium bisulfate
!     at 273 K and sea salt at 273 K.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EZWATEREQM( Mke, RHTOMAS )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN)    :: RHTOMAS
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(fp),INTENT(INOUT) :: Mke(nbins,ICOMPHARD)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP

!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
!
    INTEGER             :: k, j
    REAL(fp)            :: so4mass, naclmass, ocilmass
    REAL(fp)            :: wrso4, wrnacl, wrocil
    REAL(fp)            :: rhe

    !========================================================================
    ! EZWATEREQM begins here!
    !========================================================================

    rhe=100.e+0_fp*rhtomas
    if (lowRH == 1) THEN !JKodros RH switch
       if (rhe .gt. 90.e+0_fp) rhe=90.e+0_fp
    ELSE
       if (rhe .gt. 99.e+0_fp) rhe=99.e+0_fp
    END IF !JKodros RH switch
    if (rhe .lt. 1.e+0_fp) rhe=1.e+0_fp

    do k=1,nbins

       so4mass=Mke(k,srtso4)*1.2  !1.2 converts kg so4 to kg nh4hso4
       wrso4=waterso4(rhe)

       ! Add condition for srtnacl in case of running so4 only. (win, 5/8/06)
       if (srtnacl.gt.0) then
          naclmass=Mke(k,srtnacl) !already as kg nacl - no conv necessary
          ! wrnacl=waternacl(rhe)
          wrnacl = 1.e+0_fp
       else
          naclmass = 0.e+0_fp
          wrnacl = 1.e+0_fp
       endif

       if (srtocil.gt.0) then
          ocilmass=Mke(k,srtocil) !already as kg ocil - no conv necessary
          ! wrocil=waterocil(rhe)
          wrocil = 1.e+0_fp
       else
          ocilmass = 0.e+0_fp
          wrocil = 1.e+0_fp
       endif

       Mke(k,srth2o)=so4mass*(wrso4-1.e+0_fp)+naclmass &
                     *(wrnacl-1.e+0_fp) &
                     +ocilmass*(wrocil-1.e+0_fp)

    enddo

    RETURN

  END SUBROUTINE EZWATEREQM

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: waterso4
!
! !DESCRIPTION: Function WATERSO4 uses the current RH to calculate how much
!  water is in equilibrium with the sulfate.  Aerosol water concentrations are
!  assumed to be in equilibrium at all times and the array of concentrations is
!  updated accordingly.
!   Introduced to GEOS-CHEM by Win Trivitayanurak. 8/6/07
!   Adaptation of ezwatereqm used in size-resolved sulfate only sim
!   November, 2001
!   ezwatereqm WRITTEN BY Peter Adams, March 2000
!\\
!\\
! !INTERFACE:
!
  FUNCTION WATERSO4( RHE ) RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    REAL(fp) :: RHE ! Relative humidity (0-100 scale)
!
! !RETURN VALUE:
!
    REAL(fp) :: VALUE

! !REMARKS:
!  waterso4 is the ratio of wet mass to dry mass of a particle.  Instead
!  of calling a thermodynamic equilibrium code, this routine uses a
!  simple curve fit to estimate wr based on the current humidity.
!  The curve fit is based on ISORROPIA results for ammonium bisulfate
!  at 273 K.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! WATERSO4 begins here!
    !=================================================================

    if (rhe .gt. 99.) rhe=99.
    if (rhe .lt. 1.) rhe=1.

    if (rhe .gt. 96.) then
       value=0.7540688*rhe**3-218.5647*rhe**2+21118.19*rhe-6.801999e5
    else
       if (rhe .gt. 91.) then
          value=8.517e-2*rhe**2 -15.388*rhe +698.25
       else
          if (rhe .gt. 81.) then
             value=8.2696e-3*rhe**2 -1.3076*rhe +53.697
          else
             if (rhe .gt. 61.) then
                value=9.3562e-4*rhe**2 -0.10427*rhe +4.3155
             else
                if (rhe .gt. 41.) then
                   value=1.9149e-4*rhe**2 -8.8619e-3*rhe +1.2535
                else
                   value=5.1337e-5*rhe**2 +2.6266e-3*rhe +1.0149
                endif
             endif
          endif
       endif
    endif

    !check for error
    if (value .gt. 30.) then
       write(*,*) 'ERROR in waterso4'
       write(*,*) rhe,value
       STOP
    endif

  END FUNCTION WATERSO4
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: spinup
!
! !DESCRIPTION: Function SPINUP retuns .TRUE. or .FALSE. whether or not the
!  current time in the run have passed the spin-up period.  This would be used
!  to determine if certain errors should be fixed and let slipped or to stop a
!  run with an error message.  (win, 8/2/07)
!  ====> Be cautious that TIMEBEGIN should be changed according to
!         whatever your spin-up beginning time is
!  Example of TIMEBEGIN (in julian time)
!         2001/07/01 = 144600.0
!         2000/11/01 = 138792.0
!\\
!\\
! !INTERFACE:
!
  FUNCTION SPINUP( DAYS ) RESULT( VALUE )
!
! !USES:
!
    USE TIME_MOD!,     ONLY : GET_TAU , GET_TAUb
!
! !INPUT PARAMETERS:
!
    REAL*4,    INTENT(IN) :: DAYS   ! Spin-up duration (day)
!
! !RETURN VALUE:
!
    LOGICAL               :: VALUE
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*4                 :: TIMENOW, TIMEBEGIN, TIMEINIT, HOURS

    !========================================================================
    ! SPINUP begins here!
    !========================================================================

    TIMENOW   = GET_TAU()   ! Current time in the run (Julian time) (hrs)
    TIMEBEGIN = GET_TAUb()  ! Begin time of this run (hrs)
    TIMEINIT  = 141000. !2/1/2001    ! Start time for spin-up (hrs)
    HOURS = DAYS * 24.0     ! Period allow error to pass (hrs)
    ! Write (6, *) 'debug (BZ): In spinup(plume), TIMENOW= ', TIMENOW, ' TIMEBEGIN = ', TIMEBEGIN
    ! Criteria to let error go or to terminate the run
    !IF ( TIMENOW > MIN( TIMEBEGIN, TIMEINIT ) + HOURS  ) THEN
    IF ( TIMENOW > TIMEBEGIN + HOURS  ) THEN
       VALUE = .FALSE.
    ELSE
       VALUE = .TRUE.
    ENDIF

  END FUNCTION SPINUP
!EOC
  !------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: eznh3eqm
!
! !DESCRIPTION: Subroutine EZNH3REQM2 puts ammonia to the particle phase until
!  there is 2 moles of ammonium per mole of sulfate and the remainder
!  of ammonia is left in the gas phase. (win, 9/30/08)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EZNH3EQM( Gce, Mke )
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(fp),  INTENT(INOUT)  :: Gce(ICOMPHARD) !sfarina - fixed incorrect definition of Gc array
    REAL(fp),  INTENT(INOUT)  :: Mke(nBins,ICOMPHARD)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer       ::  ibin
    REAL(fp)        :: tot_nh3  !total kmoles of ammonia
    REAL(fp)        :: tot_so4  !total kmoles of so4
    REAL(fp)        :: sfrac    !fraction of sulfate that is in that bin

    !========================================================================
    ! EZNH3EQM begins here!
    !========================================================================

    ! get the total number of kmol nh3
    tot_nh3 = Gce(srtnh4)/17.e+0_fp
    do ibin=1,nBins
       tot_nh3 = tot_nh3 + Mke(ibin,srtnh4)/18.e+0_fp
    enddo

    ! get the total number of kmol so4
    tot_so4 = 0.e+0_fp
    do ibin=1,nBins
       tot_so4 = tot_so4 + Mke(ibin,srtso4)/96.e+0_fp
    enddo

    ! see if there is free ammonia
    if (tot_nh3/2.e+0_fp.lt.tot_so4)then  ! no free ammonia
       Gce(srtnh4) = 0.e+0_fp ! no gas phase ammonia
       do ibin=1,nBins
          sfrac = Mke(ibin,srtso4)/96.e+0_fp/tot_so4
          Mke(ibin,srtnh4) = sfrac*tot_nh3*18.e+0_fp ! put the ammonia where the sulfate is
          ! Debug
          !if ( Mke(k,srtnh4) < 0.0 ) then
          !   print *,'negative gas phase ammonia in eznh3eqm!!'
          !   print *,'bin  ', k
          !endif
       enddo
    else ! free ammonia
       do ibin=1,nBins
          Mke(ibin,srtnh4) = Mke(ibin,srtso4)/96.e+0_fp*2.e+0_fp*18.e+0_fp ! fill the particle phase
          ! Debug
          !if ( Mke(k,srtnh4) < 0.0 ) then
          !   print *,'negative gas phase ammonia in eznh3eqm!!'
          !   print *,'bin  ', k
          !endif
       enddo
       Gce(srtnh4) = (tot_nh3 - tot_so4*2.e+0_fp)*17.e+0_fp ! put whats left over in the gas phase
       ! Debug
       !if ( Gce(srtnh4) < 0.0 ) then
       !   print *,'negative gas phase ammonia in eznh3eqm!!'
       !endif

    endif

    RETURN

  END SUBROUTINE EZNH3EQM

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: storenm
!
! !DESCRIPTION: Subroutine STORENM stores values of Nk and Mk into Nkd and Mkd
!  for diagnostic purposes.  Also do gas phase concentrations. (win, 7/23/07)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE STORENM(Nk, Nkd, Mk, Mkd, Gc, Gcd )
!
! !INPUT PARAMETERS:
!
    REAL(fp),INTENT(IN)    :: Nk(nBins)
    REAL(fp),INTENT(IN)    :: Mk(nBins, ICOMPHARD)
    REAL(fp),INTENT(IN)    :: Gc(ICOMPHARD)
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),INTENT(OUT)   :: Nkd(nBins)
    REAL(fp),INTENT(OUT)   :: Mkd(nBins, ICOMPHARD)
    REAL(fp),INTENT(OUT)   :: Gcd(ICOMPHARD)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: J, K

    !sfarina
    DO J= 1, ICOMPHARD
       Gcd(J)=Gc(J)
    ENDDO
    DO K = 1, nBins
       Nkd(K)=Nk(K)
       DO J= 1, ICOMPHARD
          Mkd(K,J)=Mk(K,J)
       ENDDO
    ENDDO

    RETURN

  END SUBROUTINE STORENM
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: mnfix
!
! !DESCRIPTION: Subroutine MNFIX examines the mass and number distrubution and
!  determine if any bins have an average mass outside their normal range.  This
!  can happen because some process, e.g. advection, seems to treat the mass and
!  number species inconsistently.  If any bins are out of range, I shift some
!  mass and number to a new bin in a way that conserves both. (win, 7/23/07)
!  Originally written by Peter Adams, September 2000
!  Modified for GEOS-CHEM by Win Trivitayanurak (win@cmu.edu)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MNFIX ( NK, MK, ERRORSWITCH )
!
! !USES:
!
    USE ERROR_MOD,    ONLY : ERROR_STOP, IT_IS_NAN
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(fp),  INTENT(INOUT) :: NK(nBins),  MK(nBins, ICOMPHARD)
    LOGICAL, INTENT(INOUT) :: ERRORSWITCH
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer             :: K,J,KK !counters
    integer             :: NEWBIN !bin number into which mass is shifted
    REAL(fp)            :: XOLD, XNEW !average masses of old and new bins
    REAL(fp)            :: DRYMASS !dry mass of in a bin
    REAL(fp)            :: AVG !average dry mass of particles in bin
    REAL(fp)            :: NUM_INITIAL !number of particles initially in problem bin
    REAL(fp)            :: NSHIFT  !number to shift to new bin
    REAL(fp)            :: MSHIFT !mass to shift to new bin
    REAL(fp)            :: FJ !fraction of mass that is component j
    REAL(fp)            :: save1,save2,save3,save4,save5
    REAL(fp), PARAMETER :: EPS  = 1.e-20_fp !small number for Nk
    REAL(fp), PARAMETER :: EPS2 = 1.e-28_fp !small number for Mk
    REAL(fp), PARAMETER :: TINY = 1.e-36_fp !small number
    REAL(fp), PARAMETER :: VTINY= 1.e-50_fp !very small number

    LOGICAL             :: FIXERROR
    LOGICAL             :: PRT
    REAL(fp)            :: TOTMAS, TOTNUM !for print debug

    !=================================================================
    ! MNFIX begins here!
    !=================================================================

    FIXERROR = .TRUE.
    !PRT = .FALSE.
    PRT = ERRORSWITCH !just carrying a signal to print out value at the observed box - since mnfix does not have any information about I,J,L location. (Win, 9/27/05)
    !ERRORSWITCH = .FALSE.
    !PRT = .FALSE.             !TO AVOID THE HUGE AMOUNT OF PRINTING (JKodros 6/2/15)
    !xk(1)=xk(2)/2.e+0_fp  ! jrp for some reason xk(1) is changing?!
    save1=xk(1)
    
    ! Check for any incoming negative values or NaN
    !--------------------------------------------------------------------------
    DO K = 1, nBins
       IF ( IT_IS_NAN(NK(K)) ) THEN
          !PRINT *,'11 Found Nan in Nk at bin',K
          ERRORSWITCH = .TRUE.
          print *,'11 MNFIX(0): Found NaN in Nk(,',k,'): '
          GOTO 300
       ENDIF
       DO J = 1, ICOMPHARD
          IF ( IT_IS_NAN(MK(K,J)) ) THEN
             PRINT *,'11 Found Nan in Mk at bin',K,'component',J
             ERRORSWITCH = .TRUE.
             GOTO 300
          ENDIF
       ENDDO
       IF ( NK(K) < 0e+0_fp ) THEN
          IF ( PRT ) THEN
             PRINT *,'MNFIX[0]: FOUND NEGATIVE N'
             PRINT *, 'Bin, N', K, NK(K)
          ENDIF
          IF ( ABS(NK(K)) < 1e+0_fp .and. FIXERROR ) THEN
             NK(K) = 0e+0_fp
             IF ( PRT ) PRINT *,'Negative N < -1.0 Reset to zero'
          ELSE
             ERRORSWITCH = .TRUE.
             print *,'MNFIX(0): Found negative Nk(',k,') >-1e+0_fp, value = ', NK(K)
             GOTO 300          !exit mnfix if found negative error (win, 4/18/06)
          ENDIF
       ENDIF
       IF ( IT_IS_NAN(NK(K)) ) THEN
          !PRINT *,'Found Nan in Nk at bin',K
          ERRORSWITCH = .TRUE.
          print *,'MNFIX(0): Found NaN in Nk(,',k,')'
          GOTO 300
       ENDIF
       DO J = 1, ICOMPHARD
          IF ( MK(K,J) < 0e+0_fp ) THEN
             IF ( PRT ) THEN
                PRINT *,'MNFIX[0]: FOUND NEGATIVE M'
                PRINT *,'Bin, Comp, Mk', K, J, MK(K,J)
             ENDIF
             IF( ABS(MK(K,J)) < 1e-5_fp .and. FIXERROR ) THEN
                MK(K,J) = 0e+0_fp
                IF ( PRT ) PRINT *,'Negative M < -1.d-5 Reset to zero'
             ELSE
                ERRORSWITCH =.TRUE.
                print *,'MNFIX(0): Found negative Mk(',k,',comp',j,'), value = ', MK(K,J)
                GOTO 300       !exit mnfix if found negative error (win, 4/18/06)
             ENDIF
          ENDIF
          IF ( IT_IS_NAN(MK(K,J)) ) THEN
             PRINT *,'Found Nan in Mk at bin',K,'component',J
             ERRORSWITCH = .TRUE.
             GOTO 300
          ENDIF
       ENDDO                   !icomp
    ENDDO                     !ibins
    save2=xk(1)

    ! JRP check for neg numbers
    !DO K = 1,IBINS
    !  IF (NK(K) < 0.e+0_fp) THEN
    !     print*,'1 NK < 0 in MNFIX',K,NK(K)
    !  ENDIF
    !  DO J=1,ICOMP
    !     IF (MK(K,J) < 0.e+0_fp) THEN
    !        print*,'1 MK < 0 in MNFIX',K,J,MK(K,J)
    !     ENDIF
    !  ENDDO
    !  IF ( IT_IS_NAN(NK(K)) ) THEN
    !     PRINT *,'11 Found Nan in Nk at bin',K
    !     ERRORSWITCH = .TRUE.
    !     print *,'11 MNFIX(0): Found NaN in Nk(,',k,')'
    !     GOTO 300
    !  ENDIF
    !  DO J = 1, ICOMP
    !     IF ( IT_IS_NAN(MK(K,J)) ) THEN
    !        PRINT *,'11 Found Nan in Mk at bin',K,'component',J
    !        ERRORSWITCH = .TRUE.
    !        GOTO 300
    !     ENDIF
    !  ENDDO
    !ENDDO

    ! Check if both number and mass are zero, if yes then exit mnfix.
    !----------------------------------------------------------------
    TOTNUM = 0e+0_fp
    TOTMAS = 0e+0_fp
    DO K = 1,nBins
       TOTNUM = TOTNUM + NK(K)
       DO J=1,ICOMPHARD-2
          TOTMAS = TOTMAS + MK(K,J)
       ENDDO
    ENDDO
    IF ( TOTNUM == 0e+0_fp .AND. TOTMAS == 0e+0_fp ) THEN
       IF ( PRT ) PRINT *,'MNFIX: Nk=Mk=0. Exit now'
       GOTO 300
    ENDIF

    ! If number is tiny ( < EPS) then set it to zero
    !DO K = 1,IBINS
    !   IF ( NK(K) <= EPS ) THEN
    !      NK(K) = 0e+0_fp
    !      DO J= 1, ICOMP-1
    !         MK(K,J) = 0e+0_fp
    !      ENDDO               !STOP  !original (win, 9/1/05)
    !   ENDIF
    !ENDDO

    ! If N is tiny and M is tiny, set both to zeroes
    !--------------------------------------------------------
    DO K = 1, nBins
       IF ( IT_IS_NAN(NK(K)) ) THEN
          ! PRINT *,'22 Found Nan in Nk at bin',K
          ERRORSWITCH = .TRUE.
          print *,'22 MNFIX(0): Found NaN in Nk(,',k,')'
          GOTO 300
       ENDIF
       DO J = 1, ICOMPHARD
          IF ( IT_IS_NAN(MK(K,J)) ) THEN
             PRINT *,'22 Found Nan in Mk at bin',K,'component',J
             ERRORSWITCH = .TRUE.
             GOTO 300
          ENDIF
       ENDDO
       IF ( NK(K) <= EPS .AND. NK(K)>= 0e+0_fp ) THEN
          !print*,'1111'
          !print*,k,EPS,xk(K),xk(K+1)
          !print*,'word up'
          NK(K) = EPS
          !NK(K) = 0.e+0_fp
          !DO J = 1, ICOMP-IDIAG
          DO J = 1, ICOMPHARD
             if (J .eq. 1) then
                !MK(K,J) = EPS*sqrt(xk(K)*xk(K+1))
                MK(K,J) = EPS*AVGMASS(k)
                !MK(K,J) = 0.e+0_fp
             else
                MK(K,J) = VTINY
             endif
          enddo
          !print*,'allbins',MK(:,1)
       ENDIF ! If tiny number
       TOTMAS = SUM(MK(K,1:ICOMPHARD-2))
       if (TOTMAS.lt.eps2) then
          !print*,'2222'
          NK(K) = EPS
          !NK(K) = 0.e+0_fp
          DO J = 1, ICOMPHARD
             !DO J = 1, ICOMP-IDIAG
             if (J .eq. 1) then
                !MK(K,J) = EPS*sqrt(xk(K)*xk(K+1))
                MK(K,J) = EPS*AVGMASS(k)
                !MK(K,J) = 0.e+0_fp
             else
                MK(K,J) = VTINY
             endif
          enddo
       endif
    ENDDO
    save3=xk(1)

    ! JRP check for neg numbers
    DO K = 1,nBins
       IF (NK(K) < 0.e+0_fp) THEN
          print*,'2 NK < 0 in MNFIX',K,NK(K)
       ENDIF
       DO J=1,ICOMPHARD
          IF (MK(K,J) < 0.e+0_fp) THEN
             print*,'2 MK < 0 in MNFIX',K,J,MK(K,J)
          ENDIF
       ENDDO
       IF ( IT_IS_NAN(NK(K)) ) THEN
          PRINT *,'2 Found Nan in Nk at bin',K
          ERRORSWITCH = .TRUE.
          print *,'2 MNFIX(0): Found NaN in Nk(,',k,')'
          GOTO 300
       ENDIF
       DO J = 1, ICOMPHARD
          IF ( IT_IS_NAN(MK(K,J)) ) THEN
             PRINT *,'2 Found Nan in Mk at bin',K,'component',J
             ERRORSWITCH = .TRUE.
             GOTO 300
          ENDIF
       ENDDO
    ENDDO

    ! Check to see if any bins are completely out of bounds for min or max bin
    !-------------------------------------------------------------------------
    DO K = 1, nBins
       DRYMASS = 0.e+0_fp
       DO J = 1, ICOMPHARD-2
          DRYMASS = DRYMASS + MK(K,J)
       ENDDO

       IF ( NK(k) == 0e+0_fp ) THEN
          !AVG = SQRT( xk(K)* xk(K+1) )
          AVG = SQRT( AVGMASS(k) )
       ELSE
          AVG = DRYMASS/ NK(K)
       ENDIF

       IF ( AVG >  xk(nBins+1) ) THEN
          IF ( PRT ) PRINT *, 'MNFIX [1]: AVG > Xk(ibins+1) at bin',K
          IF ( FIXERROR ) THEN
             !out of bin range - remove some mass
             MSHIFT = NK(k)* xk(nBins+1)/ 1.2
             DO J= 1, ICOMPHARD
                MK(K,J) = MK(K,J)* MSHIFT/ (DRYMASS+EPS2)
             ENDDO
          ELSE
             ERRORSWITCH = .TRUE.
             print *,'MNFIX(1): AVG>Xk(ibins+1) at bin',K
             GOTO 300
          ENDIF
       ENDIF
       IF ( AVG < xk(1)) THEN
          IF( PRT ) PRINT *,'MNFIX [2]: AVG < Xk(1)'
          IF( FIXERROR ) THEN
             !out of bin range - remove some number
             NK(K) = DRYMASS/ ( xk(1)* 1.2 )
          ELSE
             ERRORSWITCH = .TRUE.
             print *,'MNFIX(1): AVG < Xk(1) at bin',K
             GOTO 300
          ENDIF
       ENDIF
    ENDDO

    ! JRP check for neg numbers
    DO K = 1,nBins
       IF (NK(K) < 0.e+0_fp) THEN
          print*,'3 NK < 0 in MNFIX',K,NK(K)
       ENDIF
       DO J=1,ICOMPHARD
          IF (MK(K,J) < 0.e+0_fp) THEN
             print*,'3 MK < 0 in MNFIX',K,J,MK(K,J)
          ENDIF
       ENDDO
       IF ( IT_IS_NAN(NK(K)) ) THEN
          PRINT *,'3 Found Nan in Nk at bin',K
          ERRORSWITCH = .TRUE.
          print *,'3 MNFIX(0): Found NaN in Nk(,',k,')'
          GOTO 300
       ENDIF
       DO J = 1, ICOMPHARD
          IF ( IT_IS_NAN(MK(K,J)) ) THEN
             PRINT *,'3 Found Nan in Mk at bin',K,'component',J
             ERRORSWITCH = .TRUE.
             GOTO 300
          ENDIF
       ENDDO
    ENDDO
    save4=xk(1)

    !if (PRT) then !<step5.1-temp>
    !   print *,'After_Check2 ---------------------'
    !   do k=1,ibins
    !      totmas = sum(MK(k,1:icomp-1))
    !      print *, totmas,NK(k), totmas/NK(k)
    !   enddo
    !endif

    !print*,1,NK(1),NK(2)
    !print*,1,MK(1,:)
    !print*,1,MK(2,:)

    ! Check to see if any bins are out of bounds
    !-------------------------------------------------------------------
    DO K = 1, nBins
       !if (PRT) print *,'Now at bin',k !<step4.4>tmp (win, 9/28/05)

       DRYMASS = 0.e+0_fp
       DO J = 1, ICOMPHARD-2
          DRYMASS = DRYMASS + MK(K,J)
       ENDDO

       IF ( NK(K) == 0e+0_fp ) THEN
          !AVG = SQRT(xk(K)*xk(K+1)) !set to mid-range value
          AVG = AVGMASS(k) !set to mid-range value
       ELSE
          AVG = DRYMASS/NK(K)
       ENDIF

       !if (PRT) then     !<step5.1-temp>
       !   print *,'After_Check3---------------------'
       !   totmas = sum(MK(k,1:icomp-1))
       !   print *, totmas,NK(k), totmas/NK(k)
       !endif

       ! If over boundary of the current bin
       IF ( AVG >  xk(K+1) ) THEN
          IF ( PRT ) PRINT *, 'MNFIX [3]: AVG>Xk(',K+1,')'
          !IF ( PRT ) CALL DEBUGPRINT(NK,MK,0,0,0,'inside MNFIX')
          IF ( FIXERROR ) THEN
             !Average mass is too high - shift to higher bin
             !KK = K + 1 ! jrp, this was causing errors
             !ERRORSWITCH=.TRUE.
             KK = K
             XNEW = xk(KK+1)/ 1.1
             if ( PRT ) PRINT *, 'k',k,'AVG',AVG,' XNEW ',XNEW
100          IF ( XNEW <= AVG ) THEN
                IF ( KK < nBins ) THEN
                   KK = KK + 1
                   XNEW = xk(KK+1)/ 1.1
                   if (PRT) PRINT *, '..move up to bin ',KK,' XNEW ',XNEW
                   GOTO 100
                ELSE
                   ! Already reach highest bin - must remove some mass (win, 8/1/07)
                   ! Updated by jrp 3/1/2012
                   MSHIFT = NK(k)* xk(k+1)/ 1.1
                   if( PRT ) PRINT*,' Mass being discarded: '
                   DO J= 1, ICOMPHARD
                      !if (PRT)
                      !print*,'Removing mass in MNFIX',MSHIFT, DRYMASS
                      MK(K,J) = MK(K,J)* MSHIFT/ (DRYMASS)
                   ENDDO
                   ! and recalculate dry mass (win, 8/1/07)
                   DRYMASS = 0.e+0_fp ! jrp fix 2/29/12
                   DO J = 1, ICOMPHARD-2
                      DRYMASS = DRYMASS + MK(K,J)
                   ENDDO
                   GOTO 111
                ENDIF
             ENDIF

             if(PRT)print*,'Old NK',NK(k),'Old DRYMASS',DRYMASS,'bin',k

             !XOLD = SQRT( xk(K)* xk(K+1) )
             XOLD = AVGMASS(k)
             NUM_INITIAL = NK(K)
             NSHIFT = ( DRYMASS - XOLD * NUM_INITIAL )/ ( XNEW - XOLD )
             MSHIFT = XNEW * NSHIFT
             NK(K) = NK(K) - NSHIFT
             NK(KK) =NK(KK) + NSHIFT

             if(prt) then
                print*,'NSHIFT',NSHIFT, 'MSHIFT',MSHIFT
                print*,'New NK',k,NK(k),' Nk(kk)',kk,NK(kk)
                print*,'Total mass bin',k,sum(MK(k,1:ICOMPHARD-2))
                print*,'SO4 mass bin  ',k,(MK(k,srtso4))
                print*,'Total mass bin',kk,sum(MK(kk,1:ICOMPHARD-2))
                print*,'SO4 mass bin  ' ,kk,(MK(kk,srtso4))
             endif

             DO J = 1, ICOMPHARD-2
                FJ = MK(K,J)/ DRYMASS
                MK(K,J) = XOLD * NK(K) * FJ
                MK(KK,J) = MK(KK,J) + MSHIFT * FJ
             ENDDO

             if(prt) then
                print*,'After shift mass'
                print*,'Total mass bin',k,sum(MK(k,1:ICOMPHARD-2))
                print*,'SO4 mass bin  ',k,(MK(k,srtso4))
                print*,'Total mass bin',kk,sum(MK(kk,1:ICOMPHARD-2))
                print*,'SO4 mass bin  ',kk,(MK(kk,srtso4))
             endif

          ELSE
             ERRORSWITCH = .TRUE.
             PRINT *, 'MNFIX(3) : AVG>Xk(',K+1,')'
             GOTO 300
          ENDIF    ! Fixerror
       ENDIF       ! AVG > Xk(k+1)

       !if (PRT) then     !<step5.1-temp>
       !   print *,'After_Check4---------------------'
       !   totmas = sum(MK(k,1:icomp-1))
       !   print *, totmas,NK(k), totmas/NK(k)
       !endif

       ! If under boundary of the current bin
111    IF ( AVG <  xk(K) ) THEN
          IF ( PRT ) PRINT *,'MNFIX [4]: AVG<Xk(',K,')'
          IF ( FIXERROR ) THEN
             !average mass is too low - shift number to lower bin
             !KK = K - 1 ! jrp potential for errors here
             KK = K
             XNEW = xk(KK)* 1.1
200          IF ( XNEW >= AVG ) THEN
                IF ( KK > 1 ) THEN
                   KK = KK - 1
                   XNEW = xk(KK)* 1.1
                   GOTO 200
                ELSE
                   ! Already reach lowest bin - must remove some number (win, 8/1/07)
                   NK(K) = DRYMASS/ ( xk(1)* 1.2 )
                   GOTO 222
                ENDIF
             ENDIF
             !XOLD = SQRT(xk(K)* xk(K+1))
             XOLD = AVGMASS(k)
             NUM_INITIAL = NK(K)
             NSHIFT = NUM_INITIAL - DRYMASS/XOLD !(win, 10/20/08)
             !Prior to 10/20/08 (win)
             !NSHIFT = (DRYMASS - XOLD * NUMBER)/ ( XNEW - XOLD )
             MSHIFT = XNEW * NSHIFT
             NK(K) = NK(K) - NSHIFT
             NK(KK) = NK(KK) + NSHIFT
             DO J=1,ICOMPHARD
                FJ = MK(K,J)/ DRYMASS
                MK(K,J) = XOLD * NK(K) * FJ
                MK(KK,J) = MK(KK,J) + MSHIFT * FJ
             ENDDO

          ELSE
             ERRORSWITCH = .TRUE.
             PRINT *, 'MNFIX(4): AVG < Xk(',k,')'
             GOTO 300
          ENDIF
222    ENDIF

       !if (PRT) then     !<step5.1-temp>
       !   print *,'After_Check5---------------------'
       !   totmas = sum(MK(k,1:icomp-1))
       !   print *, totmas,NK(k), totmas/NK(k)
       !endif
       !if (PRT) print *,MK(k,1),NK(k), MK(k,1)/NK(k),'Check5'!<step4.4>tmp (win, 9/28/05)

    ENDDO ! loop bin
    save5=xk(1)

    !print*,2,NK(1),NK(2)
    !print*,2,MK(1,:)
    !print*,2,MK(2,:)

    ! JRP check for neg numbers
    DO K = 1,nBins
       IF (NK(K) < 0.e+0_fp) THEN
          print*,'4 NK < 0 in MNFIX',K,NK(K)
       ENDIF
       DO J=1,ICOMPHARD
          IF (MK(K,J) < 0.e+0_fp) THEN
             print*,'4 MK < 0 in MNFIX',K,J,MK(K,J)
             print*,'saved xk1s',save1,save2,save3,save4,save5
             print*,'xk',xk
          ENDIF
       ENDDO
       IF ( IT_IS_NAN(NK(K)) ) THEN
          PRINT *,'4 Found Nan in Nk at bin',K
          ERRORSWITCH = .TRUE.
          print *,'4 MNFIX(0): Found NaN in Nk(,',k,')'
          GOTO 300
       ENDIF
       DO J = 1, ICOMPHARD
          IF ( IT_IS_NAN(MK(K,J)) ) THEN
             PRINT *,'4 Found Nan in Mk at bin',K,'component',J
             ERRORSWITCH = .TRUE.
             GOTO 300
          ENDIF
       ENDDO
    ENDDO

    !if (PRT) then !<step5.1-temp>
    ! Catch any small negative values resulting from fixing
    !--------------------------------------------------------------------------
    DO K = 1, nBins
       IF ( NK(K) < 0e+0_fp ) THEN
          IF ( PRT ) THEN
             PRINT *,'MNFIX[5]: FOUND NEGATIVE N'
             PRINT *, 'Bin, N', K, NK(K)
          ENDIF
          IF ( ABS(NK(K)) < 1e+0_fp .and. FIXERROR ) THEN
             NK(K) = 0e+0_fp
             IF ( PRT ) PRINT *,'Negative N > -1.0 Reset to zero'
          ELSE
             ERRORSWITCH = .TRUE.
             PRINT *, 'MNFIX(5): Negative N after fixing at bin',k
             GOTO 300          !exit mnfix if found negative error (win, 4/18/06)
          ENDIF
       ENDIF
       DO J = 1, ICOMPHARD
          IF ( MK(K,J) < 0e+0_fp ) THEN
             IF ( PRT ) THEN
                PRINT *,'MNFIX[6]: FOUND NEGATIVE M'
                PRINT *,'Bin, Comp, Mk', K, J, MK(K,J)
             ENDIF
             IF( ABS(MK(K,J)) < 1D-5 .and. FIXERROR ) THEN
                MK(K,J) = 0e+0_fp
                IF ( PRT ) PRINT *,'Negative M > -1.d-5 Reset to zero'
             ELSE
                ERRORSWITCH =.TRUE.
                PRINT *, 'MNFIX(6): Negative M after fixing at bin',k
                GOTO 300       !exit mnfix if found negative error (win, 4/18/06)
             ENDIF
          ENDIF
       ENDDO                   !icomp
    ENDDO                     !ibins

    ! JRP check for neg numbers
    DO K = 1,nBins
       IF (NK(K) < 0.e+0_fp) THEN
          print*,'5 NK < 0 in MNFIX',K,NK(K)
       ENDIF
       DO J=1,ICOMPHARD
          IF (MK(K,J) < 0.e+0_fp) THEN
             print*,'5 MK < 0 in MNFIX',K,J,MK(K,J)
          ENDIF
       ENDDO
       IF ( IT_IS_NAN(NK(K)) ) THEN
          PRINT *,'5 Found Nan in Nk at bin',K
          ERRORSWITCH = .TRUE.
          print *,'5 MNFIX(0): Found NaN in Nk(,',k,')'
          GOTO 300
       ENDIF
       DO J = 1, ICOMPHARD
          IF ( IT_IS_NAN(MK(K,J)) ) THEN
             PRINT *,'5 Found Nan in Mk at bin',K,'component',J
             ERRORSWITCH = .TRUE.
             GOTO 300
          ENDIF
       ENDDO
    ENDDO

    ! Check any last inconsistent M=0 or N=0
    !--------------------------------------------------------
    DO K = 1, nBins
       DRYMASS = 0.e+0_fp
       DO J = 1, ICOMPHARD-2
          DRYMASS = DRYMASS + MK(K,J)
       ENDDO
       IF ( NK(K) /= 0e+0_fp .AND. DRYMASS == 0e+0_fp .or. &
            NK(K) == 0e+0_fp .AND. DRYMASS /= 0e+0_fp     ) THEN
          PRINT *, '5.5 set nk, mk to ZERO for all bins'
          DO J = 1, ICOMPHARD
             MK(K,J)=0.e+0_fp
             NK(K) = 0.e+0_fp
          ENDDO
          MK(K,ICOMPHARD) = 0.e+0_fp !Set aerosol water to zero too
       ENDIF                  ! If tiny number
    ENDDO

    ! JRP check for neg numbers
    DO K = 1,nBins
       IF (NK(K) < 0.e+0_fp) THEN
          print*,'6 NK < 0 in MNFIX',K,NK(K)
          STOP
       ENDIF
       DO J=1,ICOMPHARD
          IF (MK(K,J) < 0.e+0_fp) THEN
             print*,'6 MK < 0 in MNFIX',K,J,MK(K,J)
             STOP
          ENDIF
       ENDDO
    ENDDO

300 CONTINUE

    IF (ERRORSWITCH) THEN
555    FORMAT (3E15.5E2)
       WRITE(6,*)'END OF MNFIX ( WHERE? )'
       WRITE(6,*)'DRYMAS-excl-NH4  NK      DRYMASS/NK'
       DO K = 1,nBins
          TOTMAS = SUM(MK(K,1:ICOMPHARD-1))
          !PRINT *, TOTMAS,NK(K), TOTMAS/NK(K)
       ENDDO

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !%%% NOTE: NK will be IBINS+1 upon exiting the loop, which will cause an
       !%%% out-of-bounds error.  Comment this out for now, unless it should be
       !%%% inserted into the DOloop
       !WRITE(6,555)
       !        TOTMAS, NK(K),
       !        TOTMAS/ NK(K)
       !print*,'-----------'
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !call debugprint( NK, MK, 0,0,0,'End of MNFIX')

       !write(*,*)'Nk'
       !write(*,*) NK(1:30)
       !write(*,*)'Mk(srtso4)'
       !write(*,*) MK(1:30,srtso4)
       !write(*,*)'Mk(srth2o)'
       !write(*,*) MK(1:30,srth2o)
       !STOP 'Negative Nk or Mk at after mnfix'  !comment out this to make it stop outside mnfix so that I can print out the i,j,l (location) of the error (win, 9/1/05)
    ENDIF

    RETURN

  END SUBROUTINE MNFIX
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cond_nuc
!
! !DESCRIPTION: This subroutine calculates the change in the aerosol size
!  distribution due to so4 condensation and binary/ternary nucleation during
!  the overal microphysics timestep.
!  WRITTEN BY Jeff Pierce, May 2007 for GISS GCM-II'
!  Put in GEOS-Chem by Win T. 9/30/08
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE COND_NUC(Nki,Mki,Gci,Nkf,Mkf,Gcf,fnavg,fn1avg, &
                      H2SO4rate,dti,num_iter,Nknuc,Mknuc,Nkcond,Mkcond, &
                      ionrate, surf_area, BOXVOL, BOXMASS, TEMPTMS, PRES, &
                      RHTOMAS, errswitch, lev)
!
! !INPUT PARAMETERS:
!
    ! Nki(nBins)            - number of particles per size bin in grid cell
    ! Nnuci                 - number of nucleation size particles per size bin in
    !                         grid cell
    ! Mnuci                 - mass of given species in nucleation pseudo-bin
    !                         (kg/grid cell)
    ! Mki(nBins, ICOMPHARD) - mass of a given species per size bin/grid cell
    ! Gci(icomp)            - amount (kg/grid cell) of all species present in the
    !                         gas phase except water
    ! H2SO4rate             - rate of H2SO4 chemical production [kg s^-1]
    ! dt                    - total model time step to be taken (s)
    REAL(fp) Nki(nBins), Mki(nBins, ICOMPHARD), Gci(ICOMPHARD)
    double precision H2SO4rate
    real(fp)             dti
!
! !OUTPUT PARAMETERS:
!
    ! Nkf, Mkf, Gcf  - same as above, but final values
    ! Nknuc,  Mknuc  - same as above, final values from just nucleation
    ! Nkcond, Mkcond - same as above, but final values from just condensation
    ! fn, fn1
    REAL(fp) Nkf(nBins), Mkf(nBins, ICOMPHARD), Gcf(ICOMPHARD)
    REAL(fp) Nknuc(nBins), Mknuc(nBins, ICOMPHARD)
    REAL(fp) Nkcond(nBins),Mkcond(nBins,ICOMPHARD)
    double precision fnavg        ! nucleation rate of clusters cm-3 s-1
    double precision fn1avg       ! formation rate of particles to first size bin cm-3 s-1
    REAL(fp)           BOXVOL, BOXMASS, TEMPTMS, RHTOMAS, PRES
    logical          errswitch    ! signal for error
    integer          lev          ! layer of the model
    REAL(fp)   surf_area
    REAL(fp)   ionrate
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    double precision dti_db
    integer          i,j,k,c      ! counters
    double precision fn           ! nucleation rate of clusters cm-3 s-1
    double precision fn1          ! formation rate of particles to first size bin cm-3 s-1
    double precision pi, R        ! pi and gas constant (J/mol K)
    double precision CSi,CSa      ! intial and average condensation sinks
    double precision CS1,CS2      ! guesses for condensation sink [s^-1]
    double precision CStest       ! guess for condensation sink
    REAL(fp)           Nk1(nBins), Mk1(nBins, ICOMPHARD), Gc1(ICOMPHARD)
    REAL(fp)           Nk2(nBins), Mk2(nBins, ICOMPHARD), Gc2(ICOMPHARD)
    REAL(fp)           Nk3(nBins), Mk3(nBins, ICOMPHARD), Gc3(ICOMPHARD)
    logical          nflg         ! returned from nucleation, says whether nucleation occurred or not
    double precision mcond,mcond1 ! mass to condense [kg]
    double precision tol          ! tolerance
    double precision eps          ! small number
    double precision sinkfrac(nBins) ! fraction of condensation sink coming from bin k
    double precision totmass      ! the total mass of H2SO4 generated during the timestep
    double precision tmass
    double precision CSch         ! fractional change in condensation sink
    double precision CSch_tol     ! tolerance in change in condensation sink
    double precision addt         ! adaptive timestep time
    double precision time_rem     ! time remaining
    integer          num_iter     ! number of iteration
    double precision sumH2SO4     ! used for finding average H2SO4 conc over timestep
    integer          iter         ! number of iteration
    double precision rnuc         ! critical radius [nm]
    double precision gasConc      ! gas concentration [kg]
    double precision mass_change  ! change in mass during nucleation
    double precision total_nh4_1,total_nh4_2
    double precision min_tstep    ! minimum timestep [s]
    integer          nuc_bin      ! the nucleation bin
    double precision sumfn, sumfn1 ! used for getting average nucleation rates
    logical          tempvar,  pdbg
    real(fp)           tnumb
!
! !DEFINED PARAMETERS:
!
    parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
    parameter(eps=1E-40)
    parameter(CSch_tol=0.01)
    parameter(min_tstep=1.0e+0_fp)

    !=================================================================
    ! COND_NUC begins here
    !=================================================================

    pdbg      = errswitch ! transfer the signal to print debug from outside
    errswitch = .false.   ! flag error to outide to terminate program.

    dti_db = dble(dti)

    ! Initialize values of Nkf, Mkf, Gcf, and time
    do j=1,ICOMPHARD
       Gc1(j)=Gci(j)
       Gcf(j)=Gci(j)
    enddo
    do k=1,nBins
       Nk1(k)=Nki(k)
       Nknuc(k)=Nki(k)
       Nkcond(k)=Nki(k)
       do j=1,ICOMPHARD
          Mk1(k,j)=Mki(k,j)
          Mknuc(k,j)=Mki(k,j)
          Mkcond(k,j)=Mki(k,j)
       enddo
    enddo

    ! Get initial condensation sink
    CS1 = 0.e+0_fp
    call getCondSink(Nk1,Mk1,srtso4,CS1,sinkfrac,surf_area,BOXVOL,TEMPTMS,PRES)
    if( pdbg) print*,'CS1', CS1
    !CS1 = max(CS1,eps)

    !Get initial H2SO4 concentration guess (assuming no nucleation)
    !Make sure that H2SO4 concentration doesn't exceed the amount generated
    !during that timestep (this will happen when the condensation sink is very low)

    ! get the steady state H2SO4 concentration
    call getH2SO4conc(Nk1, Mk1, H2SO4rate, CS1, Gc1(srtnh4), &
                      gasConc, ionrate, surf_area, &
                      BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
    if( pdbg) print*,'gasConc',gasConc
    Gc1(srtso4) = gasConc
    addt = min_tstep
    !addt = 3600.e+0_fp
    totmass = H2SO4rate*addt*96.e+0_fp/98.e+0_fp

    tempvar = pdbg

    !Get change size distribution due to nucleation with initial guess
    call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,fn,fn1,totmass,nuc_bin, &
                    addt, ionrate, surf_area, BOXVOL, BOXMASS, TEMPTMS, &
                    PRES, RHTOMAS, PDBG, lev)

    if(pdbg) then
       print*,'COND_NUC: Found an error at nucleation --> TERMINATE'
       errswitch = .true.
       return
    endif
    pdbg = tempvar !put the print debug switch back to pdbg
    !if(pdbg) call debugprint(Nk2, Mk2, 0,0,0,'After nucleation[1]')

    !print*,'after nucleation'
    !print*,'Nnuc1',Nnuc1
    !print*,'Nnuc2',Nnuc2
    !print*,'Mnuc1',Mnuc1
    !print*,'Mnuc2',Mnuc2

    mass_change = 0.e+0_fp

    do k=1,nBins
       mass_change = mass_change + (Mk2(k,srtso4)-Mk1(k,srtso4))
    enddo
    if( pdbg)  print*,'mass_change',mass_change

    mcond = totmass-mass_change ! mass of h2so4 to condense

    if( pdbg) print*,'after nucleation'
    if( pdbg)  print*,'totmass',totmass,'mass_change1',mass_change,'mcond',mcond
    if( pdbg)  print*,'cs1',CS1, Gc1(srtso4)

    if (mcond.lt.0.e+0_fp)then
       tmass = 0.e+0_fp
       do k=1,nBins
          do j=1,ICOMPHARD-2
             tmass = tmass + Mk2(k,j)
          enddo
       enddo
       !if (abs(mcond).gt.tmass*1.0D-8) then
       if (abs(mcond).gt.totmass*1.0e-8_fp) then
          if (-mcond.lt.Mk2(nuc_bin,srtso4)) then
             !if (CS1.gt.1.0D-5)then
             !   print*,'budget fudge 1 in cond_nuc'
             !endif
             tmass = 0.e+0_fp
             do j=1,ICOMPHARD-2
                tmass = tmass + Mk2(nuc_bin,j)
             enddo
             Nk2(nuc_bin) = Nk2(nuc_bin)*(tmass+mcond)/tmass
             Mk2(nuc_bin,srtso4) = Mk2(nuc_bin,srtso4) + mcond
             mcond = 0.e+0_fp
          else
             print*,'budget fudge 2 in cond_nuc'
             do k=2,nBins
                Nk2(k) = Nk1(k)
                Mk2(k,srtso4) = Mk1(k,srtso4)
             enddo
             Nk2(1) = Nk1(1)+totmass/sqrt(xk(1)*xk(2))
             Mk2(1,srtso4) = Mk1(1,srtso4) + totmass
             mcond = 0.e+0_fp
             !print*,'mcond < 0 in cond_nuc', mcond, totmass
             !stop
          endif
       else
          mcond = 0.e+0_fp
       endif
    endif

    !if (mcond.lt.0.e+0_fp)then
    !   print*,'mcond < 0 in cond_nuc', mcond
    !   stop
    !endif
    tmass = 0.e+0_fp
    do k=1,nBins
       do j=1,ICOMPHARD-2
          tmass = tmass + Mk2(k,j)
       enddo
    enddo
    if( pdbg)  print*, 'mcond',mcond,'tmass',tmass,'nuc',Nk2(1)-Nk1(1)
    tempvar = pdbg

    ! Get guess for condensation
    call ezcond(Nk2,Mk2,mcond,srtso4,Nk3,Mk3,surf_area, &
                BOXVOL, TEMPTMS, PRES, pdbg )

    if(pdbg) then
       print*,'COND_NUC: Found an error at EZCOND --> TERMINATE'
       errswitch = .true.
       return
    endif
    pdbg = tempvar
    ! if(pdbg) call debugprint(Nk3, Mk3, 0,0,0,'After EZCOND[1]')
    !print*,'after ezcond',Nk2,Nk3
    !jrp mcond1 = 0.e+0_fp
    !jrp do k=1,ibins
    !jrp    do j=1,icomp
    !jrp       mcond1 = mcond1 + (Mk3(k,j)-Mk2(k,j))
    !jrp    enddo
    !jrp enddo
    !print*,'mcond',mcond,'mcond1',mcond1

    Gc3(srtnh4) = Gc1(srtnh4)

    call eznh3eqm(Gc3,Mk3)
    call ezwatereqm(Mk3, RHTOMAS)

    ! check to see how much condensation sink changed
    call getCondSink(Nk3,Mk3,srtso4,CS2,sinkfrac,surf_area, &
                     BOXVOL,TEMPTMS, PRES)
    CSch = abs(CS2 - CS1)/CS1

    !if (CSch.gt.CSch_tol) then ! condensation sink didn't change much use whole timesteps
    ! get starting adaptive timestep to not allow condensationk sink
    ! to change that much
    ! Avoid div-by-zero (bmy, 1/28/14)
    IF ( ABS( CSch ) > 0e+0_fp ) THEN
       addt = addt*CSch_tol/CSch/2e+0_fp
    ELSE
       addt = 0e+0_fp
    ENDIF
    addt = min(addt,dti_db)
    addt = max(addt,min_tstep)

    time_rem = dti_db ! time remaining
    if( pdbg)    print*,'addt',addt,time_rem
    num_iter = 0
    sumH2SO4=0.e+0_fp
    sumfn = 0.e+0_fp
    sumfn1 = 0.e+0_fp
    ! do adaptive timesteps
    do while (time_rem .gt. 0.e+0_fp)
       num_iter = num_iter + 1
       if( pdbg) print*, 'iter', num_iter, ' addt', addt, 'time_rem', time_rem
       ! get the steady state H2SO4 concentration
       if (num_iter.gt.1)then ! no need to recalculate for first step
          call getH2SO4conc(Nk1, Mk1, H2SO4rate, CS1, Gc1(srtnh4), &
                            gasConc, ionrate, surf_area, &
                            BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
          Gc1(srtso4) = gasConc
       endif
       if( pdbg)    print*,'gasConc',gasConc

       sumH2SO4 = sumH2SO4 + Gc1(srtso4)*addt
       totmass = H2SO4rate*addt*96.e+0_fp/98.e+0_fp
       !call nucleation(Nk1,Mk1,Gc1,Nnuc1,Mnuc1,totmass,addt,Nk2, &
       !                Mk2,Gc2,Nnuc2,Mnuc2,nflg,lev)

       !Debug to see what goes in nucleation (win, 10/3/08)
       if(pdbg) then
          print*,'Temperature',TEMPTMS,'RH',RHTOMAS
          print*,'H2SO4',Gc1(srtso4)/boxvol*1000.e+0_fp/98.e+0_fp*6.022e+23_fp
          print*,'NH3ppt',Gc1(srtnh4)/17.e+0_fp/(boxmass/29.e+0_fp)*1e+12_fp
       endif

       tempvar = pdbg
       call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,fn,fn1,totmass, &
                       nuc_bin,addt, ionrate, surf_area, BOXVOL, BOXMASS, &
                       TEMPTMS, PRES, RHTOMAS, PDBG, lev)

       if(pdbg) then
          print*,'COND_NUC: Error at nucleation[2] --> TERMINATE'
          errswitch=.true.
          return
       endif
       pdbg = tempvar
       ! if(pdbg) call debugprint(Nk2, Mk2, 0,0,0, 'After nucleation[2]')
       !print*,'after nucleation iter'
       sumfn = sumfn + fn*addt
       sumfn1 = sumfn1 + fn1*addt

       !total_nh4_1 = Mnuc1(srtnh4)
       !total_nh4_2 = Mnuc2(srtnh4)
       !do i=1,ibins
       !   total_nh4_1 = total_nh4_1 + Mk1(i,srtnh4)
       !   total_nh4_2 = total_nh4_2 + Mk2(i,srtnh4)
       !enddo
       !print*,'total_nh4',total_nh4_1,total_nh4_2

       mass_change = 0.e+0_fp

       do k=1,nBins
          mass_change = mass_change + (Mk2(k,srtso4)-Mk1(k,srtso4))
       enddo
       if( pdbg)    print*,'mass_change2',mass_change

       mcond = totmass-mass_change ! mass of h2so4 to condense

       !print*,'after nucleation'
       !print*,'totmass',totmass,'mass_change',mass_change,'mcond',mcond

       !print*,'2 mass_change',mass_change,mcond,totmass
       !print*,'2 cs1',CS1, Gc1(srtso4)

       if (mcond.lt.0.e+0_fp)then
          tmass = 0.e+0_fp
          do k=1,nBins
             do j=1,ICOMPHARD-2
                tmass = tmass + Mk2(k,j)
             enddo
          enddo
          !if (abs(mcond).gt.tmass*1.0D-8) then
          if (abs(mcond).gt.totmass*1.0e-8_fp) then
             if (-mcond.lt.Mk2(nuc_bin,srtso4)) then
                !if (CS1.gt.1.0D-5)then
                !   print*,'budget fudge 1 in cond_nuc'
                !endif
                tmass = 0.e+0_fp
                do j=1,ICOMPHARD-2
                   tmass = tmass + Mk2(nuc_bin,j)
                enddo
                Nk2(nuc_bin) = Nk2(nuc_bin)*(tmass+mcond)/tmass
                Mk2(nuc_bin,srtso4) = Mk2(nuc_bin,srtso4) + mcond
                mcond = 0.e+0_fp
             else
                print*,'budget fudge 2 in cond_nuc'
                do k=2,nBins
                   Nk2(k) = Nk1(k)
                   Mk2(k,srtso4) = Mk1(k,srtso4)
                enddo
                Nk2(1) = Nk1(1)+totmass/sqrt(xk(1)*xk(2))
                Mk2(1,srtso4) = Mk1(1,srtso4) + totmass
                print*,'mcond < 0 in cond_nuc', mcond, totmass
                mcond = 0.e+0_fp
                ! should I stop or not?? (win, 10/4/08)
                !stop
                ! change from stop here to stop outside with more info (win, 10/4/08)
                print*,'COND_NUC: --> TERMINATE'
                !10/4/08 errswitch = .true.
                !10/4/08 return
             endif
          else
             mcond = 0.e+0_fp
          endif
       endif

       do k=1,nBins
          Nknuc(k) = Nknuc(k)+Nk2(k)-Nk1(k)
          do j=1,ICOMPHARD-2
             Mknuc(k,j)=Mknuc(k,j)+Mk2(k,j)-Mk1(k,j)
          enddo
       enddo

       !Gc2(srtnh4) = Gc1(srtnh4)
       !call eznh3eqm(Gc2,Mk2,Mnuc2)
       !call ezwatereqm(Mk2,Mnuc2)

       !call getCondSink(Nk2,Mk2,Nnuc2,Mnuc2,srtso4,CStest,sinkfrac)

       ! Before entering ezcond, check if there's enough aerosol to
       ! condense onto. After several iteration in the case with high
       ! H2SO4 amount but little existing aerosol and also lack the conditions
       ! for nucleation, the whole size distribution is grown out of our
       ! tracked size bins, so let's exit the loop if there is no aerosol
       ! to condense onto anymore. (win, 10/4/08)
       tmass = 0.e+0_fp
       tnumb = 0.e+0_fp
       do k=1,nBins
          tnumb = tnumb + Nk2(k)
          do j=1,ICOMPHARD-2
             tmass = tmass + Mk2(k,j)
          enddo
       enddo

       if( (tmass+mcond)/tnumb  > Xk(nBins) ) then
          if( .not. SPINUP(10.0) ) then
             print*,'Not enough aerosol for condensation!'
             print*,'  Exiting COND_NUC iteration with '
             print*,time_rem,'sec remaining time'
          endif

          Gc3(srtnh4)=Gc2(srtnh4)
          do k=1,nBins
             Nk3(k)=Nk2(k)
             do j=1,ICOMPHARD
                Mk3(k,j)=Mk2(k,j)
             enddo
          enddo
          goto 100
       endif

       tempvar = pdbg

       call ezcond(Nk2,Mk2,mcond,srtso4,Nk3,Mk3,surf_area, &
                   BOXVOL, TEMPTMS, PRES, pdbg)
       do k=1,nBins
          Nkcond(k) = Nkcond(k)+Nk3(k)-Nk2(k)
          do j=1,ICOMPHARD-2
             Mkcond(k,j)=Mkcond(k,j)+Mk3(k,j)-Mk2(k,j)
          enddo
       enddo
       Gc3(srtnh4) = Gc1(srtnh4)

       if(pdbg) then
          print*,'COND_NUC: Error at EZCOND[2] --> TERMINATE'
          errswitch=.true.
          return
       endif
       pdbg = tempvar

       !if(pdbg) call debugprint(Nk3, Mk3, 0,0,0,'After EZCOND[2]')

       if( pdbg)    print*,'after ezcond iter'
       call eznh3eqm(Gc3,Mk3)
       call ezwatereqm(Mk3, RHTOMAS)

       ! check to see how much condensation sink changed
       call getCondSink(Nk3,Mk3,srtso4,CS2,sinkfrac,surf_area, &
                        BOXVOL,TEMPTMS, PRES)

       time_rem = time_rem - addt
       if (time_rem .gt. 0.e+0_fp) then
          CSch = abs(CS2 - CS1)/CS1
          !jrp if (CSch.lt.0.e+0_fp) then
          !jrp    print*,''
          !jrp    print*,'CSch LESS THAN ZERO!!!!!', CS1,CStest,CS2
          !jrp    print*,'Nnuc',Nnuc1,Nnuc2
          !jrp    print*,''
          !jrp
          !jrp    addt = min(addt,time_rem)
          !jrp else

          ! Allow adaptive timestep to change
          ! Avoid div-by-zero error
          IF ( ABS( CSch ) > 0e+0_fp ) THEN
             addt = min(addt*CSch_tol/CSch,addt*1.5e+0_fp)
          ELSE
             addt = 0e+0_fp
          ENDIF

          ! allow adaptive timestep to change again
          addt = min(addt,time_rem)
          addt = max(addt,min_tstep)
          !jrp endif
          if( pdbg)     print*,'CS1',CS1,'CS2',CS2
          CS1 = CS2
          Gc1(srtnh4)=Gc3(srtnh4)
          do k=1,nBins
             Nk1(k)=Nk3(k)
             do j=1,ICOMPHARD
                Mk1(k,j)=Mk3(k,j)
             enddo
          enddo
       endif
    enddo ! while loop

100 continue

    Gcf(srtso4)=sumH2SO4/dti_db
    fnavg = sumfn/dti_db
    fn1avg = sumfn1/dti_db
    if( pdbg)    print*,'AVERAGE GAS CONC',Gcf(srtso4)

    !jrp else
    !jrp    num_iter = 1
    !jrp    Gcf(srtso4)=Gc1(srtso4)
    !jrp endif

    if( pdbg) print*, 'cond_nuc num_iter =', num_iter
    !T0M(1,1,1,3) = double(num_iter) ! store iterations here

    ! if(pdbg) call debugprint(Nk3, Mk3, 0,0,0,'End of COND_NUC')

    do k=1,nBins
       Nkf(k)=Nk3(k)
       do j=1,ICOMPHARD
          Mkf(k,j)=Mk3(k,j)
       enddo
    enddo
    Gcf(srtnh4)=Gc3(srtnh4)

    return

  END SUBROUTINE COND_NUC
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nucleation
!
! !DESCRIPTION: This subroutine calls the Vehkamaki 2002 and Napari 2002
!  nucleation parameterizations and gets the binary and ternary nucleation
!  rates. The number of particles added to the first size bin is calculated
!  by comparing the growth rate of the new particles to the coagulation sink.
!  WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II'
!  Introduce to GEOS-Chem by Win Trivitayanurak (win, 9/30/08)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NUCLEATION(Nki,Mki,Gci,Nkf,Mkf,Gcf,fn,fn1,totsulf, &
                        nuc_bin,dt,ionrate, surf_area, BOXVOL, BOXMASS, &
                        TEMPTMS, PRES, RHTOMAS, pdbg,lev)
!
! !USES:
!
    USE ERROR_MOD,      ONLY : ERROR_STOP, IT_IS_NAN
!
! !INPUT PARAMETERS:
!
    !Initial values of
    !=================
    !Nki(ibins) - number of particles per size bin in grid cell
    !Mki(ibins, icomp) - mass of a given species per size bin/grid cell
    !Gci(icomp-1) - amount (kg/grid cell) of all species present in the
    !               gas phase except water
    !dt - total model time step to be taken (s)
    double precision Nki(nBins), Mki(nBins, ICOMPHARD), Gci(ICOMPHARD-1)
    REAL(fp), INTENT(IN)       :: BOXVOL,  BOXMASS, TEMPTMS
    REAL(fp), INTENT(IN)       :: PRES,    RHTOMAS
!
! !OUTPUT PARAMETERS:
!
    !Nkf, Mkf, Gcf - same as above, but final values
    !fn, fn1
    double precision Nkf(nBins), Mkf(nBins, ICOMPHARD), Gcf(ICOMPHARD-1)
    integer j,i,k
    double precision totsulf
    integer nuc_bin
    double precision dt
    double precision fn       ! nucleation rate of clusters cm-3 s-1
    double precision fn1      ! formation rate of particles to first size bin cm-3 s-1

    LOGICAL  PDBG             ! Signal print for debug
    integer lev ! layer of model

    REAL(fp)                     ionrate
    REAL(fp)                     surf_area
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    double precision nh3ppt   ! gas phase ammonia in pptv
    double precision h2so4    ! gas phase h2so4 in molec cc-1
    double precision rnuc     ! critical nucleation radius [nm]
    double precision gtime    ! time to grow to first size bin [s]
    double precision ltc, ltc1, ltc2 ! coagulation loss rates [s-1]
    double precision Mktot    ! total mass in bin
    double precision neps
    double precision meps
    double precision density  ! density of particle [kg/m3]
    double precision pi
    double precision frac     ! fraction of particles growing into first size bin
    double precision d1,d2    ! diameters of particles [m]
    double precision mp       ! mass of particle [kg]
    double precision mold     ! saved mass in first bin
    double precision mnuc     !mass of nucleation
    double precision sinkfrac(nBins) ! fraction of loss to different size bins
    double precision nadd     ! number to add
    double precision CS       ! kerminan condensation sink [m-2]
    double precision Dpmean   ! the number wet mean diameter of the existing aerosol
    double precision Dp1      ! the wet diameter of bin 1
    double precision dens1    ! density in bin 1 [kg m-3]
    double precision GR       ! growth rate [nm hr-1]
    double precision gamma,eta ! used in kerminen 2004 parameterzation
    double precision drymass,wetmass,WR
    double precision fn_c     ! barrierless nucleation rate
    double precision h1,h2,h3,h4,h5,h6
    double precision dum1,dum2,dum3,dum4   ! dummy variables
    double precision rhin,tempin ! rel hum in

    LOGICAL ERRORSWITCH
!
! !DEFINED PARAMETERS:
!
    parameter (neps=1E8, meps=1E-8)
    parameter (pi=3.14159)

    !=================================================================
    ! NUCLEATION begins here
    !=================================================================

    errorswitch = .false.

    h2so4 = Gci(srtso4)/boxvol*1000.e+0_fp/98.e+0_fp*6.022e+23_fp
    nh3ppt = Gci(srtnh4)/17.e+0_fp/(boxmass/29.e+0_fp)*1e+12_fp* &
             PRES/101325.*273./TEMPTMS ! corrected for pressure (because this should be concentration)

    fn = 0.e+0_fp
    fn1 = 0.e+0_fp
    rnuc = 0.e+0_fp
    gtime = 0.e+0_fp
    nuc_bin = 1 ! added by Pengfei Liu,initialize  nuc_bin value
    ! if requirements for nucleation are met, call nucleation subroutines
    ! and get the nucleation rate and critical cluster size
    if (h2so4.gt.1.e+4_fp) then
       if (nh3ppt.gt.0.1.and.tern_nuc.eq.1) then
          call napa_nucl(TEMPTMS,RHTOMAS,h2so4,nh3ppt,fn,rnuc) !ternary nuc
          if (ion_nuc.eq.1.and.ionrate.ge.1.e+0_fp) then
             !call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS, &
             !              h1,h2,h3,h4,h5,h6)
          else
             h1=0.e+0_fp
          endif
          if (h1.gt.fn)then
             fn=h1
             rnuc=h5
          endif
       elseif (bin_nuc.eq.1) then
          call vehk_nucl(TEMPTMS,RHTOMAS,h2so4,fn,rnuc) !binary nuc
          if ((ion_nuc.eq.1).and.(ionrate.ge.1.e+0_fp)) then
          !   call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS, &
          !                 h1,h2,h3,h4,h5,h6)
          else
             h1=0.e+0_fp
          endif
          if (h1.gt.fn)then
             fn=h1
             rnuc=h5
          endif
          if (fn.lt.1.0e-6_fp)then
             fn = 0.e+0_fp
          endif
       elseif ((ion_nuc.eq.1).and.(ionrate.ge.1.e+0_fp)) then
          !call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS, &
          !              h1,h2,h3,h4,h5,h6)
          fn=h1
          rnuc=h5
       elseif(ion_nuc.eq.2) then
          ! Yu Ion nucleation
          !! first we need to calculate the available surface area
          !surf_area = 0.e+0_fp
          !do k=1, ibins
          !   if (Nki(k) .gt. Neps) then
          !      Mktot=0.e+0_fp
          !      do j=1,icomp
          !         Mktot=Mktot+Mki(k,j)
          !      enddo
          !      mp=Mktot/Nki(k)
          !      density=aerodens(Mki(k,srtso4),0.e+0_fp, &
          !                       Mki(k,srtnh4),0.e+0_fp,Mki(k,srth2o))  ! assume bisulfate
          !      ! diameter = ((mass/density)*(6/pi))**(1/3)
          !      d2 = 1.D6*((mp/density)*(6.D0/pi))**(1.D0/3.D0) ! (micrometers)
          !      ! surface area per particle = pi*diameter**2
          !      surf_area = surf_area + 1.D-6*(Nki(k)/boxvol)* &
          !                              pi*(d2**2.D0) ! (um2 cm-2)
          !   endif
          !enddo
          rhin=dble(RHTOMAS*100.e+0_fp)
          tempin=dble(TEMPTMS)

          !call YUJIMN(h2so4, rhin, tempin, ionrate, surf_area, &
          !            fn, dum1, rnuc, dum2)
          fn=0.
          rnuc=1E-9
       endif
       !if((act_nuc.eq.1).and.(lev.le.7))then
          !call bl_nucl(h2so4,fn,rnuc)
       !endif
       call cf_nucl(TEMPTMS,RHTOMAS,h2so4,nh3ppt,fn_c) ! use barrierless nucleation as a max
       fn = min(fn,fn_c)
       !if (fn.gt.1.0)then
       !   print*, 'fn',fn
       !   print*, 'Yu Yes!'
       !   print*, 'ionrate',ionrate
       !   print*, 'surf_area',surf_area
       !endif
    endif

    if (pdbg) then
       if( bin_nuc == 1 ) then
          print *, 'BINARY cluster form rate : fn',fn
       else
          print *, 'TERNARY cluster form rate: fn',fn
       endif
    endif

    ! if nucleation occured, see how many particles grow to join the first size
    ! section
    if (fn.gt.0.e+0_fp) then

       if(pdbg) print*,'Nki',Nki
       if(pdbg) print*,'Mki',Mki

       call getCondSink_kerm(Nki,Mki,CS,Dpmean,Dp1,dens1,BOXVOL,TEMPTMS,PRES)

       if(pdbg) print*,'CS',CS,'Dpmean',Dpmean,'Dp1',Dp1,'dens1',dens1

       d1 = rnuc*2.e+0_fp*1e-9_fp
       drymass = 0.e+0_fp
       do j=1,ICOMPHARD-2
          drymass = drymass + Mki(1,j)
       enddo
       wetmass = 0.e+0_fp
       do j=1,ICOMPHARD
          wetmass = wetmass + Mki(1,j)
       enddo

       ! to prevent division by zero (win, 10/1/08)
       if(drymass == 0.e+0_fp) then
          WR = 1.e+0_fp
       else
          WR = wetmass/drymass
       endif

       if(pdbg) print*,'rnuc',rnuc,'WR',WR
       if(pdbg) print*,'d1',d1,'Gci(srtso4)',Gci(srtso4),&
                       'TEMP',temptms,'boxvol',boxvol

       if( IT_IS_NAN( Gci(srtso4) )) then
          print*,'rnuc',rnuc,'WR',WR
          print*,'d1',d1,'Gci(srtso4)',Gci(srtso4)
          call ERROR_STOP('Found NaN in Gci','nucleation')
       endif
       ! print*,'[nucleation] Gci',Gci
       call getGrowthTime(d1,Dp1,Gci(srtso4)*WR,TEMPTMS, &
                          boxvol,dens1,gtime)
       if (pdbg) print*,'gtime',gtime

       GR = (Dp1-d1)*1e+9_fp/gtime*3600.e+0_fp ! growth rate, nm hr-1

       gamma = 0.23e+0_fp*(d1*1.0e+9_fp)**(0.2e+0_fp)* &
               (Dp1*1.0d9/3.e+0_fp)**0.075e+0_fp* &
               (Dpmean*1.0e+9_fp/150.e+0_fp)** &
               0.048e+0_fp*(dens1*1.0e-3_fp)** &
               (-0.33e+0_fp)*(TEMPTMS/293.e+0_fp) ! equation 5 in kerminen
       eta = gamma*CS/GR

       if (Dp1.gt.d1)then
          fn1 = fn*exp(eta/(Dp1*1.0e+9_fp)-eta/(d1*1.0e+9_fp))
       else
          fn1 = fn
       endif

       if (pdbg) print*,'eta',eta,'Dp1',Dp1,'d1',d1,'fn1',fn1

       mnuc = sqrt(xk(1)*xk(2))

       nadd = fn1

       nuc_bin = 1

       mold = Mki(nuc_bin,srtso4)
       Mkf(nuc_bin,srtso4) = Mki(nuc_bin,srtso4)+nadd*mnuc*boxvol*dt
       Nkf(nuc_bin) = Nki(nuc_bin)+nadd*boxvol*dt

       Gcf(srtso4) = Gci(srtso4) ! - (Mkf(nuc_bin,srtso4)-mold)
       Gcf(srtnh4) = Gci(srtnh4)

       if (pdbg) then
          print*, 'nadd',nadd
          print *,'Mass add to bin',nuc_bin,'=',nadd*mnuc*boxvol*dt
          print *,'Number added',nadd*boxvol*dt
          print *,'Gcf(srtso4)',Gcf(srtso4)
          print *,'Gcf(srtnh4)',Gcf(srtnh4)
       endif

       do k=1,nBins
          if (k .ne. nuc_bin)then
             Nkf(k) = Nki(k)
             do i=1,ICOMPHARD
                Mkf(k,i) = Mki(k,i)
             enddo
          else
             do i=1,ICOMPHARD
                if (i.ne.srtso4) then
                   Mkf(k,i) = Mki(k,i)
                endif
             enddo
          endif
       enddo

       do k=1,nBins
          if (Nkf(k).lt.1.e+0_fp) then
             Nkf(k) = 0.e+0_fp
             do j=1,ICOMPHARD
                Mkf(k,j) = 0.e+0_fp
             enddo
          endif
       enddo
       !print *, 'mnfix in tomas_mod:2679'

       call mnfix(Nkf,Mkf, ERRORSWITCH)
       pdbg = errorswitch ! carry the error signal from mnfix to outside
       if (errorswitch) print*,'NUCLEATION: Error after mnfix'

       ! there is a chance that Gcf will go less than zero because we are
       ! artificially growing particles into the first size bin.
       ! don't let it go less than zero.

    else

       do k=1,nBins
          Nkf(k) = Nki(k)
          do i=1,ICOMPHARD
             Mkf(k,i) = Mki(k,i)
          enddo
       enddo

    endif

    pdbg = errorswitch        ! carry the error signal from mnfix to outside

    RETURN

  END SUBROUTINE NUCLEATION
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: multicoag
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MULTICOAG( DT, Nk, Mk, BOXVOL, PRES, TEMPTMS, PDBG )
!
! !INPUT PARAMETERS:
!
    REAL(fp),    INTENT(IN)     :: DT                ! Time step (s)
    REAL(fp),    INTENT(IN)     :: PRES
    REAL(fp),    INTENT(IN)     :: TEMPTMS
    REAL(fp),    INTENT(IN)     :: BOXVOL
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(fp),  INTENT(INOUT)  :: Nk(nBins)
    REAL(fp),  INTENT(INOUT)  :: Mk(nBins, ICOMPHARD)
    LOGICAL,   INTENT(INOUT)  :: PDBG              ! For signalling print debug
!
! !REMARKS:
!  Some key variables
!  kij represents the coagulation coefficient (cm3/s) normalized by the
!      volume of the GCM grid cell (boxvol, cm3) such that its units are (s-1)
!  dNdt and dMdt are the rates of change of Nk and Mk.  xk contains
!     the mass boundaries of the size bins.  xbar is the average mass
!     of a given size bin (it varies with time in this algorithm).  phi
!     and eff are defined in the reference, equations 13a and b.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER     :: K, J, I, JJ
    REAL(fp)    :: dNdt(nBins), dMdt(nBins,ICOMPHARD-2)
    REAL(fp)    :: xbar(nBins), phi(nBins), eff(nBins)
    REAL*4      :: kij(nBins,nBins)
    REAL*4      :: Dpk(nBins)             !diameter (m) of particles in bin k
    REAL*4      :: Dk(nBins)              !Diffusivity (m2/s) of bin k particles
    REAL*4      :: ck(nBins)              !Mean velocity (m/2) of bin k particles
    REAL*4      :: olddiff                !used to iterate to find diffusivity
    REAL*4      :: density                !density (kg/m3) of particles
    REAL*4      :: mu                     !viscosity of air (kg/m s)
    REAL*4      :: mfp                    !mean free path of air molecule (m)
    REAL*4      :: Kn                     !Knudsen number of particle
    REAL(fp)    :: mp                     !particle mass (kg)
    REAL*4      :: beta                   !correction for coagulation coeff.
    !      real(fp), external ::   aerodens  !<tmp> try change to double precision (win, 1/4/06)

    !temporary summation variables
    REAL(fp)    :: k1m(ICOMPHARD-2),k1mx(ICOMPHARD-2)
    REAL(fp)    :: k1mx2(ICOMPHARD-2)
    REAL(fp)    :: k1mtot,k1mxtot
    REAL(fp)    :: sk2mtot, sk2mxtot
    REAL(fp)    :: sk2m(ICOMPHARD-2), sk2mx(ICOMPHARD-2)
    REAL(fp)    :: sk2mx2(ICOMPHARD-2)
    REAL(fp)    :: High_in
    REAL(fp)    :: mtotal, mktot

    REAL*4      :: zeta                      !see reference, eqn 6
    REAL*4      :: tlimit, dtlimit, itlimit  !fractional change in M/N allowed in one time step
    REAL*4      :: dts  !internal time step (<dt for stability)
    REAL*4      :: tsum !time so far
    REAL(fp)    :: Neps !minimum value for Nk
!dbg
    character*12 limit        !description of what limits time step

    REAL(fp)    :: mi, mf   !initial and final masses

#if defined(TOMAS12) || defined(TOMAS15)
    parameter(zeta=1.28125 , dtlimit=0.25, itlimit=10.)
#else
    parameter(zeta=1.0625, dtlimit=0.25, itlimit=10.)
#endif
    REAL*4      ::pi, kB  !kB is Boltzmann constant (J/K)
    REAL*4      ::R       !gas constant (J/ mol K)
    parameter (pi=3.141592654, kB=1.38e-23, R=8.314, Neps=1.0e-3)

    REAL(fp)      :: M_NH4

    LOGICAL     :: ERRSPOT

    !sfarina
1   format(16E15.3)

    !=================================================================
    ! MULTICOAG begins here!
    !=================================================================
    tsum = 0.0

    ! If any Nk are zero, then set them to a small value to avoid division by zero
    do k=1,nBins
       if (Nk(k) .lt. Neps) then
          Nk(k)=Neps
#if defined(TOMAS12) || defined(TOMAS15)
          Mk(k,srtso4)=Neps*sqrt( xk(k)*xk(k+1) ) !make the added particles SO4
#else
          Mk(k,srtso4)=Neps*1.4e+0_fp*xk(k) !make the added particles SO4
#endif
       endif
    enddo

    ! Calculate air viscosity and mean free path
    mu=2.5277e-7*temptms**0.75302
    mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*R*temptms)))  !S&P eqn 8.6

    !<temp>
    !write(6,*)'+++ Nk(1:30)    =',Nk(1:30)
    !write(6,*)'+++ Mk(1:30,SO4)=',Mk(1:30,srtso4)
    !write(6,*)'+++ Mk(1:30,H2O)=',Mk(1:30,srth2o)
    !if (pdbg) call debugprint(Nk,Mk,0,0,0,'Inside MULTICOAG')
    ! Calculate particle sizes and diffusivities
    do k=1,nBins

       IF ( SRTNH4 > 0 ) THEN
          M_NH4 = Mk(k,SRTNH4)
       ELSE
          M_NH4 = 0.1875e+0_fp*Mk(k,srtso4)  !assume bisulfate
       ENDIF
       !tmp write(6,*)'+++ multicoag:  Mk(',k,'srtso4)=',Mk(k,srtso4)
       !density=aerodens(Mk(k,srtso4),0.e+0_fp, M_NH4,        &
       !        Mk(k,srtnacl), Mk(k,srtecil), Mk(k,srtecob),  &
       !        Mk(k,srtocil), Mk(k,srtocob), Mk(k,srtdust),  &
       !        Mk(k,srth2o))     !use Mk for sea salt mass(win, 4/18/06)
      density=aerodens(Mk(k,srtso4),0.e+0_fp, M_NH4,        &
               0.e+0_fp,0.e+0_fp,0.e+0_fp,  &
               0.e+0_fp,0.e+0_fp,0.e+0_fp, &
               Mk(k,srth2o))
       !Update mp calculation to include all species (win, 4/18/06)

       !prior to 9/26/08 (win)
       !Mktot=0.1875e+0_fp*Mk(k,srtso4) !start with NH4 mass

       Mktot = M_NH4         ! start with ammonium (win, 9/26/08)
       Mktot = Mktot + Mk(k,srth2o) ! then include water

       do j=1, ICOMPHARD-2
          Mktot=Mktot+Mk(k,j)
       enddo
       mp=Mktot/Nk(k)
       Dpk(k)=((mp/density)*(6./pi))**(0.333)
       Kn=2.0*mfp/Dpk(k)                            !S&P Table 12.1
       Dk(k)=kB*temptms/(3.0*pi*mu*Dpk(k)) &        !S&P Table 12.1
         *((5.0+4.0*Kn+6.0*Kn**2+18.0*Kn**3)/(5.0-Kn+(8.0+pi)*Kn**2))
       ck(k)=sqrt(8.0*kB*temptms/(pi*mp))           !S&P Table 12.1
    enddo

    ! Calculate coagulation coefficients
    do i=1,nBins
       do j=1,nBins
          Kn=4.0*(Dk(i)+Dk(j)) &
             /(sqrt(ck(i)**2+ck(j)**2)*(Dpk(i)+Dpk(j))) !S&P eqn 12.51
          beta=(1.0+Kn)/(1.0+2.0*Kn*(1.0+Kn))          !S&P eqn 12.50
          !This is S&P eqn 12.46 with non-continuum correction, beta
          kij(i,j)=2.0*pi*(Dpk(i)+Dpk(j))*(Dk(i)+Dk(j))*beta
          kij(i,j)=kij(i,j)*1.0e+6_fp/boxvol  !normalize by grid cell volume
       enddo
    enddo

10  continue     !repeat process here if multiple time steps are needed

    if(pdbg) print*,'In the time steps loop +++++++++++++'

    ! Calculate xbar, phi and eff
#if defined(TOMAS12) || defined(TOMAS15)
    do k=1,nBins

       xbar(k)=0.0
       do j=1,ICOMPHARD-2
          xbar(k)=xbar(k)+Mk(k,j)/Nk(k)            !eqn 8b
       enddo
       if(k.lt.nBins-1)then !from 1 to 10 bins

          eff(k)=2./9.*Nk(k)/xk(k) *(4.-xbar(k)/xk(k)) !eqn 4 in tzivion 1999
          phi(k)=2./9.*Nk(k)/xk(k) *(xbar(k)/xk(k)-1.) !eqn 4 in tzivion 1999

          !Constraints in equation 15
          if (xbar(k) .lt. xk(k)) then
             eff(k)=2./3.*Nk(k)/xk(k)
             phi(k)=0.0

          else if (xbar(k) .gt. xk(k+1)) then
             phi(k)=2./3.*Nk(k)/xk(k)
             eff(k)=0.0
          endif
       else                      ! from 11 bins to 12 bins
          eff(k)=2./31./31.*Nk(k)/xk(k) &
                 *(32.-xbar(k)/xk(k)) !eqn 4 in tzivion 1999
          phi(k)=2./31./31.*Nk(k)/xk(k) &
                 *(xbar(k)/xk(k)-1.) !eqn 4 in tzivion 1999

          !Constraints in equation 15
          if (xbar(k) .lt. xk(k)) then
             eff(k)=2./31.*Nk(k)/xk(k)
             phi(k)=0.0

          else if (xbar(k) .gt. xk(k+1)) then
             phi(k)=2./31.*Nk(k)/xk(k)
             eff(k)=0.0
          endif
       endif

    enddo

#else
    do k=1,nBins

       xbar(k)=0.0
       do j=1,ICOMPHARD-2
          xbar(k)=xbar(k)+Mk(k,j)/Nk(k)            !eqn 8b
       enddo

       eff(k)=2.*Nk(k)/xk(k)*(2.-xbar(k)/xk(k))    !eqn 13a
       phi(k)=2.*Nk(k)/xk(k)*(xbar(k)/xk(k)-1.)    !eqn 13b

       !Constraints in equation 15
       if (xbar(k) .lt. xk(k)) then
          eff(k)=2.*Nk(k)/xk(k)
          phi(k)=0.0
       else if (xbar(k) .gt. xk(k+1)) then
          phi(k)=2.*Nk(k)/xk(k)
          eff(k)=0.0
       endif
    enddo
#endif

    ! Necessary initializations
    sk2mtot=0.0
    sk2mxtot=0.0
    do j=1,ICOMPHARD-2
       sk2m(j)=0.0
       sk2mx(j)=0.0
       sk2mx2(j)=0.0
    enddo

    ! Calculate rates of change for Nk and Mk
    do k=1,nBins

       !Initialize to zero
       do j=1,ICOMPHARD-2
          k1m(j)=0.0
          k1mx(j)=0.0
          k1mx2(j)=0.0
       enddo
       High_in=0.0
       k1mtot=0.0
       k1mxtot=0.0

       !Calculate sums
#if defined(TOMAS12) || defined(TOMAS15)
       do j=1,ICOMPHARD-2
          if (k .gt. 1.and.k.lt.nBins) then
             do i=1,k-1
                k1m(j)=k1m(j)+kij(k,i)*Mk(i,j)
                k1mx(j)=k1mx(j)+kij(k,i)*Mk(i,j)*xbar(i)*zeta
                k1mx2(j)=k1mx2(j)+kij(k,i)*Mk(i,j)*xbar(i)**2.*zeta**3.
             enddo
          elseif(k.eq.nBins)then
             k1m(j)= sk2m(j)+kij(k,k-1)*Mk(k-1,j)
             k1mx(j)=sk2mx(j)+kij(k,k-1)*Mk(k-1,j)*xbar(k-1)*4.754
             k1mx2(j)=sk2mx2(j)+kij(k,k-1)*Mk(k-1,j)*xbar(k-1)**2.*107.4365
          endif
          k1mtot=k1mtot+k1m(j)
          k1mxtot=k1mxtot+k1mx(j)
       enddo
#else
       do j=1,ICOMPHARD-2
          if (k .gt. 1) then
             do i=1,k-1
                k1m(j)=k1m(j)+kij(k,i)*Mk(i,j)
                k1mx(j)=k1mx(j)+kij(k,i)*Mk(i,j)*xbar(i)
                k1mx2(j)=k1mx2(j)+kij(k,i)*Mk(i,j)*xbar(i)**2
             enddo
          endif
          k1mtot=k1mtot+k1m(j)
          k1mxtot=k1mxtot+k1mx(j)
       enddo
#endif

       if (k .lt. nBins) then
          do i=k+1,nBins
             High_in=High_in+Nk(i)*kij(k,i)
          enddo
       endif

       !Calculate rates of change
#if defined(TOMAS12) || defined(TOMAS15)
       if(k.lt.nBins-1)then

          dNdt(k)= -Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.125 &
                   -(phi(k)*k1mtot+(eff(k)-phi(k))/6./xk(k)*k1mxtot) &
                   -kij(k,k)*(phi(k)/3.*xbar(k)*Nk(k)+(eff(k)-phi(k))/18. &
                   /xk(k)*zeta*xbar(k)*xbar(k)*Nk(k))

          if (k .gt. 1) then
             !yhl Nk*low_in changes to -0.5*Kij*Nk**2.
             dNdt(k)=dNdt(k)+0.625*kij(k-1,k-1)*Nk(k-1)**2 &
                     +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/6./xk(k-1) &
                     *sk2mxtot) &
                     +kij(k-1,k-1)*(phi(k-1)/3.*xbar(k-1)*Nk(k-1)+(eff(k-1) &
                     -phi(k-1))/18./xk(k-1)*zeta*xbar(k-1)*xbar(k-1) &
                     *Nk(k-1))

          endif

          do j=1,ICOMPHARD-2

             dMdt(k,j)= Nk(k)*k1m(j)-Mk(k,j)*High_in & ! !term5,term6
                        -(phi(k)*xk(k+1)*k1m(j)+ &
                        (eff(k)+2.*phi(k))/6.*k1mx(j) &
                        +(eff(k)-phi(k))/6./xk(k)*k1mx2(j)) & ! term3
                        - kij(k,k)*Nk(k)*Mk(k,j)/3. & ! I assume 1/2Nk and 2/3Mk for half bin
                        - kij(k,k)*(phi(k)*xk(k+1)*Mk(k,j)/3. &
                        +(eff(k)+2.*phi(k))/6.*zeta*xbar(k)*Mk(k,j)/3. &
                        +(eff(k)-phi(k))/6./xk(k)*zeta**3.*xbar(k)**2. &
                        *Mk(k,j)/3.)

             !yhl  Term9(-kij(k,k)*Nk(k)*Mk(k,j)) is cancled out by term6 (k)
             if (k .gt. 1) then
                dMdt(k,j)=dMdt(k,j) &
                          +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1) &
                          +2.*phi(k-1))/6.*sk2mx(j) &
                          +(eff(k-1)-phi(k-1))/6./xk(k-1)*sk2mx2(j)) & !term1
                          +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)/3. &
                          +kij(k-1,k-1)*(phi(k-1)*xk(k)*Mk(k-1,j)/3. &
                          +(eff(k-1)+2.*phi(k-1))/6.*zeta &
                          *xbar(k-1)*Mk(k-1,j)/3.+(eff(k-1)-phi(k-1))/6. &
                          /xk(k-1)*zeta**3.*xbar(k-1)**2.*Mk(k-1,j)/3.)
             endif
          enddo
       else if (k.eq.nBins-1)then

          dNdt(k)=0.625*kij(k-1,k-1)*Nk(k-1)**2 &
                  +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/6./xk(k-1) &
                  *sk2mxtot) &
                  +kij(k-1,k-1)*xbar(k-1)*Nk(k-1)/3.*(phi(k-1) &
                  +(eff(k-1) &
                  -phi(k-1))/6./xk(k-1)*zeta*xbar(k-1))

          !yhl updated the following
          dNdt(k)=dNdt(k)-Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.02 &
                  -(phi(k)*k1mtot+(eff(k)-phi(k))/62./xk(k)*k1mxtot) &
                  -kij(k,k)*xbar(k)*Nk(k)*0.484*(phi(k)+(eff(k) &
                  -phi(k))/62./xk(k)*4.754*xbar(k))

          !yhl I am not sure how it bring 0.5*kij(k-1,k-1)*Nk(k-1)**2 here. But
          !yhl It results in much closer result as 30 bins. Apr.27.08

          do j=1,ICOMPHARD-2
             dMdt(k,j)= &
                        +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1) &
                        +2.*phi(k-1))/6.*sk2mx(j) &
                        +(eff(k-1)-phi(k-1))/6./xk(k-1)*sk2mx2(j)) & !term1
                        +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)/3. &
                        +kij(k-1,k-1)*(phi(k-1)*xk(k)*Mk(k-1,j)/3. &
                        +(eff(k-1)+2.*phi(k-1))/6.*zeta &
                        *xbar(k-1)*Mk(k-1,j)/3.+(eff(k-1)-phi(k-1))/6. &
                        /xk(k-1)*zeta**3.*xbar(k-1)**2.*Mk(k-1,j)/3.)

             !yhl updated the following
             dMdt(k,j)= dMdt(k,j)+Nk(k)*k1m(j)-Mk(k,j)*High_in & ! !term5,term6
                        -(phi(k)*xk(k+1)*k1m(j)+(eff(k)/62.+0.484*phi(k)) &
                        *k1mx(j)+(eff(k)-phi(k))/62./xk(k)*k1mx2(j)) & ! term3
                        -kij(k,k)*Nk(k)*Mk(k,j)*0.103226 & ! I assume 1/2Nk and 2/3Mk for half bin
                        -kij(k,k)*Mk(k,j)*0.484*(phi(k)*xk(k+1)+(eff(k)/62. &
                        +0.484*phi(k))*4.754*xbar(k) &
                        +(eff(k)-phi(k))/62./xk(k)*107.4365*xbar(k)**2.)
          enddo

       else if (k.eq.nBins)then
          dNdt(k)=-Nk(k)*High_in-kij(k,k)*Nk(k)**2.*1.103226 &
                  -(phi(k)*k1mtot+(eff(k)-phi(k))/62./xk(k)*k1mxtot) &
                  -kij(k,k)*0.484*xbar(k)*Nk(k)*(phi(k)+(eff(k)-phi(k)) &
                  /62./xk(k)*4.754*xbar(k)) &
                  +0.52*kij(k-1,k-1)*Nk(k-1)**2 &
                  +(phi(k-1)*sk2mtot+(eff(k-1)-phi(k-1))/62./xk(k-1) &
                  *sk2mxtot) &
                  +kij(k-1,k-1)*xbar(k-1)*Nk(k-1)*0.484*(phi(k-1) &
                  +(eff(k-1) &
                  -phi(k-1))/62./xk(k-1)*4.754*xbar(k-1))

          do j=1,ICOMPHARD-2
             dMdt(k,j)= Nk(k)*k1m(j)-Mk(k,j)*High_in & ! !term5,term6
                        -(phi(k)*xk(k+1)*k1m(j)+(eff(k)/62.+0.484 &
                        *phi(k))*k1mx(j) &
                        +(eff(k)-phi(k))/62./xk(k)*k1mx2(j)) & ! term3
                        -kij(k,k)*Nk(k)*Mk(k,j)*0.103226 & ! I assume 1/2Nk and 2/3Mk for half bin
                        -kij(k,k)*Mk(k,j)*0.484*(phi(k)*xk(k+1)+(eff(k)/62. &
                        +0.484*phi(k))*4.754*xbar(k) &
                        +(eff(k)-phi(k))/62./xk(k)*107.4365*xbar(k)**2.) &
                        +(phi(k-1)*xk(k)*sk2m(j)+(eff(k-1)/62.+0.484 &
                        *phi(k-1)) &
                        *sk2mx(j)+(eff(k-1)-phi(k-1))/62./xk(k-1)*sk2mx2(j)) & !term1
                        +kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j)*0.103226 &
                        +kij(k-1,k-1)*Mk(k-1,j)*0.484*(phi(k-1)*xk(k) &
                        +(eff(k-1)/62.+0.484*phi(k-1))*4.754*xbar(k-1) &
                        +(eff(k-1)-phi(k-1))/62./xk(k-1)*107.4365 &
                        *xbar(k-1)**2.)
          enddo
       endif

#else
       dNdt(k)= &
                -kij(k,k)*Nk(k)**2 &
                -phi(k)*k1mtot &
                -zeta*(eff(k)-phi(k))/(2*xk(k))*k1mxtot &
                -Nk(k)*High_in
       if (k .gt. 1) then
          dNdt(k)=dNdt(k)+ &
                  0.5*kij(k-1,k-1)*Nk(k-1)**2 &
                  +phi(k-1)*sk2mtot &
                  +zeta*(eff(k-1)-phi(k-1))/(2*xk(k-1))*sk2mxtot
       endif

       do j=1,ICOMPHARD-2
          dMdt(k,j)= &
                     +Nk(k)*k1m(j) &
                     -kij(k,k)*Nk(k)*Mk(k,j) &
                     -Mk(k,j)*High_in &
                     -phi(k)*xk(k+1)*k1m(j) &
                     -0.5*zeta*eff(k)*k1mx(j) &
                     +zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j)
          if (k .gt. 1) then
             dMdt(k,j)=dMdt(k,j)+ &
                       kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j) &
                       +phi(k-1)*xk(k)*sk2m(j) &
                       +0.5*zeta*eff(k-1)*sk2mx(j) &
                       -zeta**3*(phi(k-1)-eff(k-1))/(2*xk(k-1))*sk2mx2(j)
          endif
          !dbg if (j. eq. srtso4) then
          !dbg    if (k. gt. 1) then
          !dbg       write(*,1) Nk(k)*k1m(j), kij(k,k)*Nk(k)*Mk(k,j), &
          !dbg          Mk(k,j)*in, phi(k)*xk(k+1)*k1m(j), &
          !dbg          0.5*zeta*eff(k)*k1mx(j), &
          !dbg          zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j), &
          !dbg          kij(k-1,k-1)*Nk(k-1)*Mk(k-1,j), &
          !dbg          phi(k-1)*xk(k)*sk2m(j), &
          !dbg          0.5*zeta*eff(k-1)*sk2mx(j), &
          !dbg          zeta**3*(phi(k-1)-eff(k-1))/(2*xk(k-1))*sk2mx2(j)
          !dbg    else
          !dbg       write(*,1) Nk(k)*k1m(j), kij(k,k)*Nk(k)*Mk(k,j), &
          !dbg          Mk(k,j)*in, phi(k)*xk(k+1)*k1m(j), &
          !dbg          0.5*zeta*eff(k)*k1mx(j), &
          !dbg          zeta**3*(phi(k)-eff(k))/(2*xk(k))*k1mx2(j)
          !dbg    endif
          !dbg endif
       enddo
#endif

       !dbg
       if(pdbg) write(*,*) 'k,dNdt,dMdt: ', k, dNdt(k), dMdt(k,srtso4)

       !Save the summations that are needed for the next size bin
       sk2mtot=k1mtot
       sk2mxtot=k1mxtot
       do j=1,ICOMPHARD-2
          sk2m(j)=k1m(j)
          sk2mx(j)=k1mx(j)
          sk2mx2(j)=k1mx2(j)
       enddo

    enddo  !end of main k loop

    ! Update Nk and Mk according to rates of change and time step

    !If any Mkj are zero, add a small amount to achieve finite
    !time steps
    do k=1,nBins
       do j=1,ICOMPHARD-2
          if (Mk(k,j) .eq. 0.e+0_fp) then
             !add a small amount of mass
             mtotal=0.e+0_fp
             do jj=1,ICOMPHARD-2
                mtotal=mtotal+Mk(k,jj)
             enddo
             Mk(k,j)=1.e-10_fp*mtotal
          endif
       enddo
    enddo

    call mnfix(NK, MK, PDBG)

    !Choose time step
    dts=dt-tsum      !try to take entire remaining time step
    limit='comp'
    do k=1,nBins
       if(pdbg) print*,'At bin ',k
       if (Nk(k) .gt. Neps) then
          !limit rates of change for this bin
          if (dNdt(k) .lt. 0.0) tlimit=dtlimit
          if (dNdt(k) .gt. 0.0) tlimit=itlimit
          if (abs(dNdt(k)*dts) .gt. Nk(k)*tlimit) then
             dts=Nk(k)*tlimit/abs(dNdt(k))
             if(pdbg) print*,'tlimit',tlimit,'Nk(',k,')',Nk(k), &
                             'dNdt',dNdt(k), ' == dts ',dts
             limit='number'
             if(pdbg) write(limit(8:9),'(I2)') k
             if(pdbg) write(*,*) Nk(k), dNdt(k)
          endif
          do j=1,ICOMPHARD-2
             !limit rates of change x(win, 4/22/06)
             if (dMdt(k,j) .lt. 0.0) tlimit=dtlimit
             if (dMdt(k,j) .gt. 0.0) tlimit=itlimit
             if (abs(dMdt(k,j)*dts) .gt. Mk(k,j)*tlimit) then
                mtotal=0.e+0_fp
                do jj=1,ICOMPHARD-2
                   mtotal=mtotal+Mk(k,jj)
                enddo
                !only use this criteria if this species is significant
                if ((Mk(k,j)/mtotal) .gt. 1.e-5_fp) then
                   dts=Mk(k,j)*tlimit/abs(dMdt(k,j))
                   if(pdbg) print*,'tlimit',tlimit,'Mk(',k,j,')',Mk(k,j), &
                                   'dMdt',dMdt(k,j), ' == dts ',dts
                else
                   if (dMdt(k,j) .lt. 0.0) then
                      !set dmdt to 0 to avoid very small mk going negative
                      dMdt(k,j)=0.0
                      if(pdbg) print*,' dMdt(k,j) < 0 '
                   endif
                endif
                limit='mass'
                if(pdbg) write(limit(6:7),'(I2)') k
                if(pdbg) write(limit(9:9),'(I1)') j
                if(pdbg) write(*,*) Mk(k,j), dMdt(k,j)
             endif
          enddo
       else
          !nothing in this bin - don't let it affect time step
          Nk(k)=Neps
#if defined(TOMAS12) || defined(TOMAS15)
          Mk(k,srtso4)=Neps*sqrt(xk(k)*xk(k+1)) !make the added particles SO4
#else
          Mk(k,srtso4)=Neps*1.4e+0_fp*xk(k) !make the added particles SO4
#endif
          !make sure mass/number don't go negative
          if (dNdt(k) .lt. 0.0) dNdt(k)=0.0
          if (pdbg) print*,' dNdt(k) < 0 '
          do j=1,ICOMPHARD-2
             if (dMdt(k,j) .lt. 0.0) dMdt(k,j)=0.0
          enddo
       endif
    enddo  !loop bin

    if (pdbg .and. dts .lt. 1. ) then
       write(*,*), dts, 'dts < 1. in multicoag'
    endif

    if (dts .eq. 0.) then
       write(*,*) 'time step is 0 in multicoag - inf/nan/tiny error'
       !pause
       do k = 1,nBins
          print *, 'dNdt(k)', dNdt(k)
          print *, 'dMdt(k,j)'
          do j = 1,ICOMPHARD-2
             print *, dMdt(k,j)
          end do
       end do

       !call debugprint(nk, mk, 0,0,0,'MULTICOAG before terminate: dts=0')
       PDBG = .true.
       return
       !stop
       !go to 20
    endif

    !Change Nk and Mk
    !dbg
    if(pdbg) write(*,*) 'tsum=',tsum+dts,' ',limit
    do k=1,nBins
       Nk(k)=Nk(k)+dNdt(k)*dts
       do j=1,ICOMPHARD-2
          Mk(k,j)=Mk(k,j)+dMdt(k,j)*dts
       enddo
    enddo

    !Update time and repeat process if necessary
    tsum=tsum+dts
    if (tsum .lt. dt) then
       !print*,'tsum',tsum, 'less than 3600. loop again'
       goto 10
    endif

    RETURN

  END SUBROUTINE MULTICOAG
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getcondsink
!
! !DESCRIPTION: This subroutine calculates the condensation sink (first order
!  loss rate of condensing gases) from the aerosol size distribution.
!  WRITTEN BY Jeff Pierce, May 2007 for GISS GCM-II
!  Put in GEOS-Chem by Win T. (9/30/08)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE getCondSink(Nko, Mko, spec, CS, sinkfrac, surf_area, &
            BOXVOL, TEMPTMS, PRES)
!
! !INPUT PARAMETERS:
!
    !Initial values of
    !=================
    !Nk(nBins) - number of particles per size bin in grid cell
    !Nnuc - number of particles per size bin in grid cell
    !Mnuc - mass of given species in nucleation pseudo-bin (kg/grid cell)
    !Mk(nBins, ICOMPHARD) - mass of a given species per size bin/grid cell
    !spec - number of the species we are finding the condensation sink for
    double precision Nko(nBins), Mko(nBins, ICOMPHARD)
    REAL(fp), INTENT(IN)       :: BOXVOL, TEMPTMS, PRES
    integer spec
!
! !OUTPUT PARAMETERS:
!
    !CS - condensation sink [s^-1]
    !sinkfrac(nBins) - fraction of condensation sink from a bin
    double precision CS, sinkfrac(nBins)
    REAL(fp), INTENT(OUT)    :: surf_area
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer i,j,k,c           ! counters
    double precision pi, R    ! pi and gas constant (J/mol K)
    double precision mu                  !viscosity of air (kg/m s)
    double precision mfp                 !mean free path of air molecule (m)
    double precision l_ab                !mean free path of h2so4 molecule (m)
    real Di       !diffusivity of gas in air (m2/s)
    double precision Neps     !tolerance for number
    real density  !density [kg m^-3]
    double precision mp       !mass per particle [kg]
    double precision Dpk(nBins) !diameter of particle [m]
    double precision Kn       !Knudson number
    double precision beta(nBins) !non-continuum correction factor
    double precision Mktot    !total mass in bin [kg]
    double precision c_a      !average speed of a, h2so4 molecule
!
! !DEFINED PARAMETERS:
!
    parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
    parameter(Neps=1.0e+10_fp)
    double precision alpha(ICOMPHARD) ! accomodation coef
    !data alpha/0.65,0.,0.,0.,0.,0.,0.,0.,0./
    real(fp) Sv(ICOMPHARD)         !parameter used for estimating diffusivity
    !data Sv /42.88,42.88,42.88,42.88,42.88,42.88,42.88, &
    !         42.88,42.88/

    !=================================================================
    ! getCondSink begins here
    !=================================================================

    ! have to find a better way to simply assign contants to these array
    ! The problem is I declare the array with ICOMP - its value will be
    ! determined at time of run, so I can't use DATA statement
    DO J=1,ICOMPHARD
       !IF ( J == SRTSO4 ) THEN
       alpha(J) = 0.65
       !ELSE
       !   alpha(J) = 0.
       !ENDIF
       Sv(J) = 42.88
    ENDDO


    ! get some parameters

    !mu=2.5277e-7 * TEMPTMS**0.75302
    !mfp=2.0*mu / ( pres*sqrt( 8.0 * 0.6589 / (pi*R*TEMPTMS) ) )  !S&P eqn 8.6

    !mfp=2.0*mu / ( pres*sqrt( 8.0 * 0.0289 / (pi*R*TEMPTMS) ) )  !S&P eqn 8.6

    Di=gasdiff(TEMPTMS,pres,98.0_fp,Sv(spec))

    c_a  = sqrt(8.0 * TEMPTMS * R / 0.098)
    l_ab = 2.0 * Di / c_a

    ! get size dependent values
    do k=1,nBins
       if (Nko(k) .gt. Neps) then
          Mktot=0.e+0_fp
          do j=1,ICOMPHARD
             Mktot=Mktot+Mko(k,j)
          enddo
          !kpc  Density should be changed due to more species involed.
          !density=aerodens(Mko(k,srtso4),0.e+0_fp, &
          !        Mko(k,srtnh4),Mko(k,srtnacl),Mko(k,srtecil), &
          !        Mko(k,srtecob),Mko(k,srtocil),Mko(k,srtocob), &
          !        Mko(k,srtdust),Mko(k,srth2o)) !assume bisulfate
          density=aerodens(Mko(k,srtso4),0.e+0_fp, &
                  Mko(k,srtnh4),0.e+0_fp,0.e+0_fp,&
                  0.e+0_fp,0.e+0_fp,0.e+0_fp, &
                  0.e+0_fp,Mko(k,srth2o)) !assume bisulfate
          mp=Mktot/Nko(k)
       else
          !nothing in this bin - set to "typical value"
          density=1500.
#if defined(TOMAS12) || defined(TOMAS15)
          mp=sqrt(xk(k)*xk(k+1))
#else
          mp=1.4*xk(k)
#endif
       endif
       Dpk(k)  = ( (mp/density)*(6./pi) )**(0.333)
       !Kn     = 2.0 * mfp  / Dpk(k)     !S&P eqn 11.35 (text)
       Kn      = 2.0 * l_ab / Dpk(k)     !S&Pv2 chapter 12 - Kn for Dahneke correction factor
       beta(k) = ( 1.+Kn )  / ( 1.+2.*Kn*(1.+Kn)/alpha(spec) )   !S&P eqn 11.35
    enddo
    
    ! get condensation sink
    CS = 0.e+0_fp
    surf_area = 0.e+0_fp
    do k=1,nBins
       CS = CS + Dpk(k)*Nko(k)*beta(k)
       surf_area = surf_area+Nko(k)*pi*(Dpk(k)*1.0e+6_fp)**2
    enddo
    !bc 21/01/2022 - check if divide by zero below -added 2 if 
    do k=1,nBins
       sinkfrac(k) = 0.e-0_fp
       if (CS > 0.e-0_fp) then
          sinkfrac(k) = Dpk(k)*Nko(k)*beta(k)/CS
       endif
    enddo
    CS = 2.e+0_fp*pi*dble(Di)*CS/(dble(boxvol)*1.e-6_fp)
    surf_area = 0.e-0_fp
    if (CS  > 0.e-0_fp) then
       surf_area = surf_area/(dble(boxvol))
    endif
    
    return

  end subroutine getcondsink
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getH2SO2conc
!
! !DESCRIPTION: This subroutine uses newtons method to solve for the steady
!  state H2SO4 concentration when nucleation is occuring.
!  It solves for H2SO4 in 0 = P - CS*H2SO4 - M(H2SO4)
!  where P is the production rate of H2SO4, CS is the condensation sink
!  and M(H2SO4) is the loss of mass towards making new particles.
!  WRITTEN BY Jeff Pierce, May 2007 for GISS GCM-II
!  Put in GEOS-CHEM by Win T. (9/30/08)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE getH2SO4conc(Nk, Mk, H2SO4rate, CS, NH3conc, gasConc, &
                          ionrate, surf_area, BOXVOL, BOXMASS, &
                          TEMPTMS, PRES, RHTOMAS, lev)
!
! !USES:
!
    USE ERROR_MOD,      ONLY : ERROR_STOP, IT_IS_NAN
!
! !INPUT PARAMETERS:
!
    !Initial values of
    !=================
    ! H2SO4rate - H2SO4 generation rate [kg box-1 s-1]
    ! CS - condensation sink [s-1]
    ! NH3conc - ammonium in box [kg box-1]
    REAL(fp)            :: Nk(nBins)
    REAL(fp)            :: Mk(nBins, ICOMPHARD)
    double precision       H2SO4rate
    double precision       CS
    double precision       NH3conc
    REAL(fp), INTENT(IN)  :: BOXVOL,  BOXMASS, TEMPTMS
    REAL(fp), INTENT(IN)  :: PRES,    RHTOMAS
    integer                lev
!
! !OUTPUT PARAMETERS:
!
    ! gasConc - gas H2SO4 [kg/box]
    double precision       gasConc
    REAL(fp)            :: ionrate
    REAL(fp)            :: surf_area
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer i,j,k,c           ! counters
    double precision fn, rnuc ! nucleation rate [# cm-3 s-1] and critical radius [nm]
    double precision mnuc, mnuc1 ! mass of nucleated particle [kg]
    double precision fn1, rnuc1 ! nucleation rate [# cm-3 s-1] and critical radius [nm]
    double precision res      ! d[H2SO4]/dt, need to find the solution where res = 0
    double precision massnuc     ! mass being removed by nucleation [kg s-1 box-1]
    double precision gasConc1 ! perturbed gasConc
    double precision gasConc_hi, gasConc_lo
    double precision res1     ! perturbed res
    double precision res_new  ! new guess for res
    double precision dresdgasConc ! derivative for newtons method
    double precision Gci(ICOMPHARD)      !array to carry gas concentrations
    logical nflg              !says if nucleation occured
    double precision H2SO4min !minimum H2SO4 concentration in parameterizations (molec/cm3)
    double precision pi
    integer iter,iter1
    double precision CSeps    ! low limit for CS
    double precision max_H2SO4conc !maximum H2SO4 concentration in parameterizations (kg/box)
    double precision nh3ppt   !ammonia concentration in ppt
!
! !DEFINED PARAMETERS:
!
    parameter(pi=3.141592654)
    !parameter(H2SO4min=1.D4) !molecules cm-3
    parameter(CSeps=1.0e-20_fp)

    !=================================================================
    ! getH2SO4conc begins here
    !=================================================================

    do i=1,ICOMPHARD
       Gci(i)=0.e+0_fp
    enddo
    Gci(srtnh4)=NH3conc

    ! make sure CS doesn't equal zero
    !CS = max(CS,CSeps)

    ! some specific stuff for napari vs. vehk
    if (ion_nuc.eq.1) then
       H2SO4min=1.0e+5_fp
    elseif (ion_nuc.eq.2) then
       H2SO4min=5.0e+5_fp
    else
       H2SO4min=1.0e+4_fp
    endif

    if ((bin_nuc.eq.1).or.(tern_nuc.eq.1).or.(ion_nuc.le.2))then
       nh3ppt = Gci(srtnh4)/17.e+0_fp/(boxmass/29.e+0_fp)*1e+12_fp* &
  &             PRES/101325.*273./TEMPTMS ! corrected for pressure (because this should be concentration)
       if (ion_nuc.eq.1)then
          max_H2SO4conc=1.0e+8_fp*boxvol/1000.e+0_fp*98.e+0_fp/6.022e+23_fp
       elseif (ion_nuc.eq.2)then
          max_H2SO4conc=5.0e+8_fp*boxvol/1000.e+0_fp*98.e+0_fp/6.022e+23_fp
       elseif ((nh3ppt.gt.1.0e+0_fp).and.(tern_nuc.eq.1))then
          max_H2SO4conc=1.0e+9_fp*boxvol/1000.e+0_fp*98.e+0_fp/6.022e+23_fp
       elseif (bin_nuc.eq.1)then
          max_H2SO4conc=1.0e+11_fp*boxvol/1000.e+0_fp*98.e+0_fp/6.022e+23_fp
       else
          max_H2SO4conc = 1.0e+100_fp
       endif
    else
       max_H2SO4conc = 1.0e+100_fp
    endif

    ! Checks for when condensation sink is very small
    if (CS.gt.CSeps) then
       gasConc = H2SO4rate/CS
    else
       if((bin_nuc.gt.0).or.(tern_nuc.gt.0).or. (ion_nuc.gt.0))then
          gasConc = max_H2SO4conc
       else
          print*,'condesation sink too small in getH2SO4conc'
          STOP
       endif
    endif

    gasConc = min(gasConc,max_H2SO4conc)
    Gci(srtso4) = gasConc
    call getNucRate(Nk, Mk, Gci,fn,mnuc,nflg,ionrate, surf_area, &
                    BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)

    if (fn.gt.0.e+0_fp) then      ! nucleation occured
       !convert to kg/box
       gasConc_lo = H2SO4min*boxvol/(1000.e+0_fp/98.e+0_fp*6.022e+23_fp)

       ! Test to see if gasConc_lo gives a res < 0
       ! (this means ANY nucleation is too high)
       Gci(srtso4) = gasConc_lo*1.000001e+0_fp
       call getNucRate(Nk,Mk,Gci,fn1,mnuc1,nflg,ionrate,surf_area, &
                       BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
       if (fn1.gt.0.e+0_fp) then
          massnuc = mnuc1*fn1*boxvol*98.e+0_fp/96.e+0_fp
          !massnuc = 4.e+0_fp/3.e+0_fp*pi*(rnuc1*1.e-9_fp)**3*1350.*fn1*boxvol*
          !massnuc = 4.e+0_fp/3.e+0_fp*pi*(rnuc1*1.e-9_fp)**3*1800.*fn1*boxvol*%
          !          98.e+0_fp/96.e+0_fp
          !jrp print*,'res',res
          !jrp print*,'H2SO4rate',H2SO4rate
          !jrp print*,'CS*gasConc_lo',CS*gasConc_lo
          !jrp print*,'mnuc',mnuc
          res = H2SO4rate - CS*gasConc_lo - massnuc
          if (res.lt.0.e+0_fp) then ! any nucleation too high
             ! if (.not. spinup(14.0)) print*,'nucleation cuttoff'
             ! have nucleation occur and fix mass balance after
             gasConc = gasConc_lo*1.000001
             return
          endif
       endif

       ! we know this must be the upper limit (since no nucleation)
       gasConc_hi = gasConc
       !take density of nucleated particle to be 1350 kg/m3
       massnuc = mnuc*fn*boxvol*98.e+0_fp/96.e+0_fp
       !print*,'H2SO4rate',H2SO4rate,'CS*gasConc',CS*gasConc,'mnuc',mnuc
       res = H2SO4rate - CS*gasConc - massnuc

       ! check to make sure that we can get solution
       if (res.gt.H2SO4rate*1.e-10_fp) then
          print*,'gas production rate too high in getH2SO4conc'
          print*,H2SO4rate,CS,gasConc,massnuc,res
          return
          !STOP
       endif

       iter = 0
       !jrp print*, 'iter',iter
       !jrp print*,'gasConc_lo',gasConc_lo,'gasConc_hi',gasConc_hi
       !jrp print*,'res',res
       do while ((abs(res/H2SO4rate).gt.1.e-4_fp).and.(iter.lt.40))
          iter = iter+1
          if (res .lt. 0.e+0_fp) then ! H2SO4 concentration too high, must reduce
             gasConc_hi = gasConc ! old guess is new upper bound
          elseif (res .gt. 0.e+0_fp) then ! H2SO4 concentration too low, must increase
             gasConc_lo = gasConc ! old guess is new lower bound
          endif
          !print*, 'iter',iter
          !print*,'gasConc_lo',gasConc_lo,'gasConc_hi',gasConc_hi
          gasConc = sqrt(gasConc_hi*gasConc_lo) ! take new guess as logmean
          Gci(srtso4) = gasConc
          call getNucRate(Nk, Mk,Gci,fn,mnuc,nflg,ionrate,surf_area, &
                          BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
          massnuc = mnuc*fn*boxvol*98.e+0_fp/96.e+0_fp
          res = H2SO4rate - CS*gasConc - massnuc
          !print*,'res',res
          !print*,'H2SO4rate',H2SO4rate,'CS',CS,'gasConc',gasConc
          if (iter.eq.40.and.CS.gt.1.0e-4_fp)then
             print*,'getH2SO4conc iter break'
             print*,'H2SO4rate',H2SO4rate,'CS',CS
             print*,'gasConc',gasConc,'massnuc',massnuc
             print*,'max_H2SO4conc',max_H2SO4conc
             print*,'fn',fn
             print*,'res/H2SO4rate',res/H2SO4rate
          endif
       enddo

       !print*,'IN getH2SO4conc'
       !print*,'fn',fn
       !print*,'H2SO4rate',H2SO4rate
       !print*,'massnuc',massnuc,'CS*gasConc',CS*gasConc

    else
       ! nucleation didn't occur
    endif

    return

  end SUBROUTINE GETH2SO4CONC
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getnucrate
!
! !DESCRIPTION: This subroutine calls the Vehkamaki 2002 and Napari 2002
!  nucleation parameterizations and gets the binary and ternary nucleation
!  rates.
!  WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II
!  Put in GEOS-Chem by win T. 9/30/08
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE getNucRate(Nk, Mk, Gci,fn,mnuc,nflg, ionrate,surf_area, &
                        BOXVOL, BOXMASS, TEMPTMS, PRES, RHTOMAS, lev)
!
! !USES:
!
    USE ERROR_MOD,      ONLY : ERROR_STOP, IT_IS_NAN
!
! !INPUT PARAMETERS:
!
    !Initial values of
    !=================
    ! Gci(icomp-1) - amount (kg/grid cell) of all species present in the
    !                gas phase except water
    REAL(fp),   INTENT(IN)       :: BOXVOL,  BOXMASS, TEMPTMS
    REAL(fp),   INTENT(IN)       :: PRES,    RHTOMAS
    REAL(fp), INTENT(IN)       :: Gci(ICOMPHARD)
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(INOUT)    :: Nk(nBins)
    REAL(fp), INTENT(INOUT)    :: Mk(nBins, ICOMPHARD)
!
! !OUTPUT PARAMETERS:
!
    ! fn - nucleation rate [# cm-3 s-1]
    ! rnuc - radius of nuclei [nm]
    ! nflg - says if nucleation happend
    REAL(fp)                   :: surf_area
    REAL(fp)                   :: ionrate

    integer j,i,k
    double precision fn       ! nucleation rate to first bin cm-3 s-1
    double precision mnuc     !mass of nucleating particle [kg]
    logical nflg
    integer lev
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    double precision nh3ppt   ! gas phase ammonia in pptv
    double precision h2so4    ! gas phase h2so4 in molec cc-1
    double precision gtime    ! time to grow to first size bin [s]
    double precision ltc, ltc1, ltc2 ! coagulation loss rates [s-1]
    double precision Mktot    ! total mass in bin
    double precision neps
    double precision meps
    double precision density  ! density of particle [kg/m3]
    double precision pi
    double precision frac     ! fraction of particles growing into first size bin
    double precision d1,d2    ! diameters of particles [m]
    double precision mp       ! mass of particle [kg]
    double precision mold     ! saved mass in first bin
    double precision rnuc     ! critical nucleation radius [nm]
    double precision sinkfrac(nBins) ! fraction of loss to different size bins
    double precision nadd     ! number to add
    double precision CS       ! kerminan condensation sink [m-2]
    double precision Dpmean   ! the number wet mean diameter of the existing aerosol
    double precision Dp1      ! the wet diameter of bin 1
    double precision dens1    ! density in bin 1 [kg m-3]
    double precision GR       ! growth rate [nm hr-1]
    double precision gamma,eta ! used in kerminen 2004 parameterzation
    double precision drymass,wetmass,WR
    double precision fn_c     ! barrierless nucleation rate
    double precision h1,h2,h3,h4,h5,h6
    double precision dum1,dum2,dum3,dum4   ! dummy variables
    double precision rhin,tempin ! rel hum in

    real(fp)    mydummy
!
! !DEFINED PARAMETERS:
!
    parameter (neps=1E8, meps=1E-8)
    parameter (pi=3.14159)

    !=================================================================
    ! getNucRate begins here
    !=================================================================

    h2so4 = Gci(srtso4)/boxvol*1000.e+0_fp/98.e+0_fp*6.022e+23_fp
    nh3ppt = Gci(srtnh4)/17.e+0_fp/(boxmass/29.e+0_fp)*1e+12_fp* &
             PRES/101325.*273./TEMPTMS ! corrected for pressure (because this should be concentration)

    fn = 0.e+0_fp
    rnuc = 0.e+0_fp

    !print*,'h2so4',h2so4,'nh3ppt',nh3ppt

    ! if requirements for nucleation are met, call nucleation subroutines
    ! and get the nucleation rate and critical cluster size
    if (h2so4.gt.1.e+4_fp) then
       if ((nh3ppt.gt.0.1).and.(tern_nuc.eq.1)) then
          ! print*, 'napari'
          call napa_nucl(TEMPTMS,RHTOMAS,h2so4,nh3ppt,fn,rnuc) !ternary nuc
          if (ion_nuc.eq.1.and.ionrate.ge.1.e+0_fp) then
             !call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS, &
             !              h1,h2,h3,h4,h5,h6)
          else
             h1=0.e+0_fp
          endif
          if (h1.gt.fn)then
             fn=h1
             rnuc=h5
          endif
          nflg=.true.
       elseif (bin_nuc.eq.1) then
          ! print*, 'vehk'
          call vehk_nucl(TEMPTMS,RHTOMAS,h2so4,fn,rnuc) !binary nuc
          if ((ion_nuc.eq.1).and.(ionrate.ge.1.e+0_fp)) then
             !call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS, &
             !              h1,h2,h3,h4,h5,h6)
          else
             h1=0.e+0_fp
          endif
          if (h1.gt.fn)then
             fn=h1
             rnuc=h5
          endif
          if (fn.gt.1.0e-6_fp)then
             nflg=.true.
          else
             fn = 0.e+0_fp
             nflg=.false.
          endif
       elseif ((ion_nuc.eq.1).and.(ionrate.ge.1.e+0_fp)) then
          !call ion_nucl(h2so4,surf_area,TEMPTMS,ionrate,RHTOMAS, &
          !              h1,h2,h3,h4,h5,h6)
          fn=h1
          rnuc=h5
          nflg=.true.
       elseif(ion_nuc.eq.2) then
          ! Yu Ion nucleation
          !! first we need to calculate the available surface area
          !surf_area = 0.e+0_fp
          !do k=1, ibins
          !   if (Nki(k) .gt. Neps) then
          !      Mktot=0.e+0_fp
          !      do j=1,icomp
          !         Mktot=Mktot+Mki(k,j)
          !      enddo
          !      mp=Mktot/Nki(k)
          !      density=aerodens(Mki(k,srtso4),0.e+0_fp, &
          !                       Mki(k,srtnh4),0.e+0_fp,Mki(k,srth2o))  ! assume bisulfate
          !      ! diameter = ((mass/density)*(6/pi))**(1/3)
          !      d2 = 1.D6*((mp/density)*(6.D0/pi))**(1.D0/3.D0) ! (micrometers)
          !      ! surface area per particle = pi*diameter**2
          !      surf_area = surf_area + 1.D-6*(Nki(k)/boxvol)* &
          !                  pi*(d2**2.D0) ! (um2 cm-2)
          !   endif
          !enddo
          rhin=dble(RHTOMAS*100.e+0_fp)
          tempin=dble(TEMPTMS)
          !call YUJIMN(h2so4, rhin, tempin, ionrate, surf_area, &
          !            fn, dum1, rnuc, dum2)
          fn=0.
          rnuc=1E-9
          nflg=.true.
       else
          nflg=.false.
       endif
       if((act_nuc.eq.1).and.(lev.le.7))then
          ! call bl_nucl(h2so4,fn,rnuc)
          nflg=.true.
       endif
       call cf_nucl(TEMPTMS,RHTOMAS,h2so4,nh3ppt,fn_c) ! use barrierless nucleation as a max for ternary
       fn = min(fn,fn_c)
    else
       nflg=.false.
    endif

    if (fn.gt.0.e+0_fp) then
       call getCondSink_kerm(Nk,Mk,CS,Dpmean,Dp1,dens1, &
                             BOXVOL, TEMPTMS, PRES)
       d1 = rnuc*2.e+0_fp*1e-9_fp
       drymass = 0.e+0_fp
       do j=1,ICOMPHARD-2
          drymass = drymass + Mk(1,j)
       enddo
       wetmass = 0.e+0_fp
       do j=1,ICOMPHARD
          wetmass = wetmass + Mk(1,j)
       enddo
       !prior 10/15/08
       !WR = wetmass/drymass

       ! prevent division by zero (win, 10/15/08)
       if( drymass == 0.e+0_fp ) then
          WR = 1.e+0_fp
       else
          WR = wetmass/drymass
       endif

       !print*,'[getnucrate] Gci',Gci
       !print*,'WR',WR, 'drymass',drymass, 'wetmass',wetmass
       call getGrowthTime(d1,Dp1,Gci(srtso4)*WR,TEMPTMS, &
                          boxvol,dens1,gtime)
       GR = (Dp1-d1)*1e+9_fp/gtime*3600.e+0_fp ! growth rate, nm hr-1

       gamma = 0.23e+0_fp*(d1*1.0e+9_fp)**(0.2e+0_fp)* &
               (Dp1*1.0e+9_fp/3.e+0_fp)**0.075e+0_fp* &
               (Dpmean*1.0e+9_fp/150.e+0_fp)** &
               0.048e+0_fp*(dens1*1.0e-3_fp)** &
               (-0.33e+0_fp)*(TEMPTMS/293.e+0_fp) ! equation 5 in kerminen
       eta = gamma*CS/GR
       !print*,'fn1',fn
       if (Dp1.gt.d1)then
          fn = fn*exp(eta/(Dp1*1.0e+9_fp)-eta/(d1*1.0e+9_fp))
       endif
       !print*,'fn2',fn
       if( IT_IS_NAN( fn ) ) then
          print*, '---------------->>> Found NAN in GETNUCRATE'
          print*,'fn',fn
          print*,'eta',eta, 'Dp1',Dp1,'d1',d1
          print*,'gamma',gamma,'CS',CS,'GR',GR,'gtime',gtime
          call ERROR_STOP('Found NaN in fn','getnucrate')
       endif

       mnuc = sqrt(xk(1)*xk(2))
    endif

    return

  end SUBROUTINE GETNUCRATE
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: vehk_nucl
!
! !DESCRIPTION: Subroutine vehk_nucl calculates the binary nucleation rate and
!  radius of the critical nucleation cluster using the parameterization of...
!  .
!    Vehkamaki, H., M. Kulmala, I. Napari, K. E. J. Lehtinen, C. Timmreck,
!    M. Noppel, and A. Laaksonen. "An Improved Parameterization for Sulfuric
!    Acid-Water Nucleation Rates for Tropospheric and Stratospheric Conditions."
!    Journal of Geophysical Research-Atmospheres 107, no. D22 (2002).
!  .
!  WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II'
!  Introduce to GEOS-Chem by Win Trivitayanurak Sep 29,2008
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE VEHK_NUCL (tempi,rhi,cnai,fn,rnuc)
!
! !INPUT PARAMETERS:
!
    real(fp),   intent(in)   :: tempi ! temperature of air [K]
    real(fp),   intent(in)   :: rhi ! relative humidity of air as a fraction
    real(fp), intent(in)   :: cnai ! concentration of gas phase sulfuric acid [molec cm-3]
!
! !OUTPUT PARAMETERS:
!
    real(fp), intent(out)  :: fn ! nucleation rate [cm-3 s-1]
    real(fp), intent(out)  :: rnuc ! critical cluster radius [nm]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)  :: fb0(10),fb1(10),fb2(10),fb3(10),fb4(10),fb(10)
    REAL(fp)  :: gb0(10),gb1(10),gb2(10),gb3(10),gb4(10),gb(10) ! set parameters
    REAL(fp)  :: temp    ! temperature of air [K]
    REAL(fp)  :: rh      ! relative humidity of air as a fraction
    REAL(fp)  :: cna     ! concentration of gas phase sulfuric acid [molec cm-3]
    REAL(fp)  :: xstar   ! mole fraction sulfuric acid in cluster
    REAL(fp)  :: ntot    ! total number of molecules in cluster
    integer   :: i       ! counter

    ! Nucleation Rate Coefficients
    data fb0 /0.14309, 0.117489, -0.215554, -3.58856, 1.14598, &
              2.15855, 1.6241, 9.71682, -1.05611, -0.148712        /
    data fb1 /2.21956, 0.462532, -0.0810269, 0.049508, -0.600796, &
              0.0808121, -0.0160106, -0.115048, 0.00903378, 0.00283508/
    data fb2 /-0.0273911, -0.0118059, 0.00143581, -0.00021382, &
               0.00864245, -0.000407382, 0.0000377124, 0.000157098, &
              -0.0000198417, -9.24619e-6_fp /
    data fb3 /0.0000722811, 0.0000404196, &
             -4.7758e-6_fp, 3.10801e-7_fp, &
             -0.0000228947, -4.01957e-7_fp, &
              3.21794e-8_fp, 4.00914e-7_fp, &
              2.46048e-8_fp, 5.00427e-9_fp /
    data fb4 /5.91822, 15.7963, -2.91297, -0.0293333, -8.44985, &
              0.721326, -0.0113255, 0.71186, -0.0579087, -0.0127081  /

    ! Coefficients of total number of molecules in cluster
    data gb0 /-0.00295413, -0.00205064, 0.00322308, 0.0474323, &
              -0.0125211, -0.038546, -0.0183749, -0.0619974, &
               0.0121827, 0.000320184 /
    data gb1 /-0.0976834, -0.00758504, 0.000852637, -0.000625104, &
               0.00580655, -0.000672316, 0.000172072, 0.000906958, &
              -0.00010665, -0.0000174762 /
    data gb2 /0.00102485, 0.000192654, &
             -0.0000154757, 2.65066e-6_fp, &
             -0.000101674, 2.60288e-6_fp, &
             -3.71766e-7_fp, -9.11728e-7_fp, &
             2.5346e-7_fp, 6.06504e-8_fp /
    data gb3 /-2.18646e-6_fp, -6.7043e-7_fp, &
               5.66661e-8_fp, -3.67471e-9_fp, &
               2.88195e-7_fp, 1.19416e-8_fp, &
              -5.14875e-10_fp, -5.36796e-9_fp, &
              -3.63519e-10_fp, -1.42177e-11_fp /
    data gb4 /-0.101717, -0.255774, 0.0338444, -0.000267251, &
               0.0942243, -0.00851515, 0.00026866, -0.00774234, &
               0.000610065, 0.000135751 /

    !=================================================================
    ! VEHK_NUCL begins here!
    !=================================================================
    temp=dble(tempi)
    rh=dble(rhi)
    cna=cnai

    ! Respect the limits of the parameterization
    if (cna .lt. 1.e4_fp) then ! limit sulf acid conc
       fn = 0.
       rnuc = 1.
       !print*,'cna < 1D4', cna
       goto 10
    endif
    if (cna .gt. 1.0e+11_fp) cna=1.0e11 ! limit sulfuric acid conc
    if (temp .lt. 230.15) temp=230.15 ! limit temp
    if (temp .gt. 305.15) temp=305.15 ! limit temp
    if (rh .lt. 1e-4_fp) rh=1e-4_fp ! limit rh
    if (rh .gt. 1.) rh=1. ! limit rh

    ! Mole fraction of sulfuric acid
    xstar=0.740997-0.00266379*temp-0.00349998*log(cna) &
         +0.0000504022*temp*log(cna)+0.00201048*log(rh) &
         -0.000183289*temp*log(rh)+0.00157407*(log(rh))**2. &
         -0.0000179059*temp*(log(rh))**2. &
         +0.000184403*(log(rh))**3. &
         -1.50345e-6_fp*temp*(log(rh))**3.

    ! Nucleation rate coefficients
    do i=1, 10
       fb(i) = fb0(i)+fb1(i)*temp+fb2(i)*temp**2. &
              +fb3(i)*temp**3.+fb4(i)/xstar
    enddo

    ! Nucleation rate (1/cm3-s)
    fn = exp(fb(1)+fb(2)*log(rh)+fb(3)*(log(rh))**2. &
         +fb(4)*(log(rh))**3.+fb(5)*log(cna) &
         +fb(6)*log(rh)*log(cna)+fb(7)*(log(rh))**2.*log(cna) &
         +fb(8)*(log(cna))**2.+fb(9)*log(rh)*(log(cna))**2. &
         +fb(10)*(log(cna))**3.)

    !print*,'in vehk_nuc, fn',fn
    !print*,'cna',cna,'rh',rh,'temp',temp
    !print*,'xstar',xstar

    ! Cap at 10^6 particles/s, limit for parameterization
    if (fn.gt.1.0e+6_fp) then
       fn=1.0e+6_fp
    endif

    ! Coefficients of total number of molecules in cluster
    do i=1, 10
       gb(i) = gb0(i)+gb1(i)*temp+gb2(i)*temp**2. &
              +gb3(i)*temp**3.+gb4(i)/xstar
    enddo
    ! Total number of molecules in cluster
    ntot=exp(gb(1)+gb(2)*log(rh)+gb(3)*(log(rh))**2. &
         +gb(4)*(log(rh))**3.+gb(5)*log(cna) &
         +gb(6)*log(rh)*log(cna)+gb(7)*log(rh)**2.*log(cna) &
         +gb(8)*(log(cna))**2.+gb(9)*log(rh)*(log(cna))**2. &
         +gb(10)*(log(cna))**3.)

    ! cluster radius
    rnuc=exp(-1.6524245+0.42316402*xstar+0.3346648*log(ntot)) ! [nm]

10  return

  end SUBROUTINE VEHK_NUCL
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: napa_nucl
!
! !DESCRIPTION:  Subroutine NAPA_NUCL calculates the ternary nucleation rate
!  and radius of the critical nucleation cluster using the parameterization of
!  .
!     Napari, I., M. Noppel, H. Vehkamaki, and M. Kulmala. "Parametrization of
!     Ternary Nucleation Rates for H2so4-Nh3-H2o Vapors." Journal of Geophysical
!     Research-Atmospheres 107, no. D19 (2002).
!  .
!  WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II'
!  Introduce to GEOS-Chem by Win Trivitayanurak Sep 29, 2008
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE napa_nucl(tempi,rhi,cnai,nh3ppti,fn,rnuc)
!
! !INPUT PARAMETERS:
!
    real(fp),   intent(in) :: tempi ! temperature of air [K]
    real(fp),   intent(in) :: rhi ! relative humidity of air as a fraction
    real(fp), intent(in) :: cnai ! concentration of gas phase sulfuric acid [molec cm-3]
    real(fp), intent(in) :: nh3ppti ! concentration of gas phase ammonia
!
! !OUTPUT PARAMETERS:
!
    real(fp), intent(out):: fn  ! nucleation rate [cm-3 s-1]
    real(fp), intent(out):: rnuc ! critical cluster radius [nm]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(fp)    ::  aa0(20),a1(20),a2(20),a3(20),fa(20) ! set parameters
    real(fp)    ::  fnl     ! natural log of nucleation rate
    real(fp)    ::  temp    ! temperature of air [K]
    real(fp)    ::  rh      ! relative humidity of air as a fraction
    real(fp)    ::  cna     ! concentration of gas phase sulfuric acid [molec cm-3]
    real(fp)    ::  nh3ppt  ! concentration of gas phase ammonia
    integer     ::  i       ! counter

    ! Adjustable parameters
    data aa0 /-0.355297, 3.13735, 19.0359, 1.07605, 6.0916, &
               0.31176, -0.0200738, 0.165536, &
               6.52645, 3.68024, -0.066514, 0.65874, &
               0.0599321, -0.732731, 0.728429, 41.3016, &
               -0.160336, 8.57868, 0.0530167, -2.32736        /

    data a1 /-33.8449, -0.772861, -0.170957, 1.48932, -1.25378, &
               1.64009, -0.752115, 3.26623, -0.258002, -0.204098, &
              -7.82382, 0.190542, 5.96475, -0.0184179, 3.64736, &
              -0.35752, 0.00889881, -0.112358, -1.98815, 0.0234646/

    data a2 /0.34536, 0.00561204, 0.000479808, -0.00796052, &
             0.00939836, -0.00343852, 0.00525813, -0.0489703, &
             0.00143456, 0.00106259, 0.0122938, -0.00165718, &
            -0.0362432, 0.000147186, -0.027422, 0.000904383, &
            -5.39514d-05, 0.000472626, 0.0157827, -0.000076519/

    data a3 /-0.000824007, -9.74576e-06_fp, &
             -4.14699e-07_fp, 7.61229e-06_fp, &
             -1.74927e-05_fp, -1.09753e-05_fp, &
             -8.98038e-06_fp, 0.000146967, &
             -2.02036e-06_fp, -1.2656e-06_fp, &
              6.18554e-05_fp, 3.41744e-06_fp, &
              4.93337e-05_fp, -2.37711e-07_fp, &
              4.93478e-05_fp, -5.73788e-07_fp, &
              8.39522e-08_fp, -6.48365e-07_fp, &
             -2.93564e-05_fp, 8.0459e-08_fp   /

    !=================================================================
    ! NAPA_NUCL begins here!
    !=================================================================
    temp=dble(tempi)
    rh=dble(rhi)
    cna=cnai
    nh3ppt=nh3ppti

    ! Napari's parameterization is only valid within limited area
    if ((cna .lt. 1.e+4_fp).or.(nh3ppt.lt.0.1)) then ! limit sulf acid and nh3 conc
       fn = 0.
       rnuc = 1
       goto 10
    endif
    if (cna .gt. 1.0e+9_fp) cna=1.0e+9_fp ! limit sulfuric acid conc
    if (nh3ppt .gt. 100.) nh3ppt=100. ! limit temp
    if (temp .lt. 240.) temp=240. ! limit temp
    if (temp .gt. 300.) temp=300. ! limit temp
    if (rh .lt. 0.05) rh=0.05 ! limit rh
    if (rh .gt. 0.95) rh=0.95 ! limit rh

    do i=1,20
       fa(i)=aa0(i)+a1(i)*temp+a2(i)*temp**2.+a3(i)*temp**3.
    enddo

    fnl=-84.7551+fa(1)/log(cna)+fa(2)*log(cna)+fa(3)*(log(cna))**2. &
       +fa(4)*log(nh3ppt)+fa(5)*(log(nh3ppt))**2.+fa(6)*rh &
       +fa(7)*log(rh)+fa(8)*log(nh3ppt)/log(cna)+fa(9)*log(nh3ppt) &
       *log(cna)+fa(10)*rh*log(cna)+fa(11)*rh/log(cna) &
       +fa(12)*rh &
       *log(nh3ppt)+fa(13)*log(rh)/log(cna)+fa(14)*log(rh) &
       *log(nh3ppt)+fa(15)*(log(nh3ppt))**2./log(cna)+fa(16)*log(cna) &
       *(log(nh3ppt))**2.+fa(17)*(log(cna))**2.*log(nh3ppt) &
       +fa(18)*rh &
       *(log(nh3ppt))**2.+fa(19)*rh*log(nh3ppt)/log(cna)+fa(20) &
       *(log(cna))**2.*(log(nh3ppt))**2.

    fn=exp(fnl)

    ! Try scaling down the rate by 1e-5 to see how the param is
    ! doing on the false positive nucleation (win, 12/18/08)
    !sensitivity simulation, change scaling factor down to 1e-4
    fn = fn * 1.e-5

    ! Cap at 10^6 particles/cm3-s, limit for parameterization
    if (fn.gt.1.0e+6_fp) then
       fn=1.0e+6_fp
       fnl=log(fn)
    endif

    rnuc=0.141027-0.00122625*fnl-7.82211e-6_fp*fnl**2. &
        -0.00156727*temp-0.00003076*temp*fnl &
        +0.0000108375*temp**2.

10  return

  end subroutine napa_nucl
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cf_nucl
!
! !DESCRIPTION: This subroutine calculates the barrierless nucleation rate and
!  radius of the critical nucleation cluster using the parameterization of...
!     Clement and Ford (1999) Atmos. Environ. 33:489-499
!     WRITTEN BY Jeff Pierce, April 2007
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE cf_nucl(tempi,rhi,cna,nh3ppt,fn)
!
! !INPUT PARAMETERS:
!
    real(fp) tempi                ! temperature of air [K]
    real(fp) rhi                  ! relative humidity of air as a fraction
    double precision cna      ! concentration of gas phase sulfuric acid [molec cm-3]
    double precision nh3ppt   ! mixing ratio of ammonia in ppt
!
! !OUTPUT PARAMETERS:
!
    double precision fn                   ! nucleation rate [cm-3 s-1]
    double precision rnuc                 ! critical cluster radius [nm]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    double precision temp                 ! temperature of air [K]
    double precision rh                   ! relative humidity of air as a fraction
    double precision alpha1

    temp=dble(tempi)
    rh=dble(rhi)

    if (nh3ppt .lt. 0.1) then
       alpha1=4.276e-10*sqrt(temp/293.15) ! For sulfuric acid
    else
       alpha1=3.684e-10*sqrt(temp/293.15) ! For ammonium sulfate
    endif
    fn = alpha1*cna**2*3600.
    ! sensitivity       fn = 1.e-3 * fn ! 10^-3 tuner
    if (fn.gt.1.0e9) fn=1.0e9 ! For numerical conversion

10  return

  end subroutine cf_nucl
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getcondsink_kerm
!
! !DESCRIPTION: Subroutine GETCONDSINK\_KERM calculates the condensation sink
!  (first order loss rate of condensing gases) from the aerosol size
!  distribution.
!  .
!  This is the cond sink in kerminen et al 2004 Parameterization for
!  new particle formation AS&T Eqn 6.
!  .
!  Written by Jeff Pierce, May 2007 for GISS GCM-II'
!  Introduced to GEOS-Chem by Win Trivitayanurak, Sep 29, 2008
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE getCondSink_kerm(Nko,Mko,CS,Dpmean,Dp1,dens1, &
                              BOXVOL, TEMPTMS, PRES)
!
! !INPUT PARAMETERS:
!
    ! Nk(nBins) - number of particles per size bin in grid cell
    ! Mk(nBins, ICOMPHARD) - mass of a given species per size bin/grid cell
    REAL(fp), INTENT(IN)        :: Nko(nBins), Mko(nBins, ICOMPHARD)
    REAL(fp),   INTENT(IN)        :: BOXVOL, TEMPTMS, PRES
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT)       :: CS       ! CS - condensation sink [s^-1]
    REAL(fp), INTENT(OUT)       :: Dpmean   ! the number mean diameter [m]
    REAL(fp), INTENT(OUT)       :: Dp1      ! the size of the first size bin [m]
    REAL(fp), INTENT(OUT)       :: dens1    ! the density of the first size bin [kg/m3]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Nnuc - number of particles per size bin in grid cell
    ! Mnuc - mass of given species in nucleation pseudo-bin (kg/grid cell)
    ! spec - number of the species we are finding the condensation sink for
    ! sinkfrac(nBins) - fraction of condensation sink from a bin
    integer        :: i,j,k,c           ! counters
    REAL(fp)       :: pi, R    ! pi and gas constant (J/mol K)
    REAL(fp)       :: mu                  !viscosity of air (kg/m s)
    REAL(fp)       :: mfp                 !mean free path of air molecule (m)
    REAL*4         :: Di       !diffusivity of gas in air (m2/s)
    REAL(fp)       :: Neps     !tolerance for number
    REAL*4         :: density  !density [kg m^-3]
    REAL(fp)       :: mp       !mass per particle [kg]
    REAL(fp)       :: Dpk(nBins) !diameter of particle [m]
    REAL(fp)       :: Kn       !Knudson number
    REAL(fp)       :: beta(nBins) !non-continuum correction factor
    REAL(fp)       :: Mktot    !total mass in bin [kg]
    REAL(fp)       :: Dtot,Ntot ! used on getting the number mean diameter
!
! !DEFINED PARAMETERS:
!
    parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
    parameter(Neps=1.0e+10_fp)

    !=================================================================
    ! GETCONDSINK_KERM  begins here!
    !=================================================================

    ! get some parameters
    mu=2.5277e-7*TEMPTMS**0.75302
    !mfp=2.0*mu / ( pres*sqrt( 8.0 * 0.6589 / (pi*R*TEMPTMS) ) )  !S&P eqn 8.6
    mfp=2.0*mu / ( pres*sqrt( 8.0 * 0.0289 / (pi*R*TEMPTMS) ) )  !S&P eqn 8.6
    !Di=gasdiff(temp,pres,98.0,Sv(srtso4))
    !print*,'Di',Di

    ! get size dependent values
    CS = 0.e+0_fp
    Ntot = 0.e+0_fp
    Dtot = 0.e+0_fp
    do k=1,nBins
       if (Nko(k) .gt. Neps) then
          Mktot=0.e+0_fp
          do j=1,ICOMPHARD
             Mktot=Mktot+Mko(k,j)
          enddo
          !kpc Density should be changed due to more species involed.
          !density=aerodens(Mko(k,srtso4),0.e+0_fp, &
          !        Mko(k,srtnh4),Mko(k,srtnacl),Mko(k,srtecil), &
          !        Mko(k,srtecob),Mko(k,srtocil),Mko(k,srtocob), &
          !        Mko(k,srtdust),Mko(k,srth2o))
          density=aerodens(Mko(k,srtso4),0.e+0_fp, &
                  Mko(k,srtnh4),0.e+0_fp,0.e+0_fp, &
                  0.e+0_fp,0.e+0_fp,0.e+0_fp, &
                  0.e+0_fp,Mko(k,srth2o))
          mp=Mktot/Nko(k)
       else
          !nothing in this bin - set to "typical value"
          density=1500.
          mp=1.4*xk(k)
       endif
       Dpk(k)=((mp/density)*(6./pi))**(0.333)
       Kn=2.0*mfp/Dpk(k)      !S&P eqn 11.35 (text)
       CS=CS+0.5e+0_fp*(Dpk(k)*Nko(k)/(dble(boxvol)*1.0e-6_fp)*(1+Kn)) &
            /(1.e+0_fp+0.377e+0_fp*Kn+1.33e+0_fp*Kn*(1+Kn))
       Ntot = Ntot + Nko(k)
       Dtot = Dtot + Nko(k)*Dpk(k)
       if (k.eq.1)then
          Dp1=Dpk(k)
          dens1 = density
       endif
    enddo

    if (Ntot.gt.1e+15_fp)then
       Dpmean = Dtot/Ntot
    else
       Dpmean = 150.e+0_fp
    endif

    return

  END SUBROUTINE GETCONDSINK_KERM
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getgrowthtime
!
! !DESCRIPTION: This subroutine calculates the time it takes for a particle to
!  grow from one size to the next by condensation of sulfuric acid (and
!  associated NH3 and water) onto particles.
!  .
!  This subroutine assumes that the growth happens entirely in the kinetic
!  regine such that the dDp/dt is not size dependent.  The time for growth
!  to the first size bin may then be approximated by the time for growth via
!  sulfuric acid (not including nh4 and water) to the size of the first size bin
!  (not including nh4 and water).
!  WRITTEN BY Jeff Pierce, April 2007 for GISS GCM-II'
!  Introduce to GEOS-Chem by Win Trivitayanurak (win, 9/29/08)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE getGrowthTime (d1,d2,h2so4,temp,boxvol,density,gtime)
!
! !USES:
!
    USE ERROR_MOD,      ONLY : ERROR_STOP, IT_IS_NAN
!
! !INPUT PARAMETERS:
!
    ! d1: intial diameter [m]
    ! d2: final diameter [m]
    ! h2so4: h2so4 ammount [kg]
    ! temp: temperature [K]
    ! boxvol: box volume [cm3]
    REAL(fp), INTENT(IN)  ::  d1,d2    ! initial and final diameters [m]
    REAL(fp), INTENT(IN)  ::  h2so4    ! h2so4 amount [kg]
    real(fp),   INTENT(IN)  ::  temp     ! temperature [K]
    real(fp),   INTENT(IN)  ::  boxvol  ! box volume [cm3]
    REAL(fp), INTENT(IN)  ::  density  ! density of particles in first bin [kg/m3]
!
! !OUTPUT PARAMETERS:
!
    ! gtime: the time it takes the particle to grow to first size bin [s]
    REAL(fp), INTENT(OUT) ::  gtime    ! the time it will take the particle to
                                       ! grow to first size bin [s]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)     ::  pi, R, MW
    REAL(fp)     ::  csulf    ! concentration of sulf acid [kmol/m3]
    REAL(fp)     ::  mspeed   ! mean speed of molecules [m/s]
    REAL(fp)     ::  alpha    ! accomidation coef
!
! !DEFINED PARAMETERS:
!
    parameter(pi=3.141592654e+0_fp, R=8.314e+0_fp) !pi and gas constant (J/mol K)
    parameter(MW=98.e+0_fp) ! density [kg/m3], mol wgt sulf [kg/kmol]
    parameter(alpha=0.65)

    !=================================================================
    ! GETGROWTHTIME begins here!
    !=================================================================
    !print *,'h2so4',h2so4,'MW',MW,'boxvol',boxvol,dble(boxvol)

    csulf = h2so4/MW/(dble(boxvol)*1e-6_fp) ! SA conc. [kmol/m3]
    mspeed = sqrt(8.e+0_fp*R*dble(temp)*1000.e+0_fp/(pi*MW))

    ! Kinetic regime expression (S&P 11.25) solved for T
    gtime = (d2-d1)/(4.e+0_fp*MW/density*mspeed*alpha*csulf)

    if ( IT_IS_NAN(gtime) ) then
       !jrp
       print*,'IN GET GROWTH TIME'
       print*,'d1',d1,'d2',d2
       print*,'h2so4',h2so4
       print*,'boxvol',boxvol
       print*,'csulf',csulf,'mspeed',mspeed
       print*,'density',density,'gtime',gtime
       call ERROR_STOP('Found NaN in fn','getnucrate')
    endif

    RETURN

  END SUBROUTINE GETGROWTHTIME
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ezcond
!
! !DESCRIPTION: This subroutine takes a given amount of mass and condenses it
!     across the bins accordingly.
!     WRITTEN BY Jeff Pierce, May 2007 for GISS GCM-II'
!     Put in GEOS-Chem by Win T. 9/30/08
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EZCOND (Nki,Mki,mcondi,spec,Nkf,Mkf,surf_area, &
                     BOXVOL, TEMPTMS, PRES, errswitch)
!
! !INPUT PARAMETERS:
!
    !Initial values of
    !=================
    !Nki(nBins) - number of particles per size bin in grid cell
    !Mki(nBins, ICOMPHARD) - mass of a given species per size bin/grid cell [kg]
    !mcond - mass of species to condense [kg/grid cell]
    !spec - the number of the species to condense
    double precision Nki(nBins), Mki(nBins, ICOMPHARD)
    double precision mcondi
    REAL(fp), INTENT(IN)       :: BOXVOL, TEMPTMS, PRES
    LOGICAL ERRSWITCH   ! signal error to outside

!
! !OUTPUT PARAMETERS:
!
    !Nkf, Mkf - same as above, but final values
    double precision Nkf(nBins), Mkf(nBins, ICOMPHARD)
    REAL(fp)           surf_area
    integer spec
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer i,j,k,c           ! counters
    double precision mcond
    double precision pi, R    ! pi and gas constant (J/mol K)
    double precision CS       ! condensation sink [s^-1]
    double precision sinkfrac(nBins+1) ! fraction of CS in size bin
    double precision Nk1(nBins), Mk1(nBins, ICOMPHARD)
    double precision Nk2(nBins), Mk2(nBins, ICOMPHARD)
    double precision madd     ! mass to add to each bin [kg]
    double precision maddp(nBins)    ! mass to add per particle [kg]
    double precision mconds ! mass to add per step [kg]
    integer          nsteps            ! number of condensation steps necessary
    integer          my_floor, my_ceil       ! the floor and ceiling (temporary)
    double precision eps     ! small number
    double precision tdt      !the value 2/3
    double precision mpo,mpw  !dry and "wet" mass of particle
    double precision WR       !wet ratio
    double precision tau(nBins) !driving force for condensation
    double precision totsinkfrac ! total sink fraction not including nuc bin
    double precision CSeps    ! lower limit for condensation sink
    double precision tot_m,tot_s    !total mass, total sulfate mass
    double precision ratio    ! used in mass correction
    double precision fracch(nBins,ICOMPHARD)
    double precision totch

    double precision tot_i,tot_f,tot_fa ! used for conservation of mass check
    LOGICAL          PDBG,  ERRORSWITCH
    real(fp)         zeros(nBins)
!
! !DEFINED PARAMETERS:
!
    parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
    parameter(eps=1.e-40_fp)
    parameter(CSeps=1.e-20_fp)

    !=================================================================
    ! EZCOND begins here
    !=================================================================

    pdbg = errswitch ! take the signal for print debug from outside
    errswitch = .false. !signal to terminate with error. Initialize with .false.

    tdt=2.e+0_fp/3.e+0_fp

    mcond=mcondi

    ! initialize variables
    do k=1,nBins
       Nk1(k)=Nki(k)
       do j=1,ICOMPHARD
          Mk1(k,j)=Mki(k,j)
       enddo
    enddo

    !print *, 'mnfix in tomas_mod:2804'
    call mnfix(Nk1,Mk1, errorswitch)
    if(errorswitch) then
       print *, 'EZCOND: MNFIX[1] found error --> TERMINATE'
       errswitch=.true.
       return
    endif

    ! get the sink fractions
    ! set Nnuc to zero for this calc
    call getCondSink(Nk1,Mk1,spec,CS,sinkfrac,surf_area, &
                     BOXVOL,TEMPTMS, PRES)

    ! make sure that condensation sink isn't too small
    if (CS.lt.CSeps) then     ! just make particles in first bin
       Mkf(1,spec) = Mk1(1,spec) + mcond
       Nkf(1) = Nk1(1) + mcond/sqrt(xk(1)*xk(2))
       do j=1,ICOMPHARD
          if (j.ne.spec) then ! Bug confirmed by Jeff, (ICOMPHARD.ne.spec) -> (j.ne.spec)
             Mkf(1,j) = Mk1(1,j)
          endif
       enddo
       do k=2,nBins
          Nkf(k) = Nk1(k)
          do j=1,ICOMPHARD
             Mkf(k,j) = Mk1(k,j)
          enddo
       enddo
       return
    endif

    if (pdbg) then
       print*,'CS',CS
       print*,'sinkfrac',sinkfrac
       print*,'mcond',mcond
    endif

    ! determine how much mass to add to each size bin
    ! also determine how many condensation steps we need
    totsinkfrac = 0.e+0_fp
    do k=1,nBins
       totsinkfrac = totsinkfrac + sinkfrac(k) ! get sink frac total not including nuc bin
    enddo
    nsteps = 1
    do k=1,nBins
       if (sinkfrac(k).lt.1.0e-20_fp)then
          madd = 0.e+0_fp
       else
          madd = mcond*sinkfrac(k)/totsinkfrac
       endif
       mpo=0.0
       do j=1,ICOMPHARD-2
          mpo=mpo + Mk1(k,j)
       enddo
       if(mpo == 0.0 ) then  ! prevent division by zero (win, 10/16/08)
          my_floor = 0
       else
          my_floor = int(madd*0.00001/mpo)
       endif
       my_ceil = my_floor + 1
       nsteps = max(nsteps,my_ceil) ! don't let the mass increase by more than 10%
    enddo

    if(pdbg) print*,'nsteps',nsteps

    ! mass to condense each step
    mconds = mcond/nsteps

    ! do steps of condensation
    do i=1,nsteps
       if (i.ne.1) then
          ! set Nnuc to zero for this calculation
          call getCondSink(Nk1,Mk1,spec,CS,sinkfrac,surf_area, &
                           BOXVOL,TEMPTMS, PRES)
          totsinkfrac = 0.e+0_fp
          do k=1,nBins
             totsinkfrac = totsinkfrac + sinkfrac(k) ! get sink frac total not including nuc bin
          enddo
       endif

       tot_m=0.e+0_fp
       tot_s=0.e+0_fp
       do k=1,nBins
          do j=1,ICOMPHARD-2
             tot_m = tot_m + Mk1(k,j)
             if (j.eq.srtso4) then
                tot_s = tot_s + Mk1(k,j)
             endif
          enddo
       enddo

       if (pdbg) print *,'tot_s ',tot_s,' tot_m ',tot_m

       ! change criteria to bigger amount (win, 9/30/08)
       if (mcond.gt.tot_m*5.0e-2_fp) then
          if (pdbg) print *,'Entering TMCOND '

          do k=1,nBins
             mpo=0.0
             mpw=0.0
             !WIN'S CODE MODIFICATION 6/19/06
             !THIS MUST CHANGED WITH THE NEW dmdt_int
             do j=1,ICOMPHARD-2
                mpo = mpo+Mk1(k,j) !accumulate dry mass
             enddo
             do j=1,ICOMPHARD
                mpw = mpw+Mk1(k,j) ! have wet mass include amso4
             enddo
             if( mpo > 0.0 ) then    ! prevent division by zero (win, 10/16/08)
                WR = mpw/mpo  !WR = wet ratio = total mass/dry mass
             else
                WR = 1.0
             endif
             if (Nk1(k) .gt. 0.e+0_fp) then
                !Change maddp(k) from mass/no. to be just mass (win,10/3/08)
                ! this is because in tmcond here, the moxd argument takes
                ! mass to add for each bin array, not mass/no. array.
                maddp(k) = mconds*sinkfrac(k)/totsinkfrac
                !Prior to 10/3/08 (win)
                !maddp(k) = mconds*sinkfrac(k)/totsinkfrac/Nk1(k)
                mpw=mpw/Nk1(k)

                if(pdbg) print*,'mpw',mpw,'maddp',maddp(k),'WR',WR
                !Change the maddp(k) to accordingly -- adding the /Nk1(k) (win, 10/3/08)
                tau(k)=1.5e+0_fp*((mpw+maddp(k)/Nk1(k)*WR)**tdt-mpw**tdt)
                ! Prior to 10/3/08 (win)
                !tau(k)=1.5e+0_fp*((mpw+maddp(k)*WR)**tdt-mpw**tdt) !added WR to moxid term (win, 5/15/06)
                !     tau(k)=0.e+0_fp
                !     maddp(k)=0.e+0_fp
             else
                !nothing in this bin - set tau to zero
                tau(k)=0.e+0_fp
                maddp(k) = 0.e+0_fp
             endif
          enddo
          !print*,'tau',tau
          !print *, 'mnfix in tomas_mod:2942'
          call mnfix(Nk1,Mk1, errorswitch)
          if (errorswitch) then
             print *, 'EZCOND: MNFIX[2] found error --> TERMINATE'
             errswitch=.true.
             return
          endif
          ! do condensation
          errorswitch = pdbg
          !prior to 9/30/08 from Jeff's version
          call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec,errorswitch,maddp)

          ! For SO4 condensation, the last argument should be zeroes (win, 9/30/08)
          !zeros(:) = 0.e+0_fp
          !call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec,errorswitch,zeros)

          if( errorswitch) then
             errswitch=.true.
             print *,'EZCOND: error after TMCOND --> TERMINATE'
             return
          endif
          errorswitch = pdbg

          !call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec)
          !jrp totch=0.0
          !jrp do k=1,ibins
          !jrp    do j=1,icomp
          !jrp       fracch(k,j)=(Mk2(k,j)-Mk1(k,j))
          !jrp       totch = totch + (Mk2(k,j)-Mk1(k,j))
          !jrp    enddo
          !jrp enddo
          !print*,'fracch',fracch,'totch',totch

       elseif (mcond.gt.tot_s*1.0e-12_fp) then
          if (pdbg) print *,'Small mcond: distrib w/ sinkfrac '
          if (pdbg) print *, 'maddp(bin) to add to SO4'
          do k=1,nBins
             if (Nk1(k) .gt. 0.e+0_fp) then
                maddp(k) = mconds*sinkfrac(k)/totsinkfrac
             else
                maddp(k) = 0.e+0_fp
             endif
             if(pdbg) print *, maddp(k)
             Mk2(k,srtso4)=Mk1(k,srtso4)+maddp(k)
             do j=1,ICOMPHARD
                if (j.ne.srtso4) then
                   Mk2(k,j)=Mk1(k,j)
                endif
             enddo
             Nk2(k)=Nk1(k)
          enddo
          if(pdbg) errorswitch = .true.

          !print *, 'mnfix in tomas_mod:2999'
          call mnfix(Nk2,Mk2, errorswitch)
          if(errorswitch) then
             print *, 'EZCOND: MNFIX[3] found error --> TERMINATE'
             errswitch=.true.
             return
          endif
       else ! do nothing
          if (pdbg) print *,'Very small mcond: do nothing!'
          mcond = 0.e+0_fp
          do k=1,nBins
             Nk2(k)=Nk1(k)
             do j=1,ICOMPHARD
                Mk2(k,j)=Mk1(k,j)
             enddo
          enddo
       endif
       if (i.ne.nsteps)then
          do k=1,nBins
             Nk1(k)=Nk2(k)
             do j=1,ICOMPHARD
                Mk1(k,j)=Mk2(k,j)
             enddo
          enddo
       endif

    enddo

    do k=1,nBins
       Nkf(k)=Nk2(k)
       do j=1,ICOMPHARD
          Mkf(k,j)=Mk2(k,j)
       enddo
    enddo

    ! check for conservation of mass
    tot_i = 0.e+0_fp
    tot_fa = mcond
    tot_f = 0.e+0_fp
    do k=1,nBins
       tot_i=tot_i+Mki(k,srtso4)
       tot_f=tot_f+Mkf(k,srtso4)
       tot_fa=tot_fa+Mki(k,srtso4)
    enddo

    if(pdbg) then
       print *,'Check conserv of mass after mcond is distrib'
       print *,' Initial total so4 ',tot_i
       print *,' Final total so4   ',tot_f
       print *,'Percent error=',abs((mcond-(tot_f-tot_i))/mcond)*1e2
    endif

    if ( mcond > 0.0_fp ) then
       if ( abs((mcond-(tot_f-tot_i))/mcond).gt.0.e+0_fp) then
          IF(mcond > 1.e-8_fp .and. tot_i > 5.e-2_fp)  THEN
             !Add a check to check error if mcond is significant (win, 10/2/08)

             IF (abs((mcond-(tot_f-tot_i))/mcond).lt.1.e+0_fp .OR. &
                  spinup(31.0) ) THEN
                !Prior to 10/2/08 (win)   .. original was Jeff's fix
                !! do correction of mass
                !ratio = (tot_f-tot_i)/mcond
                !if(pdbg) print *,'Mk at mass correction '
                !if(pdbg) print *,'  ratio',ratio
                !do k=1,ibins
                !   Mkf(k,srtso4)=Mki(k,srtso4)+
                !   &              (Mkf(k,srtso4)-Mki(k,srtso4))/ratio
                !   if(pdbg) print *,Mkf(k,srtso4)
                !enddo

                ! Do mass correction (win, 10/2/08)
                ratio = (tot_i+mcond)/tot_f
                if(pdbg) print *,'Mk at mass correction apply ratio= ',ratio
                do k=1,nBins
                   Mkf(k,srtso4)=Mkf(k,srtso4) * ratio
                   if(pdbg) print *,Mkf(k,srtso4)
                enddo

                if(pdbg) errorswitch=.true.
                !print *, 'mnfix in tomas_mod:3079'
                call mnfix(Nkf,Mkf, errorswitch)
                if(errorswitch) then
                   print *, 'EZCOND: MNFIX[4] found error --> TERMINATE'
                   errswitch=.true.
                   return
                endif
             else
                print*,'ERROR in ezcond'
                print*,'Condensation error',(mcond-(tot_f-tot_i))/mcond
                print*,'mcond',mcond,'change',tot_f-tot_i
                print*,'tot_i',tot_i,'tot_fa',tot_fa,'tot_f',tot_f
                print*,'Nki',Nki
                print*,'Nkf',Nkf
                print*,'Mki',Mki
                print*,'Mkf',Mkf
                !Prior to 10/2/08 (win)
                !STOP
                ! Send error signal to outside and terminate with more info
                ! (win, 10/2/08)
                !!as of 10/27/08, try comment out this signal to stop the run
                ! (win, 10/27/08)
                !! the problem is that maybe or mostly the mass conservation is
                ! ruined becuase of the fudging inside mnfix.
                !ERRSWITCH=.TRUE.
                !RETURN
             ENDIF
          ENDIF
       endif
    endif

    !jrp if (abs(tot_f-tot_fa)/tot_i.gt.1.0D-8)then
    !jrp    print*,'No S conservation in ezcond'
    !jrp    print*,'initial',tot_fa
    !jrp    print*,'final',tot_f
    !jrp    print*,'mcond',mcond,'change',tot_f-tot_i
    !jrp    print*,'ERROR',(mcond-(tot_f-tot_i))/mcond
    !jrp endif

    ! check for conservation of mass
    tot_i = 0.e+0_fp
    tot_f = 0.e+0_fp
    do k=1,nBins
       tot_i=tot_i+Mki(k,srtnh4)
       tot_f=tot_f+Mkf(k,srtnh4)
    enddo
    if (.not. spinup(14.0)) then
       if (abs(tot_f-tot_i)/tot_i.gt.1.0e-8_fp)then
          if ( tot_i > 1.0e-20_fp ) then
             print*,'No N conservation in ezcond'
             print*,'initial',tot_i
             print*,'final  ',tot_f
          endif
       endif
    endif

    return

  end SUBROUTINE EZCOND
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tmcond
!
! !DESCRIPTION: Subroutine TMCOND do condensation calculation.
!  Original code from Peter Adams
!  Modified for GEOS-CHEM by Win Trivitayaurak (win@cmu.edu)
!  CONDENSATION
!   Based on Tzivion, Feingold, Levin, JAS 1989 and
!   Stevens, Feingold, Cotton, JAS 1996
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TMCOND(TAU,X,AMKD,ANKD,AMK,ANK,CSPECIES,pdbug,moxd)
!
! !INPUT PARAMETERS:
!
    ! TAU(k) ...... Forcing for diffusion = (2/3)*CPT*ETA_BAR*DELTA_T
    ! X(K) ........ Array of bin limits in mass space
    ! AMKD(K,J) ... Input array of mass moments
    ! ANKD(K) ..... Input array of number moments
    ! CSPECIES .... Index of chemical species that is condensing
    REAL(fp)       :: TAU(nBins)
    REAL(fp)       :: X(nBins+1),AMKD(nBins,ICOMPHARD),ANKD(nBins)
    INTEGER        :: CSPECIES
    LOGICAL        :: pdbug !(win, 4/10/06)
    REAL(fp)       :: moxd(nBins) ! condensing mass distributed to size bins
                       ! according to the selected absorbing media (win, 3/5/08)
!
! !OUTPUT PARAMETERS:
!
    ! AMK(K,J) .... Output array of mass moments
    ! ANK(K) ...... Output array of number moments
    REAL(fp)       :: AMK(nBins,ICOMPHARD),ANK(nBins)
!
! !REMARKS:
! The supersaturation is calculated outside of the routine and assumed
! to be constant at its average value over the timestep.
! .
! The method has three basic components:
! (1) first a top hat representation of the distribution is construced
!     in each bin and these are translated according to the analytic
!     solutions
! (2) The translated tophats are then remapped to bins.  Here if a
!     top hat entirely or in part lies below the lowest bin it is
!     not counted.
!     .
! Additional notes (Peter Adams)
!     .
!     I have changed the routine to handle multicomponent aerosols.  The
!     arrays of mass moments are now two dimensional (size and species).
!     Only a single component (CSPECIES) is allowed to condense during
!     a given call to this routine.  Multicomponent condensation/evaporation
!     is accomplished via multiple calls.  Variables YLC and YUC are
!     similar to YL and YU except that they refer to the mass of the
!     condensing species, rather than total aerosol mass.
!     .
!     I have removed ventilation variables (VSW/VNTF) from the subroutine
!     call.  They still exist internally within this subroutine, but
!     are initialized such that they do nothing.
!     .
!     I have created a new variable, AMKDRY, which is the total mass in
!     a size bin (sum of all chemical components excluding water).  I
!     have also created WR, which is the ratio of total wet mass to
!     total dry mass in a size bin.
!     .
!     AMKC(k,j) is the total amount of mass after condensation of species
!     j in particles that BEGAN in bin k.  It is used as a diagnostic
!     for tracking down numerical errors.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER        :: L,I,J,K,IMN
    REAL(fp)       :: DN,DM,DYI,XL,XU,YL,YLC,YU,YUC
    REAL(fp)       :: TEPS,NEPS,NEPS2,EX2,ZERO
    REAL(fp)       :: XI,XX,XP,YM,WTH,W1,W2,WW,AVG
    REAL(fp)       :: VSW,VNTF(nBins)
    REAL(fp)       :: TAU_L, maxtau

    REAL(fp)       :: AMKDRY(nBins), WR(nBins), AMKWET(nBins)
    REAL(fp)       :: AMKDRYSOL(nBins)

    LOGICAL        :: errspot !(win, 4/12/06)

    REAL(fp)       :: c1, c2 !correction factor (win, 5/25/06)
    REAL(fp)       :: madd(nBins) !condensing mass to be added by aqoxid
                                  !or SOAcond. For error fixing (win, 9/27/07)
    REAL(fp)       :: xadd(nBins) !mass per particle to be added by aqoxid
                                  ! or SOAcond. For error fixing (win, 9/27/07)
    REAL(fp)       :: macc !accumulating the condensing mass (win, 7/24/06)
    REAL(fp)       :: delt1,delt2 !the delta = mass not conserved (win, 7/24/06)
    REAL(fp)       :: dummy, xtra,maddtot ! for mass conserv fixing (win, 9/27/07)
    integer        :: kk !counter (wint, 7/24/06)
    REAL(fp)       :: AMKD_tot

    PARAMETER (TEPS=1.0e-40_fp,NEPS=1.0e-20_fp)
    PARAMETER (EX2=2.e+0_fp/3.e+0_fp,ZERO=0.0e+0_fp)
    PARAMETER (NEPS2=1.0e-10_fp)

    !=================================================================
    ! TMCOND begins here!
    !=================================================================

3   format(I4,200E20.11)

    !<step4.5> This first check cause the error of 'number not conserved'
    ! though only with the small amounts because when ANKD(k) = 0.e+0_fp from start,
    ! the original check just give it a value NEPS = 1.d-20, and then undergo
    ! tmcond calculation.   I'm changing the check to if ANKD(k)= 0.e+0_fp,
    ! then keep it that way and make the following calculations skip when
    ! ANKD(k) is zero (win, 10/18/05)

    ! If any ANKD are zero, set them to a small value to avoid division by zero
    !do k=1,ibins
    !   if (ANKD(k) .lt. NEPS) then
    !      ANKD(k)=NEPS
    !      AMKD(k,srtso4)=NEPS*1.4*xk(k) !make the added particles SO4
    !   endif
    !enddo

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !<step5.1> Add print for debugging (win, 4/10/06)
    if (pdbug) then
       ! call debugprint(ANKD, AMKD, 0,0,0,'Entering TMCOND')
       ! print *, 'TMCOND:entering*************************'
       ! print *,'Nk(1:30)'
       ! print *, ANKD(1:30)
       ! print *,'Mk(1:30,comp)'
       ! do j=1,icomp
       ! print *,'comp',j
       ! print *, AMKD(1:30,j)
       ! enddo
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    errspot = .false. !initialize error signal as false (Win, 4/12/06)

    !pja Sometimes, after repeated condensation calls, the average bin mass
    !pja can be just above the bin boundary - in that case, transfer a some
    !pja to the next highest bin
    !sfarina this is also true when small particles are growing really fast?
    !sfarina SOACOND throws thousands of errors for XI < 1
    !sfarina what this really means is AVG particle massfor bin k > XK(k+1)
    !sfarina through the debugger I found that mostly the difference is small
    do k=1,nBins-1
       if ( ANKD(k) .lt. NEPS2) goto 300 !<step4.5> (win, 10/18/05)
       ! Modify the check to include all dry mass (win, 10/3/08)
       AMKD_tot = 0.e+0_fp
       do kk=1,ICOMPHARD-2
          AMKD_tot = AMKD_tot + AMKD(k,kk)
       enddo
       if ((AMKD_tot)/ANKD(k).gt.xk(k+1)) then
          !Prior to 10/3/08 (win)
          !if ((AMKD(k,srtso4))/ANKD(k).gt.xk(k+1)) then
          !sfarina: this does noting to help our avg mass per particle
          !         falling outside of bin boundaries:
          !         amkd_tot / ankd(k) = (amkd_tot * 0.9) / (ankd(k) * 0.9)
          !         we need to shift more mass than number.
          !         assuming we have some kind of distributionof particle sizes in bin K
          !         the largest ones will have more mass than average, so we can safly move
          !         more mass than number.
          !         that or we redistribute mass before SOAcond
          !
          do j=1,ICOMPHARD-2
             AMKD(k+1,j)=AMKD(k+1,j)+0.1e+0_fp*AMKD(k,j)
             AMKD(k,j)=AMKD(k,j)*0.9e+0_fp
          enddo
          ANKD(k+1)=ANKD(k+1)+0.1e+0_fp*ANKD(k)
          ANKD(k)=ANKD(k)*0.9e+0_fp
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !<step5.1> Add print for debugging (win, 4/10/06)
          if (pdbug) then
             print *, 'Modified at checkpoint1: BIN',k
             print *,'ANKD(k)',ANKD(k),'ANKD(k+1)',ANKD(k+1)
             print *,'Mk(k,comp)       Mk(k+1,comp)'
             do j=1,ICOMPHARD
                print *,'comp',j
                print *, AMKD(k,j), AMKD(k+1,j)
             enddo
          endif
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       endif
300    continue   !<step4.5> If aerosol number is zero (win, 10/18/05)

    enddo

    !pja Initialize ventilation variables so they don't do anything
    VSW=0.0e+0_fp
    DO L=1,nBins
       VNTF(L)=0.0e+0_fp
    ENDDO

    !pja Initialize AMKDRY and WR
    DO L=1,nBins
       AMKDRY(L)=0.e+0_fp
       AMKWET(L)=0.e+0_fp
       AMKDRYSOL(L) = 0.e+0_fp
       DO J=1,ICOMPHARD-2     ! dry mass excl. nh4 (win, 9/26/08)
          AMKDRY(L)=AMKDRY(L)+AMKD(L,J)
          ! Accumulate the absorbing media (win, 3/5/08)
          IF ( J == SRTOCIL  ) &
               AMKDRYSOL(L) = AMKDRYSOL(L) + AMKD(L,J)
       ENDDO
       DO J=1,ICOMPHARD
          AMKWET(L) = AMKWET(L) + AMKD(L,J)
       ENDDO
       if (AMKDRY(L) .gt. 0.e+0_fp) &   !<step4.5> In case there is no mass, then just skip (win, 10/18/05)
            WR(L)= AMKWET(L) / AMKDRY(L)
    ENDDO

    !debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(pdbug)then
       print*,'AMKDRY(1:nBins)'
       print *,AMKDRY(1:nBins)
       print *,'WR(1:nBins)'
       print *,WR(1:nBins)
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !pja Initialize X() array of particle masses based on xk()
    DO L=1,nBins
       X(L)=xk(L)
    ENDDO

    !
    ! Only solve when significant forcing is available
    !
    maxtau=0.0e+0_fp
    do l=1,nBins
       maxtau=max(maxtau,abs(TAU(l)))
    enddo

    !debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(pdbug) then
       print*,'tau(1:nBins)'
       print *,tau(1:nBins)
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    IF(ABS(maxtau).LT.TEPS)THEN
       DO L=1,nBins
          DO J=1,ICOMPHARD
             AMK(L,J)=AMKD(L,J)
          ENDDO
          ANK(L)=ANKD(L)
       ENDDO
    ELSE
       !<step5.3> Try to fix the error of mass conservation
       ! during aqueous oxidation. Too little mass is used up
       ! (win, 7/24/06)
       IF ( MAXVAL(MOXD(:)) >  0e+0_fp ) THEN
          IF( PDBUG ) PRINT *,'Mass_to_add_by_aqoxid_or_SOAcond'
          maddtot = 0e+0_fp
          DO L = 1, nBins
             IF(TAU(L) >  0e+0_fp ) THEN
                MADD(L) = MOXD(L)
                XADD(L) = MOXD(L) / ANKD(L)
                !IF( CSPECIES == SRTSO4 ) THEN
                !   MADD(L) = MOXD * ANKD(L)  ! absolute condensing mass
                !   XADD(L) = MOXD            ! mass per particle
                !ELSE IF ( CSPECIES == SRTOCIL ) THEN
                !   MADD(L) = MOXD * AMKDRYSOL(L)
                !   XADD(L) = MADD(L) / ANKD(L)
                !ELSE
                !   PRINT *,'TMCOND ERROR : mass fixing not supported'
                !ENDIF
             ELSE
                MADD(L) = 0e+0_fp
                XADD(L) = 0e+0_fp
             ENDIF
             IF ( PDBUG ) PRINT *,L,madd(L), xadd(L)
             maddtot = maddtot + madd(L)
          ENDDO
       ENDIF

       DO L=1,nBins
          DO J=1,ICOMPHARD
             AMK(L,J)=0.e+0_fp
          ENDDO
          ANK(L)=0.e+0_fp
       ENDDO
       WW=0.5e+0_fp
       ! IF(TAU.LT.0.)WW=.5e+0_fp
       !
       ! identify tophats and do lagrangian growth
       !
       DO L=1,nBins
          IF(ANKD(L) .LT. NEPS2)GOTO 200 !skip if Number is effectively zero

          !if tau is zero, leave everything in same bin
          IF (TAU(L) .EQ. 0.) THEN
             ANK(L)=ANK(L)+ANKD(L)
             DO J=1,ICOMPHARD
                AMK(L,J)=AMK(L,J)+AMKD(L,J)
             ENDDO
          ENDIF
          IF (TAU(L) .EQ. 0.) GOTO 200

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !<step5.1> Add print for debugging (win, 4/10/06)
          if (pdbug) then
             print *, 'Identify_tophat_and_grow-BIN',L
             print *,'Starting_Nk(1:nBins)'
             print *, ANK(1:nBins)
             print *,'Starting_Mk(1:nBins,comp)'
             do j=1,ICOMPHARD-1
                print *,'comp',j
                print *, AMK(1:nBins,j)
             enddo
          endif
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          !pja Limiting AVG, the average particle size to lie within the size
          !pja bounds causes particles to grow or shrink arbitrarily and is
          !pja wreacking havoc with choosing condensational timesteps and
          !pja conserving mass.  I have turned them off.
          !AVG=MAX(X(L),MIN(X(L+1),AMKDRY(L)/(NEPS+ANKD(L))))
          !try bring the above line back, win 4/10/06
          !win 4/10/06

          AVG=AMKDRY(L)/ANKD(L)
          XX=X(L)/AVG

#if defined(TOMAS12) || defined(TOMAS15)
          if(l.lt.nBins-1)then ! bin quadrupuling
             XI=.5e+0_fp + XX*(2.5e+0_fp - 2.0e+0_fp*XX)
             !XI<1 means the AVG falls out of bin bounds
             if (XI .LT. 1.e+0_fp) then
                !W1 will have sqrt of negative number
                write(*,*)'ERROR: tmcond - XI<1 for bin: ',L
                write(*,*)'AVG is ',AVG
                write(*,*)'Nk is ', ANKD(L)
                write(*,*)'Mk are ', (AMKD(L,j),j=1,ICOMPHARD)
                write(*,*)'Initial N and M are: ',ANKD(L),AMKDRY(L)
                errspot = .true.
                RETURN
             endif
             W1 =SQRT(12.e+0_fp*(XI-1.e+0_fp))*AVG/4.0e+0_fp ! cyhl 4.0=xk(k+1)/xk(k)
             W2 =(MIN(X(L+1)-AVG,AVG-X(L)))*2.0e+0_fp
          else ! final 2 bins mass*32
             XI=.5e+0_fp + XX*(16.5e+0_fp - 16.0e+0_fp*XX)
             if (XI .LT. 1.e+0_fp) then
                !W1 will have sqrt of negative number
                write(*,*)'ERROR: tmcond - XI<1 for bin: ',L
                write(*,*)'lower limit is',X(L)
                write(*,*)'AVG is ',AVG
                write(*,*)'Nk is ', ANKD(L)
                write(*,*)'Mk are ', (AMKD(L,j),j=1,ICOMPHARD)
                write(*,*)'Initial N and M are: ',ANKD(L),AMKDRY(L)
                errspot = .true.
                RETURN
             endif
             W1 =SQRT(12.e+0_fp*(XI-1.e+0_fp))*AVG/32.0e+0_fp ! cyhl 32.0=xk(k+1)/xk(k)
             W2 =(MIN(X(L+1)-AVG,AVG-X(L)))*2.0e+0_fp
          endif
#else
          XI=.5e+0_fp + XX*(1.5e+0_fp - XX)
          !XI<1 means the AVG falls out of bin bounds

          if (XI .LT. 1.e+0_fp) then
             !W1 will have sqrt of negative number
             write(*,*)'ERROR: tmcond - XI<1 for bin: ',L
             write(*,*)'AVG is ',AVG
             write(*,*)'Nk is ', ANKD(L)
             write(*,*)'Mk are ', (AMKD(L,j),j=1,ICOMPHARD)
             write(*,*)'Initial N and M are: ',ANKD(L),AMKDRY(L)
             errspot = .true.
             RETURN
          endif
          W1 =SQRT(12.e+0_fp*(XI-1.e+0_fp))*AVG
          W2 =MIN(X(L+1)-AVG,AVG-X(L))
#endif

          WTH=W1*WW+W2*(1.e+0_fp-WW)
          IF(WTH.GT.1.) then
             write(*,*)'WTH>1 in cond, bin #',L
             errspot = .true.
             RETURN
          ENDIF

          XU=AVG+WTH*.5e+0_fp
          XL=AVG-WTH*.5e+0_fp
          ! Ventilation added bin-by-bin
          TAU_L=TAU(l)*MAX(1.e+0_fp,VNTF(L)*VSW)
          IF(TAU_L/TAU(l).GT. 6.) THEN
             PRINT *,'TAU..>6.',TAU(l),TAU_L,VSW,L
          ENDIF
          IF(TAU_L.GT.TAU(l)) THEN
             PRINT *,'TAU...',TAU(l),TAU_L,VSW,L
          ENDIF
          ! prior to 5/25/06 (win)
          !YU=DMDT_INT(XU,TAU_L,WR(L))
          !YUC=XU*AMKD(L,CSPECIES)/AMKDRY(L)+YU-XU
          !IF (YU .GT. X(ibins+1) ) THEN
          !   YUC=YUC*X(ibins+1)/YU
          !   YU=X(ibins+1)
          !ENDIF
          !YL=DMDT_INT(XL,TAU_L,WR(L))
          !YLC=XL*AMKD(L,CSPECIES)/AMKDRY(L)+YL-XL
          !add new correction factor to YU and YL (win, 5/25/06)
          YU=DMDT_INT(XU,TAU_L,WR(L))
          YL=DMDT_INT(XL,TAU_L,WR(L))

          ! change to check MOXD of current bin (win, 10/3/08)
          IF( MOXD(L) == 0e+0_fp) THEN
             !Prior to 10/3/08 (win)
             !IF( MAXVAL(MOXD(:)) == 0e+0_fp ) THEN
             C1=1.e+0_fp          !for so4cond call, without correction factor.
          ELSE
             C1 = XADD(L)*2.e+0_fp/(YU+YL-XU-XL)
          ENDIF
          C2 = C1 - ( C1 - 1.e+0_fp ) * ( XU + XL )/( YU + YL )
          !prior to 10/2/08 (win)
          YU = YU * C2
          YL = YL * C2
          ! Run into a problem that YU < XU creating YUC<0
          ! So let's limit the application of C2 to only if
          ! it does not result in YU < XU and YL < XL (win, 10/2/08)
          !IF(TAU_L > 0.e+0_fp) YU = max( YU*C2, XU )
          !IF(TAU_L > 0.e+0_fp) YL = max( YL*C2, XL )

          !end part for fudging to get higher AVG

          YUC=XU*AMKD(L,CSPECIES)/AMKDRY(L)+YU-XU
          IF (YU .GT. X(nBins+1) ) THEN
             !IF(.not.SPINUP(60.)) write(116,*) &
             !     'YU > Xk(30+1) ++++++++++++' !debug (win, 7/17/06)
             YUC=YUC*X(nBins+1)/YU
             YU=X(nBins+1)
             !errspot=.true.  !just try temp (win, 7/30/07)
          ENDIF
          YLC=XL*AMKD(L,CSPECIES)/AMKDRY(L)+YL-XL
          DYI=1.e+0_fp/(YU-YL)

          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !<step5.2> Debug why there is extra mass added when called
          ! by aqoxid. (win, 5/10/06)
          if (pdbug) then
             print *, 'XU',XU,'YU',YU,'YUC',YUC,'c2',c2
             print *, 'XL',XL,'YL',YL,'YLC',YLC
          endif
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          !deal with tiny negative (win, 5/28/06)
          if(YUC.lt.0e+0_fp .or. YLC.lt.0e+0_fp)then
             if(YLC.lt.0e+0_fp) YLC=0e+0_fp
             if(YUC.lt.0e+0_fp) then
                YUC = 0e+0_fp
                YLC = 0e+0_fp
             endif
             if(pdbug) print *,'Fudge negative YUC, YLC to zero'
          endif
          !
          ! deal with portion of distribution that lies below lowest gridpoint
          !
          IF(YL.LT.X(1))THEN

             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             !<step5.2> Debug step-by-step (win, 5/10/06)
             if (pdbug) print *,'YL<X(1)_Just_condensing_to_current_bin'
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             !pja Instead of the following, I will just add all new condensed
             !pja mass to the same size bin
             !if ((YL/XL-1.e+0_fp) .LT. 1.e-3_fp) then
             !   !insignificant growth - leave alone
             !   ANK(L)=ANK(L)+ANKD(L)
             !   DO J=1,icomp-1
             !      AMK(L,J)=AMK(L,J)+AMKD(L,J)
             !   ENDDO
             !   GOTO 200
             !else
             !   !subtract out lower portion
             !   write(*,*)'ERROR in cond - low portion subtracted'
             !   write(*,*) 'Nk,Mk: ',ANKD(L),AMKD(L,1),AMKD(L,2)
             !   write(*,*) 'TAU: ', TAU_L
             !   write(*,*) 'XL, YL, YLC: ',XL,YL,YLC
             !   write(*,*) 'XU, YU, YUC: ',XU,YU,YUC
             !   ANKD(L)=ANKD(L)*MAX(ZERO,(YU-X(1)))*DYI
             !   YL=X(1)
             !   YLC=X(1)*AMKD(1,CSPECIES)/AMKDRY(1)
             !   DYI=1.e+0_fp/(YU-YL)
             !endif
             ANK(L)=ANK(L)+ANKD(L)
             do j=1,ICOMPHARD
                if (J.EQ.CSPECIES) then
                   AMK(L,J)=AMK(L,J)+(YUC+YLC)*.5e+0_fp*ANKD(L)
                else
                   AMK(L,J)=AMK(L,J)+AMKD(L,J)
                endif
             enddo
             GOTO 200
          ENDIF
          IF(YU.LT.X(1))GOTO 200
          !
          ! Begin remapping (start search at present location if condensation)
          !
          IMN=1
          IF(TAU(l).GT.0.)IMN=L
          DO I=IMN,nBins
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             !<step5.2> Debug step-by-step (win, 5/10/06)
             if(pdbug) print *,'Now_remapping_in_bin',I
             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             IF(YL.LT.X(I+1))THEN
                ![1] lower bound of new tophat in the current I bin
                IF(YU.LE.X(I+1))THEN
                   ![2] upper bound of new tophat also in the current I bin
                   DN=ANKD(L)      ! DN = number from the bin L being remapped
                   do j=1,ICOMPHARD
                      DM=AMKD(L,J)
                      IF (J.EQ.CSPECIES) THEN
                         !Add mass from new tophat to the existing mass of bin I
                         AMK(I,J)=(YUC+YLC)*.5e+0_fp*DN+AMK(I,J)
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         !<step5.2> Debug step-by-step (win, 5/10/06)
                         if (pdbug) then
                            print *,'CASE_1:_New_Tophat_in_a_single_bin'
                            print *,'SO4_from_tophat=',(YUC+YLC)*.5e+0_fp*DN
                         endif
                         !<step5.3> Check mass conservation (win, 7/24/06)
                         if(MAXVAL(moxd(:)).gt.0e+0_fp)then
                            delt1 = (YUC+YLC)*.5e+0_fp*DN-AMKD(L,J)-madd(L)
                            if( abs(delt1)/madd(L).gt.1e-6_fp .and. &
                                 madd(L).gt.1e-4_fp)then
                               ! Just print out this for debugging
                               IF(.not.SPINUP(60.) .and. pdbug ) then
                                  !write(116,*)'CASE1_mass_conserv_fix'
                                  write(116,13) L, madd(L), delt1
13                                FORMAT('CASE_1 Bin ',I2,' moxid ', &
                                         E13.5,' delta ',E13.5 )
                                  !errspot=.true. !just try temp (win, 7/30/07)
                               ENDIF
                               AMK(I,J) = AMK(I,J)-delt1 !fix the error
                            endif
                         endif
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      ELSE
                         !For non-condensing, migrate the mass to bin I
                         AMK(I,J)=AMK(I,J)+DM
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         !<step5.2> Debug step-by-step (win, 5/10/06)
                         if (pdbug) then
                            !print *,' Migrating_mass(',j,')',DM   !use this debugging line if there are more than seasalt+so4
                            print *,'Migrating_mass',DM
                         endif
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      ENDIF
                   enddo
                   !Add number of old bin to ANK (which is blank for the first loop of bin I)
                   ANK(I)=ANK(I)+DN
                   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   !<step5.2> (win, 5/10/06)
                   if(pdbug) print*,'Migrating_number',DN
                   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ELSE
                   ![3] upper bound of new tophat grow beyond the upper bound of bin I
                   DN=ANKD(L)*(X(I+1)-YL)*DYI !DN= proportion of the number from tophat that still stays in the bin I
                   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   !<step5.2> (win, 5/10/06)
                   if ( pdbug) then
                      print*,'Case_2:_Tophat_cross_bin_boundary'
                      print *,'Number_that_remain_in_low_bin',DN
                   endif
                   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                   !<step5.3> For fixing mass conserv problem (win, 7/24/06)
                   macc=0e+0_fp

                   do j=1,ICOMPHARD
                      !DM= proporation of the mass that is still in bin I
                      DM=AMKD(L,J)*(X(I+1)-YL)*DYI
                      IF (J.EQ.CSPECIES) THEN
                         !XP= what would have grown to be X(I+1)
                         XP=DMDT_INT(X(I+1),-1.0e+0_fp*TAU_L,WR(L))
                         YM=XP*AMKD(L,J)/AMKDRY(L)+X(I+1)-XP
                         !add the condensing mass to the existing sulfate of bin I
                         AMK(I,J)=DN*(YM+YLC)*0.5e+0_fp+AMK(I,J)
                         !<step5.3>Accumulating the condensing mass for error check (win, 7/24/06)
                         macc = macc + DN*(YM+YLC)*0.5e+0_fp
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         !<step5.2> (win, 5/10/06)
                         if(pdbug)then
                            print *,'XP',XP,'YM',YM
                            print *,'Cond_TophatLowEnd',DN*(YM+YLC)*0.5e+0_fp
                         endif
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      ELSE
                         !Add DM to AMK (which is blank for the first loop of bin I)
                         AMK(I,J)=AMK(I,J)+DM
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         if(pdbug) print*,'Other___in_low_end',DM
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      ENDIF
                   enddo
                   ANK(I)=ANK(I)+DN ! Add DN number to ANK (which is blank for the first loop of bin I)
                   ! Remapping loop from bin I+1 to bin30
                   DO K=I+1,nBins
                      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      if(pdbug) print *,'Spreading_to_bin',K
                      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      IF(YU.LE.X(K+1))GOTO 100
                      ![4] Found the bin where the high end of the tophat is in --> do the final loop

                      ![5.1] This part for distributing to the bins in between
                      !      the original and the furthest bin that growing occurs

                      !Use width of bin K to proportionate number from old bin wrt. to the top hat (YU-YL)
                      DN=ANKD(L)*(X(K+1)-X(K))*DYI

                      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      if(pdbug) then
                         print *,'Number_migrated',DN
                      endif
                      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      do j=1,ICOMPHARD
                         !Proportion of old-bin mass that falls in this current bin K
                         DM=AMKD(L,J)*(X(K+1)-X(K))*DYI
                         IF (J.EQ.CSPECIES) THEN
                            XP=DMDT_INT(X(K),-1.0e+0_fp*TAU_L,WR(L)) !what would have grown to be X(k)
                            YM=XP*AMKD(L,J)/AMKDRY(L)+X(K)-XP !what would have grown to be X(k) but just for sulfate
                            AMK(K,J)=DN*1.5e+0_fp*YM+AMK(K,J)    ! A factor of 1.5 is from averaging (YM+2*YM)
                            !<step5.3> Accumulating condensing mass for error check (win, 7/24/06)
                            macc = macc+DN*1.5e+0_fp*YM
                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            !<step5.2> (win, 5/10/06)
                            if(pdbug)then
                               print *,'XP',XP,'YM',YM
                               print *,'Cond_mass_spread',DN*1.5e+0_fp*YM
                            endif
                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         ELSE
                            AMK(K,J)=AMK(K,J)+DM    !Add migrating mass of non-condensing species
                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if(pdbug) print*,'No-cond_mass_migrate',DM
                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         ENDIF
                      enddo
                      ANK(K)=ANK(K)+DN  !Add migrating number to the exising number of bin K
                   ENDDO
                   !This STOP is for when there's excessive growth over bin30
                   STOP 'Trying to put stuff in bin ibins+1'

100                CONTINUE
                   ![5.2] Final section that the tophat grows to.
                   DN=ANKD(L)*(YU-X(K))*DYI
                   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   if(pdbug) then
                      print *,'Found_right_edge_for_tophat'
                      print *,'Number_migrated',DN
                   endif
                   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   do j=1,ICOMPHARD
                      DM=AMKD(L,J)*(YU-X(K))*DYI  ! proportion of old mass that gets to this furthest bin.
                      IF (J.EQ.CSPECIES) THEN
                         XP=DMDT_INT(X(K),-1.0e+0_fp*TAU_L,WR(L))   !what would have grown to be X(k)
                         YM=XP*AMKD(L,J)/AMKDRY(L)+X(K)-XP !=XP for just sulfate
                         AMK(K,J)=DN*(YUC+YM)*0.5e+0_fp+AMK(K,J) !add condensing mass to existing sulfate of bin K
                         !<step5.3>Accumulating condensing mass for error check (win, 7/24/06)
                         macc = macc+DN*(YUC+YM)*0.5e+0_fp
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         !<step5.2> (win, 5/10/06)
                         if(pdbug)then
                            print *,'XP',XP,'YM',YM
                            print *,'Cond_mass_spread_final', &
                                     DN*(YM+YUC)*0.5e+0_fp
                         endif
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      ELSE
                         AMK(K,J)=AMK(K,J)+DM  !This adds the migrating mass to the exising mass of non-condensing species
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                         if(pdbug) print*,'No-cond_mass_migrated',DM
                         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      ENDIF
                   enddo
                   ANK(K)=ANK(K)+DN   !This adds the migrating number to the existing number of bin K

                   !<step5.3> Check mass conservation (win, 7/24/06)
                   if(MAXVAL(moxd(:)).gt.0e+0_fp)then
                      delt2 = 0e+0_fp
                      delt2 = macc-AMKD(L,CSPECIES)-madd(L)
                      if(abs(delt2)/ madd(L) > 1e-6)then
                         if( madd(L) > 10.e+0_fp .and. &
                             abs(delt2)/ madd(L) > 15e-2_fp ) then
                            !print *,'TMCOND ERROR: mass condensation', &
                            !  'discrep >15% during aqoxid or SOAcond'
                            IF(.not.SPINUP(60.))  THEN
14                             FORMAT('CASE_2 Bin',I2,' moxid',F7.1, &
                                      ' delta',F7.1 )
                               write(116,14) L, madd(L),delt2
                               !write(116,*)'CASE_2_mass_not_conserve'
                               !write(116,*)'For_bin',L,'moxid',madd(L) &
                               !     ,'delta',delt2
                            ENDIF
                            errspot=.true. !just try temp (win, 7/30/07)
                         endif !significant mass add (10 kg) - then print error.
                         !<step5.3> Fix the problem of mass not conserved
                         !in case of aqueous oxidation by find the missing mass
                         !and spread them equally into the bins that the final
                         !tophat has grown to. (win, 7/24/06)
                         xtra  = 0e+0_fp
                         dummy = 0e+0_fp
                         do kk = I,K
                            !AMK(kk,CSPECIES) = AMK(kk,CSPECIES)-delt2/(K-I+1)
                            dummy = AMK(kk,CSPECIES) - &
                                    ( delt2/(K-I+1) + xtra )
                            if(dummy < 0.e+0_fp )then
                               xtra = xtra + delt2/(K-I+1)
                            else
                               AMK(kk,CSPECIES) = dummy
                               xtra = 0.e+0_fp
                            endif
                         enddo
                      endif   !error>treshold
                   endif      !moxd>0

                ENDIF  !YU.LE.X(I+1)
                GOTO 200
             ELSE    !YL > X(I+1)
                IF(I == nBins .and.(madd(L)/maddtot)> 1.5e-1_fp) THEN
11                 FORMAT( 'Tophat>Xk(31) at bin ',I3,' loosing ', &
                           E13.5,' kg = ',F5.1,'%')
                   if(MAXVAL(moxd(:)) > 0e+0_fp) then
                      print 11, L, madd(L),(madd(L)/maddtot)*1.e+2_fp
                      !write(116,11) L, madd(L),(madd(L)/maddtot)*1.e+2_fp
                      !write(117,*) madd(L)  !for accumulating mass loss
                      !PRINT *,'Tophat > Xk(31): growth over bin30,Loss%'
                      !if(moxd >0e+0_fp)print *,madd(L),(madd(L)/maddtot)*1.d2
                      !errspot = .true.
                   endif
                ENDIF
             ENDIF   !YL.LT.X(I+1)
          ENDDO !I loop
200       CONTINUE
       ENDDO    !L loop
    ENDIF

    !Signal error out to so4cond so the run can stop in aerophys and show i,j,l (win, 4/12/06)
    pdbug = errspot

    RETURN

  END SUBROUTINE TMCOND
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gasdiff
!
! !DESCRIPTION: Function GASDIFF returns the diffusion constant of a species in
!  air (m2/s). It uses the method of Fuller, Schettler, and Giddings as
!  described in Perry's Handbook for Chemical Engineers.
!  WRITTEN BY Peter Adams, May 2000
!\\
!\\
! !INTERFACE:
!
  FUNCTION GASDIFF( TEMP, PRES, MW, SV ) RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    real(fp) temp, pres  !temperature (K) and pressure (Pa) of air
    real(fp) mw          !molecular weight (g/mol) of diffusing species
    real(fp) Sv          !sum of atomic diffusion volumes of diffusing species
!
! !RETURN VALUE:
!
    real(fp) VALUE
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(fp) mwair, Svair   !same as above, but for air
    real(fp) mwf, Svf
    parameter(mwair=28.9, Svair=20.1)

    !========================================================================
    ! GASDIFF begins here!
    !========================================================================

    mwf=sqrt((mw+mwair)/(mw*mwair))
    Svf=(Sv**(1./3.)+Svair**(1./3.))**2.
    VALUE =1.0e-7*temp**1.75*mwf/pres*1.0e5/Svf

  END FUNCTION GASDIFF
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aerodens
!
! !DESCRIPTION: Function AERODENS calculates the density (kg/m3) of a sulfate-
!  nitrate-ammonium-nacl-OC-EC-dust-water mixture.  Inorganic mass (sulfate-
!  nitrate-ammonium-nacl-water) is assumed to be internally mixed.  Then the
!  density of inorg and EC, OC, and dust is combined weighted by mass.
!  WRITTEN BY Peter Adams, May 1999 in GISS GCM-II' and extened to include
!  carbonaceous aerosol in Jan, 2002.
!\\
!\\
! !INTERFACE:
!
  FUNCTION AERODENS( MSO4, MNO3, MNH4, MNACL, MECIL, MECOB, MOCIL, &
                     MOCOB, MDUST, MH2O )  RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN)  ::  MSO4, MNO3, MNH4, MNACL, MH2O
    REAL(fp),  INTENT(IN)  ::  MECIL, MECOB, MOCIL, MOCOB, MDUST
!
! !RETURN VALUE:
!
    REAL(fp)                  :: VALUE
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(fp)                  :: IDENSITY, DEC, DOC, DDUST, MTOT
    parameter(dec=2200., doc=1400., ddust=2650.)

    !=================================================================
    ! AERODENS begins here!
    !=================================================================

    IDENSITY = INODENS( MSO4, MNO3, MNH4, MNACL, MH2O )
    MTOT = MSO4+MNO3+MNH4+MNACL+MH2O+MECIL+MECOB+MOCIL+MDUST+MOCOB
    IF ( MTOT > 0.e+0_fp ) THEN
       VALUE = ( IDENSITY*(MSO4+MNO3+MNH4+MNACL+MH2O) + &
                 DEC*(MECIL+MECOB) + DOC*(MOCIL+MOCOB)+ &
                 DDUST*MDUST                            )/MTOT
    ELSE
       VALUE = 1400.
    ENDIF

  END FUNCTION AERODENS
!EOC
  !------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: inodens
!
! !DESCRIPTION: Function INODENS calculates the density (kg/m3) of a sulfate-
!  nitrate-ammonium-nacl-water mixture that is assumed to be internally mixed.
!  WRITTEN BY Peter Adams, May 1999 in GISS GCM-II'
!  Introduced to GEOS-CHEM by Win Trivitayanurak (win@cmu.edu) 8/6/07 first
!  as AERODENS, then change to INODENS on 9/3/07
!\\
!\\
! !INTERFACE:
!
  FUNCTION INODENS( MSO4_, MNO3_, MNH4_, MNACL_, MH2O_ ) &
       RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    ! mso4, mno3, mnh4, mh2o, mnacl - These are the masses of each aerosol
    ! component.  Since the density is an intensive property,
    ! these may be input in a variety of units (ug/m3, mass/cell, etc.).
    REAL(fp),  INTENT(IN)  ::  MSO4_, MNO3_, MNH4_, MNACL_, MH2O_
!
! !RETURN VALUE:
!
    REAL(fp)               :: VALUE
!
! !REMARKS:
! ----Literature cited----
!     I. N. Tang and H. R. Munkelwitz, Water activities, densities, and
!       refractive indices of aqueous sulfates and sodium nitrate droplets
!       of atmospheric importance, JGR, 99, 18,801-18,808, 1994
!     Ignatius N. Tang, Chemical and size effects of hygroscopic aerosols
!       on light scattering coefficients, JGR, 101, 19,245-19,250, 1996
!     Ignatius N. Tang, Thermodynamic and optical properties of mixed-salt
!       aerosols of atmospheric importance, JGR, 102, 1883-1893, 1997
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(fp) MSO4, MNO3, MNH4, MNACL, MH2O
    !real(fp) so4temp, no3temp, nh4temp, nacltemp, h2otemp
    real(fp) mwso4, mwno3, mwnh4, mwnacl, mwh2o            !molecular weights
    real(fp) ntot, mtot                          !total number of moles, mass
    real(fp) nso4, nno3, nnh4, nnacl, nh2o       !moles of each species
    real(fp) xso4, xno3, xnh4, xnacl, xh2o       !mole fractions
    real(fp) rso4, rno3, rnh4, rnacl, rh2o       !partial molar refractions
    real(fp) ran, rs0, rs1, rs15, rs2       !same, but for solute species
    real(fp) asr                            !ammonium/sulfate molar ratio
    real(fp) nan, ns0, ns1, ns15, ns2, nss  !moles of dry solutes (nss = sea salt)
    real(fp) xan, xs0, xs1, xs15, xs2, xss  !mass % of dry solutes - Tang (1997) eq. 10
    real(fp) dan, ds0, ds1, ds15, ds2, dss  !binary solution densities - Tang (1997) eq. 10
    real(fp) mwan, mws0, mws1, mws15, mws2  !molecular weights
    real(fp) yan, ys0, ys1, ys15, ys2, yss  !mole fractions of dry solutes
    real(fp) yh2o
    real(fp) d                              !mixture density
    real(fp) xtot

    ! In the lines above, "an" refers to ammonium nitrate, "s0" to
    ! sulfuric acid, "s1" to ammonium bisulfate, and "s2" to ammonium sulfate.
    ! "nacl" or "ss" is sea salt.
    parameter(mwso4=96.e+0_fp, &
              mwno3=62.e+0_fp, &
              mwnh4=18.e+0_fp, &
              mwh2o=18.e+0_fp, &
              mwnacl=58.45e+0_fp)
    parameter(mwan=mwnh4+mwno3,          &
              mws0=mwso4+2.e+0_fp,       &
              mws1=mwso4+1.e+0_fp+mwnh4, &
              mws2=2.e+0_fp*mwnh4+mwso4)

    !=================================================================
    ! INODENS begins here!
    !=================================================================

    ! Pass initial component masses to local variables
    mso4=mso4_
    mno3=mno3_
    mnh4=mnh4_
    mnacl=mnacl_
    mh2o=mh2o_

    !so4temp=mso4
    !no3temp=mno3
    !nh4temp=mnh4
    !h2otemp=mh2o
    !nacltemp=mnacl

    ! [Pengfei Liu, avoid equality test with floating-point real numbers
    !<step4.7> if the aerosol mass is zero - then just return the
    !typical density = 1500 kg/m3 (win, 1/4/06)
    !if (mso4 .eq. 0.e+0_fp .and. mno3 .eq.0.e+0_fp &
    !    .and. mnh4.eq.0.e+0_fp .and. mnacl .eq. 0.e+0_fp ) then
    !   VALUE = 1500.e+0_fp !kg/m3
    !   goto 10
    !endif
    if ((mso4+mno3+mnh4+mnacl) .gt. 0.e+0_fp) then
       CONTINUE
    else
       VALUE = 1500.e+0_fp !kg/m3
       RETURN
    endif
    ! Pengfei Liu, 2018/02/07]

    ! Calculate mole fractions
    mtot  = mso4+mno3+mnh4+mnacl+mh2o
    nso4  = mso4/mwso4
    nno3  = mno3/mwno3
    nnh4  = mnh4/mwnh4
    nnacl = mnacl/mwnacl
    nh2o  = mh2o/mwh2o
    ntot  = nso4+nno3+nnh4+nnacl+nh2o
    xso4  = nso4/ntot
    xno3  = nno3/ntot
    xnh4  = nnh4/ntot
    xnacl = nnacl/ntot
    xh2o  = nh2o/ntot

    ! If there are more moles of nitrate than ammonium, treat unneutralized
    ! HNO3 as H2SO4
    if (nno3 .gt. nnh4) then
       !make the switch
       nso4=nso4+(nno3-nnh4)
       nno3=nnh4
       mso4=nso4*mwso4
       mno3=nno3*mwno3

       !recalculate quantities
       mtot = mso4+mno3+mnh4+mnacl+mh2o
       nso4 = mso4/mwso4
       nno3 = mno3/mwno3
       nnh4 = mnh4/mwnh4
       nnacl = mnacl/mwnacl
       nh2o = mh2o/mwh2o
       ntot = nso4+nno3+nnh4+nnacl+nh2o
       xso4 = nso4/ntot
       xno3 = nno3/ntot
       xnh4 = nnh4/ntot
       xnacl = nnacl/ntot
       xh2o = nh2o/ntot

    endif

    ! Calculate the mixture density
    ! Assume that nitrate exists as ammonium nitrate and that other ammonium
    ! contributes to neutralizing sulfate
    nan=nno3
    if (nnh4 .gt. nno3) then
       !extra ammonium
       asr=(nnh4-nno3)/nso4
    else
       !less ammonium than nitrate - all sulfate is sulfuric acid
       asr=0.e+0_fp
    endif
    if (asr .ge. 2.e+0_fp) asr=2.e+0_fp
    if (asr .ge. 1.e+0_fp) then
       !assume NH4HSO4 and (NH4)2(SO4) mixture
       !NH4HSO4
       ns1=nso4*(2.e+0_fp-asr)
       !(NH4)2SO4
       ns2=nso4*(asr-1.e+0_fp)
       ns0=0.e+0_fp
    else
       !assume H2SO4 and NH4HSO4 mixture
       !NH4HSO4
       ns1=nso4*asr
       !H2SO4
       ns0=nso4*(1.e+0_fp-asr)
       ns2=0.e+0_fp
    endif

    !Calculate weight percent of solutes
    xan=nan*mwan/mtot*100.e+0_fp
    xs0=ns0*mws0/mtot*100.e+0_fp
    xs1=ns1*mws1/mtot*100.e+0_fp
    xs2=ns2*mws2/mtot*100.e+0_fp
    xnacl=nnacl*mwnacl/mtot*100.e+0_fp
    xtot=xan+xs0+xs1+xs2+xnacl

    ! [Pengfei Liu, fix the polynomial issue
    !Calculate binary mixture densities (Tang, eqn 9)
    !dan=0.9971e+0_fp +4.05e-3_fp*xtot +9.0e-6_fp*xtot**2.e+0_fp
    !ds0=0.9971e+0_fp +7.367e-3_fp*xtot -4.934d-5*xtot**2.e+0_fp &
    !     +1.754e-6_fp*xtot**3.e+0_fp - 1.104d-8*xtot**4.e+0_fp
    !ds1=0.9971e+0_fp +5.87e-3_fp*xtot -1.89e-6_fp*xtot**2.e+0_fp &
    !     +1.763e-7_fp*xtot**3.e+0_fp
    !ds2=0.9971e+0_fp +5.92e-3_fp*xtot -5.036e-6_fp*xtot**2.e+0_fp &
    !     +1.024d-8*xtot**3.e+0_fp
    !dss=0.9971e+0_fp +7.41e-3_fp*xtot -3.741d-5*xtot**2.e+0_fp &
    !     +2.252e-6_fp*xtot**3.e+0_fp   -2.06d-8*xtot**4.e+0_fp
    dan=0.9971e+0_fp + xtot * (4.05e-3_fp + 9.0e-6_fp * xtot)
    ds0=0.9971e+0_fp &
        +xtot*(7.367e-3_fp &
        +xtot*(-4.934d-5 &
        +xtot*(1.754e-6_fp &
        +xtot*(-1.104d-8  ))))
    ds1=0.9971e+0_fp &
        +xtot*(5.87e-3_fp &
        +xtot*(-1.89e-6_fp &
        +xtot*(1.763e-7_fp )))
    ds2=0.9971e+0_fp &
        +xtot*(5.92e-3_fp &
        +xtot*(-5.036e-6_fp &
        +xtot*(1.024d-8    )))
    dss=0.9971e+0_fp &
        +xtot*(7.41e-3_fp &
        +xtot*(-3.741d-5 &
        +xtot*(2.252e-6_fp &
        +xtot*(-2.06d-8     ))))
    ! Pengfei Liu, 2018/02/07]

    !Convert x's (weight percent of solutes) to fraction of dry solute (scale to 1)
    xtot=xan+xs0+xs1+xs2+xnacl
    xan=xan/xtot
    xs0=xs0/xtot
    xs1=xs1/xtot
    xs2=xs2/xtot
    xnacl=xnacl/xtot

    !Calculate mixture density
    d=1.e+0_fp/(xan/dan+xs0/ds0+xs1/ds1+xs2/ds2+xnacl/dss)  !Tang, eq. 10

    if ((d .gt. 2.e+0_fp) .or. (d .lt. 0.997e+0_fp)) then
       write(*,*) 'ERROR in aerodens'
       write(*,*) mso4,mno3,mnh4,mnacl,mh2o
       print *, 'xtot',xtot
       print *, 'xs1',xs1, 'ns1',ns1,'mtot',mtot,'asr',asr
       write(*,*) 'density(g/cm3)',d
       STOP
    endif

    ! Restore masses passed
    !mso4=so4temp
    !mno3=no3temp
    !mnh4=nh4temp
    !mnacl=nacltemp
    !mh2o=h2otemp

    ! Return the density
    VALUE = 1000.e+0_fp*d    !Convert g/cm3 to kg/m3

    !<step4.7> negative value check (win, 1/4/06)
    if ( VALUE < 0e+0_fp ) then
       print *, 'ERROR :: aerodens - negative', VALUE
       STOP
    endif

10  CONTINUE

  END FUNCTION INODENS
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dmdt_int
!
! !DESCRIPTION: Function DMDT_INT apply the analytic solution to the droplet
!  growth equation in mass space for a given scale length which mimics the
!  inclusion of gas kinetic effects. (win, 7/23/07)
!  Originally written by Peter Adams
!  Modified for GEOS-CHEM by Win Trivitayanurak (win@cmu.edu)
!\\
!\\
! !INTERFACE:
!
  FUNCTION DMDT_INT ( M0, TAU, WR ) RESULT( VALUE )
!
! !INPUT PARAMETERS:
!
    ! M0  initial mass
    ! L0  length scale
    ! Tau forcing from vapor field
    REAL(fp),   INTENT(IN)  ::  M0,  TAU,  WR
!
! !RETURN VALUE:
!
    REAL(fp)                :: VALUE
!
! !REMARKS:
!  Original note from Peter Adams:
!  I have changed the length scale.  Non-continuum effects are
!  assumed to be taken into account in choice of tau (in so4cond subroutine).
!  .
!  I have also added another argument to the function call, WR.  This
!  is the ratio of wet mass to dry mass of the particle.  I use this
!  information to calculate the amount of growth of the wet particle,
!  but then return the resulting dry mass.  This is the appropriate
!  way to implement the condensation algorithm in a moving sectional
!  framework.
!  .
!  Reference: Stevens et al. 1996, Elements of the Microphysical Structure
!           of Numerically Simulated Nonprecipitating Stratocumulus,
!           J. Atmos. Sci., 53(7),980-1006.
! This calculates a solution for m(t+dt) using eqn.(A3) from the reference
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)                ::  X,  L0,  C,  ZERO,  MH2O
    PARAMETER (C=2.e+0_fp/3.e+0_fp,L0=0.0e+0_fp,ZERO=0.0e+0_fp)

    !=================================================================
    ! DMDT_INT begins here!
    !=================================================================

    MH2O = ( WR - 1.e+0_fp ) * M0
    X = ( ( M0 + MH2O ) ** C + L0 )
    X = MAX( ZERO, SQRT(MAX(ZERO,C*TAU+X))-L0 )

    !<step5.3> Do aqueous oxidation dry - so no need to select process (win, 7/14/06)
    !<step5.3> For so4cond condensation, use constant water amount.
    ! For aqueous oxidation, use constant wet ratio. (win, 7/13/06)
    !prior to 10/2/08
    !VALUE = X * X * X - MH2O
    !!DMDT_INT = X*X*X/WR    !<step5.2> change calculation to keep WR constant after condensation/evap (win, 5/14/06)

    !<step6.3> bring back the previously reverted back (win, 10/2/08)
    VALUE = X*X*X/WR
    !pja Perform some numerical checks on dmdt_int
    IF ((TAU > 0.0) .and. (VALUE < M0)) VALUE = M0
    IF ((TAU < 0.0) .and. (VALUE > M0)) VALUE = M0

  END FUNCTION DMDT_INT
!EOC
END MODULE Lagrange_singlebox_Mod
