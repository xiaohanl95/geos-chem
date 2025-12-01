#undef AER_CONSTANT_WP
#undef AER_PVH2SO4
!#define AER_CONSTANT_WP
!#define AER_PVH2SO4
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
! !MODULE: sect_aer_mod.F90
!
! !DESCRIPTION: Module Sect\_Aer_\_Mod contines arrays and routines for the
!  sectional representation of stratospheric aerosols.
!\\
!\\
! !INTERFACE: 
!
module sect_aer_mod
!
! !USES:
!
  ! For GEOS-Chem Precision (fp)
  USE Precision_Mod,        ONLY : fp, &
                                   dp=>f8

  ! Data which will likely go to State_Chm
  USE Sect_Aer_Data_Mod,    ONLY : n_aer_bin, &
                                   aer_mass,  & ! Mass per particle in each bin
                                   aer_molec, & ! Molec per particle in each bin
                                   NAStep       ! Number of aerosol substeps

  ! Arrays which should be put into State_Chm
  USE Sect_Aer_Data_Mod,    ONLY : bvp,       & ! Eq. vapor pressure of H2SO4
                                   st,        & ! Surface tension
                                   airvdold,  & ! ???
                                   ck           ! Coagulation kernel

  ! Constants
  USE PhysConstants,        ONLY : mw=>h2omw, &
                                   av=>avo,   &
                                   pi

  USE Sect_Aer_Data_Mod,    ONLY : Ma=>aer_mwh2so4, &
                                   den_h2so4,       &
                                   Rgas=>aer_Rgas,  &
                                   akb=>aer_akb,    &
                                   aer_dry_rad,     &
                                   aer_Vrat

#if defined( MDL_BOX )
  USE sect_aux_mod, only : error_stop, debug_msg
#else
  USE error_mod,    only : error_stop, debug_msg
  Use errcode_mod
#endif

#if defined( USE_TIMERS )
  USE GEOS_Timers_Mod
#endif

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Do_Sect_Aer
  PUBLIC  :: Init_Sect_Aer
  PUBLIC  :: Cleanup_Sect_Aer
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: AER_wpden
  PRIVATE :: Cal_Coag_Kernel
  PRIVATE :: Cal_Coag_Kernel_Implicit
  PRIVATE :: AER_grow
  PRIVATE :: AER_nucleation
  PRIVATE :: AER_coagulation
  PRIVATE :: AER_coagulation_Implicit
  PRIVATE :: AER_mass_check
  PRIVATE :: AER_sedimentation
  PRIVATE :: localminmax
  PRIVATE :: Vehkamaeki
  PRIVATE :: rnucnum_f
  PRIVATE :: rnucmas_f
  PRIVATE :: so4fracnuc_f
  PRIVATE :: density
  PRIVATE :: surftension
  PRIVATE :: pvolume_h2o
  PRIVATE :: pvolume_h2so4
  PRIVATE :: sath2o
  PRIVATE :: sath2so4
  PRIVATE :: solh2o
  PRIVATE :: solh2so4

  ! Temporary
  PUBLIC  :: ID_Bins
!
! !PRIVATE TYPES:
  INTEGER, ALLOCATABLE        :: ID_Bins(:,:)        ! IDs of the aerosols
  REAL(dp),ALLOCATABLE        :: aero_grow_d(:,:,:)  ! Aerosol growth rate
  REAL(dp),ALLOCATABLE        :: aero_nrate_d(:,:,:) ! Aerosol nucleation rate
  INTEGER                     :: id_H2O              ! Index of the H2O tracer
  INTEGER                     :: id_H2SO4            ! Index of gas-phase H2SO4
  
! !REVISION HISTORY:
!  28 Nov 2018 - S.D.Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
!  LOGICAL               :: Module_Logical_Example
!  REAL(f4), ALLOCATABLE         :: JvSumMon   (:,:,:,:)  ! For monthly avg of J-values
!
! !DEFINED PARAMTERS:
!
   Real(dp), PARAMETER :: so4gmin=1.0e+04_dp
   Real(dp), PARAMETER :: so4gmax=1.0e+11_dp
   INTEGER,  PARAMETER :: n_aer_type=1   ! Only SUL for now

CONTAINS
!EOC
#if !defined( MDL_BOX )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_sect_aer
!
! !DESCRIPTION: Subroutine Do\_Sect\_Aer is the driver subroutine for
!  sectional stratospheric aerosols.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Sect_Aer( Input_Opt,  State_Met,             &
                          State_Chm, State_Diag, State_Grid, &
                          Run_Micro, RC                      )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,        ONLY : OptInput
    USE PhysConstants,        ONLY : AVO
    USE PRESSURE_MOD        
    USE Species_Mod,          ONLY : Species
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Chm_Mod,        ONLY : Ind_
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Met_Mod,        ONLY : MetState
    USE TIME_MOD,             ONLY : GET_TS_CHEM
    Use UnitConv_Mod,         ONLY : Convert_Spc_Units
!
! !INPUT VARIABLES:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
    LOGICAL,        INTENT(IN)    :: Run_Micro  ! Run full microphysics calc?
    TYPE(GrdState), INTENT(IN)    :: State_Grid ! Grid description object
!
! !INPUT/OUTPUT VARIABLES:
!
    TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag ! Diagnostics State object
!
! !OUTPUT VARIABLES:
!
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure
! 
! !REVISION HISTORY:
!  28 Nov 2018 - S.D.Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    Integer                :: n_lagrangian

    ! Loop indices
    Integer                :: I,J,K,L

    ! Local properties
    Real(fp)               :: T_K, P_hPa, air_dens
    Real(fp)               :: H2O_vv, H2SO4_vv, H2O_vv_old, H2SO4_vv_old
    Real(fp)               :: Sul_vv(n_aer_bin), Sul_vv_old(n_aer_bin)
    Real(fp), Allocatable  :: Sul_col(:,:)

    ! For implicit coagulation
    Real(fp)               :: fijk_Box(n_aer_bin,n_aer_bin,n_aer_bin)

    ! Substep length
    Real(fp)               :: dt_substep, dt_step

    ! Diagnostics
    Real(fp)               :: box_grow_d, box_nrate_d, box_wp_d
    Real(fP)               :: s_start, s_end
    !Real(fp)               :: aero_vt_d(NZ)

    ! Strings
    Character(Len=255)     :: Src_Unit

    ! Pointers to properties for a single box
    Real(fp), Pointer      :: aWP_Box(:)
    Real(fp), Pointer      :: aDen_Box(:)
    Real(fp), Pointer      :: ST_Box(:)
    Real(fp), Pointer      :: BVP_Box(:)
    Real(fp), Pointer      :: CK_Box(:,:)
    Real(fp), Pointer      :: Spc_vv(:)

    ! Pointers to larger sets
    Real(fp), Pointer      :: Spc_col(:,:)
    Real(fp), Pointer      :: T_Col(:)
    Real(fp), Pointer      :: P_Col(:)
    Real(fp), Pointer      :: dz_Col(:)
    Real(fp), Pointer      :: air_dens_Col(:)
    Real(fp), Pointer      :: aWP_Col(:,:)
    Real(fp), Pointer      :: aDen_Col(:,:)

    ! Temporary
    Logical                :: LAER_nuc
    Logical                :: LAER_grow
    Logical                :: LAER_coag
    Logical                :: LAER_coag_imp
    Logical                :: LAER_sed

    ! For diagnostics
    Logical                :: Archive_RWet
    Logical                :: Archive_VGrav

    Real(fp)               :: rWet( n_aer_bin)

    Integer                :: NX, NY, NZ

    ! Internal return codes
    Integer                :: RCOMP

    !=======================================================================
    ! Do_Sect_Aer begins here!
    !=======================================================================

#if defined( USE_TIMERS )
    CALL GEOS_Timer_Start( "==> Sectional aero.", RC )
#endif

    ! Assume success
    RC = 0

    ! Grid dimensions
    NX = State_Grid%NX
    NY = State_Grid%NY
    NZ = State_Grid%NZ

    ! Debug options
    LAER_nuc  = Input_Opt%LAER_nuc
    LAER_grow = Input_Opt%LAER_grow
    LAER_coag = Input_Opt%LAER_coag
    LAER_sed  = Input_Opt%LAER_sed

    ! Implicit coagulation?
    LAER_coag_imp = .False.

    If (Input_Opt%amIRoot) Then
       write(*,*) 'Running Do_Sect_Aer'
    End If

    ! Get the length of the microphysics substeps
    dt_step    = Real(GET_TS_CHEM())
    dt_substep = dt_step/Real(NAStep)

    ! Convert Eulerian data to VMR
    Call Convert_Spc_Units(Input_Opt,  State_Chm, &
                           State_Grid, State_Met, 'v/v dry', &
                           RC, Src_Unit)

    ! Are diagnostics allocated?
    ! TODO: Reenable once diagnostics are coded
    !Archive_RWet  = Associated(State_Diag%Strat_RWet)
    !Archive_VGrav = Associated(State_Diag%Strat_VGrav)
    Archive_RWet  = .False.
    Archive_VGrav = .False.

    ! First, loop over the standard Eulerian grid cells
    !$OMP PARALLEL  DO                                               &
    !$OMP DEFAULT(  SHARED                                         ) &
    !$OMP PRIVATE(  I,          J,        K,          L            ) &
    !$OMP PRIVATE(  aWP_Box,    aDen_Box, Sul_vv,     Spc_vv       ) &
    !$OMP PRIVATE(  H2O_vv,     H2SO4_vv, H2O_vv_old, H2SO4_vv_old ) &
    !$OMP PRIVATE(  Sul_vv_old, T_K,      P_hPa,      ST_Box       ) &
    !$OMP PRIVATE(  air_dens,   CK_Box,   BVP_Box,    rWet         ) &
    !$OMP PRIVATE(  box_grow_d, box_wp_d, box_nrate_d,fijk_Box     ) &
    !$OMP PRIVATE(  s_start,    s_end,    RCOMP                    ) &
    !$OMP SCHEDULE( DYNAMIC, 1                                     )
    Do I=1,NX
    Do J=1,NY
    Do L=1,NZ
       ! Point to properties of the target box
       aWP_Box  => State_Chm%Strat_WPcg(I,J,L,:)
       aDen_Box => State_Chm%Strat_Den(I,J,L,:)
       Spc_vv   => State_Chm%Species(I,J,L,:)
       
       ! Get properties of the box
       H2O_vv   = Spc_vv(id_H2O)
       H2SO4_vv = Spc_vv(id_H2SO4)
       Do k=1,n_aer_bin
          Sul_vv(k) = Spc_vv(ID_Bins(k,1))
       End Do

       ! Store the starting values
       H2O_vv_old    = H2O_vv
       H2SO4_vv_old  = H2SO4_vv
       Sul_vv_old(:) = Sul_vv(:)
       If (State_Met%InTroposphere(I,J,L)) Then
          ! Apply some feasible stand-in values
          aWP_Box(:) = 60.0e+0_fp
          aDen_Box(:) = 1.5e+0_fp
          ! Dump all sulfate into the "H2SO4" container
          ! This represents tropospheric aerosol
          Spc_vv(id_H2SO4) = H2SO4_vv + Sum(Sul_vv)
          Do k=1,n_aer_bin
             Spc_vv(ID_Bins(k,1)) = 0.0e+0_fp
          End Do
          aero_grow_d(I,J,L)  = 0.0e+0_fp
          aero_nrate_d(I,J,L) = 0.0e+0_fp
          rWet(:)             = aer_dry_rad(:)
       Else
          ! Get properties of location
          T_K      = State_Met%T(I,J,L)
          P_hPa    = State_Met%PMID(I,J,L)

          ! Get other aerosol properties
          ! Surface tension (ST) and equilibrium vapor pressure (BVP)
          ST_Box  => st(I,J,L,:)
          BVP_Box => bvp(I,J,L,:)

          ! Calculate updated weight % H2SO4 and density
          Call AER_wpden(aWP_Box, aDen_Box, ST_Box, BVP_Box,&
                         T_K,     P_hPa,    H2O_vv, rWet)

          ! If we aren't running a microphysics step, stop here
          If (Run_Micro) Then
             ! Get more met properties
             air_dens = State_Met%AIRNUMDEN(I,J,L)

             ! Calculate coagulation kernel
             CK_Box => ck(I,J,L,:,:)
             If (LAER_coag) &
             Call Cal_Coag_Kernel(CK_Box,aWP_Box,aDen_Box,T_K,P_hPa)          

             ! Calculate volume fractions
             fijk_Box = 0.0e+0_fp
             If (LAER_coag.and.LAER_Coag_Imp) &
             Call Cal_Coag_Kernel_Implicit(fijk_Box,aWP_Box,aDen_Box)

             ! Run microphysical processes
             box_grow_d  = 0.0e+0_fp
             box_nrate_d = 0.0e+0_fp

             Do K=1,NAStep
                If (LAER_grow) &
                Call AER_grow(T_K,      air_dens,   dt_substep, aWP_Box, &
                              aDen_box, ST_Box,     BVP_Box,    H2SO4_vv, &
                              Sul_vv,   box_grow_d)
                If (LAER_nuc) &
                Call AER_nucleation(T_K,     P_hPa,       air_dens, dt_substep, &
                                    aWP_Box, aDen_Box,    H2O_vv,   H2SO4_vv,   &
                                    Sul_vv,  box_nrate_d, box_wp_d              )
                If (LAER_coag.and.(.not.LAER_Coag_Imp)) &
                Call AER_coagulation(Sul_vv,CK_box,air_dens,dt_substep)
                If (LAER_coag.and.LAER_Coag_Imp) &
                Call AER_coagulation_Implicit(Sul_vv,fijk_box,CK_box,aWP_Box,aDen_Box,air_dens,dt_substep)
             End Do
             aero_grow_d(I,J,L)  = box_grow_d/real(NAStep)
             aero_nrate_d(I,J,L) = box_nrate_d/real(NAStep)
         
             ! Make sure that the mass hasn't drifted too far 
             Call AER_mass_check(H2SO4_vv, Sul_vv,  H2SO4_vv_old, Sul_vv_old, &
                                 air_dens, s_start, s_end,        RCOMP       )
             If (RCOMP.ne.0) RC = RCOMP
             !Call AER_Washout(?,?,?,?)
          End If

          ! Update properties that WEREN'T held in pointers
          Spc_vv(id_H2O)   = H2O_vv
          Spc_vv(id_H2SO4) = H2SO4_vv
          Do k=1,n_aer_bin
             Spc_vv(ID_Bins(k,1)) = Sul_vv(k)
          End Do

       End If

       ! Update diagnostics
       ! TODO: Reenable
       If (Archive_RWet) Then
       !   Do k=1,n_aer_bin
       !      State_Diag%Strat_RWet(I,J,L,K) = rWet(K)
       !   End Do
       End If

       ! Free up pointers
       aWP_Box  => Null()
       aDen_Box => Null()
       BVP_Box  => Null()
       ST_Box   => Null()
       CK_Box   => Null()
       Spc_vv   => Null()
    End Do ! L
    End Do ! J
    End Do ! I
    !$OMP END PARALLEL DO

    ! Everything still going OK?
    If (RC.ne.0) Then
       Call Error_Stop('Failure in sectional strat aerosol', &
                       'sect_aer_mod.F90')
    End If

    ! Do sedimentation separately - has to be performed in columns
    If (Run_Micro.and.LAER_sed) Then
       !$OMP PARALLEL                                                              &
       !$OMP DEFAULT(  SHARED                                                    ) &
       !$OMP PRIVATE(  I,      J,     K,            aWP_Col, aDen_Col, Spc_Col   ) &
       !$OMP PRIVATE(  T_Col,  P_Col, air_dens_Col, dz_col,  Sul_Col,  RCOMP     )
       Allocate(Sul_Col(NZ,n_aer_bin),Stat=RCOMP)
       If (RCOMP.ne.0) RC = RCOMP
       !$OMP DO 
       Do I=1,NX
       Do J=1,NY
          ! Don't bother if we had an allocation issue
          If (RCOMP.ne.0) Cycle
          ! Point to properties of the target column (chemical)
          aWP_Col  => State_Chm%Strat_WPcg(I,J,:,:)
          aDen_Col => State_Chm%Strat_Den(I,J,:,:)
          Spc_Col  => State_Chm%Species(I,J,:,:)
          Do k=1,n_aer_bin
             Sul_col(:,k) = Spc_col(:,ID_Bins(k,1))
          End Do

          ! Now met properties
          T_col        => State_Met%T(I,J,:)
          P_col        => State_Met%PMID(I,J,:)
          air_dens_col => State_Met%AIRNUMDEN(I,J,:)
          dz_col       => State_Met%BxHeight(I,J,:)

          ! Run sedimentation in the column
          Call AER_sedimentation(Sul_Col, T_Col, P_Col, aWP_Col, aDen_Col,   & 
                                 dt_step, air_dens_col, dz_col,  NZ)

          ! Store updated VMRs
          Do k=1,n_aer_bin
             Spc_col(:,ID_Bins(k,1)) = Sul_col(:,k)
          End Do

          ! Diagnostics
          ! TODO: Reenable
          If (Archive_VGrav) Then
             !State_Diag%Strat_VGrav(I,J,:) = aero_vt_d
          End If

          ! Free up pointers
          aWP_Col      => Null()
          aDen_Col     => Null()
          Spc_Col      => Null()
          T_Col        => Null()
          P_Col        => Null()
          air_dens_col => Null()
          dz_col       => Null()
       End Do ! J
       End Do ! I
       !$OMP END DO
       If (Allocated(Sul_Col)) Deallocate(Sul_Col)
       !$OMP END PARALLEL
    End If ! LAER_sed and Run_Micro

    ! Everything still going OK?
    If (RC.ne.0) Then
       Call Error_Stop('Failure in sectional strat aerosol settling', &
                       'sect_aer_mod.F90')
    End If

    ! Convert back to whatever we had before
    Call Convert_Spc_Units(Input_Opt,  State_Chm, &
                           State_Grid, State_Met, Src_Unit,  &
                           RC)

    ! Now loop over the Lagrangian grid cells
    !If (n_lagrangian > 0) Then
    !   Do I = 1, n_lagrangian
    !      ! Retrieve box properties
    !      Box => Box_List(I)
    !      T = Box%T(I)
    !      P = Box%P(I)
    !      H2O_vv  = Box%Species(id_h2o)
    !      aWP_Box => Box%aWP
    !   End Do
    !End If

#if defined( USE_TIMERS )
    CALL GEOS_Timer_End(   "==> Sectional aero.", RC )
#endif

    ! Success
    RC = 0

  end subroutine do_sect_aer
!EOC
#endif
#if defined( MDL_BOX )
  subroutine do_sect_aer(n_boxes,aWP_Arr,aDen_Arr,&
                         vvSO4_Arr,Sfc_Ten_Arr,vvH2O_Vec,&
                         vvH2SO4_Vec,rWet_Arr,T_K_Vec,p_hPa_Vec,&
                         ndens_Vec,ts_sec,ts_coag,&
                         LAER_Nuc,LAER_Grow,LAER_Coag,&
                         LAER_Coag_Imp,RC)
    
    integer, intent(in)     :: n_boxes
    real(fp), intent(inout) :: aWP_Arr     (n_boxes,n_aer_bin)
    real(fp), intent(inout) :: aDen_Arr    (n_boxes,n_aer_bin)
    real(fp), intent(inout) :: vvSO4_Arr   (n_boxes,n_aer_bin)
    real(fp), intent(inout) :: rWet_Arr    (n_boxes,n_aer_bin)
    real(fp), intent(inout) :: Sfc_Ten_Arr (n_boxes,n_aer_bin)
    real(fp), intent(inout) :: vvH2O_Vec   (n_boxes)
    real(fp), intent(inout) :: vvH2SO4_Vec (n_boxes)
    real(fp), intent(in)    :: T_K_Vec     (n_boxes)
    real(fp), intent(in)    :: p_hPa_Vec   (n_boxes)
    real(fp), intent(in)    :: ndens_Vec   (n_boxes)
    integer, intent(in)     :: ts_sec
    integer, intent(in)     :: ts_coag
    logical, intent(in)     :: LAer_Coag_Imp
    integer, intent(out)    :: RC

    ! Atmospheric conditions
    real(fp)                :: T_K, p_hPa, air_dens
    
    ! Aerosol properties taken from input
    real(fp)                :: aWP_Box      (n_aer_bin)
    real(fp)                :: aDen_Box     (n_aer_bin)
    real(fp)                :: vvSO4_Box    (n_aer_bin)
    real(fp)                :: ST_Box       (n_aer_bin)
    real(fp)                :: vvSO4_Box_0  (n_aer_bin)
    
    real(fp)                :: vvH2O
    real(fp)                :: vvH2SO4
    real(fp)                :: vvH2SO4_0

    ! Aerosol properties calculated during operation
    real(fp)                :: BVP_Box      (n_aer_bin)
    real(fp)                :: rWet         (n_aer_bin)

    ! Coagulation kernel
    real(fp), pointer       :: CK_Box(:,:)
    real(fp)                :: fijk_Box(n_aer_bin,n_aer_bin,n_aer_bin)
    
    ! Loop indices 
    integer                 :: I,J,K,L,N
    integer                 :: I_Bin, I_Box

    ! For time step control
    Integer                 :: NAStep
    Real(fp)                ::dt_substep

    ! Diagnostics
    Real(fp)                :: box_grow_d, box_nrate_d, box_wp_d
    Real(fp)                :: S_Start, S_End

    ! Logicals
    logical, intent(in)     :: LAer_Nuc
    logical, intent(in)     :: LAer_Grow
    logical, intent(in)     :: LAer_Coag
    logical                 :: Run_Micro

    ! Fake locations
    I = 1
    J = 1
    L = 1

    Run_Micro = (LAer_Nuc.or.LAer_Grow.or.LAer_Coag)

    ! How many substeps?
    NAStep = NINT(dble(ts_sec)/dble(ts_coag))
    dt_substep = dble(ts_sec)/dble(NAStep)

!$OMP PARALLEL  DO                                               &
!$OMP DEFAULT(  SHARED                                         ) &
!$OMP PRIVATE(  I_Box, vvH2O, vvH2SO4, ST_Box, aWP_Box         ) &
!$OMP PRIVATE(  vvSO4_Box_0, vvH2SO4_0, T_K, p_hPa             ) &
!$OMP PRIVATE(  air_dens,   CK_Box,   BVP_Box,    rWet         ) &
!$OMP PRIVATE(  box_grow_d, box_wp_d, box_nrate_d,fijk_Box     ) &
!$OMP PRIVATE(  s_start,    s_end,    RC, K, N                 ) &
!$OMP SCHEDULE( DYNAMIC, 1                                     )
    Do I_Box=1,N_Boxes
       ! Copy properties into temporary variables
       vvH2O        = vvH2O_Vec  (I_Box)
       vvH2SO4      = vvH2SO4_Vec(I_Box)
       ST_Box(:)    = Sfc_Ten_Arr(I_Box,:)
       aWP_Box(:)   = aWP_Arr    (I_Box,:)
       aDen_Box(:)  = aDen_Arr   (I_Box,:)
       vvSO4_Box(:) = vvSO4_Arr  (I_Box,:)

       ! Store for later mass check
       vvSO4_Box_0(:) = vvSO4_Box(:)
       vvH2SO4_0      = vvH2SO4

       ! Get properties of location
       T_K      = T_K_Vec(I_Box)
       P_hPa    = p_hPa_Vec(I_Box)

       ! For safety's sake (these will be recalculated)
       BVP_Box(:) = 0.0e+0_fp
       rWet(:)    = 0.0e+0_fp

       ! Calculate updated weight % H2SO4 and density
       Call AER_wpden(aWP_Box, aDen_Box, ST_Box, BVP_Box,&
                      T_K,     P_hPa,    vvH2O,  rWet)

       ! If we aren't running a microphysics step, stop here
       If (Run_Micro) Then
          ! Get more met properties
          air_dens = ndens_vec(I_Box)

          ! Calculate coagulation kernel
          CK_Box => ck(I_Box,1,1,:,:)
          If (LAER_coag) &
          Call Cal_Coag_Kernel(CK_Box,aWP_Box,aDen_Box,T_K,P_hPa)          

          ! If using implicit coagulation..
          fijk_Box = 0.0e+0_fp
          If (LAER_coag.and.LAER_Coag_Imp) &
          Call Cal_Coag_Kernel_Implicit(fijk_Box,aWP_Box,aDen_Box)

          ! Run microphysical processes
          box_grow_d  = 0.0e+0_fp
          box_nrate_d = 0.0e+0_fp

          Do K=1,NAStep
             If (LAER_grow) &
             Call AER_grow(T_K,       air_dens,   dt_substep, aWP_Box, &
                           aDen_box,  ST_Box,     BVP_Box,    vvH2So4, &
                           vvSO4_Box, box_grow_d)
             If (LAER_nuc) &
             Call AER_nucleation(T_K,      P_hPa,       air_dens, dt_substep, &
                                 aWP_Box,  aDen_Box,    vvH2O,    vvH2SO4,    &
                                 vvSO4_Box,box_nrate_d, box_wp_d              )
             If (LAER_coag.and.(.not.LAER_Coag_Imp)) &
             Call AER_coagulation(vvSO4_Box,CK_box,air_dens,dt_substep)
             If (LAER_coag.and.LAER_Coag_Imp) &
             Call AER_Coagulation_Implicit(vvSO4_Box,fijk_Box,CK_box,aWP_Box,aDen_Box,air_dens,dt_substep)
          End Do
          !aero_grow_d(I,J,L)  = box_grow_d/real(NAStep)
          !aero_nrate_d(I,J,L) = box_nrate_d/real(NAStep)
          CK_Box => Null()
       
          ! Make sure that the mass hasn't drifted too far 
          Call AER_mass_check(vvH2SO4, vvSO4_Box,  vvH2SO4_0, vvSO4_Box_0, &
                              air_dens, s_start, s_end,        RC          )
          If (RC.ne.0) Then
             Call Error_Stop('Failure in strat aerosol mass balance', &
                             'sect_aer_mod.F90')
          End If
          !Call AER_Washout(?,?,?,?)
       End If

       ! Copy results back out
       vvH2O_Vec  (I_Box) = vvH2O
       vvH2SO4_Vec(I_Box) = vvH2SO4
       Sfc_Ten_Arr(I_Box,:) = ST_Box(:)
       aWP_Arr    (I_Box,:) = aWP_Box(:)
       aDen_Arr   (I_Box,:) = aDen_Box(:)
       vvSO4_Arr  (I_Box,:) = vvSO4_Box(:)
       rWet_Arr   (I_Box,:) = rWet(:)
    End Do
!$OMP END PARALLEL DO

    RC = 0
  end subroutine do_sect_aer
#endif
#if !defined( MDL_BOX )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_sect_aer
!
! !DESCRIPTION: Subroutine Init\_Sect\_Aer initializes module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Sect_Aer( Input_Opt,  State_Chm, State_Grid, RC )
!
! !USES:
!
    USE sect_aer_data_mod,    Only : aer_allocate_ini
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Chm_Mod,        ONLY : Ind_
    USE Species_Mod,          ONLY : Species
!
! !INPUT VARIABLES:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid ! Grid description object
!
! !INPUT/OUTPUT VARIABLES:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !OUTPUT VARIABLES:
!
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure
! 
! !REVISION HISTORY:
!  28 Nov 2018 - S.D.Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
!    LOGICAL                :: prtDebug,  IsLocNoon
    TYPE(Species), POINTER  :: ThisSpc
    Integer                 :: AS, K, N, N_Bins
    Character(Len=12)       :: SpcName
    Character(Len=255)      :: Err_Msg

    Integer                 :: NX, NY, NZ

    !=======================================================================
    ! Init_Sect_Aer begins here!
    !=======================================================================

    ! Assume success
    RC = 0

    ! Get grid information
    NX = State_Grid%NX
    NY = State_Grid%NY
    NZ = State_Grid%NZ

    ! Figure out how many bins we have
    ! For now we only have one aerosol type (!)
    ThisSpc => Null()
    n_bins = 0
    Do k=1,State_Chm%nAeroSpc
       n = State_Chm%Map_Aero(k)
       ThisSpc => State_Chm%SpcData(N)%Info
       If (ThisSpc%MP_BinNumber > n_bins ) Then
          n_bins = ThisSpc%MP_BinNumber
       End If
       ThisSpc => Null()
    End Do
    ThisSpc => Null()

    ! Save the count
    n_aer_bin = n_bins

    ! Run allocation routine
    Call AER_allocate_ini(Input_Opt,State_Chm,State_Grid,RC)

    ! If showing debug information, print out nominal radius information
    If (Input_Opt%LPRT.and.Input_Opt%amIRoot) Then
       write(*,'(a,I3,a)') 'Initializing sectional aerosols with ',N_Bins, ' bins'
       write(*,'(a)') 'Radius information:'
       Do k=1,n_aer_bin
          Write(*,'(a,I3,a,F12.6,a)') ' --> Bin ', k, ': ', aer_dry_rad(k), ' um'
       End Do
    End If

    ! Allocate module variables
    Allocate(ID_Bins(n_aer_bin,n_aer_type),STAT=AS)
    If (AS.ne.0) Then
       Call Debug_Msg('Failed to allocate ID_Bins')
       RC = -1
       Return
    End If
    ID_Bins = 0

    Allocate(aero_grow_d(NX,NY,NZ),STAT=AS)
    If (AS.ne.0) Then
       Call Debug_Msg('Failed to allocate aero_grow_d')
       RC = -1
       Return
    End If
    aero_grow_d = 0.0e+0_dp

    Allocate(aero_nrate_d(NX,NY,NZ),STAT=AS)
    If (AS.ne.0) Then
       Call Debug_Msg('Failed to allocate aero_nrate_d')
       RC = -1
       Return
    End If
    aero_nrate_d = 0.0e+0_dp

    ! Set all tracer IDs
    Do k=1,n_bins
       Write(SpcName,'(a,a,I0.3)') 'AERSCT','SUL',k
       N = Ind_(Trim(SpcName))
       ! Throw error if a bin is missing
       If (N.le.0) Then
          Write(err_msg,'(a,I3,a,I3)') 'Missing tracer for bin ',k,' of ', n_bins
          Call Debug_Msg(Trim(err_msg))
          RC = -1
       End If
       ID_Bins(k,1) = N
    End Do

    ! Check the gas phase species too
    id_H2SO4 = Ind_('SO4')
    If (id_H2SO4.le.0) Then
       Write(err_msg,'(a)') 'ID not found for SO4'
       Call Debug_Msg(Trim(err_msg))
       RC = -1
    End If
    id_H2O   = Ind_('H2O')
    If (id_H2O  .le.0) Then
       Write(err_msg,'(a)') 'ID not found for H2O'
       Call Debug_Msg(Trim(err_msg))
       RC = -1
    End If

  end subroutine init_sect_aer
#endif
#if defined( MDL_BOX )
  subroutine init_sect_aer( n_boxes, n_bins, RC )

    USE sect_aer_data_mod,    Only : aer_allocate_ini
    integer, intent(in) :: n_bins, n_boxes
    integer, optional, intent(out) :: RC
    integer             :: k, RC_temp, AS
    character(len=255)  :: err_msg
 
    n_aer_bin = n_bins 
    Call AER_allocate_ini(n_boxes,RC_temp)
    If (RC_temp.ne.0) Then
       If (present(RC)) RC = rc_temp
       Write(err_msg,'(a,I4)') 'Bad allocation exit code: ', RC_temp
       Call error_stop(err_msg,'init_sect_aer (sect_aer_mod.F90)')
    End If

    ! Allocate module variables
    Allocate(ID_Bins(n_aer_bin,n_aer_type),STAT=AS)
    If (AS.ne.0) Then
       Call Debug_Msg('Failed to allocate ID_Bins')
       RC = -1
       Return
    End If
    ID_Bins = 0

    Do k=1,n_bins
       ID_Bins(k,1) = k
    End Do
    id_H2O = n_bins + 1
    id_H2SO4 = n_bins + 2
    RC_temp = 0

    if (present(RC)) RC = RC_temp

  end subroutine
#endif
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cleanup_sect_aer
!
! !DESCRIPTION: Subroutine Cleanup\_Sect\_Aer cleans up any leftover module
!  variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Sect_Aer( RC )
!
! !USES:
!
    USE sect_aer_data_mod,    Only : aer_cleanup
!
! !OUTPUT VARIABLES:
!
    INTEGER,        INTENT(OUT)   :: RC                  ! Return code
! 
! !REVISION HISTORY:
!  28 Nov 2018 - S.D.Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
!    LOGICAL                :: prtDebug,  IsLocNoon

    !=======================================================================
    ! Cleanup_Sect_Aer begins here!
    !=======================================================================

    ! Run cleanup routine
    Call AER_cleanup()

    If (Allocated(ID_Bins      )) Deallocate(ID_Bins      )
    If (Allocated(aero_grow_d  )) Deallocate(aero_grow_d  )
    If (Allocated(aero_nrate_d )) Deallocate(aero_nrate_d )

#if defined( MDL_BOX )
    RC = 0
#else
    RC = GC_SUCCESS
#endif

  end subroutine cleanup_sect_aer
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: aer_wpden
!
! !DESCRIPTION: Subroutine AER\_WPDen calculates the weight percentage and
! density of all aerosol bins.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AER_wpden(aWP_Box, aDen_Box, ST_Box, BVP_Box,&
                       T_K,     P_hPa,    H2O_vv, rWet)
!
! !USES:
!
! -- Nothing --
!
! !INPUT VARIABLES:
!
    REAL(fp),       INTENT(IN)    :: T_K         ! Temperature in K
    REAL(fp),       INTENT(IN)    :: P_hPa       ! Pressure in hPa
    REAL(fp),       INTENT(IN)    :: H2O_vv      ! VMR of water vapor
! 
! !OUTPUT VARIABLES:
!
    REAL(fp),       INTENT(INOUT) :: aWP_Box(n_aer_bin)  ! Weight percentage SO4
    REAL(fp),       INTENT(INOUT) :: aDen_Box(n_aer_bin) ! Aerosol density (g/cm3)
    REAL(fp),       INTENT(INOUT) :: ST_Box(n_aer_bin)   ! Surface tension (N/m?)
    REAL(fp),       INTENT(INOUT) :: BVP_Box(n_aer_bin)  ! Eq. vap. pressure (?)
    REAL(fp),       INTENT(OUT)   :: rWet(n_aer_bin)     ! Wet radius (um)
! 
! !REVISION HISTORY:
!  28 Nov 2018 - S.D.Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    real(dp) :: pph2o, pph2osat, aw, h2o, &
                sulfmolal, wpp, xa,  &
                wp1, wp2, den1, den2, fac, change, &
                y1, y2, temp, vpice, pvh2o

    real(dp) :: aer_vol_wet(n_aer_bin), aer_r_wet(n_aer_bin) !eth_af_dryS

    real(dp) :: pph2o_ref

    integer :: i,j,k,nsize

    !=======================================================================
    ! AER_WPDEN begins here!
    !=======================================================================

    nsize   = n_aer_bin

    if (sum(aWP_Box) .le. 0. .or. sum(aDen_Box) .le. 0.) then 
       ! For the first model time step, no values set for aWP or aDEN
       aer_vol_wet(:)=aer_mass(:)/den_h2so4
    else 
       ! Calculate wet volume of aerosol size bins !eth_af_dryS            
       aer_vol_wet(:)=aer_mass(:)/aWP_Box(:)/.01/aDen_Box(:)*1.E12
    end if
    aer_r_wet(:)=(aer_vol_wet(:)*3./4./PI)**(1./3.)!wet radius !eth_af_dryS

    ! Store for output
    rWet(:) = aer_r_wet(:)

    ! Get VMR of H2O in v/v
    h2o=max(epsilon(1._dp),H2O_vv)
    temp=max(190.15_dp,T_K)
    temp=min(300._dp,temp)

    ! Calculate weight % h2so4 in aerosols as fn of tmp
    ! and h2o partial pressure using parameterization of
    ! Tabazadeh, Toon, Clegg, and Hamill, GRL, 24, 1931, 1997
    vpice = exp(9.550426_dp - 5723.265_dp/temp +  &
             3.53068_dp*log(temp) - 0.00728332*temp)
    vpice=vpice/100._dp  
    ! Partial pressure of ambient h2o in hPa/mbar
    pph2o=h2o*p_hPa
    pph2o=min(vpice,pph2o)
    pph2o_ref = pph2o

    ! Saturation water vapor partial pressure in hPa
    pph2osat=18.452406985_dp-3505.1578807_dp/temp &
             -330918.55082_dp/(temp**2)             &
             +12725068.262/(temp**3)
    pph2osat=exp(pph2osat)
    call pvolume_h2o(temp,pvh2o)

    do k=nsize,1,-1
#if !defined( AER_CONSTANT_WP )
       pph2o = pph2o_ref
       if (k<nsize.and.ST_Box(k)==0.) then
          pph2o=pph2o/exp(2.*ST_Box(k+1)*pvh2o/(aer_r_wet(k)*1E-4*8.314E7*temp))
       else
          pph2o=pph2o/exp(2.*ST_Box(k)*pvh2o/(aer_r_wet(k)*1E-4*8.314E7*temp))
       endif
#endif
!          water activity
       aw=min(0.999_dp,max(epsilon(1._dp),pph2o/pph2osat))

       if(aw .le. 0.05_dp .and. aw > 0._dp) then
          y1=12.372089320_dp*aw**(-0.16125516114_dp) &
              -30.490657554_dp*aw -2.1133114241_dp
          y2=13.455394705_dp*aw**(-0.19213122550_dp) &
              -34.285174607_dp*aw -1.7620073078_dp
       else if(aw .le. 0.85_dp .and. aw > 0.05_dp) then
          y1=11.820654354_dp*aw**(-0.20786404244_dp) &
              -4.8073063730_dp*aw -5.1727540348_dp
          y2=12.891938068_dp*aw**(-0.23233847708_dp) &
              -6.4261237757_dp*aw -4.9005471319_dp
       else
          y1=-180.06541028_dp*aw**(-0.38601102592_dp) &
              -93.317846778_dp*aw +273.88132245_dp
          y2=-176.95814097_dp*aw**(-0.36257048154_dp) &
              -90.469744201_dp*aw +267.45509988_dp
       end if
       !sulfmolal=y1+((T_K-190._dp)*(y2-y1)/70._dp)
       sulfmolal=y1+((temp-190._dp)*(y2-y1)/70._dp)
       aWP_Box(k)=9800.*sulfmolal/(98.*sulfmolal+1000._dp)

       if (aWP_Box(k) .lt. 15._dp) then
          aWP_Box(k)=15._dp
       end if
       if(aWP_Box(k) .gt. 100._dp) then
          aWP_Box(k)=100._dp
       end if
!           mass fraction
       wpp=aWP_Box(k)*0.01_dp
       
!           mole fraction
       xa=18.*wpp/(18.*wpp+98.*(1.-wpp))

!          CALCULATE DENSITY OF AEROSOLS (GM/CC) AS FN OF WT %
!          AND TEMPERATURE
       aDen_Box(k)=density(temp,wpp)
!          CALCULATE SURFACE TENSION OF AEROSOLS AGAINST AIR

       ST_Box(k)=surftension(temp,xa)

!      ! Calculate equilibrium vapor pressure of h2so4
       BVP_Box(k)=solh2so4(temp,xa)
    end do

  END SUBROUTINE AER_wpden
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: density
!
! !DESCRIPTION: Function Density calculates the density of SO4 aerosol bins
! as a function of local temperature and acid mass fraction. Based on
! Vehkamaeki et al (2002).
!\\
!\\
! !INTERFACE:
!
  function density(temp,so4mfrac)
!
! !USES:
!
!    -- NOTHING --
!
! !INPUT VARIABLES:
!
    real(dp), intent(in) :: temp, so4mfrac
! 
! !OUTPUT VARIABLES:
!
    REAL(fp)             :: density ! Density (g/cm3)
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(dp), parameter :: &
               a1= 0.7681724_dp,&
               a2= 2.184714_dp, &
               a3= 7.163002_dp, &
               a4=-44.31447_dp, &
               a5= 88.74606_dp, &
               a6=-75.73729_dp, &
               a7= 23.43228_dp
    real(dp), parameter :: &
               b1= 1.808225e-3_dp,&
               b2=-9.294656e-3_dp,&
               b3=-3.742148e-2_dp,&
               b4= 2.565321e-1_dp,&
               b5=-5.362872e-1_dp,&
               b6= 4.857736e-1_dp,&
               b7=-1.629592e-1_dp
    real(dp), parameter :: &
               c1=-3.478524e-6_dp,&
               c2= 1.335867e-5_dp,&
               c3= 5.195706e-5_dp,&
               c4=-3.717636e-4_dp,&
               c5= 7.990811e-4_dp,&
               c6=-7.458060e-4_dp,&
               c7= 2.581390e-4_dp
    real(dp) :: a,b,c,so4m2,so4m3,so4m4,so4m5,so4m6

    !=======================================================================
    ! DENSITY begins here!
    !=======================================================================

    so4m2=so4mfrac*so4mfrac
    so4m3=so4mfrac*so4m2
    so4m4=so4mfrac*so4m3
    so4m5=so4mfrac*so4m4
    so4m6=so4mfrac*so4m5

    a=+a1+a2*so4mfrac+a3*so4m2+a4*so4m3 &
            +a5*so4m4+a6*so4m5+a7*so4m6
    b=+b1+b2*so4mfrac+b3*so4m2+b4*so4m3 &
            +b5*so4m4+b6*so4m5+b7*so4m6
    c=+c1+c2*so4mfrac+c3*so4m2+c4*so4m3 &
            +c5*so4m4+c6*so4m5+c7*so4m6
    density=(a+b*temp+c*temp*temp) ! units are gm/cm**3
  end function density

!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: surftension
!
! !DESCRIPTION: Function SurfTension calculates the surface tension of
! the aerosol in air. Formulae from Vehkamaeki et al (2002).
!\\
!\\
! !INTERFACE:
!
      function surftension(temp,so4frac)
!
! !USES:
!
!
! !INPUT VARIABLES:
!
    real(dp),intent(in) :: temp, so4frac
! 
! !OUTPUT VARIABLES:
!
    REAL(fp)            :: surftension ! Surface tension (N/m?)
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      real(dp) :: a,b,so4mfrac,so4m2,so4m3,so4m4,so4m5,so4sig
      real(dp), parameter :: &
                 a1= 0.11864_dp,&
                 a2=-0.11651_dp,&
                 a3= 0.76852_dp,&
                 a4=-2.40909_dp,&
                 a5= 2.95434_dp,&
                 a6=-1.25852_dp
      real(dp), parameter :: &
                 b1=-1.5709e-4_dp,&
                 b2= 4.0102e-4_dp,&
                 b3=-2.3995e-3_dp,&
                 b4= 7.611235e-3_dp,&
                 b5=-9.37386e-3_dp,&
                 b6= 3.89722e-3_dp
      real(dp), parameter :: convfac=1.e3_dp  ! convert from newton/m to dyne/cm

!        so4 mass fraction
      so4mfrac=Ma*so4frac/(Ma*so4frac+Mw*(1._dp-so4frac))
      so4m2=so4mfrac*so4mfrac
      so4m3=so4mfrac*so4m2
      so4m4=so4mfrac*so4m3
      so4m5=so4mfrac*so4m4

      a=+a1+a2*so4mfrac+a3*so4m2+a4*so4m3+a5*so4m4+a6*so4m5
      b=+b1+b2*so4mfrac+b3*so4m2+b4*so4m3+b5*so4m4+b6*so4m5
      so4sig=a+b*temp
      surftension=so4sig*convfac

      end function surftension
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: solh2so4
!
! !DESCRIPTION: Function SolH2SO4 calculates the activity of acid.
!\\
!\\
! !INTERFACE:
!
  function solh2so4(T,xa)
!
! !USES:
!
! -- Nothing --
!
! !INPUT VARIABLES:
!
    real(dp),intent(in) :: T, xa
! 
! !OUTPUT VARIABLES:
!
    REAL(fp)            :: solh2so4 ! ???
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(dp) :: xw, a12,b12, cacta

    !=======================================================================
    ! SOLH2SO4 begins here!
    !=======================================================================

    xw=1.0-xa
    ! Compute activity of acid
    a12=5.672E3_dp -4.074E6_dp/T +4.421E8_dp/(T*T)
    b12=1._dp/0.527_dp
    cacta=10._dp**(a12*xw*xw/(xw+b12*xa)**2/T)
    solh2so4=cacta*xa*sath2so4(T)
  end function solh2so4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: solh2o
!
! !DESCRIPTION: Function SolH2O calculates the activity of water
!\\
!\\
! !INTERFACE:
!
  function solh2o(T,xa)
!
! !USES:
!
! -- Nothing --
!
! !INPUT VARIABLES:
!
    real(dp),intent(in) :: T, xa
! 
! !OUTPUT VARIABLES:
!
    REAL(fp)            :: solh2o ! ???
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(dp) :: xw, a11,b11, cactw

    !=======================================================================
    ! SOLH2O begins here!
    !=======================================================================

    xw=1.0-xa
    ! Compute activity of acid
    a11=2.989E3 -2.147E6_dp/T + 2.33E8_dp/(T*T)
    b11=0.527_dp
    cactw=10.**(a11*xa*xa/(xa+b11*xw)**2/T)
    solh2o=cactw*xw*sath2o(T)
  end function solh2o
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sath2o
!
! !DESCRIPTION: Function SatH2O calculates saturation water vapor partial
! pressure in hPa (mbar) based on different approximations for different
! temperature ranges.
!\\
!\\
! !INTERFACE:
!
  function sath2o(T)
!
! !USES:
!
! -- Nothing --
!
! !INPUT VARIABLES:
!
    real(dp),intent(in) :: T
! 
! !OUTPUT VARIABLES:
!
    REAL(dp)            :: sath2o ! Saturation water vapor pressure in hPa
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(dp), parameter :: akb=1.38054e-16_dp
    real(dp) :: ppwater, pph2osat

    !=======================================================================
    ! SATH2O begins here!
    !=======================================================================

    if(T.gt.229.) then
!       Preining et al., 1981 (from Kulmala et al., 1998)
!       saturation vapor pressure (N/m**2)
       ppwater=exp(77.34491296_dp - 7235.424651_dp/T &
                -8.2_dp*log(T) + 5.7133e-3*T)
!       saturation number densities (cm**-3)
       sath2o=10._dp/(akb*T)*ppwater
    else
!       Tabazadeh et al., 1997, parameterization for 185<T<260
!             saturation water vapor partial pressure (mb)
       pph2osat=18.452406985_dp-3505.1578807_dp/T &
                -330918.55082_dp/T/T             &
                +12725068.262_dp/(T**3)
       pph2osat=exp(pph2osat)
       sath2o=1.e3_dp/(akb*T)*pph2osat
    end if

  end function sath2o
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: sath2so4
!
! !DESCRIPTION: Function sath2so4 calculates the saturation number density of 
! H2SO4 in air.
!\\
!\\
! !INTERFACE:
!
  function sath2so4(T)
!
! !USES:
!
! -- Nothing --
!
! !INPUT VARIABLES:
!
    real(dp),intent(in) :: T ! Temperature (K)
! 
! !OUTPUT VARIABLES:
!
    REAL(dp)            :: sath2so4 ! Saturation vapor # density of H2SO4
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(dp) :: ppacid
    !real(dp), parameter ::      &
    !          b1=1.01325e5_dp, &
    !          b2=11.5_dp,  &
    !          b3=1.0156e4_dp,  &
    !          b4=0.38_dp/545._dp, &
    !          tref=360.15_dp
    real(dp), parameter ::      &
              b1=1.01325e5_dp, &
              b2=11.695_dp,  &
              b3=1.0156e4_dp,  &
              b4=0.38_dp/545._dp, &
              tref=360.15_dp
    real(dp), parameter :: akb=1.38054e-16_dp

    !=======================================================================
    ! SATH2SO4 begins here!
    !=======================================================================

    ! Saturation vapor pressure (N/m**2)
    ppacid=b1*exp(-b2+b3*(1._dp/tref-1._dp/T &
           +b4*(1._dp+log(tref/T)-tref/T)))
    ! Saturation number density (cm**-3)
    sath2so4=ppacid*10._dp/(akb*T)

  end function sath2so4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pvolume_h2o
!
! !DESCRIPTION: Subroutine PVolume\_H2O calculates ???
!\\
!\\
! !INTERFACE:
!
  subroutine pvolume_h2o(t1,pvh2o)
!
! !USES:
!
! -- Nothing --
!
! !INPUT VARIABLES:
!
    real(dp), intent(in) :: t1 
! 
! !OUTPUT VARIABLES:
!
    real(dp), intent(out) :: pvh2o
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(dp), parameter :: x1=2.393284E-02,x2=-4.359335E-05,x3=7.961181E-08

    !=======================================================================
    ! PVOLUME_H2O begins here!
    !=======================================================================

    ! 1000: L/mole -> cm3/mole
     pvh2o=(x1+x2*T1+x3*T1**2)*1000. 
       
  end subroutine pvolume_h2o
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pvolume_h2so4
!
! !DESCRIPTION: Subroutine PVolume\_H2SO4 calculates ???
!\\
!\\
! !INTERFACE:
!
  subroutine pvolume_h2so4(t1,ws,pvh2so4)
!
! !USES:
!
! -- Nothing --
!
! !INPUT VARIABLES:
!
      real(dp), intent(in) :: t1, ws
! 
! !OUTPUT VARIABLES:
!
      real(dp), intent(out) :: pvh2so4
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      real(dp), dimension(22),parameter :: x=(/  &
       2.393284E-02,-4.359335E-05,7.961181E-08,0.0,-0.198716351, &
       1.39564574E-03,-2.020633E-06,0.51684706,-3.0539E-03,4.505475E-06, &
      -0.30119511,1.840408E-03,-2.7221253742E-06,-0.11331674116, &
       8.47763E-04,-1.22336185E-06,0.3455282,-2.2111E-03,3.503768245E-06, &
      -0.2315332,1.60074E-03,-2.5827835E-06/)
      real(dp) :: w
      
    !=======================================================================
    ! PVOLUME_H2SO4 begins here!
    !=======================================================================

    w=ws*0.01
    pvh2so4=x(5)+x(6)*T1+x(7)*T1**2+(x(8)+x(9)*T1+x(10)*T1**2)*w &
         +(x(11)+x(12)*T1+x(13)*T1**2)*w*w
    pvh2so4=pvh2so4*1000.
 
  end subroutine pvolume_h2so4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cal_coag_kernel_implicit
!
! !DESCRIPTION: Subroutine Cal\_Coag\_Kernel\_Implicit calculates the
! coagulation kernel for the implicit scheme
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE cal_coag_kernel_implicit(fijk_Box,aWP_Box,aDen_Box)
!
! !USES:
!
! -- Nothing --
!
! !INPUT VARIABLES:
!
      real(dp), intent(in)    :: aWP_Box(n_aer_bin)
      real(dp), intent(in)    :: aDen_Box(n_aer_bin)
!
! !OUPUT PARAMETERS:
!
      real(dp), intent(inout) :: fijk_Box(n_aer_bin,n_aer_bin,n_aer_bin)
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer :: i,j,k,nsize
    real(dp) :: aer_vol_wet(n_aer_bin)
    real(dp) :: vij
      
    !=======================================================================
    ! CAL_COAG_KERNEL_IMPLICIT begins here!
    !=======================================================================

    ! Calculate wet volume of aerosol size bins
    aer_vol_wet(:)=aer_mass(:)/aWP_Box(:)/.01/aDen_Box(:)*1.E12 
    
    ! Copy this over for expediency
    nsize=n_aer_bin
    do i=1,nsize
    do j=1,nsize
       vij = aer_vol_wet(i) + aer_vol_wet(j)
       do k=1,nsize
          fijk_box(i,j,k) = 0._dp
          if (k.eq.nsize) then
             if (vij.ge.aer_vol_wet(k)) fijk_box(i,j,k) = 1.0_dp
          end if
          if (k.lt.nsize.and.k.gt.1) then
             if (vij.lt.aer_vol_wet(k+1).and.vij.ge.aer_vol_wet(k)) then
                fijk_box(i,j,k) = &
                   (aer_vol_wet(k+1)-vij)/(aer_vol_wet(k+1)-aer_vol_wet(k))*aer_vol_wet(k)/vij
             end if
             if (vij.lt.aer_vol_wet(k).and.vij.gt.aer_vol_wet(k-1)) then
                fijk_box(i,j,k) = 1._dp - fijk_box(i,j,k-1)
             end if
          end if
       end do
    end do
    end do

  END SUBROUTINE cal_coag_kernel_implicit
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cal_coag_kernel
!
! !DESCRIPTION: Subroutine Cal\_Coag\_Kernel calculates the coagulation
! kernel
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE cal_coag_kernel(CK_Box,aWP_Box,aDen_Box,T_K,P_hPa)
!
! !USES:
!
! -- Nothing --
!
! !INPUT VARIABLES:
!
      real(dp), intent(in)    :: aWP_Box(n_aer_bin)
      real(dp), intent(in)    :: aDen_Box(n_aer_bin)
      real(dp), intent(in)    :: T_K, P_hPa
!
! !OUPUT PARAMETERS:
!
      real(dp), intent(inout) :: CK_Box(n_aer_bin,n_aer_bin)
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer :: j,k,nsize

    real(dp), dimension(n_aer_bin) :: D,G,DEL,AR

    real(dp) :: LB, FP, AMASS, XKN, B, AAA, BBB, RR, DD, GG, DDEL, &
                FAC1, FAC2, FAC3, vis

    real(dp), parameter :: FP0=6.6E-6_dp, P0=1013._dp, T0=293.15_dp
    
    real(dp) :: aer_vol_wet(n_aer_bin), aer_r_wet(n_aer_bin) !eth_af_dryS
      
    !=======================================================================
    ! CAL_COAG_KERNEL begins here!
    !=======================================================================

    ! Calculate wet volume of aerosol size bins
    aer_vol_wet(:)=aer_mass(:)/aWP_Box(:)/.01/aDen_Box(:)*1.E12 
    ! Wet radius
    aer_r_wet(:)=(aer_vol_wet(:)*3./4./PI)**(1./3.)
    
    AR=aer_r_wet*1.e-4_dp
    
    if (T_K<150. .or. P_hPa<1.) then
        CK_Box=0.
        return
    endif

!   Viscosity
    vis=1.8325E-4_dp*(416.16_dp/(T_K+120._dp))*(T_K/296.16_dp)**1.5_dp

!         MEAN FREE PATH OF AIR MOLECULES
    FP=FP0*(P0/P_hPa)*(T_K/T0)

    ! Copy this over for expediency
    nsize=n_aer_bin
    do  J=1,NSIZE
       AMASS=aer_mass(J)/aWP_Box(J)/0.01
       XKN=FP/AR(J)
       ! Kasten (1968) give 1.249, not 1.246 (SDE 2019-04-16)
       B=(1._dp/(6._dp*pi*VIS*AR(J)))* &
         (1._dp+1.249_dp*XKN+0.42_dp*XKN*EXP(-0.87_dp/XKN))
       D(J)=akb*T_K*B
       G(J)=SQRT(8._dp*akb*T_K/(pi*AMASS))
       !LB=2.546_dp*D(J)/G(J)   ! Seinfeld and Pandis (2006)
       LB=(2._dp/pi)*D(J)/G(J)  ! Jacobson (1999)
       AAA=(2._dp*AR(J)+LB)**3
       BBB=(4._dp*AR(J)*AR(J)+LB*LB)**(3./2.)
       DEL(J)=(1._dp/(6._dp*AR(J)*LB))*(AAA-BBB)-2._dp*AR(J)
    enddo
    do J=1,NSIZE
       do K=1,NSIZE
          RR=AR(J)+AR(K)
          DD=D(J)+D(K)
          GG=SQRT(G(J)*G(J)+G(K)*G(K))
          DDEL=SQRT(DEL(J)*DEL(J)+DEL(K)*DEL(K))
          FAC1=RR/(RR+DDEL)
          FAC2=4._dp*DD/(GG*RR)
          FAC3=1._dp/(FAC1+FAC2)
          CK_Box(J,K)=4._dp*pi*RR*DD*FAC3
       enddo
    enddo

  END SUBROUTINE cal_coag_kernel
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: aer_grow
!
! !DESCRIPTION: Subroutine AER\_grow alculates and applies particle growth
! and evaporation 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AER_grow (T_K,dens,dt,aWP_Box,aDen_Box,ST_Box,  &
                       BVP_Box,H2SO4_vv,Sul_vv,box_grow_d)
!
! !USES:
!
! -- Nothing --
!
! !INPUT VARIABLES:
!
    real(dp), intent(in) :: T_K                 ! Temperature (K)
    real(dp), intent(in) :: dens                ! Air molecular # density (molec/cm3)
    real(dp), intent(in) :: dt                  ! Timestep length (s)
    real(dp), intent(in) :: aWP_Box(n_aer_bin)  ! SO4 weight percentage
    real(dp), intent(in) :: aDen_Box(n_aer_bin) ! Aerosol density (g/cm3)
    real(dp), intent(in) :: ST_Box(n_aer_bin)   ! Surface tension (N/m?)
    real(dp), intent(in) :: BVP_Box(n_aer_bin)  ! Eq vapor pressure (#/cm3)
!
! !OUTPUT VARIABLES:
!
    REAL(dp), INTENT(inout) :: box_grow_d
    REAL(dp), INTENT(inout) :: H2SO4_vv
    REAL(dp), INTENT(inout) :: Sul_vv(n_aer_bin)
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer  :: k, nsize
    real(dp) :: DEN1, AIR, SULFATE, HAMU, SULF, FRA, avgmolwt, &
                volume, diff, AD, per, ave, rad, vknud, corr,  &
                curvature(n_aer_bin), prediff, growth,&
                ratio, EFF, D, ratio_aer !eth_af_ratio

    real(dp) :: dn(n_aer_bin), dgas, change(n_aer_bin)
    real(dp) :: h2so4, h2so40, an(n_aer_bin), an0(n_aer_bin), &
                total1,total2, temp, pvh2so4
    
    real(dp) :: aer_vol_wet(n_aer_bin), aer_r_wet(n_aer_bin)  !eth_af_dryS

    real(dp) :: pv_old, wpp

    !=======================================================================
    ! AER_GROW begins here!
    !=======================================================================

    DEN1=0.2*1.14+0.8*0.808
!C  AIR IS DIAMETER OF AN AIR MOLECULE, SULFATE IS DIAMETER OF AN H2SO4
    AIR=2.*(28.9/((4./3.)*PI*DEN1*AV))**0.333333
    SULFATE=2.*(98./((4./3.)*PI*1.8409*AV))**0.333333
    D=0.5*(AIR+SULFATE)
!C  HAMU IS AMU TO THE POWER OF 0.5
    HAMU=(28.9/(28.9+98.))**0.5
    SULF=98./av
      
    aer_vol_wet(:)=aer_mass(:)/aWP_Box(:)/.01/aDen_Box(:)*1.E12 !calculate wet volume of aerosol size bins !eth_af_dryS            
    aer_r_wet(:)=(aer_vol_wet(:)*3./4./PI)**(1./3.)!wet radius !eth_af_dryS
   
    h2so4 = h2so4_vv*dens 
    an(:)    = Sul_vv(:)*dens/aer_molec !particle number density
    h2so40=h2so4
    an0=an
    
    total1=h2so4 + sum(an*aer_molec)
    temp = T_K
    if (sum(an) .eq. 0.) return
    
    DN=0._dp
    DGAS=0._dp
    !This scheme is based of Hamill et al. 1977 -> "Microphysical processes affecting stratospheric aerosol particles"    
    AD=dens
    EFF=1./((PI*AD*D**2)*HAMU)
!   DIFF IS THE DIFFUSION COEFFICIENT
    DIFF=5./(16.*D**2*AD)*SQRT(AKB*Temp/(2.*PI*SULF)* &
         ((28.9+98.)/28.9))

    nsize = n_aer_bin
    DO K=1,n_aer_bin
#if defined( AER_PVH2SO4 )
       ! Calculation in original code
       wpp = 0.01 * aWP_Box(k)
       fra = 18.*wpp/(18.*wpp+98.*(1.-wpp))
       avgmolwt = fra*98.+(1.-fra)*18.
       pv_old = avgmolwt/aDen_box(k)
       pvh2so4 = pv_old
#else
       ! Calculation in SOCOL (Krieger et al 2000, 
       ! "Measurement of the refractive indices of 
       ! H2SO4–HNO3–H2O solutions to stratospheric 
       ! temperatures")
       call pvolume_h2so4(temp,aWP_Box(k),pvh2so4)
#endif

       RAD=aer_r_wet(K)*1.E-4
       VKNUD=EFF/RAD
       CORR=(1.3333+0.71/VKNUD)/(1.+1./VKNUD)
!   CURVATURE CORRECTION (KELVIN EFFECT)
       curvature(k)=exp(2.*pvh2so4*ST_Box(k)/(rad*Rgas*T_K))
       !curvature(k)=exp(2.*avgmolwt*st(i,j,k)/(rad*Rgas*t(i,j)*aDEN(i,j,k)))
       PREDIFF=H2SO4-BVP_Box(k)*curvature(k)
!   CHANGE IN number of condensed H2SO4 molecules PER SECOND !eth_af_dryS 
       GROWTH=4.*PI*RAD*DIFF*PREDIFF/(1.+CORR*VKNUD) !d[molec]/dt
       GROWTH=GROWTH*dt !d[molec]
       CHANGE(k)=(GROWTH/aer_molec(k))*an(K) !(d[molec])/molec*particle number density
       if (an(k) <= 0.) cycle
       IF (GROWTH.LT.0.) THEN
          IF (K.GT.1) then
             IF (aer_Vrat/(aer_Vrat-1.)*abs(CHANGE(k)).GT.an(K)) &
                 CHANGE(k)=-0.95*an(K)*(aer_Vrat-1.)/aer_Vrat
             DN(K)=DN(K)+CHANGE(k)*aer_Vrat/(aer_Vrat-1.)
             DN(K-1)=DN(K-1)-CHANGE(k)*aer_Vrat/(aer_Vrat-1.)
          else
             if(abs(change(k)).gt.an(K)) change=-0.95*an(K)
             DN(k)=DN(k)+change(k)
          end if
       ELSE IF (GROWTH.GT.0.) THEN
          IF (K.LT.NSIZE) THEN
             IF (CHANGE(k)/(aer_Vrat-1.).GT.an(K)) &
                 CHANGE(k)=0.95*an(K)*(aer_Vrat-1.)
             DN(K)=DN(K)-CHANGE(k)/(aer_Vrat-1.)
             DN(K+1)=DN(K+1)+CHANGE(k)/(aer_Vrat-1.)
          else
             DN(k)=DN(k)+change(k)
          END IF
       END IF
    enddo
       
    dgas=-min(h2so4,sum(dn(:)*aer_molec(:)))    
    ratio=1.
    ratio_aer=1.

    !ensure that h2so4 concentration is at least equal to saturation vapour pressure
    if(dgas .lt. 0._dp) then
       if( (h2so4+dgas) .lt. bvp_box(n_aer_bin)) then
           ratio=(h2so4-bvp_box(n_aer_bin))/abs(dgas)
       end if
    end if 
    dgas=dgas*ratio
    ! ensure that the flux from/into the gas phase is equal to the flux into/from
    ! the particle phase
    ratio_aer=ratio
    if(abs(dgas) .gt. epsilon(1._dp)) then
       ratio_aer=abs(dgas)/abs(sum(dn(:)*aer_molec(:)))
    end if
    dn=dn*ratio_aer      
    !-eth_af_ratio

    !update h2so4 and particle number densities
    H2SO4=min(total1,max(0.,H2SO4+dgas))
    an=max(0.,an+dn)
    
    total2=h2so4 + sum(an(:)*aer_molec(:))
    
!    if ( abs(sum(an0)-sum(an)) > 1E-4) then     
!      write(*,*) '===== evaporation' 
!      write(*,*) krow,i,j
!      write(*,*) h2so40, bvp(i,j)
!      do k=1,nsize
!         write(*,*) an0(k),an(k),change(k)
!      enddo
!      write(*,*) '=====', abs((total2-total1)/total1)
!      stop
!    endif
    
    if (abs((total2-total1)/total1) > 1E-3 ) then 
       !write(*,*) 'mass issue in grow'
       !write(*,*) krow,i,j
       !write(*,*) total1,total2, (total1-total2)/total1
       h2so4=h2so4*total1/total2
       an=an*total1/total2
       total2=h2so4 + sum(an*aer_molec)
       !stop
    endif
    ! positive: condensation, h2so4(g)-> aerosol; negtaive: evaporate aerosol->h2so4(g)
    box_grow_d = box_grow_d + (h2so4_vv - h2so4/dens)/dt 
    h2so4_vv  = h2so4/dens
    sul_vv(:) = an(1:n_aer_bin)*aer_molec/dens

  END SUBROUTINE AER_grow
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: aer_nucleation
!
! !DESCRIPTION: Subroutine AER\_nucleation calculates the total nucleation
! rate for stratospheric sulfates
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AER_nucleation (T_K,P_hPa,dens,dt,aWP_Box,aDen_Box,  &
                             H2O_vv,H2SO4_vv,Sul_vv,box_nrate_d,  &
                             box_wp_d)
!
! !USES:
!
! -- Nothing --
!
!
! !INPUT VARIABLES:
!
    real(dp), intent(in) :: T_K                 ! Temperature (K)
    real(dp), intent(in) :: P_hPa               ! Pressure (hPa)
    real(dp), intent(in) :: dens                ! Air molecular # density (molec/cm3)
    real(dp), intent(in) :: dt                  ! Timestep length (s)
    real(dp), intent(in) :: aWP_Box(n_aer_bin)  ! SO4 weight percentage
    real(dp), intent(in) :: aDen_Box(n_aer_bin) ! Aerosol density (g/cm3)
!
! !OUTPUT VARIABLES:
!
    REAL(dp), INTENT(inout) :: box_nrate_d
    REAL(dp), INTENT(out)   :: box_wp_d    !  aero_wp_d?
    REAL(dp), INTENT(inout) :: H2O_vv
    REAL(dp), INTENT(inout) :: H2SO4_vv
    REAL(dp), INTENT(inout) :: Sul_vv(n_aer_bin)
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
       real(dp) :: n_rate, cluster_r, rupper, dn, dh2so4, h2o,h2so4
       real(dp) :: h2osat, rh,h2so4_threshold, vpice, pph2o, an(n_aer_bin)
       real(dp) :: awp1, aden1, voln(n_aer_bin)
       real(dp) :: aer_vol_wet(n_aer_bin),aer_r_wet(n_aer_bin)  !eth_af_dryS
       
       integer :: itmp, ih2o, jl, jk, k 
       
       real, dimension(7,9),parameter :: h2so4min=reshape( (/ &
      1.00E+03,3.16E+02,1.00E+02,3.16E+00,1.00E-02,1.00E-02,1.00E-02,&
      1.00E+04,3.16E+03,3.16E+02,3.16E+01,3.16E+00,1.00E+00,1.00E+00,&
      3.16E+04,1.00E+04,3.16E+03,1.00E+03,3.16E+01,3.16E+00,3.16E+00,&
      3.16E+05,1.00E+05,3.16E+04,3.16E+03,3.16E+02,3.16E+01,1.00E+01,&
      1.00E+07,3.16E+05,1.00E+05,3.16E+04,3.16E+03,3.16E+02,3.16E+01,&
      1.00E+07,1.00E+07,1.00E+07,1.00E+05,3.16E+04,3.16E+03,1.00E+02,&
      1.00E+07,1.00E+07,1.00E+07,1.00E+07,1.00E+05,3.16E+04,1.00E+03,&
      1.00E+07,1.00E+07,1.00E+07,1.00E+07,1.00E+07,1.00E+05,1.00E+04,&
      1.00E+07,1.00E+07,1.00E+07,1.00E+07,1.00E+07,1.00E+07,1.00E+05/),&
       shape(h2so4min))

    !=======================================================================
    ! AER_NUCLEATION begins here!
    !=======================================================================

      ! Need to do this here rather than later, to enable voln to be calculated
      h2so4=h2so4_vv*dens
      an = Sul_vv(:)*dens

       n_rate=0._dp
       cluster_r=0._dp
       aer_vol_wet(:)=aer_mass(:)/aWP_Box(:)/.01/aDEN_Box(:)*1.E12 !calculate wet volume of aerosol size bins in um^3 !eth_af_dryS            
      aer_r_wet(:)=(aer_vol_wet(:)*3./4./PI)**(1./3.)!wet radius in um !eth_af_dryS
      voln = aer_vol_wet(:)*an(:)/aer_molec(:) !particle volume mixing ratio
      if (sum(voln)>0.) then
         awp1 = sum(voln*awp_Box(:))/sum(voln)
         aden1= sum(voln*aden_Box(:))/sum(voln)
         box_wp_d = sum(voln*awp_box(:))/sum(voln)
      else
         awp1=sum(aer_vol_wet*awp_box(:))/sum(aer_vol_wet)
         aden1=sum(aer_vol_wet*aden_box(:))/sum(aer_vol_wet)
         box_wp_d = sum(aer_vol_wet*awp_box(:))/sum(aer_vol_wet)
      endif
      
      if (t_k>=252._dp .or. p_hPa<=1._dp) return

      h2o=max(1e-6_dp,H2O_vv)
      vpice = exp(9.550426_dp - 5723.265_dp/t_k +  &
               3.53068_dp*log(t_k) - 0.00728332*t_k)
      vpice=vpice/100._dp  
!         partial pressure of ambient h2o
      pph2o=h2o*P_hPa
      pph2o=min(vpice,pph2o)
      h2o = pph2o/P_hPa*dens

      itmp=(t_k-190.)/5. + 1.
      if(itmp.lt.1) itmp=1
      if(itmp.gt.9) itmp=9
      ih2o=(log10(h2o)-12.)/0.5 + 1.95
      if(ih2o.lt.1) ih2o=1
      if(ih2o.gt.7) ih2o=7
      
      if(h2so4.le.h2so4min(ih2o,itmp)) return
      
      call VEHKAMAEKI(t_k,h2o,h2so4,awp1,aden1,n_rate,cluster_r)
      
      if (n_rate<=epsilon(1._dp)) then
         n_rate=0.
         return
      endif
      do k=1, n_aer_bin-1
         !rupper=aer_r_wet(k+1)!assume no nucleation in largest bin
         rupper = aer_r_wet(k) * (2.*aer_vrat/(1.+aer_vrat))**(1./3.)
         if (cluster_r*1.e4_dp .le. rupper) then
            dh2so4=min(n_rate*dt,h2so4)
            !dn=dh2so4*98._dp/(aer_vol_wet(k)*1e-12_dp*aden(i,j,k,krow)*awp(i,j,k,krow)*0.01_dp*Av)                
            h2so4_vv=(h2so4-dh2so4)/dens
            sul_vv(k)=sul_vv(k)+dh2so4/dens
            box_nrate_d = box_nrate_d+dh2so4/dens/dt
            !write(*,*) '==nuclea: dn',krow,i,j,n_rate,dn,dh2so4,aden(i,j),awp(i,j),Av,h2so4,dt,aer_vol_wet(k)
            exit
         endif
      enddo
     
  END SUBROUTINE AER_nucleation
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: vehkamaeki
!
! !DESCRIPTION: Subroutine Vehkamaeki estimates nucleation rates based on
! Vehkamaeki et al (2002)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE VEHKAMAEKI(TEMP,WATER,H2SO4C,WPBACK,DENBACK,DRATE,RADIUS)
!
! !USES:
!
! -- Nothing --
!
! !INPUT VARIABLES:
!
      real(dp), intent(in) :: TEMP,WATER,H2SO4C,WPBACK,DENBACK
!
! !OUTPUT VARIABLES:
!
      real(dp), intent(out) :: drate,radius
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!---->define limit values for nucleation rate and mass of nucleus

      real(dp), parameter :: rnucmax=1.E12_dp
      real(dp), parameter :: rmasmin=2._dp,rmasmax=300._dp
      real(dp), parameter :: tmin=190._dp

      real(dp) :: h2so4sat, h2osat, RH, RA, tempnew, so4new,    &
                  so4nuc, so4mnuc, rnucnum, rnuctot, rnucso4,   &
                  rnucmas,rnucrad,densnuc,volnuc,radnuc

    !=======================================================================
    ! VEHKAMAEKI begins here!
    !=======================================================================

      drate=0._dp
      radius=0._dp

      if(h2so4c.lt.so4gmin) return
!         saturation number densities (cm**-3)

      h2so4sat=sath2so4(temp)
      h2osat=sath2o(temp)
!         RH=relative humidity, RA=relative acidity
      RH=min(1.,water/h2osat)
      RA=min(1.,h2so4c/h2so4sat)
!      write(*,*) 'In Vehkamaeki, T=',temp,'  h2so4=',h2so4c, &
!                  '  h2o=',water,'  h2so4sat=',h2so4sat,     &
!                  '  h2osat=',h2osat,'  rh=',rh,'  ra=',ra


!---->calculate h2so4 mole/mass fraction in critical nucleus
      
      tempnew=max(tmin,temp)
      so4new=min(so4gmax,h2so4c)
!     write(*,*) 'tempnew=',tempnew,'  so4new=',so4new,'  rh=',rh

!---->calculate h2so4 mole fraction of critical nucleus
      so4nuc=so4fracnuc_f(tempnew,rh,so4new)
!      write(*,*) 'so4nuc=so4fracnuc_f(tempnew,rh,so4new)', so4nuc
      so4nuc=min(1.,max(0.,so4nuc))
      so4mnuc=so4nuc*Ma/(so4nuc*Ma+(1.-so4nuc)*Mw)
!      write(*,*) 'so4nuc=',so4nuc,'  so4mnuc=',so4mnuc

!---- >calculate nucleation rate: rnucnum [#particles/cm3]
      rnucnum=rnucnum_f(tempnew,rh,so4new,so4nuc)
!      rnucnum=min(rnucmax,rnucnum)
!      write(*,*) 'rnucnum=',rnucnum
      if (rnucnum.le.0.) return

!---->calculate the number of total and h2so4 molecules in the critical cluster
!         write(*,*) '=tempnew,rh,so4new,so4nuc', tempnew,rh,so4new,so4nuc
         rnuctot=rnucmas_f(tempnew,rh,so4new,so4nuc)
         rnucso4=so4nuc*rnuctot   ! H2SO4 molecules per particle
!           H2SO4 molecules per second nucleated
         rnucso4=min(rmasmax,max(rmasmin,rnucso4))
         rnucmas=rnucnum*rnucso4
!         write(*,*) 'rnuctot=',rnuctot,'  rnucso4=',rnucso4,'  rnucmas=',rnucmas
!!$         if(rnucmas.lt.0.) then
!!$            write(*,*) 'Trouble in Vehkamaeki'
!!$            write(*,*) 'Inputs: T=',TEMP,'  H2O=',WATER,'  H2SO4=', H2SO4C
!!$            write(*,*) 'so4nuc=',so4nuc,so4mnuc,'  rnucnum=',rnucnum
!!$            write(*,*) 'rnuctot=',rnuctot,rnucso4,'  rnucmas=',rnucmas
!!$            return
!!$         end if
!---->calculate the radius of the critical cluster (centimeters)

         rnucrad=1.E-7_dp*exp(-1.6524245_dp + 0.42316402_dp*so4nuc + 0.3346648_dp*log(rnuctot))
         densnuc=density(tempnew,so4mnuc)
         volnuc=rnucso4/Av*Ma/so4mnuc/densnuc
         radnuc=(volnuc*3./4._dp/pi)**(1./3.)
!    adjust particle radius for ambient density and wp
         radius=radnuc*(so4mnuc*densnuc/(wpback*0.01_dp)/denback)**(1./3.)
!         write(*,*) 'rnucrad=',rnucrad,'  volnuc=',volnuc, '  wpback=', wpback, '  denback=',denback, &
!             '  radnuc=',radnuc,'  densnuc=', densnuc,' radius=',radius,'  Av=',Av,'  pi=',pi, '  so4mnuc=',so4mnuc
         drate=rnucmas


  END subroutine VEHKAMAEKI
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: rnucnum_f
!
! !DESCRIPTION: Function rNucNum\_f calculates the nucleation rate in
! particles/cm3 based on Vehkamaeki et al (2002)
!\\
!\\
! !INTERFACE:
!
  function rnucnum_f(temp,relhum,so4gas,so4frac)
!
! !USES:
!
! -- Nothing --
!
! !INPUT VARIABLES:
!
      real(dp), intent(in) :: temp    ! Temperature, K
      real(dp), intent(in) :: relhum  ! RH, fraction (0 to 1)
      real(dp), intent(in) :: so4gas  ! ???
      real(dp), intent(in) :: so4frac ! ???
!
! !OUTPUT VARIABLES:
!
      real(dp) :: rnucnum_f
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      real(dp), parameter ::       &
                 a1=0.14309_dp, &
                 a2=2.21956_dp, &
                 a3=2.73911e-2_dp,&
                 a4=7.22811e-5_dp,&
                 a5=5.91822_dp
      real(dp), parameter ::       &
                 b1=0.117489_dp,&
                 b2=0.462532_dp,&
                 b3=1.18059e-2_dp,&
                 b4=4.04196e-5_dp,&
                 b5=15.7963_dp
      real(dp), parameter ::       &
                 c1=0.215554_dp,&
                 c2=8.10269e-2_dp,&
                 c3=1.43581e-3_dp,&
                 c4=4.7758e-6_dp, &
                 c5=2.91297_dp
      real(dp), parameter ::       &
                 d1=3.58856_dp, &
                 d2=4.9508e-2_dp, &
                 d3=2.1382e-4_dp, &
                 d4=3.10801e-7_dp,&
                 d5=2.93333e-2_dp
      real(dp), parameter ::       &
                 e1=1.14598_dp, &
                 e2=6.00796e-1_dp,&
                 e3=8.64245e-3_dp,&
                 e4=2.28947e-5_dp,&
                 e5=8.44985
      real(dp), parameter ::       &
                 f1=2.15855_dp, &
                 f2=8.08121e-2_dp,&
                 f3=4.07382e-4_dp,&
                 f4=4.01947e-7_dp,&
                 f5=7.21326e-1_dp
      real(dp), parameter ::       &
                 g1=1.6241_dp,  &
                 g2=1.60106e-2_dp,&
                 g3=3.77124e-5_dp,&
                 g4=3.21794e-8_dp,&
                 g5=1.13255e-2_dp
      real(dp), parameter ::       &
                 h1=9.71682_dp, &
                 h2=1.15048e-1_dp,&
                 h3=1.57098e-4_dp,&
                 h4=4.00914e-7_dp,&
                 h5=0.71186_dp
      real(dp), parameter ::       &
                 p1=1.05611_dp, &
                 p2=9.03378e-3_dp,&
                 p3=1.98417e-5_dp,&
                 p4=2.46048e-8_dp,&
                 p5=5.79087e-2_dp
      real(dp), parameter ::       &
                 q1=0.148712_dp,&
                 q2=2.83508e-3_dp,&
                 q3=9.24619e-6_dp,&
                 q4=5.00427e-9_dp,&
                 q5=1.27081e-2_dp
      real(dp), parameter :: expmax=46._dp

      real(dp) :: sfracinv, temp2, temp3
      real(dp) :: a,b,c,d,e,f,g,h,p,q
      real(dp) :: rhln,rhln2,rhln3
      real(dp) :: so4ln,so4ln2,so4ln3,expon

    !=======================================================================
    ! RNUCNUM_F begins here!
    !=======================================================================

      sfracinv=1./so4frac
      temp2=temp*temp
      temp3=temp*temp2

      a=+a1+a2*temp-a3*temp2+a4*temp3+a5*sfracinv
      b=+b1+b2*temp-b3*temp2+b4*temp3+b5*sfracinv
      c=-c1-c2*temp+c3*temp2-c4*temp3-c5*sfracinv
      d=-d1+d2*temp-d3*temp2+d4*temp3-d5*sfracinv
      e=+e1-e2*temp+e3*temp2-e4*temp3-e5*sfracinv
      f=+f1+f2*temp-f3*temp2-f4*temp3+f5*sfracinv
      g=+g1-g2*temp+g3*temp2+g4*temp3-g5*sfracinv
      h=+h1-h2*temp+h3*temp2+h4*temp3+h5*sfracinv
      p=-p1+p2*temp-p3*temp2+p4*temp3-p5*sfracinv
      q=-q1+q2*temp-q3*temp2+q4*temp3-q5*sfracinv


      rhln=log(relhum)
      rhln2=rhln*rhln
      rhln3=rhln*rhln2
      so4ln=log(min(so4gmax,max(so4gmin,so4gas)))
      so4ln2=so4ln*so4ln
      so4ln3=so4ln*so4ln2

      expon=min(expmax,a+b*rhln+c*rhln2+d*rhln3           &
                      +e*so4ln+f*rhln*so4ln+g*rhln2*so4ln &
                      +h*so4ln2+p*rhln*so4ln2+q*so4ln3)
      rnucnum_f=exp(expon)

  end function rnucnum_f
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: so4fracnuc_f
!
! !DESCRIPTION: Function SO4FracNuc\_f calculates the H2SO4 molar fraction
! in a critical nucleus based on Vehkamaeki et al (2002)
!\\
!\\
! !INTERFACE:
!
  function so4fracnuc_f(temp,relhum,so4gas)
!
! !USES:
!
! -- Nothing --
!
! !INPUT VARIABLES:
!
      real(dp), intent(in) :: temp   ! Temperature, K
      real(dp), intent(in) :: relhum ! Relative humidity, fraction (0-1)
      real(dp), intent(in) :: so4gas ! ???
!
! !OUTPUT VARIABLES:
!
      real(dp) :: so4fracnuc_f
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      real(dp), parameter ::        &
                 a1=7.40997e-1_dp, &
                 a2=2.66379e-3_dp, &
                 a3=3.49998e-3_dp, &
                 a4=5.04022e-5_dp, &
                 a5=2.01048e-3_dp, &
                 a6=1.83289e-4_dp, &
                 a7=1.57407e-3_dp, &
                 a8=1.79059e-5_dp, &
                 a9=1.84403e-4_dp, &
                 a10=1.50345e-6_dp


      real(dp) :: rhln,rhln2,rhln3
      real(dp) :: so4ln,so4ln2,so4ln3

    !=======================================================================
    ! SO4FRACNUC_F begins here!
    !=======================================================================

      rhln=log(relhum)
      rhln2=rhln*rhln
      rhln3=rhln*rhln2
      so4ln=log(min(so4gmax,max(so4gmin,so4gas)))
      so4fracnuc_f=a1       -a2*temp           &
                  -a3*so4ln +a4*temp*so4ln     &
                  +a5*rhln  -a6*temp*rhln      &
                  +a7*rhln2 -a8*temp*rhln2     &
                  +a9*rhln3-a10*temp*rhln3

  end function so4fracnuc_f
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: rnucmas_f
!
! !DESCRIPTION: Function rNucMas\_f calculates the total number of molecules
! in a critical cluster (# molec/particle) based on Vehkamaeki et al (2002)
!\\
!\\
! !INTERFACE:
!
  function rnucmas_f(temp,relhum,so4gas,so4frac)
!
! !USES:
!
! -- Nothing --
!
!
! !INPUT VARIABLES:
!
      real(dp), intent(in) :: temp    ! Temperature, K
      real(dp), intent(in) :: relhum  ! Relative humidity, fraction (0-1)
      real(dp), intent(in) :: so4gas  ! ???
      real(dp), intent(in) :: so4frac ! ???
!
! !OUTPUT VARIABLES:
!
      real(dp) :: rnucmas_f
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      real(dp), parameter ::        &
                 a1=2.95413e-3_dp, &
                 a2=9.76834e-2_dp, &
                 a3=1.02485e-3_dp, &
                 a4=2.18646e-6_dp, &
                 a5=1.01717e-1_dp
      real(dp), parameter ::        &
                 b1=2.05064e-3_dp, &
                 b2=7.58504e-3_dp, &
                 b3=1.92654e-4_dp, &
                 b4=6.70430e-7_dp, &
                 b5=2.55774e-1_dp
      real(dp), parameter ::        &
                 c1=3.22308e-3_dp, &
                 c2=8.52637e-4_dp, &
                 c3=1.54757e-5_dp, &
                 c4=5.66661e-8_dp, &
                 c5=3.38444e-2_dp
      real(dp), parameter ::        &
                 d1=4.74323e-2_dp, &
                 d2=6.25104e-4_dp, &
                 d3=2.65066e-6_dp, &
                 d4=3.67471e-9_dp, &
                 d5=2.67251e-4_dp
      real(dp), parameter ::        &
                 e1=1.25211e-2_dp, &
                 e2=5.80655e-3_dp, &
                 e3=1.01674e-4_dp, &
                 e4=2.88195e-7_dp, &
                 e5=9.42243e-2_dp
      real(dp), parameter ::        &
                 f1=3.85460e-2_dp, &
                 f2=6.72316e-4_dp, &
                 f3=2.60288e-6_dp, &
                 f4=1.19416e-8_dp, &
                 f5=8.51515e-3_dp
      real(dp), parameter ::        &
                 g1=1.83749e-2_dp, &
                 g2=1.72072e-4_dp, &
                 g3=3.71766e-7_dp, &
                 g4=5.14875e-10_dp,&
                 g5=2.68660e-4_dp
      real(dp), parameter ::        &
                 h1=6.19974e-2_dp, &
                 h2=9.06958e-4_dp, &
                 h3=9.11728e-7_dp, &
                 h4=5.36796e-9_dp, &
                 h5=7.74234e-3_dp
      real(dp), parameter ::        &
                 p1=1.21827e-2_dp, &
                 p2=1.06650e-4_dp, &
                 p3=2.53460e-7_dp, &
                 p4=3.63519e-10_dp,&
                 p5=6.10065e-4_dp
      real(dp), parameter ::        &
                 q1=3.20184e-4_dp, &
                 q2=1.74762e-5_dp, &
                 q3=6.06504e-8_dp, &
                 q4=1.42177e-11_dp,&
                 q5=1.35751e-4_dp

      real(dp) :: sfracinv, temp2, temp3
      real(dp) :: a,b,c,d,e,f,g,h,p,q
      real(dp) :: rhln,rhln2,rhln3
      real(dp) :: so4ln,so4ln2,so4ln3

    !=======================================================================
    ! RNUCMAS_F begins here!
    !=======================================================================

      sfracinv=1./so4frac
      temp2=temp*temp
      temp3=temp*temp2

      a=-a1-a2*temp+a3*temp2-a4*temp3-a5*sfracinv
      b=-b1-b2*temp+b3*temp2-b4*temp3-b5*sfracinv
      c=+c1+c2*temp-c3*temp2+c4*temp3+c5*sfracinv
      d=+d1-d2*temp+d3*temp2-d4*temp3-d5*sfracinv
      e=-e1+e2*temp-e3*temp2+e4*temp3+e5*sfracinv
      f=-f1-f2*temp+f3*temp2+f4*temp3-f5*sfracinv
      g=-g1+g2*temp-g3*temp2-g4*temp3+g5*sfracinv
      h=-h1+h2*temp-h3*temp2-h4*temp3-h5*sfracinv
      p=+p1-p2*temp+p3*temp2-p4*temp3+p5*sfracinv
      q=+q1-q2*temp+q3*temp2-q4*temp3+q5*sfracinv


      rhln=log(relhum)
      rhln2=rhln*rhln
      rhln3=rhln*rhln2
      so4ln=log(min(so4gmax,max(so4gmin,so4gas)))
      so4ln2=so4ln*so4ln
      so4ln3=so4ln*so4ln2

      rnucmas_f=exp(a+b*rhln+c*rhln2+d*rhln3           &
                   +e*so4ln+f*rhln*so4ln+g*rhln2*so4ln &
                   +h*so4ln2+p*rhln*so4ln2+q*so4ln3)

  end function rnucmas_f
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: aer_coagulation_implicit
!
! !DESCRIPTION: Function AER\_coagulation\_implicit calculates the effect of
! coagulation, redistributing mass between the bins accordingly. This 
! uses an implicit, rather than explicit, numerical integration scheme.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AER_coagulation_implicit(Sul_vv,fijk_box,CK_box,aWP_Box,aDen_Box,dens,dt)
!
! !USES:
!
! -- Nothing --
!
!
! !INPUT VARIABLES:
!
    real(dp), intent(in) :: fijk_box(n_aer_bin,n_aer_bin,n_aer_bin) ! Volume fractions
    real(dp), intent(in) :: CK_box(n_aer_bin,n_aer_bin)             ! Coagulation kernel
    real(dp), intent(in) :: aWP_Box(n_aer_bin)                      ! Wt pcg (%)
    real(dp), intent(in) :: aDen_Box(n_aer_bin)                     ! Density (g/cm3)
    real(dp), intent(in) :: dens                                    ! Air dens (#/cm3)
    real(dp), intent(in) :: dt                                      ! Substep length (s)
!
! !OUTPUT VARIABLES:
!
    REAL(dp), INTENT(inout) :: Sul_vv(n_aer_bin)        ! SO4 in each bin (v/v)
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer :: nlt,nha, jd, nsize, ik, jk, k, I, J
    real(dp) :: dn(n_aer_bin), an0(n_aer_bin),an(n_aer_bin)
    real(dp) :: s,f,al,total1,total2,tmin,mass_error

    real(dp) :: aer_vol_wet(n_aer_bin), aer_r_wet(n_aer_bin) !eth_af_dryS
    integer :: i_limit
   
    real(dp) :: nm, dm, dm1(n_aer_bin) 
    real(dp) :: vkold(n_aer_bin), vknew(n_aer_bin)

    !=======================================================================
    ! AER_COAGULATION_IMPLICIT begins here!
    !=======================================================================

    nsize=n_aer_bin

    !particles/cm^3
    an0=sul_vv(:)*dens/aer_molec

    an = an0
    if (sum(an)<=0.) return

    ! Calculate wet volume of aerosol size bins in um^3 !eth_af_dryS            
    aer_vol_wet(:)=aer_mass(:)/aWP_Box(:)/.01/aDEN_Box(:)*1.E12 

    ! Get the volumetric mixing ratio of wet aerosol
    vkold(:) = an(:) * aer_vol_wet(:)
    vknew(:) = 0.0e+0_dp

    total1=sum(an0*aer_molec)

    !write(*,'(a)') 'Before coag:'
    !do k=1,nsize
    !   write(*,'("--> ",I3,3(x,E16.5E4))') k,an0(k),sul_vv(k),aer_vol_wet(k)
    !end do

    do k=1,nsize
       nm  = 0.0e+0_dp 
       dm1(:) = (1.0e+0_dp - fijk_box(k,:,k))*ck_box(k,:)*an(:)
       dm = 1.0e+0_dp + dt * sum(dm1)
       if (k.eq.1) then
          vknew(k) = vkold(k)/dm
       else
          do jk=1,k
          do ik=1,k-1
             if (fijk_box(ik,jk,k).ne.0.0e+0_dp) then
                nm = nm + fijk_box(ik,jk,k) * ck_box(ik,jk) * vknew(ik) * an(jk)
             end if
          end do
          end do
          vknew(k) = (vkold(k) + dt*nm)/dm
       end if
       an(k) = vknew(k)/aer_vol_wet(k)
    end do
    
    !write(*,'(a)') 'After coag:'
    !do k=1,nsize
    !   write(*,'("--> ",I3,3(x,E16.5E4))') k,an(k),an(k)*aer_molec(k)/dens,aer_vol_wet(k)
    !end do

    total2=sum(an*aer_molec)

    If (total2.ne.total2) then
      call Error_Stop('NaN in coagulation','AER_coagulation_implicit')
    End If

    mass_error = abs((total2-total1)/total1)
    if (mass_error > 1E-6 ) then 
      ! Problem - "mass issue"?
      an=an*total1/total2
      if (mass_error > 0.01e+0 ) then
        write(*,'(a,F10.4,a)') 'MASS ERROR IN I-COAG: ',100.0*mass_error,'%'
      endif
    endif
    Sul_vv(:) = an*aer_molec/dens

  END SUBROUTINE AER_coagulation_implicit
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: aer_coagulation
!
! !DESCRIPTION: Function AER\_coagulation calculates the effect of particle
! coagulation, redistributing mass between the bins accordingly.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AER_coagulation(Sul_vv,CK_box,dens,dt)
!
! !USES:
!
! -- Nothing --
!
!
! !INPUT VARIABLES:
!
    real(dp), intent(in) :: CK_box(n_aer_bin,n_aer_bin) ! Coagulation kernel
    real(dp), intent(in) :: dens                        ! Air dens (#/cm3)
    real(dp), intent(in) :: dt                          ! Substep length (s)
!
! !OUTPUT VARIABLES:
!
    REAL(dp), INTENT(inout) :: Sul_vv(n_aer_bin)        ! SO4 in each bin (v/v)
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    integer :: nlt,nha, jd, nsize, ik, jk, k, I, J
    real(dp) :: dn(n_aer_bin), an0(n_aer_bin),an(n_aer_bin)
    real(dp) :: s,f,al,total1,total2,tmin,mass_error

    integer :: i_limit
    real(dp) :: t_rem
    
    !=======================================================================
    ! AER_COAGULATION begins here!
    !=======================================================================

    nsize=n_aer_bin

    an0=sul_vv(:)*dens/aer_molec !particles/cm^3
    an =an0
    if (sum(an)<=0.) return
    total1=sum(an0*aer_molec)

    t_rem = dt
    do while (t_rem > 1.0e-3_dp)
       dn = 0._dp 
       DO I=1,n_aer_bin-1
          DO J=1,I
             S=1.
             IF (I.EQ.J) S=0.5 
             F=aer_mass(J)/aer_mass(I) !eth_af_dryS
             AL=an(I)*an(J)*S*ck_box(I,J)
             DN(J)=DN(J)-AL
             DN(I)=DN(I)-AL*F/(aer_Vrat-1.)
             DN(I+1)=DN(I+1)+AL*F/(aer_Vrat-1.)
          ENDDO
       ENDDO
       tmin = t_rem
       !i_limit = 0
       do k=1,n_aer_bin
          if(-dn(k)*t_rem > 0.9*an(k)) tmin=min(tmin,-0.9*an(k)/dn(k))
          !if(-dn(k)*dt > 0.9*an(k)) then
          !   tmin=min(tmin,-0.9*an(k)/dn(k))
          !   i_limit = k
          !end if
       enddo
       t_rem = t_rem - tmin

       ! HACK
       !t_rem = 0.0d0
       !if (i_limit > 0) then
       !  write(*,'(a,I3,a,F12.7)') 'LIMITED by bin ', i_limit, ' to ', tmin
       !end if
       an=max(0.,an+dn*tmin)
    end do
    total2=sum(an*aer_molec)
    mass_error = abs((total2-total1)/total1)
    if (mass_error > 1E-3 ) then 
      ! Problem - "mass issue"?
      an=an*total1/total2
      write(*,'(a,F10.4,a)') 'MASS ERROR IN X-COAG: ',100.0*mass_error,'%'
    endif
    
    Sul_vv(:) = an*aer_molec/dens

  END SUBROUTINE AER_coagulation
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: aer_mass_check
!
! !DESCRIPTION: Subroutine AER\_mass\_check enforces mass conservation after
! aerosol microphysics is completed.
!\\
!\\
! !INTERFACE:
!
  subroutine AER_mass_check(h2so4_vv,sul_vv,h2so4_vv_old,sul_vv_old, &
                            dens,s1,s2,rc)
!
! !USES:
!
!
!
! !INPUT VARIABLES:
!
    real(dp), intent(in)    :: h2so4_vv_old
    real(dp), intent(in)    :: sul_vv_old(n_aer_bin)
    real(dp), intent(in)    :: dens
!
! !OUTPUT VARIABLES:
!
    real(dp), intent(inout) :: h2so4_vv
    real(dp), intent(inout) :: sul_vv(n_aer_bin)
    real(dp), intent(out)   :: s1,s2
    integer,  intent(out)   :: rc
! 
! !REVISION HISTORY:
!  04 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(dp) :: h1,h2,t1,t2
    integer :: i,j,k
    character(len=255) :: err_msg

    !=======================================================================
    ! AER_MASS_CHECK begins here!
    !=======================================================================

    ! Assume success
    RC = 0

    ! s1 and s2 are the total sulfur burdens in molec/cm3
    s1=0.
    s2=0.
    h1=h2so4_vv_old*dens
    t1=h1 + sum(sul_vv_old(:)*dens)
    s1=s1+t1
    h2=h2so4_vv*dens
    t2=h2 + sum(sul_vv(:)*dens)
    s2=s2+t2
    if (t1 <= tiny(1.0e+0_dp) .and. t2 <= tiny(1.0e+0_dp)) then
       ! Started and ended with no sulfur - this is fine
       rc = 0
       return
    elseif (t1 <= tiny(1.0e+0_dp) .and. t2 > tiny(1.0e+0_dp)) then
       ! This won't end well
       Call Debug_msg('Mass balance failure in stratospheric aerosol scheme')
       Write(Err_msg,'(a,E16.5e4,a,E16.5e4,a)') 'Started with ', t1, ' molecSO4/cm3, ended with ', t2, ' molecSO4/cm3'
       Call Debug_Msg(Trim(Err_Msg)) 
       RC = -1
    elseif (abs(t2-t1)/t1 > 0.001) then
       ! Attempt a scaling correction
       h2so4_vv=h2so4_vv*t1/t2
       sul_vv(:) = sul_vv(:)*t1/t2
       h2=h2so4_vv*dens
       t2=h2 + sum(sul_vv(:)*dens)

       if (abs(t2-t1)/t1 > 0.1) then
          Call Debug_Msg('Mass balance failure in stratospheric aerosol scheme')
          write(Err_msg,'(a,E16.5E4,E16.5E4,F0.3,a)') 'Totals: prev, curr, error: ', t1,t2,100.0*(t2-t1)/t1,'%'
          Call Debug_Msg(Trim(Err_Msg))
          write(Err_msg,'(a7,E16.5E4,a,E16.5E4)') 'SO4aq  : ', sum(sul_vv_old),' --> ',sum(sul_vv)
          Call Debug_Msg(Trim(Err_Msg))
          write(Err_msg,'(a7,E16.5E4,a,E16.5E4)') 'H2SO4  : ', h2so4_vv_old,' --> ',h2so4_vv
          Call Debug_Msg(Trim(Err_Msg))
          do k=1,n_aer_bin
             write(Err_msg,'(a4,I3,a,E16.5E4,a,E16.5E4)') 'Bin ',k,' : ',sul_vv_old(k),' --> ',sul_vv(k)
             Call Debug_Msg(Trim(Err_Msg))
          enddo
          ! Attempt a correction anyway - but this should really be causing an error
          h2so4_vv = h2so4_vv*t1/t2
          sul_vv(:) = sul_vv(:)*t1/t2
          h2=h2so4_vv*dens
          t2=h2 + sum(sul_vv(:)*dens)
          !write(*,*) 'aer mass correction'
          !write(*,*) t1,t2,(t2-t1)/t1
          RC = -1
       endif
    endif

  end subroutine AER_mass_check
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: aer_sedimentation
!
! !DESCRIPTION: Subroutine AER\_sedimentation simulates gravitational settling
! of sectional stratospheric aerosols.
!\\
!\\
! !INTERFACE:
!
#if defined( MDL_BOX ) 
  subroutine aer_sedimentation()
  write(*,*) 'SEDIMENTATION NOT VALID IN BOX MODEL'
  end subroutine aer_sedimentation
#else
  SUBROUTINE AER_sedimentation(Sul_vv,T_K,P_hPa,aWP_Col,aDen_Col,dt,dens,dz,NZ)
!
! !USES:
!
    USE PhysConstants,   Only : g=>g0
!
!
! !INPUT VARIABLES:
!           
    integer,  intent(in) :: NZ                       ! Number of levels
    real(dp), intent(in) :: T_K(NZ)                  ! Temperature (K)
    real(dp), intent(in) :: P_hPa(NZ)                ! Pressure (hPa)
    !integer,  intent(in) :: LTROP
    real(dp), intent(in) :: dt                          ! Timestep length (s)
    real(dp), intent(in) :: dens(NZ)                 ! Air molecular number density (#/cm3)
    real(dp), intent(in) :: dz(NZ)                   ! Layer thickness (m)
    real(dp), intent(in) :: aWP_Col(NZ,n_aer_bin)    ! H2SO4 wt pcg (%)
    real(dp), intent(in) :: aDen_Col(NZ,n_aer_bin)   ! Aerosol density (g/cm3)
!    real(dp), intent(in) :: dt, xtm1(nbdim,nlev,n_aer_bin), pqm1(nbdim,nlev), pqte(nbdim,nlev)
!
! !OUTPUT VARIABLES:
!           
    real(dp), intent(inout)  :: Sul_vv(NZ,n_aer_bin)
    !real(dp), intent(out)    :: aero_vt_d(NZ)
! 
! !REVISION HISTORY:
!  05 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(dp) :: vtte(NZ,n_aer_bin)
    integer  :: NZ1, jl, jk
    
    integer  :: ii,jj !counters !eth_af_dryS
    
    real(dp) ::  tc(NZ), vis(NZ), &
                 pq(NZ), fpair(NZ), dz2(NZ), &
                 vt(NZ,n_aer_bin), &
                 yt(NZ,n_aer_bin),yt1(NZ,n_aer_bin)
    logical :: ismaxmin(NZ)
    real(dp) :: v12, c12, a1, a2 
    
    real(dp) :: aer_vol_wet(NZ), aer_r_wet(NZ) !eth_af_dryS
    real(dp) :: aer_ar_wet(NZ) !eth_af_dryS
    
    integer :: minlev, nlev, k
 
    !=======================================================================
    ! AER_SEDIMENTATION begins here!
    !=======================================================================

    !aero_vt_d(:) = 0.

    ! Viscosity of air
    vis(:)=10._dp*1.8325e-5_dp*(416.16_dp/(T_K(:)+120._dp))*(T_K(:)/296.16_dp)**1.5_dp ! [g/cm s]
    
    ! Free mean path of air [cm]
    fpair(:) = 6.6E-6*(1013./P_hPa(:))*(T_K(:)/293.15) ![cm]
    
    !nlev1=nlev
    !if (nlev .eq. 39) then
    !       nlev1 = nlev-4
    !elseif (nlev .eq. 90) then
    !       nlev1 = nlev-5
    !end if
    ! Lowest level we want to consider
    minlev = 1
    ! Highest level we want to consider
    nlev = NZ
    
    do k=1,n_aer_bin
       aer_vol_wet=aer_mass(k)/aWP_Col(:,k)/.01/aDEN_Col(:,k)*1.E12 !calculate wet volume of aerosol size bins !eth_af_dryS            
       aer_r_wet=(aer_vol_wet*3./4./PI)**(1./3.) !wet radius !eth_af_dryS
       aer_ar_wet=aer_r_wet*1.e-4_dp !wet radius in cm !eth_af_dryS

       ! Calculate fall speed vt in m/s
       vt(:,k)=2._dp/9._dp*(aden_col(:,k))*(aer_ar_wet(:))**2._dp* g*100._dp &
               /vis(:)
       vt(:,k) = vt(:,k) *                                   &
                 (1._dp+ 1.257_dp*fpair(:)/aer_ar_wet(:) +  &
                  0.4_dp*fpair(:)/aer_ar_wet(:)*       &
                  EXP(-1.1_dp/(fpair(:)/aer_ar_wet(:))) )

       ! Force CFL <= 1
       vt(:,k) = MIN( vt(:,k)*0.01 , dz(:)/dt )!vt[cm/s] - [m/s]
       
       ! Number density of aerosol droplets in particles/cm3
       yt(:,k) = Sul_vv(:,k)*dens(:)/aer_molec(k) ![mol/mol]->[particles/cm3]
       !call localminmax(yt(1:nproma,:,k),yt1(1:nproma,:,k),nproma,nlev)
      
       ! Calculate intermediate value - # of particles/cm3 * m/s
       vtte(nlev,k) = yt(nlev,k)*vt(nlev,k)
       vtte(1,k) = 0.

       ! Unchanged from original but should (?) be symmetrical
       ismaxmin(1)=.false.
       ismaxmin(nlev)=.false.
       ismaxmin(2:nlev-1) = yt(2:nlev-1,k)>=max(yt(1:nlev-2,k),yt(3:nlev,k)) .or. &
                            yt(2:nlev-1,k)<=min(yt(1:nlev-2,k),yt(3:nlev,k))
       do jk=2,nlev-1
           v12 = vt(jk+1,k)*dz(jk+1)/sum(dz(jk:jk+1)) + vt(jk,k)*dz(jk)/sum(dz(jk:jk+1))
           c12 = min(1., v12 * dt/dz(jk))
           vtte(jk,k) = yt(jk,k)+ (yt(jk-1,k)-yt(jk+1,k))*(1-c12)/4.
           a1=1.75-0.45*c12
           a2=max(1.5, 1.2+0.6*c12)
           if (ismaxmin(jk+1)) vtte(jk,k) = yt(jk,k)+ (yt(jk-1,k)-yt(jk+1,k))*(1-c12)*a1/4.
           if (ismaxmin(jk-1)) vtte(jk,k) = yt(jk,k)+ (yt(jk-1,k)-yt(jk+1,k))*(1-c12)*a2/4.
           vtte(jk,k)=  MIN( MAX( vtte(jk,k), MIN(yt(jk,k),yt(jk-1,k)) ), MAX(yt(jk,k),yt(jk-1,k)) )
           vtte(jk,k) = vtte(jk,k)*v12 ! [kg/kg]*[m/s]
       enddo

       vtte(:,k) = min(yt(:,k)*dz/dt,vtte(:,k))
       yt(:,k) = max(0., yt(:,k) - vtte(:,k)*dt/dz(:))

       ! Diagnostic: Fall speed?
       !aero_vt_d(minlev:nlev) = aero_vt_d(minlev:nlev) + vtte(minlev:nlev,k)*dt/dz(:)

       ! Calculate new concentrations and clip to be >= 0
       yt(1:nlev-1,k) = max(0., yt(1:nlev-1,k) + vtte(2:nlev,k)*dt/dz(1:nlev-1)) !eth_af_dryS 
                               
       Sul_vv(1:nlev,k) = yt(1:nlev,k)*aer_molec(k)/dens(1:nlev) ![particles/cm3]->[mol/mol]
    end do

  END SUBROUTINE AER_sedimentation
#endif
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: localminmax
!
! !DESCRIPTION: Subroutine localminmax finds local minima and maxima. Not 
! currently used.
!\\
!\\
! !INTERFACE:
!
  subroutine localminmax(aa,aa1,nx,ny)
!
! !USES:
!
! -- Nothing --
!
!
! !INPUT VARIABLES:
!           
       integer, intent(in) :: nx,ny
       real(dp),intent(in) :: aa(nx,ny)
!
! !OUTPUT VARIABLES:
!           
       real(dp),intent(out):: aa1(nx,ny)
! 
! !REVISION HISTORY:
!  05 Dec 2018 - S.D.Eastham - Adapted from SOCOL version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!       
       integer :: ix,iy

    !=======================================================================
    ! LOCALMINMAX begins here!
    !=======================================================================

       do iy=2,ny-1
          do ix=1,nx
             if ( (aa(ix,iy)>aa(ix,iy-1).and.aa(ix,iy)>aa(ix,iy+1)) .or. & 
                  (aa(ix,iy)<aa(ix,iy-1).and.aa(ix,iy)<aa(ix,iy+1))      ) then
                aa1(ix,iy)=aa(ix,iy)
             else
                aa1(ix,iy)=aa(ix,iy+1)
             endif
          enddo
       enddo
       aa1(1:nx,1) = aa(1:nx,1)
       aa1(1:nx,ny)= aa(1:nx,ny)

  end subroutine localminmax

end module sect_aer_mod
