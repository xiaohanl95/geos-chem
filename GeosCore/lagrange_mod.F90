!-------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!--------------------------------------------------------------------------

! New version
! use linked link to store injected plumes' character(concentration, location,
! etc.)
! this version could successfully run !!!


! 1. dissolve plume once hitting tropopause

! 2. add variable 'N_split' to define the number of small parts after 
!    plume split

! 3. Oct 9th, 2020 
!  Update the horizontal stretch code, originally the horizontal plume segment is
!  calculated based on adjacent segments, but this requires all the plume
!  segments are exit in the whole simulation without dissolving. The new code
!  will calculate the horizontal stretch based on the center and one endpoint 
!  (box_lat_edge, box_lon_edge) of the segment, which will not rely on the
!   adjacent segment any more.

! 4. Oct 9th, 2020
! Assume that the box_length would not change because of PV=nRT, 
! only the plume cross-section will be used to represent the volume change because
! of PV=nRT

! 5. Oct 10th, 2020
! Add Third judge in 2D cross-section model.
! At the end of simulation, not only the 1D slab model will dissolve, 2D model
! will also dissolve.

! 6. Oct 10th, 2020
! Inn stead of using adjacent plume segment, the box_alpha is calculated based
! on the center and one endpoint inside one plume segment

! 7. Oct 12nd, 2020
! delete the output for grid area (State_Met_AREA_M2.txt), since this variable
! can be found in the formal output as "AREA".

! 8. Oct 12nd, 2020
! correct the runge-kutta method for advection process (dx/dt=u).
! since the wind speed is 3-h average value, assume wind speed is independent of
! time, only change with locations [u(z,y,x) instead of u(t,z,y,x)]
! RK_Dlev(Ki) = RK_Dt(1) * curr_omeg / 100.0 

! 9. Oct 12nd, 2020
! use runge-kutta method when calcualte box_u(i_box), box_v(i_box),
! box_omeg(i_box)

! 10. Oct 12nd, 2020
! modify the box_alpha and ensure the wind speed along plume cross-section is
! correct. range of box_alpha is [PI, -PI)

! 11. Oct 12nd, 2020
! modify the 1D slab model splitting, correct the box_lon and  box_lon_edge 
! for the new plume after splitting
! to get the correct new box_lon, box_length should times COS(box_alpha) first.

! 12. Oct 13rd, 2020
! correct the calculate for plume splitting box_lon and box_lat, as well as
! box_lon_edge and box_lat_edge.

! Oct 15th, 2020
! add forth judge to 2D cross-section model, let plume dissolve when transport
! from stratosphere into troposphere.

! Nov 1st, 2020
! Instead of setting Pdt=10s, adjust Pdt based on the CFL condition in 2D grids.
!
! Dt2 for 1D grid can only be 60s as minimum. If a less than 60s values is
! needed for Dt2, let 1D slab begin to combine 2 grid cells into 1 grid cell.

! Nov 2, 2020
! Correct the injected aerosol concentration at the beginning
! pay attention to the unit

! Nov 2, 2020
! Add Xscale and Yscale for 2D grid, which can be used to identify the box_theta
! value, when concentratino don't follow gaussian distribution after involve
! settling speed etc.
! 
! when change 2D to 1D, not use degree to judge, but use length (Xscale, Yscale) 
! containing 95% mass to judge

! Nov 3, 2020
! add a lifetime variable Lifetime(n_box) to record the plume lifetim in plume
! model. When lifetime is larger than 3 hours, begin to consider change 2D to 1D
! change Sigma2D to lifetime

! Nov 3, 2020
! redifine LenB = Pdy(Ibox) in Slab_init() function,
! set n_slab_max = 300, which is close to n_y_max;

! Nov 7, 2020
! 1. DlonN, DlatN are corrected: the circumference of the lat circle is 2*pi*r,
! instead of 2*pi*r**2 !!!
!
! 2. check the shadow plume, to ensure the total mass is conseved 
! Mr + Mb = (Mp - Ms) + Mb = (Mp + D_Mp) - (Ms + D_Ms) + [Mb - ( D_Mp - D_Ms)]
! Mr: total mass of the real plume;
! Mp: total mass of the plume model in the code;
! Ms: total mass of the shadow plume model in the code;
! Mb: total mass of the background grid cell containing plume model
! D_Mp: change mass after one time step in plume model
! D_Ms: change mass after one time step in shadow plume model

! Nov 8, 2020
! Instead of using 2nd order chemical reaction to decide the initail lenght of
! 1D grid, use the lenght containing 95% of total mass.

! Nov 11, 2020
! add following code to check cpu time for each section of code:
!       call cpu_time(start)
!       put code to test here
!       call cpu_time(finish)
!       print '("Time = ",f6.3," seconds.")',finish-start

! Nov 11, 2020
! CFL_2d is added to speed up the code in advection in 2D grid model

! Nov 14, 2020
! the concentration will change after the volume change because of idea gas law
! (PV=nRT)

! Nov 14, 2020
! track the whole lifetime of the first plume segment to check:
! (1) shadow plume model performance;
! (2) evolution of the plume is reasonable;

! Nov 16, 2020
! delete shadow plume model, replace by a exrta mass container
!
! reset the injected aerosol concentration in GCM model in lagrangain_initial
! module, to make sure the background concentration is 0 at beginning.
!
! For changing from 2D to 1D in function Interplt_2D() in module Slab_init(), 
! when the interpolate point is out of 2D grid, instead of setting as 0, set 
! as background concentration.
!       IF(Ix0<=1 .or. Ix0>=n_x_max-1)THEN
!        Interplt_2D = 0.0e+0_fp

! Nov 17, 2020
! Add two processes:
! When combine 9 to 1 grid cell in 2D grid model, extra mass will increase.
! When combine 2 to 1 grid cell in 1D grid model, extra mass will increase.

! Nov 19, 2020
! assume plume segment lenght is always along the background wind direction, so
! the box_alpha is calculated based on the background wind speed now.
!
! use Lyapunov component (Ly) to calculate the horizontal stretch of the plume
! segment.
!
! delete the box_lon_edge(i_box), box_lat_edge(i_box), box_lev_edge(i_box)

! Nov 21, 2020
! modify the next_i_lon and next_i_lat to ensure they are in the range of 
! [1, IIPAR/JJPAR].
!
!        if(next_i_lon>IIPAR) next_i_lon=next_i_lon-IIPAR
!        if(next_i_lon<1) next_i_lon=next_i_lon+IIPAR
!
!        if(next_i_lat>JJPAR) next_i_lat=JJPAR
!        if(next_i_lat<1) next_i_lat=1

! Nov 27, 2020
!
! update the 1D grid domain, make sure the number of 1D grid cell cam be divided
! by 4, to ensure n_slab_25 is an integer.
!
! add a new variable n_slab_max2=n_slab_max, Concnt1D_bdy(n_slab_max2), to add
! two more grid cell temporarily when calculate the diffusion process inside 1D
! grid.
!       Concnt1D_bdy(1)          = backgrd_concnt
!       Concnt1D_bdy(n_slab_max2) = backgrd_concnt
!       Concnt1D_bdy(2:n_slab_max2-1) = box_concnt_1D(i_box,1:n_slab_max,1)

! Nov 29, 2020
! find the error about total injected mass increasing in plume model:
! box_Ra and box_Rb are commented previously. So only box_length keeps
! increasing. Volume is not conserved.
!      length0           = box_length(i_box)
!      box_length(i_box) = EXP(Ly*Dt) * box_length(i_box)
!      IF(Judge_plume(i_box)==1) THEN
!        box_Ra(i_box) = box_Ra(i_box)*SQRT(length0/box_length(i_box))
!        box_Rb(i_box) = box_Rb(i_box)*SQRT(length0/box_length(i_box))
!
! ERROR: total mass is not conservsed between Plume model and Eulerian model
! Check the third judgement:
! ERROR solved by 
! changing  "IF ( ITS_TIME_FOR_EXIT() ) THEN" 
! to        "IF(DAY==31 .and. HOUR==23)THEN"

! Dec 5, 2020
! add bondary grid cell to 2D grid, similar like what we did for 1D grid.
! By doing so, the boundary of plume would not belong to plume mass, when 
! combine plume grid cells (9 to 1, 2 to 1), these extra mass in plume 
! boundary grid cell will be included.
!
! the dimension order of box_concnt_2D and box_concnt_1D are changed, "n_box"
! dimension are put in the last dimension.

! Dec 6, 2020
! add box_alpha(n_box), and delete box_u, box_v, box_w
!
! delete box_Ptemp(i_box)

! Dec 21, 2020
! ERROR: number of 2D plume model become super larger at 25th day
! Guess the 2D plume model could not change to 1D plume model for some reason,
! check the XY scale and time.
! add a new criteria, once the 2D plume model live for more than 6 hours, force
! the 2D plume model change to 1D plume model.
!
! add curr_Ptemp interplating in plume_module


! Jan 8, 2021 (BIG CHANGE)
! add linked list to store plume information (box_lat, box_concnt_2D, etc)

! Jan 8, 2021
! make the interaction between plume and background as one-way feedback, the
! the concentration in plume boundary is always 0.0, and the injected aerosol
! only transport from plume into background.


! Feb 5, 2021
! introduce more than 1 species
! first use for the fake 2nd order chemical reaction

! Feb 15, 2021
! add parallel computing (OMP) to the loop in code

! Feb 21, 2021
! add a plume_inject() module to decide whether use plume model or directly
! dissolve the injected plume into Eulerain model grid at the beginning.

! Mar 15, 2021
! add sink for injected aerosols
! set all the chemical species in troposphere as 0

! Mar 27, 2021
! let plume injection only for the 60 degree by adding:
! Num_inject<N_total

! Mar 30, 2021
! (1) If there is no plume left, stop using Lagrangian plume module
!       IF(Num_Plume1d+Num_Plume2d==0) use_lagrange=99
!
! (2) 
!       IF(N_stop_inject<0.or.Num_inject<N_stop_inject)THEN
!       N_stop_inject = -1, the plume will be injected in the whole simulation
!       N_stop_inject = N_total, the plume will be injected in 60 degrees
!
! (3) Make sure the last plume could be deleted successfully
!           IF(.NOT.ASSOCIATED(Plume2d%next) .and. &
!                          .NOT.ASSOCIATED(Plume2d_prev))THEN
!
!             WRITE(6,*)'         ***         '
!             WRITE(6,*)'Only left the last Plume2d', Num_Plume2d
!             WRITE(6,*)'         ***         '
!             DEALLOCATE(Plume2d_head)
!
!             i_box = i_box-1
!             Num_Plume2d = Num_Plume2d - 1
!
!             Stop_loop = 1
!             GOTO 901
!
!           ENDIF

! Apr 4, 2021
! decide when change 2D grid to 1D grid:
!          IF( Xscale/Yscale>25.0 ...


! Apr 15, 2021
! (1) debug: Num_inject<N_stop_inject only used for inject new plumes
!      IF(N_stop_inject<0.or.Num_inject<N_stop_inject)THEN
!        ! N_stop_inject<0: keep injecting plume in the whole simulation;
!        ! Num_inject<N_stop_inject: only inject N_stop_inject plumes at
!        !                           beginning;
!
!      DO i_box = Num_inject+1, Num_inject+N_parcel, 1
!
! (2) add a new variable [TROPP_sink] to turn on/off the tropopause sink for all
! species.
!       IF(TROPP_sink = 1)THEN => turn on tropopause sink
!       IF(TROPP_sink = 0)THEN => turn off tropopause sink
!
! (3) add new variable for initial injection location:
!  real, parameter       :: Inject_lon = -141.0e+0_fp
!  real, parameter       :: Inject_hPa = 52.0e+0_fp

! May 15, 2021
! Output the detailed plume concentration distribution in the plume model.
!  IF(mod(tt*Dt,24*60*60)==0)THEN   ! output once every day (24 hours)
!
!    FileConcnt   = 'Lagrange_Concnt_Volume_' // TRIM(ADJUSTL(YEAR_C)) // '-' &
!      //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) // '-' &
!      // TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) &
!      // ':' // TRIM(ADJUSTL(SECOND_C)) // '.txt'
!
!    OPEN( 485,      FILE=TRIM( FileConcnt   ), STATUS='REPLACE',  &
!        FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!    CLOSE(485)
!
!  ENDIF ! IF(mod(tt,1440)==0)THEN

! May 17, 2021
!
! (1)
! change: 
!    IF(mod(tt,144)==0)THEN   ! output once every day (24 hours)
! to:
!    IF(mod(tt*NINT(Dt),24*3600)==0)THEN   ! output once every day (24 hours)
! Once the time step is changed, hourly output still works
!
! (2)
! Move the plume concentration output from lagrange_run to plume_run
! make sure the total mass (Plume model + GOES-Chem) is conserved

! May 25, 2021
! add plume location (i_lon,i_lat,i_lev) to the output txt.

! May 26, 2021
!
! (1)
! change:         if(RK_lat<0.0)then
! to:             if(box_lat<0.0)then
!
! (2)
!    ! second interpolate vertically (Linear)
!    IF(P_mid(init_lev+1)==P_mid(init_lev))THEN
!      WRITE(6,*)"*** WARNING: two same pressure level happens! ***"
!      WRITE(6,*)"init_lev, P_mid(init_lev), P_mid(init_lev+1):"
!      WRITE(6,*)init_lev, P_mid(init_lev), P_mid(init_lev+1)
!      wind_lonlat_lev = wind_lonlat(1)
!    ELSE
!      wind_lonlat_lev = wind_lonlat(1) + (wind_lonlat(2)-wind_lonlat(1)) &
!                                 / (P_mid(init_lev+1)-P_mid(init_lev)) &
!                                     * (curr_pressure-P_mid(init_lev))
!    ENDIF

! May 31, 2021
! Put three product/time/location criterion for transferring in the sam place
!       ! ------------------------------------------------------------------
!       ! 1. Product criterion
!       ! 2. Time criterion
!       ! 3. Location critetion
!       ! ------------------------------------------------------------------
!       IF( ABS(Prod_before-Prod_after)/Prod_after<Dissolve_critiria &
!           .OR. box_life/3600/24>30.0 &
!           .OR. box_lev>State_Met%TROPP(i_lon,i_lat) )THEN

! June 3, 2021
! Count total plume volume for each Eulerian grid cell
!    ! cumulative volume of all plume in a Eulerian grid cell
!    real(fp), dimension(:), allocatable  :: SumV_Plume

! June 8, 2021
! add a new variable PASV_LA4
! which is a temporary variable containing tracer concentration that is
! transferred from the plume grid cell with a low concentration.

! June 9, 2020
! (1)
! Volume criterion: first delete the oldest plume segment
!
! (2) 
! add Entropy calculation in the code

! June 10, 2020
! add new function SortList(), which could sort the linked list based on volume
! criterion.

! June 11, 2020
!
! (1) 
!Put all the Entropy calculation in the final end
! make sure that all the interaction between plume model and GEOS-Chem are
! finished when the entropy calculation begins
!
! (2)
! Debug: use SortList() to sort Plume1d list would change both Plume1d_head and
! Plume1d_tail node, make sure update both instead of just Plume1d_head !
!
!     !====================================================================
!     ! Always make sure the Plume1d_tail refer to the last node of the list!
!     !
!     ! WARNING: the sort list after SortList() will change both Plume1d_head
!     ! and Plume1d_tail, make sure the Plume1d_tail is also updated !!! 
!     !==================================================================== 
!     Plume1d_tail => Plume1d_prev

! June 14, 2021
! 
! (1) delete PASV_LA4
!
! (2) decide to use DX or 2*DX to triger the plume split
!
!       integer, parameter    :: Split_length = 2 ! how many times of DX
!         IF( box_Ra > Split_length * DX*1.0e+5 ) THEN
!
! (3) add variable Type_transfer to record the number of transferred plume
! segments by 4 different transferring criterion
!
!       integer :: Type_transfer
!         IF( ABS(Prod_before-Prod_after)/Prod_after<Dissolve_critiria ) &
!                                                      Type_transfer = 11
!         IF( Plume1d%Is_transfer==1 )                 Type_transfer = 22
!         IF( box_lev>State_Met%TROPP(i_lon,i_lat) )   Type_transfer = 33
!         IF( box_life/3600/24>30.0 )                  Type_transfer = 44

! June 17, 2021
! revise the product criterion, make the criterion less depend on the plume
! segment size
!       ! (1) Only consider the product producted by the tracer in Plume model
!       ! (2) Consider the interaction between Plume model and Eulerian model
!       ! Reference: Karamchandani et al., (2000)
!       !            Korsakissok and Mallet, (2010)
!       backgrd_concnt = State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA)
!
!       Prod_before = SUM( ( box_concnt_1D(1:n_slab_max,1)+backgrd_concnt )**2 &
!                                               - backgrd_concnt**2 ) *V_grid_1D
!
!       Eul_concnt = ( SUM(box_concnt_1D(1:n_slab_max,1)) *V_grid_1D &
!                    + backgrd_concnt*grid_volumn ) / grid_volumn
!
!       Prod_after = (Eul_concnt**2 - backgrd_concnt**2) *grid_volumn


! June 18, 2021
! change initial plume segment length from 1km to more than 10km
! consider the host Eulerian grid resolution (2degree, 200km)

! June 21, 2021
!
! (1)
! Add plume segment length splitting
!       ! -------------------------------------
!       ! 2. Splitting for plume segment length
!       ! -------------------------------------
!
! (2)
! 2.1 put Lyaponov calculation into a function Calc_Ly()
!
! 2.2 update lyaponov calculation by Introducing the location comparison between 
! plume segment and its host Eulerian grid cell
!       box_lon>=X_mid(i_lon)
!       box_lat>=Y_mid(i_lat)
!
! 2.3 revise lyaponov equation:
!      Ly = ABS(u1-u2)/(D_x) or  Ly = ABS(v1-v2)/(D_y)
!
! 2.4 correct lyaponov equation, remove the ABS() functiono, so that Ly could be
! both positive and negative, corresponding to stretch and contract. Previous
! calculation only have stretch, NO contract.


! August 18, 2021
!
! (1)
! delete all the code about the assumed 2nd order process.
!
! (2)
! add i_inject and i_entropy to replace i_tracer
!
! (3)
! set n_species=40 to include all 40 bins of aerosol microphysical bins
!
! (4)
! correct the initial inject, making sure the Eulerian simulation and PiG
! simulation have the same initial injection.


! Aug 28, 2021
! (1)
! Add aerosol module to the code
! (2)
! edd code to excahnge unit between [molec/cm3] and [v/v] for species in plume segments
! (3)
! Use function Ind_() to find the species index in GEOS-Chem.



! double check the P_edge, whether P_edge can change at different (lon, lat)


! reconsider plume dissolve criteria, which should involve extra mass

!
! dissolve the plume when touch the top of the atmosphere
! i_lev >= LLPAR


! check: delete or release all the comment code


! Sep 30, 2021
! turn off all the interactions between plume model and GEOS-Chem, to check the
! aerosol module


MODULE Plume_list

  USE precision_mod

  TYPE :: Plume2d_list
    integer :: IsNew ! 1: the plume is new injected
    integer :: label ! injected rank

    real(fp) :: LON, LAT, LEV
    real(fp) :: LENGTH, ALPHA
    real(fp) :: LIFE

    real(fp) :: PDX, PDY
 
    real(fp), DIMENSION(:),     POINTER :: rWet
    real(fp), DIMENSION(:,:,:), POINTER :: aWP
    real(fp), DIMENSION(:,:,:), POINTER :: aDen
    real(fp), DIMENSION(:,:,:), POINTER :: Sfc_Ten

    real(fp), DIMENSION(:,:,:), POINTER :: CONCNT2d 
    ! [n_x_max, n_y_max, n_species]
    !!! Don't assign dimension value first ???


    TYPE(Plume2d_list), POINTER :: next
  END TYPE


  TYPE :: Plume1d_list
    integer :: label ! injected rank
    integer :: Is_transfer ! if =1, transfer to the host Euletian model


    real(fp) :: LON, LAT, LEV
    real(fp) :: LENGTH, ALPHA
    real(fp) :: LIFE

    real(fp) :: RA, RB, THETA

    real(fp), DIMENSION(:),   POINTER :: rWet
    real(fp), DIMENSION(:,:), POINTER :: aWP
    real(fp), DIMENSION(:,:), POINTER :: aDen
    real(fp), DIMENSION(:,:), POINTER :: Sfc_Ten

    real(fp), DIMENSION(:,:),   POINTER :: CONCNT1d ! [n_slab_max,n_species]


    TYPE(Plume1d_list), POINTER :: next
  END TYPE

END MODULE Plume_list



MODULE Lagrange_Mod

  USE Plume_list

  USE precision_mod
  USE ERROR_MOD
  USE ErrCode_Mod
  USE PhysConstants, ONLY : PI, Re, g0, AIRMW, AVO, BOLTZ
!  USE CMN_SIZE_Mod,    ONLY : IIPAR, JJPAR, LLPAR

  USE TIME_MOD,      ONLY : GET_YEAR
  USE TIME_MOD,      ONLY : GET_MONTH
  USE TIME_MOD,      ONLY : GET_DAY
  USE TIME_MOD,      ONLY : GET_HOUR
  USE TIME_MOD,      ONLY : GET_MINUTE
  USE TIME_MOD,      ONLY : GET_SECOND

  IMPLICIT NONE


  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: lagrange_init
  PUBLIC :: plume_inject
!  PUBLIC :: lagrange_run
!  PUBLIC :: plume_run
!  PUBLIC :: lagrange_write_std
  PUBLIC :: lagrange_cleanup

  ! PUBLIC VARIABLES:

  PUBLIC :: n_x_max, n_y_max
  PUBLIC :: n_x_mid, n_y_mid

  PUBLIC :: n_slab_max, n_slab_25, n_slab_50, n_slab_75

  PUBLIC :: use_lagrange

  integer               :: use_lagrange = 1
  integer               :: TROPP_sink   = 0
  integer               :: Volume_Sort  = 1 
  ! 1 = use SortList() function, transfer largest (not oldest) plume segment for
  ! volume criterion
  integer               :: Calc_entropy = 0 ! 1 = turn on entropy calculation

  integer, parameter    :: n_x_max  = 207   !number of x grids in 2D, should be (9 x odd)
  integer, parameter    :: n_y_max  = 81    !number of y grids in 2D, should be (9 x odd)

  ! the odd number of n_x_max can ensure a center grid
  integer, parameter    :: n_x_mid  = (n_x_max+1)/2 
  integer, parameter    :: n_y_mid  = (n_y_max+1)/2 

  integer, parameter    :: n_x_max2 = n_x_max+2
  integer, parameter    :: n_y_max2 = n_y_max+2

  integer, parameter    :: n_x_mid2 = (n_x_max2+1)/2 
  integer, parameter    :: n_y_mid2 = (n_y_max2+1)/2 

  ! n_slab_max should be divided by 4, to ensure n_slab_25 is an integer.
  integer, parameter    :: n_slab_max  = (INT(n_y_max/4)+1)*4 ! close to n_y_max, number of slabs in 1D
  integer, parameter    :: n_slab_max2 = n_slab_max+2
  ! add 2 more slab grid to containing background concentration

  integer               :: IIPAR, JJPAR, LLPAR

  integer               :: n_slab_25, n_slab_50, n_slab_75
  integer               :: id_PASV_LA3, id_PASV_LA2, id_PASV_LA 
  integer               :: id_PASV_EU2, id_PASV_EU


  real, parameter       :: Dx_init     = 100
  real, parameter       :: Dy_init     = 10
  real, parameter       :: Length_init = 40.0*1000.0 ! [m], 20km

  real, parameter       :: Inject_lon   = -141.0e+0_fp
  real, parameter       :: Inject_hPa   = 50.0e+0_fp
  real, parameter       :: Inject_rate  = 30.0e+0_fp ! [kg/km]
  ! 25.0e+0_fp ! [hPa] at about 25 km
  ! 50.0e+0_fp       ! [hPa] at about 20 km

  ! some parameter for sensitive test
  integer, parameter    :: N1_split            = 5     ! Cross-section splitting
  integer, parameter    :: N2_split            = 5     ! length splitting
  integer, parameter    :: Split_length        = 1     ! how many times of DX
  real, parameter       :: Dissolve_critiria   = 10*0.01
  real, parameter       :: Volume_percent      = 30*0.01
  real, parameter       :: Critical_day        = 28.0         ! [day]

  real(fp), pointer     :: X_mid(:), Y_mid(:), P_mid(:)
  real(fp), pointer     :: P_edge(:)

  real(fp)              :: DX, DY

  real(fp) :: mass_eu, mass_la, mass_la2

  integer               :: N_parcel   ! 131        
  integer               :: Num_inject, Num_Plume2d, Num_Plume1d, Num_dissolve     
  integer               :: tt     
  ! Aircraft would release 131 aerosol parcels every time step

  ! use for plume injection
  integer               :: N_total, N_stop_inject
  real(fp)              :: Length_lat

  ! SO2, SO4, AERSctSul001 to AERSctSul040
  integer, parameter    :: n_species    = 42
  integer, parameter    :: n_bins       = 40
  integer, parameter    :: i_inject     = 1     ! GEOS-Chem index of injected tracer
  integer, parameter    :: i_entropy    = 1     ! GEOS-Chem index of tracer for entropy
!  integer, parameter    :: i_product = 2

  integer               :: id_GC(n_species)     ! Species id in GEOS-Chem

!  real(fp), parameter       :: Kchem = 1.0e-20_fp ! chemical reaction rate
! use Kchem = 1.0e-20_fp for >1 year simulation

!  integer               :: Stop_inject ! 1: stop injecting; 0: keep injecting
!                           ! used for contiuing injecting scenario

!  TYPE :: Plume2d_list
!    integer :: IsNew ! 1: the plume is new injected
!    integer :: label ! injected rank
!
!    real(fp) :: LON, LAT, LEV
!    real(fp) :: LENGTH, ALPHA
!    real(fp) :: LIFE
!
!    real(fp) :: PDX, PDY
!    real(fp), DIMENSION(:,:,:), POINTER :: CONCNT2d ! [n_x_max, n_y_max, n_species]
!    !!! Don't assign dimension value first ???
!
!
!    TYPE(Plume2d_list), POINTER :: next
!  END TYPE
!
!
!  TYPE :: Plume1d_list
!    integer :: label ! injected rank
!    integer :: Is_transfer ! if =1, transfer to the host Euletian model
!
!    real(fp) :: LON, LAT, LEV
!    real(fp) :: LENGTH, ALPHA
!    real(fp) :: LIFE
!
!    real(fp) :: RA, RB, THETA
!    real(fp), DIMENSION(:,:), POINTER :: CONCNT1d ! [n_slab_max,n_species]
!
!
!    TYPE(Plume1d_list), POINTER :: next 
!  END TYPE
!

  TYPE(Plume2d_list), POINTER :: Plume2d_tail, Plume2d_head
  TYPE(Plume1d_list), POINTER :: Plume1d_tail, Plume1d_head


CONTAINS


!-----------------------------------------------------------------

  SUBROUTINE lagrange_init(am_I_root, Input_Opt, State_Chm, State_Grid, State_Met, RC)

    USE Plume_list

    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState

    USE UnitConv_Mod,   ONLY : Convert_Spc_Units

!    USE GC_GRID_MOD,   ONLY : XMID, YMID
!    USE GC_GRID_MOD,   ONLY : GET_AREA_M2 ! new

    USE State_Grid_Mod, ONLY : GrdState

    USE State_Chm_Mod,  ONLY : Ind_

    USE TIME_MOD,       ONLY : GET_TS_DYN

    USE TIME_MOD,       ONLY : GET_YEAR
    USE TIME_MOD,       ONLY : GET_MONTH
    USE TIME_MOD,       ONLY : GET_DAY
    USE TIME_MOD,       ONLY : GET_HOUR
    USE TIME_MOD,       ONLY : GET_MINUTE
    USE TIME_MOD,       ONLY : GET_SECOND



    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
    TYPE(MetState), intent(in)    :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(OptInput), intent(in)    :: Input_Opt
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State objectgg
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    INTEGER                       :: i_box, i_slab
    INTEGER                       :: ii, jj, kk
    CHARACTER(LEN=255)            :: FILENAME, FileEntropy, File996, File580
    CHARACTER(LEN=255)            :: FILENAME2, FILENAME3

    integer                       :: i_lon, i_lat, i_lev            !1:IIPAR
    real(fp)                      :: box_lon, box_lat, box_lev

    REAL(fp)                      :: lon1, lat1, lon2, lat2
    REAL(fp)                      :: box_lon_edge, box_lat_edge        

    real(fp), dimension(:,:,:), allocatable :: box_concnt_2D
    real(fp), dimension(:,:), allocatable   :: box_concnt_1D

    TYPE(Plume2d_list), POINTER :: Plume2d_new, Plume2d, Plume2d_prev
    TYPE(Plume1d_list), POINTER :: Plume1d_new, Plume1d, Plume1d_prev

    real(fp), pointer   :: X_edge(:), Y_edge(:)
    real(fp)            :: X_edge2, Y_edge2

    real(fp), pointer   :: PASV_EU
!    REAL(fp)            :: MW_g
    REAL(fp)            :: Dt

    ! Strings
    CHARACTER(LEN=255)	:: MSG, LOC
    CHARACTER(LEN=63)  	:: OrigUnit
    CHARACTER(LEN=255) 	:: ErrMsg


    Dt = GET_TS_DYN()
    N_parcel = NINT(132 *1000/Length_init /600*Dt)


    ! use for plume injection
    N_total    = 60 * 1.0e+5 / Length_init
    Length_lat = Length_init / 1.0e+5

    ! define how many plume segments will be injected
    ! -1 means keep inecting in the whole simulation
    N_stop_inject = 1 ! N_total ! -1

!    Stop_inject = 0


    ! Check 2D domain grid
    IF(MOD(n_x_max,9).ne.0) WRITE(6,*)"*** ERROR ***"
    IF(MOD(n_y_max,9).ne.0) WRITE(6,*)"*** ERROR ***"


    IIPAR = State_Grid%NX
    JJPAR = State_Grid%NY
    LLPAR = State_Grid%NZ

    DX = State_Grid%DX
    DY = State_Grid%DY


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



    File580   = 'Plume_total_concentration.txt'

    OPEN( 580,      FILE=TRIM( File580  ), STATUS='REPLACE',  &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

    CLOSE(580)



    allocate(box_concnt_2D(n_x_max, n_y_max, n_species))
    allocate(box_concnt_1D(n_slab_max,n_species))


    WRITE(6,'(a)') '--------------------------------------------------------'
    WRITE(6,'(a)') '--------------------------------------------------------'
    WRITE(6,'(a)') ' Initial Lagrnage Module (Using Dynamic time step)'
    WRITE(6,'(a)') '--------------------------------------------------------'
    WRITE(6,*)     'time step=', Dt
    WRITE(6,*)     'Injected plume every time step: ', N_parcel
    WRITE(6,'(a)') '--------------------------------------------------------'
    WRITE(6,'(a)') '--------------------------------------------------------'




    n_slab_25 = (n_slab_max-2)/4*1
    n_slab_50 = (n_slab_max-2)/4*2
    n_slab_75 = (n_slab_max-2)/4*3


    X_mid  => State_Grid%XMid(:,1) ! Grid box longitude [degrees] ! XMID(:,1,1)   ! IIPAR ! new
    Y_mid  => State_Grid%YMid(1,:) ! Grid box latitude center [degree] ! YMID(1,:,1)
    P_mid  => State_Met%PMID(1,1,:)  ! Pressure at level centers (hPa)

    P_edge => State_Met%PEDGE(1,1,:)  ! Wet air press @ level edges [hPa]

    X_edge => State_Grid%XEdge(:,1) !XEDGE(:,1,1)   ! IIPAR+1 ! new
    Y_edge => State_Grid%YEdge(1,:) !YEDGE(1,:,1)  
    ! Use second YEDGE, because sometimes YMID(2)-YMID(1) is not DLAT

    X_edge2       = X_edge(2)
    Y_edge2       = Y_edge(2)



    tt   = 0


      ! Find the 42 species id in GEOS-Chem
      id_GC(1) = Ind_('SO2')    !State_Chm%SpcData(N)%Info%Name
      id_GC(2) = Ind_('SO4')

      DO ii=3,42,1 ! From 'AERSctSul001' to 'AERSctSul040'
        id_GC(ii) = Ind_('AERSctSul001') + (ii-3)
      ENDDO

      WRITE(6,*)'shw, name:', id_GC(42), State_Chm%SpcData(id_GC(42))%Info%ModelID, &
                                         State_Chm%SpcData(id_GC(42))%Info%Name


      ! check Plume model select the correct species
    IF(TRIM(State_Chm%SpcData(id_GC(42))%Info%Name) /= 'AERSctSul040') THEN

        WRITE(6,*) ' '
        WRITE(6,*) 'The species used in Plume model include:'
        WRITE(6,*) State_Chm%SpcData(id_GC(1))%Info%Name, &
                   State_Chm%SpcData(id_GC(2))%Info%Name, &
                   State_Chm%SpcData(id_GC(3))%Info%Name, &
                   '...', &
                   State_Chm%SpcData(id_GC(42))%Info%Name
        WRITE(6,*) ' '

        MSG = 'Incorrect species used in Plume model'
        LOC = 'SUBROUTINE lagrange_init() in lagrange.F90'
        CALL GC_Error( MSG, RC, LOC )

    ENDIF




    DO ii=1,42,1
      State_Chm%Species(:,:,:,id_GC(ii)) = 0.0        ! !!! shw
    ENDDO



! !!! shw
!    IF(use_lagrange==0)THEN
    ! instantly dissolve injected plume into Eulerian grid 

      WRITE(6,'(a)') ' '
      WRITE(6,'(a)') '********************************************************'
      WRITE(6,'(a)') ' You are not using lagrange_mod now'
      WRITE(6,'(a)') ' set variable use_lagrange = 1 to turn on lagrnage_mod  '
      WRITE(6,'(a)') '********************************************************'
      WRITE(6,'(a)') ' '


      WRITE(6,'(a)') ' '
      WRITE(6,'(a)') '********************************************************'
      WRITE(6,'(a)') ' WARNING: all initial SULF concentrations are set as 0  '
      WRITE(6,'(a)') '********************************************************'
      WRITE(6,'(a)') ' '



      box_lon    = Inject_lon
      box_lat    = -29.95e+0_fp
      box_lev    = Inject_hPa


      ! check this grid box is the neast one or the one located in left-bottom
      ! ???
      i_lon = Find_iLonLat(box_lon, DX, X_edge2)
      if(i_lon>IIPAR) i_lon=i_lon-IIPAR
      if(i_lon<1) i_lon=i_lon+IIPAR

      i_lat = Find_iLonLat(box_lat, DY, Y_edge2)
      if(i_lat>JJPAR) i_lat=JJPAR
      if(i_lat<1) i_lat=1

      i_lev = Find_iPLev(box_lev,P_edge)
      if(i_lev>LLPAR) i_lev=LLPAR

      IF(i_lev==LLPAR) WRITE(6,*) 'box_lev:', box_lev, i_lev
      IF(ABS(box_lat)>40) WRITE(6,*)'*** ERROR in box_lat:', i_box, box_lat

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
      PASV_EU => State_Chm%Species(i_lon,i_lat,i_lev,id_GC(i_inject))  ! [kg/kg]

      ! MW_g = State_Chm%SpcData(nAdv)%Info%MW_g
      ! Here assume the injection rate is 30 kg/km for H2SO4: 
      PASV_EU = PASV_EU + Length_init*1.0e-3_fp*Inject_rate &
                                  /State_Met%AD(i_lon,i_lat,i_lev)

      WRITE(6,*)'shw: ', i_lon, i_lat, i_Lev

! !!! shw
!    ELSE

      WRITE(6,'(a)') ' '
      WRITE(6,'(a)') '********************************************************'
      WRITE(6,'(a)') ' You are using lagrange_mod now'
      WRITE(6,'(a)') ' set variable use_lagrange = 0 to turn off lagrnage_mod '
      WRITE(6,'(a)') ' '
      WRITE(6,'(a)') ' WARNING: all the initial concentration is set as 0 now'
      WRITE(6,'(a)') '********************************************************'
      WRITE(6,'(a)') ' '



    ! -----------------------------------------------------------
    ! create first node (head) for 2d linked list:
    ! 
    ! At the beginning, there are only one node in the 2d list, 
    ! so Plume2d_tail and Plume2d_head refer to the same node
    ! -----------------------------------------------------------
    ALLOCATE(Plume2d_tail)
    Plume2d_tail%IsNew = 1
    Plume2d_tail%label = 1

    Plume2d_tail%LON = Inject_lon
    Plume2d_tail%LAT = -29.95e+0_fp
    Plume2d_tail%LEV = Inject_hPa

    Plume2d_tail%LENGTH = Length_init ! 1000m 
    Plume2d_tail%ALPHA  = 0.0e+0_fp

    Plume2d_tail%LIFE = 0.0e+0_fp

    Plume2d_tail%PDX = Dx_init
    Plume2d_tail%PDY = Dy_init


    ALLOCATE(Plume2d_tail%rWet(n_bins))
    Plume2d_tail%rWet  = 0.0e+0_fp

    ALLOCATE(Plume2d_tail%aWP(n_x_max, n_y_max, n_bins))
    Plume2d_tail%aWP  = 0.0e+0_fp

    ALLOCATE(Plume2d_tail%aDen(n_x_max, n_y_max, n_bins))
    Plume2d_tail%aDen = 0.0e+0_fp

    ALLOCATE(Plume2d_tail%Sfc_Ten(n_x_max, n_y_max, n_bins))
    Plume2d_tail%Sfc_Ten = 0.0e+0_fp


    ! Here assume the injection rate is 30 kg/km (=30 g/m) for H2SO4
    ALLOCATE(Plume2d_tail%CONCNT2d(n_x_max, n_y_max, n_species))
    Plume2d_tail%CONCNT2d = 0.0e+0_fp
    Plume2d_tail%CONCNT2d(n_x_mid,n_y_mid,i_inject) = &
		  Plume2d_tail%LENGTH*Inject_rate &
                  /(Plume2d_tail%PDX*Plume2d_tail%PDY*Plume2d_tail%LENGTH)
    ! From [g/m3] to [molec/cm3], 98.0 g/mol for H2SO4
    Plume2d_tail%CONCNT2d(n_x_mid,n_y_mid,i_inject) = &
                              Plume2d_tail%CONCNT2d(n_x_mid,n_y_mid,i_inject) &
                                                      / 1.0e+6_fp / 98.0 * AVO


    mass_la  = Plume2d_tail%CONCNT2d(n_x_mid,n_y_mid,i_inject)* &
             (Plume2d_tail%PDX*Plume2d_tail%PDY*Plume2d_tail%LENGTH*1.0e+6_fp)
    mass_la2 = 0.0
    mass_eu  = 0.0

    NULLIFY(Plume2d_tail%next)
    Plume2d_head => Plume2d_tail


    Num_Plume2d = 1

    ! Check Num_Plume1d to determine the first 1d (changed from from 2d), and
    ! create the first node for 1d
    Num_Plume1d = 0

    Num_dissolve = 0

! !!! shw
!    ENDIF ! IF(use_lagrange==0)THEN



! !!! shw
      !======================================================================
      ! Convert species to [molec/cm3] (ewl, 8/16/16)
      !======================================================================
      CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                              'molec/cm3', RC, OrigUnit=OrigUnit )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Unit conversion error!'
         CALL GC_Error( ErrMsg, RC, 'plume_run in lagrange_mod.F90')
         RETURN
      ENDIF


	State_Chm%Species(i_lon,i_lat,i_lev,id_GC(i_inject)) = &
				Plume2d_tail%CONCNT2d(n_x_mid,n_y_mid,i_inject)

      
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
! !!! shw



    WRITE(6,*)'*** Finish lagrange_init() ***'

    ! use this value to set the initial latutude for injected plume
    Num_inject = 1


    deallocate(box_concnt_2D)
    deallocate(box_concnt_1D)

    nullify(Plume2d_new)
    nullify(Plume2d)
    nullify(Plume2d_prev)
    nullify(Plume1d_new)
    nullify(Plume1d)
    nullify(Plume1d_prev)


  END SUBROUTINE lagrange_init


!=================================================================

!=================================================================
  SUBROUTINE plume_inject(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)

    USE Plume_list

    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Met_Mod,   ONLY : MetState
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Grid_Mod,  ONLY : GrdState

    USE TIME_MOD,        ONLY : GET_TS_DYN
    USE UnitConv_Mod,    ONLY : Convert_Spc_Units


    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
    TYPE(MetState), intent(in)    :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(OptInput), intent(in)    :: Input_Opt
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State objectgg
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure


    real(fp), pointer :: X_edge(:), Y_edge(:)
    real(fp)          :: X_edge2, Y_edge2


    integer  :: i_box

    integer  :: i_lon, i_lat, i_lev
    real(fp) :: box_lon, box_lat, box_lev

    real(fp), pointer  :: PASV_EU
!    integer            :: nAdv
!    REAL(fp)           :: MW_g

    integer  :: i_advect1

    REAL(fp) :: Dt

    real(fp) :: Entropy2_Concnt, Entropy2_V, Entropy2
    real(fp) :: tracer2_mol, air2_mol, mix2_ratio

    CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: ErrMsg
    CHARACTER(LEN=255)     :: FileEntropy

    ErrMsg = ''



    !=======================================================================
    !
    ! instantly dissolve injected plume into Eulerian grid 
    ! NOT use plume model
    !
    !=======================================================================
    IF(use_lagrange==0)THEN 

      Dt = GET_TS_DYN()

      X_edge => State_Grid%XEdge(:,1) !XEDGE(:,1,1)   ! IIPAR+1 ! new
      Y_edge => State_Grid%YEdge(1,:) !YEDGE(1,:,1)  
      ! Use second YEDGE, because sometimes YMID(2)-YMID(1) is not DLAT

      X_edge2       = X_edge(2)
      Y_edge2       = Y_edge(2) 

      ! -----------------------------------------------------------
      ! add new box every time step
      ! -----------------------------------------------------------
      IF(N_stop_inject<0.or.Num_inject<N_stop_inject)THEN
        ! N_stop_inject<0: keep injecting plume in the whole simulation;
        ! Num_inject<N_stop_inject: only inject N_stop_inject plumes at
        !                           beginning;

      DO i_box = Num_inject+1, Num_inject+N_parcel, 1

        box_lon    = Inject_lon
        box_lat    = ( -30.005e+0_fp + Length_lat * MOD(i_box,N_total) ) &
                     * (-1.0)**FLOOR(1.0*i_box/N_total) ! -29.995S:29.995N:0.01
        box_lev    = Inject_hPa


        ! check this grid box is the neast one or the one located in left-bottom
        ! ???
        i_lon = Find_iLonLat(box_lon, DX, X_edge2)
        if(i_lon>IIPAR) i_lon=i_lon-IIPAR
        if(i_lon<1) i_lon=i_lon+IIPAR

        i_lat = Find_iLonLat(box_lat, DY, Y_edge2)
        if(i_lat>JJPAR) i_lat=JJPAR
        if(i_lat<1) i_lat=1

        i_lev = Find_iPLev(box_lev,P_edge)
        if(i_lev>LLPAR) i_lev=LLPAR

        IF(i_lev==LLPAR) WRITE(6,*) 'box_lev:', box_lev, i_lev
        IF(ABS(box_lat)>40) WRITE(6,*)'*** ERROR in box_lat:', i_box, box_lat

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
!        PASV_EU => State_Chm%Species(i_lon,i_lat,i_lev,id_GC(2))  ! [kg/kg]
!
!        ! MW_g = State_Chm%SpcData(nAdv)%Info%MW_g
!        ! Here assume the injection rate is 30 kg/km for H2SO4: 
!        PASV_EU = PASV_EU + Length_init*1.0e-3_fp*30.0 &
!                                    /State_Met%AD(i_lon,i_lat,i_lev)

      ENDDO ! i_box = Num_inject+1, Num_inject+N_parcel, 1

      ! use this value to set the initial latutude for injected plume
      Num_inject = Num_inject + N_parcel

      ENDIF ! IF(N_stop_inject<0.or.Num_inject<N_stop_inject)THEN


      !======================================================================
      ! Convert species to [molec/cm3] (ewl, 8/16/16)
      !======================================================================
      CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                              'molec/cm3', RC, OrigUnit=OrigUnit )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Unit conversion error!'
         CALL GC_Error( ErrMsg, RC, 'plume_run in lagrange_mod.F90')
         RETURN
      ENDIF


      !======================================================================
      ! run the fake 2nd order chemical reaction
      !======================================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( i_lev, i_lat, i_lon )
      DO i_lev = 1, LLPAR
      DO i_lat = 1, JJPAR
      DO i_lon = 1, IIPAR


        !!! shw: for a test
        IF(i_lev/=39 .or. i_lat/=16 .or. i_lon/=9)THEN
          State_Chm%Species(i_lon,i_lat,i_lev,:) = 0.0
        ENDIF


!        ! from EU to EU2
!        State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_EU2) = &
!              State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_EU2) &
!            + Dt* Kchem*State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_EU)**2


        !  For all the grid cells in the troposphere, let concentration to be zero
        IF(TROPP_sink == 1)THEN
        IF ( State_Met%PEDGE(i_lon,i_lat,i_lev) > State_Met%TROPP(i_lon,i_lat) ) THEN
          State_Chm%Species(i_lon,i_lat,i_lev,:) = 0.0    
      ENDIF
      ENDIF


      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO



    ! Calculate Entropy for Eulerian model
    IF(mod(tt*NINT(Dt),24*60*60)==0 .and. Calc_entropy==1)THEN   ! output once every day (24 hours)

      DO i_lon = 1, IIPAR
      DO i_lat = 1, JJPAR
      DO i_lev = 1, LLPAR

        Entropy2_Concnt = State_Chm%Species(i_lon, i_lat, i_lev, i_entropy)
        Entropy2_V = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp

        tracer2_mol = Entropy2_Concnt*Entropy2_V/AVO
        air2_mol    = Entropy2_V/1e+6_fp*State_Met%AIRDEN(i_lon,i_lat,i_lev)*1000/AIRMW
        mix2_ratio  = tracer2_mol/air2_mol

        IF(mix2_ratio>0) Entropy2 = Entropy2 - BOLTZ*tracer2_mol*log(mix2_ratio)

      ENDDO
      ENDDO
      ENDDO

!     WRITE(6,*) "Entropy:", Entropy2
!     WRITE(6,*) "Total Mass:", SUM( State_Chm%Species(:,:,:,id_PASV_EU) &
!                                           * State_Met%AIRVOL(:,:,:)*1e+6_fp )


      FileEntropy   = 'Plume_entropy.txt'
      OPEN( 487,      FILE=TRIM( FileEntropy   ), STATUS='OLD',  &
           POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

      WRITE(*,*) "Entropy:", Entropy2
      WRITE(487,*) Entropy2

      CLOSE(487)


    ENDIF ! IF(mod(tt*NINT(Dt),24*60*60)==0 .and. Calc_entropy==1)THEN


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

    ENDIF ! IF(use_lagrange==0)THEN




    !=======================================================================
    !
    ! use plume model
    !
    !=======================================================================
    IF(use_lagrange==1)THEN

    ! call the lagrnage_run() and plume_run() to calculate injected plume

!      State_Chm%Species(:,:,:,id_PASV_LA3) = 0.0e+0_fp  ! [kg/kg]

      CALL lagrange_run(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
      CALL plume_run(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
      CALL Plume_Aerosol(am_I_Root, State_Chm, State_Grid, &
                                        State_Met, Input_Opt, RC)
      CALL lagrange_write_std( am_I_Root, RC )

    ENDIF

    tt = tt + 1


  END SUBROUTINE plume_inject

!=================================================================

  SUBROUTINE lagrange_run(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)

    USE Plume_list

    USE Input_Opt_Mod, ONLY : OptInput

    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState

    USE TIME_MOD,      ONLY : GET_TS_DYN

!    USE GC_GRID_MOD,   ONLY : XEDGE, YEDGE
              
!    USE CMN_SIZE_Mod,  ONLY : DLAT, DLON !new
    USE State_Grid_Mod,  ONLY : GrdState

    USE UnitConv_Mod,  ONLY : Convert_Spc_Units


    logical, intent(in)           :: am_I_Root
    TYPE(MetState), intent(in)    :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State objectgg
    TYPE(OptInput), intent(in)    :: Input_Opt
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    REAL(fp) :: Dt, RK_Dt(5)                  ! = 600.0e+0_fp          

    real(fp), pointer  :: PASV_EU           
!    integer            :: nAdv        
!    REAL(fp)           :: MW_g

    integer :: i_box

    ! 1:IIPAR
    integer :: i_lon, i_lat, i_lev
    integer :: next_i_lon, next_i_lat, next_i_lev

    integer :: ii, jj, kk

    integer :: Ki

    integer :: i_species
    integer :: i_x, i_y
    integer :: i_slab


    real(fp) :: curr_lon, curr_lat, curr_pressure
    real(fp) :: RK_lon, RK_lat
    real(fp) :: X_edge2, Y_edge2

    real(fp) :: dbox_lon, dbox_lat, dbox_lev
    real(fp) :: dbox_x_PS, dbox_y_PS
    real(fp) :: box_x_PS, box_y_PS

    real(fp) :: RK_Dlon(4), RK_Dlat(4), RK_Dlev(4)
    real(fp) :: RK_x_PS, RK_y_PS
    real(fp) :: RK_Dx_PS, RK_Dy_PS
    real(fp) :: RK_u(4), RK_v(4), RK_omeg(4)

    real(fp), pointer :: u(:,:,:)
    real(fp), pointer :: v(:,:,:)
    real(fp), pointer :: omeg(:,:,:)
    real(fp), pointer :: Ptemp(:,:,:)
    real(fp), pointer :: T1(:,:,:)
    real(fp), pointer :: T2(:,:,:)

    real(fp) :: curr_T1, next_T2, ratio
    real(fp) :: V_prev, V_new

    real(fp) :: box_u, box_v, box_omeg
    real(fp) :: curr_u, curr_v, curr_omeg
    real(fp) :: curr_u_PS, curr_v_PS

    real(fp), pointer :: X_edge(:), Y_edge(:)


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

    TYPE(Plume2d_list), POINTER :: Plume2d_new, Plume2d, Plume2d_prev
    TYPE(Plume1d_list), POINTER :: Plume1d_new, Plume1d, Plume1d_prev



    CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: ErrMsg


    Dt = GET_TS_DYN()


!    IF(Stop_inject==1) GOTO 400

    allocate(box_concnt_2D(n_x_max, n_y_max, n_species))
    allocate(box_concnt_1D(n_slab_max, n_species))


    RC     =  GC_SUCCESS
    ErrMsg = ''



    RK_Dt(1) = 0.0
    RK_Dt(2) = 0.5*Dt
    RK_Dt(3) = 0.5*Dt
    RK_Dt(4) = Dt
    RK_Dt(5) = 0.0

    u     => State_Met%U   ! figure out state_met%U is based on
                           !  lat/lon or modelgrid(i,j)
    v     => State_Met%V   ! V [m s-1]
    omeg  => State_Met%OMEGA     ! Updraft velocity [Pa/s]
    Ptemp => State_Met%THETA     ! Potential temperature [K]
    T1    => State_Met%TMPU1     ! Temperature at start of
                                 !  timestep [K]
    T2    => State_Met%TMPU2     ! Temperature at end of
                                 !  timestep [K]

    X_edge => State_Grid%XEdge(:,1) !XEDGE(:,1,1)   ! IIPAR+1 ! new
    Y_edge => State_Grid%YEdge(1,:) !YEDGE(1,:,1)  
    ! Use second YEDGE, because sometimes YMID(2)-YMID(1) is not DLAT

    X_edge2       = X_edge(2)
    Y_edge2       = Y_edge(2)
     
    


    ! -----------------------------------------------------------
    ! add new box every time step
    ! -----------------------------------------------------------
    IF(N_stop_inject<0.or.Num_inject<N_stop_inject)THEN
    DO i_box = Num_inject+1, Num_inject+N_parcel, 1

      ALLOCATE(Plume2d_new)

      Plume2d_new%IsNew = 1
      Plume2d_new%label = i_box

      Plume2d_new%LON = Inject_lon
      Plume2d_new%LAT = ( -30.005e+0_fp + Length_lat * MOD(i_box,N_total) ) &
                     * (-1.0)**FLOOR(1.0*i_box/N_total) 
!      Plume2d_new%LAT = ( -30.005e+0_fp + 0.01e+0_fp * MOD(i_box,6000) ) &
!                     * (-1.0)**FLOOR(i_box/6000.0) ! -29.995S:29.995N:0.01
      Plume2d_new%LEV = Inject_hPa

      Plume2d_new%LENGTH = Length_init ! 1000m 
      Plume2d_new%ALPHA  = 0.0e+0_fp

      Plume2d_new%LIFE = 0.0e+0_fp

      Plume2d_new%PDX = Dx_init
      Plume2d_new%PDY = Dy_init

      ALLOCATE(Plume2d_tail%rWet(n_bins))
      Plume2d_tail%rWet  = 0.0e+0_fp

      ALLOCATE(Plume2d_tail%aWP(n_x_max, n_y_max, n_bins))
      Plume2d_tail%aWP  = 0.0e+0_fp

      ALLOCATE(Plume2d_tail%aDen(n_x_max, n_y_max, n_bins))
      Plume2d_tail%aDen = 0.0e+0_fp

      ALLOCATE(Plume2d_tail%Sfc_Ten(n_x_max, n_y_max, n_bins))
      Plume2d_tail%Sfc_Ten = 0.0e+0_fp

      ALLOCATE(Plume2d_new%CONCNT2d(n_x_max, n_y_max, n_species))
      Plume2d_new%CONCNT2d = 0.0e+0_fp
      Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,i_inject) = &
		   Plume2d_new%LENGTH*Inject_rate &
                   / (Plume2d_new%PDX*Plume2d_new%Pdy*Plume2d_new%LENGTH)
      ! From [g/m3] to [molec/cm3], 98.0 g/mol for H2SO4
      Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,i_inject) = &
                Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,i_inject)/1.0e+6_fp/98.0*AVO


      mass_la = mass_la + Plume2d_new%CONCNT2d(n_x_mid,n_y_mid,i_inject)* &
             (Plume2d_new%PDX*Plume2d_new%Pdy*Plume2d_new%LENGTH*1.0e+6_fp)


!      WRITE(6,*) 'mass_la: ', Plume2d_new%CONCNT2d(n_x_mid,n_y_mid)* &
!             (Plume2d_new%PDX*Plume2d_new%Pdy*Plume2d_new%LENGTH*1.0e+6_fp)

      NULLIFY(Plume2d_new%next)
      Plume2d_tail%next => Plume2d_new
      Plume2d_tail => Plume2d_tail%next


    ENDDO ! i_box = Num_inject+1, Num_inject+N_parcel, 1

      Num_Plume2d = Num_Plume2d + N_parcel

      ! use this value to set the initial latutude for injected plume
      Num_inject = Num_inject + N_parcel


    ENDIF ! IF(N_stop_inject<0.or.Num_inject<N_stop_inject)THEN


    !--------------------------------------------------------
    ! Set the initial value of plume character:
    ! box_theta, box_length, Lifetime
    !--------------------------------------------------------
!    DO i_box = n_box_prev+1, n_box-1
!      box_lon_edge = 0.5 * ( box_lon(i_box) + box_lon(i_box+1) )
!      box_lat_edge = 0.5 * ( box_lat(i_box) + box_lat(i_box+1) )
!      box_length(i_box) = 2.0 *Distance_Circle(box_lon(i_box), box_lat(i_box), &
!                                box_lon_edge, box_lat_edge)  ! [m]
!    ENDDO
!
!
!    i_box = n_box
!    box_lon_edge = box_lon(i_box) + 0.5* ( box_lon(i_box) - box_lon(i_box-1) )
!    box_lat_edge = box_lat(i_box) + 0.5* ( box_lat(i_box) - box_lat(i_box-1) )
!    box_length(i_box) = 2.0 * Distance_Circle(box_lon(i_box), box_lat(i_box), &
!                                box_lon_edge, box_lat_edge)  ! [m]



    !=======================================================================
    ! for 2D plume: Run Lagrangian trajectory-track HERE
    !=======================================================================

    IF(.NOT.ASSOCIATED(Plume2d_head)) GOTO 401

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
!      curr_T1       = box_T1(i_box)        ! K


      ! check this grid box is the neast one or the one located in left-bottom
      ! ???
      i_lon = Find_iLonLat(curr_lon, DX, X_edge2)
      if(i_lon>IIPAR) i_lon=i_lon-IIPAR
      if(i_lon<1) i_lon=i_lon+IIPAR

      i_lat = Find_iLonLat(curr_lat, DY, Y_edge2)
      if(i_lat>JJPAR) i_lat=JJPAR
      if(i_lat<1) i_lat=1

      i_lev = Find_iPLev(curr_pressure,P_edge)
      if(i_lev>LLPAR) i_lev=LLPAR

      IF(i_lev==LLPAR) WRITE(6,*) 'box_lev:', curr_pressure, i_lev

      ! For new injected plume:
!      if(Plume2d%IsNew==1)then
!
!         Plume2d%IsNew = 0
!
!        ! ====================================================================
!        ! Add concentraion of PASV into conventional Eulerian GEOS-Chem in 
!        ! corresponding with injected parcels in Lagrangian model
!        ! For conventional GEOS-Chem for comparison with Lagrangian Model:
!        !
!        !  AD(I,J,L) = grid box dry air mass [kg]
!        !  AIRMW     = dry air molecular wt [g/mol]
!        !  MW_G(N)   = species molecular wt [g/mol]
!        !     
!        ! the conversion is:
!        ! 
!        !====================================================================
!         PASV_EU => State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_EU)  ! [kg/kg]
!
!         MW_g = State_Chm%SpcData(nAdv)%Info%MW_g
!         ! Here assume the injection rate is 30 kg/km for H2SO4: 
!         PASV_EU = PASV_EU + box_length*1.0e-3_fp*30.0 &
!                                        /State_Met%AD(i_lon,i_lat,i_lev)
!
!
!          mass_eu = mass_eu + box_length*1.0e-3_fp*30.0 * 1000.0/98*AVO
!
!          WRITE(6,*)'mass_eu: ', box_length*1.0e-3_fp*30.0 * 1000.0/98*AVO

         !======================================================================
         ! Set up the new injected plumes/box
         ! identify the location of the plume box, 
         ! put the corresponding concentration of backgound grid into the initial 
         ! concentration inside plume in unit of [molec/cm3].
         !======================================================================

         !======================================================================
         ! Convert species to [molec/cm3] (ewl, 8/16/16)
         !======================================================================
!         CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
!                                 State_Chm, 'molec/cm3', RC, OrigUnit=OrigUnit )
!
!         IF ( RC /= GC_SUCCESS ) THEN
!            ErrMsg = 'Unit conversion error!'
!            CALL GC_Error( ErrMsg, RC, 'lagrange_run in lagrange_mod.F90')
!            RETURN
!         ENDIF


!         do i_species = 1, N_species
!           box_concnt_2D(:,:,i_species)=box_concnt_2D(:,:,i_species) &
!                    + State_Chm%Species(i_lon,i_lat,i_lev,State_Chm%nAdvect-1)
!        
!         ![molec]
!
!         enddo 


         !=======================================================================
         ! Convert species back to original units (ewl, 8/16/16)
         !=======================================================================
!         CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
!                                 State_Chm, OrigUnit,  RC )
!         IF ( RC /= GC_SUCCESS ) THEN
!            ErrMsg = 'Unit conversion error!'
!            CALL GC_Error( ErrMsg, RC, 'plume_mod.F90' )
!            RETURN
!         ENDIF



!      endif ! if(Plume2d%IsNew==1)then


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

       ! For 
!       dbox_lev = RK_Dt(Ki) * curr_omeg / 100.0     ! Pa => hPa
!       curr_pressure = box_lev(i_box) + dbox_lev

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



      box_u    = ( RK_u(1) + 2.0*RK_u(2) &
                         + 2.0*RK_u(3) + RK_u(4) ) / 6.0
      box_v    = ( RK_v(1) + 2.0*RK_v(2) &
                         + 2.0*RK_v(3) + RK_v(4) ) / 6.0
      box_omeg = ( RK_omeg(1) + 2.0*RK_omeg(2) &
                         + 2.0*RK_omeg(3) + RK_omeg(4) ) / 6.0




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
      if(i_lev>LLPAR)i_lev=LLPAR


      if(abs(box_lat)>Y_mid(JJPAR))then
         next_T2 = Interplt_wind_RLL_polar(T2, i_lon, i_lat, i_lev, &
                           box_lon, box_lat, box_lev)
      else
         next_T2 = Interplt_wind_RLL(T2, i_lon, i_lat, i_lev, &
                           box_lon, box_lat, box_lev)
      endif


! !!! shw
!      ! PV=nRT, V2 = T2/P2 : T1/P1 * V1
!      ratio = ( (next_T2/box_lev)/(curr_T1/curr_pressure) )**(1/2)
!
!      ! assume the volume change mainly apply to the cross-section, 
!      ! the box_length would not change 
!
!
!      V_prev = Pdx*Pdy*box_length*1.0e+6_fp
!
!      Pdx = Pdx *ratio
!      Pdy = Pdy *ratio
!
!      V_new = Pdx*Pdy*box_length*1.0e+6_fp
!
!
!      DO i_species = 1, n_species, 1
!      !$OMP PARALLEL DO           &
!      !$OMP DEFAULT( SHARED     ) &
!      !$OMP PRIVATE( i_y, i_x )
!      DO i_y = 1,n_y_max,1
!      DO i_x = 1,n_x_max,1
!         box_concnt_2D(i_x,i_y,i_species) = &
!                        box_concnt_2D(i_x,i_y,i_species)*V_prev/V_new
!      ENDDO
!      ENDDO
!      !$OMP END PARALLEL DO
!      ENDDO


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

! !!! shw
!      length0           = box_length
!      box_length = EXP(Ly*Dt) * length0
!
!
!      Pdx = Pdx*SQRT(length0/box_length)
!      Pdy = Pdy*SQRT(length0/box_length)



      ! ----------------------------------------------------------------
      ! update     
      ! ----------------------------------------------------------------

      ! !!! shw
      Plume2d%LON    = Inject_lon
      Plume2d%LAT    = -29.95e+0_fp
      Plume2d%LEV    = Inject_hPa

!      Plume2d%LON = box_lon
!      Plume2d%LAT = box_lat
!      Plume2d%LEV = box_lev
      Plume2d%LENGTH = box_length
      Plume2d%ALPHA  = box_alpha
      Plume2d%label  = box_label
      Plume2d%LIFE   = box_life


      Plume2d%PDX = Pdx
      Plume2d%PDY = Pdy

      Plume2d%CONCNT2d    = box_concnt_2D

      Plume2d => Plume2d%next

    ENDDO  ! DO WHILE(ASSOCIATED(Plume2d))

401 CONTINUE

!    WRITE(6,*)'=== loop for 2d plume in lagrange_run: ', i_box, Num_Plume2d




    !=======================================================================
    ! For 1d plume: Run Lagrangian trajectory-track HERE
    !=======================================================================

    IF(.NOT.ASSOCIATED(Plume1d_head)) GOTO 400


    Plume1d => Plume1d_head
    i_box = 0

    DO WHILE(ASSOCIATED(Plume1d))

      i_box = i_box+1


      box_lon    = Plume1d%LON
      box_lat    = Plume1d%LAT
      box_lev    = Plume1d%LEV
      box_length = Plume1d%LENGTH
      box_alpha  = Plume1d%ALPHA
      box_label  = Plume1d%label
      box_life   = Plume1d%LIFE

      box_Ra    = Plume1d%RA
      box_Rb    = Plume1d%RB

!      WRITE(6,*) i_box, SHAPE(box_concnt_1D), SHAPE(Plume1d%CONCNT1d)

      box_concnt_1D = Plume1d%CONCNT1d


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


!      call cpu_time(start)
!
!      !$OMP PARALLEL DO           &
!      !$OMP DEFAULT( SHARED     ) &
!      !$OMP PRIVATE( i_species, i_slab )
!      DO i_species = 1, n_species, 1
!      DO i_slab = 1,n_slab_max,1
!
!         box_concnt_1D(i_slab,i_species) = &
!                        box_concnt_1D(i_slab,i_species)*V_prev/V_new
!
!      ENDDO
!      ENDDO
!      !$OMP END PARALLEL DO
!
!      call cpu_time(finish)
!      WRITE(6,*)'Time2 (finish-start) for 1D:', finish-start


!      IF(box_label==1)WRITE(6,*)'mass test after idea gas law:', i_box, &
!                                    V_new * SUM(box_concnt_1D(1:n_slab_max))

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



!      IF(box_label==1)WRITE(6,*)'mass test after stretch:', i_box, &
!                                box_Ra*box_Rb*box_length*1.0e+6_fp &
!                                  * SUM(box_concnt_1D(1:n_slab_max))

!      IF(box_label==871)THEN
!        WRITE(6,*)'Length (km):', NINT(box_length/1000), Ly
!        WRITE(6,*) box_lon, box_lat, box_lev
!        WRITE(6,*) '  '
!      ENDIF


! !!! shw
      Plume1d%LON    = Inject_lon
      Plume1d%LAT    = -29.95e+0_fp
      Plume1d%LEV    = Inject_hPa

!      Plume1d%LON    = box_lon
!      Plume1d%LAT    = box_lat
!      Plume1d%LEV    = box_lev
      Plume1d%LENGTH = box_length
      Plume1d%ALPHA  = box_alpha
      Plume1d%label  = box_label
      Plume1d%LIFE   = box_life
      Plume1d%RA     = box_Ra
      Plume1d%RB     = box_Rb

      Plume1d%CONCNT1d = box_concnt_1D


      Plume1d => Plume1d%next

    ENDDO  ! DO WHILE(ASSOCIATED(Plume1d))

!    WRITE(6,*) '=== loop for 1d plume in lagrange_run: ', i_box, Num_Plume1d


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
    IF(ASSOCIATED(Plume2d)) nullify(Plume2d)
    IF(ASSOCIATED(Plume2d_prev)) nullify(Plume2d_prev)
    IF(ASSOCIATED(Plume1d_new)) nullify(Plume1d_new)
    IF(ASSOCIATED(Plume1d)) nullify(Plume1d)
    IF(ASSOCIATED(Plume1d_prev)) nullify(Plume1d_prev)

  END SUBROUTINE lagrange_run



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

    ! i,j means the four grid point value that surround around the particle
!    do i=1,2
!    do j=1,2

!      ii = i + init_lon - 1
!      jj = j + init_lat - 1
!
      ! For lon=180 deg:
!      if(ii==IIPAR+1)then
!        ii = 1
!      endif
!      if(ii==0)then
!        ii = IIPAR
!      endif

      ! Get the right ii and jj for interpolation at polar point:

      ! For South Polar Point:
 !     if(jj==0)then
 !       jj = jj+1
 !     endif

      ! For North Polar Point:
!      if(jj==JJPAR+1)then
!        jj = jj-1
!      endif

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

    ! angle between travel direction and lon: beta belongs to [-PI/2,3*PI/2)

    ! same as angle between cross-section (right when 
    ! toward travel direction) and lon: box_alpha belongs to [0, 2*PI)
!    real(fp)          :: box_alpha ! beta! [radians]

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


    TYPE(Plume2d_list), POINTER :: Plume2d_new, Plume2d, Plume2d_prev
    TYPE(Plume1d_list), POINTER :: Plume1d_new, Plume1d, Plume1d_prev


    CHARACTER(LEN=255)     :: FILENAME2
    CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: ErrMsg

    INTEGER :: YEAR
    INTEGER :: MONTH
    INTEGER :: DAY
    INTEGER :: HOUR
    INTEGER :: MINUTE
    INTEGER :: SECOND

    CHARACTER(LEN=255)     :: FileConcnt, FileXYZ, FileEntropy, File996

    CHARACTER(LEN=25) :: YEAR_C
    CHARACTER(LEN=25) :: MONTH_C
    CHARACTER(LEN=25) :: DAY_C
    CHARACTER(LEN=25) :: HOUR_C
    CHARACTER(LEN=25) :: MINUTE_C
    CHARACTER(LEN=25) :: SECOND_C


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

!    IF(mod(tt*NINT(Dt),24*60*60)==0)THEN   ! output once every day (24 hours)
!
!      FileConcnt   = 'Lagrange_Concnt_Volume_' // TRIM(ADJUSTL(YEAR_C)) // '-' &
!        //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) // '-' &
!        // TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) &
!        // ':' // TRIM(ADJUSTL(SECOND_C)) // '.txt'
!
!
!      OPEN( 485,      FILE=TRIM( FileConcnt   ), STATUS='REPLACE',  &
!          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!      CLOSE(485)
!
!
!      FileXYZ   = 'Lagrange_xyz_' // TRIM(ADJUSTL(YEAR_C)) // '-' &
!        //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) // '-' &
!        // TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) &
!        // ':' // TRIM(ADJUSTL(SECOND_C)) // '.txt'
!
!
!      OPEN( 486,      FILE=TRIM( FileXYZ   ), STATUS='REPLACE',  &
!          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!      CLOSE(486)
!
!
!
!    ENDIF ! IF(mod(tt,1440)==0)THEN




!    IF(Stop_inject==1) GOTO 999

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
                            'molec/cm3', RC, OrigUnit=OrigUnit )
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


    !!! shw: for a test
      IF(i_lev/=39 .or. i_lat/=16 .or. i_lon/=9)THEN
        State_Chm%Species(i_lon,i_lat,i_lev,:) = 0.0
      ENDIF



      ! set initial value for Entropy_V
      i_cell = (i_lev-1)*IIPAR*JJPAR+(i_lat-1)*IIPAR+i_lon
      Entropy_V(i_cell) = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp


    ! from EU to EU2
!    State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_EU2) = &
!              State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_EU2) &
!            + Dt* Kchem*State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_EU)**2

!      ! from LA to LA2
!      State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA2) = &
!                 State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA2) &
!               + Dt*Kchem*State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA)**2


      !  For all the grid cells in the troposphere, let concentration to be zero
      IF(TROPP_sink == 1)THEN
      IF ( State_Met%PEDGE(i_lon,i_lat,i_lev)>State_Met%TROPP(i_lon,i_lat) ) THEN
        State_Chm%Species(i_lon,i_lat,i_lev,:) = 0.0
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



!
!       backgrd_concnt = State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA)

!      !$OMP PARALLEL DO           &
!      !$OMP DEFAULT( SHARED     ) &
!      !$OMP PRIVATE( i_y, i_x )
!      DO i_y = 1,n_y_max,1
!      DO i_x = 1,n_x_max,1
!
!      ! run the fake 2nd order chemical reaction
!      box_concnt_2D(i_x,i_y,2) = box_concnt_2D(i_x,i_y,2) &
!         + Dt* Kchem*(box_concnt_2D(i_x,i_y,1)+backgrd_concnt)**2 &
!         - Dt* Kchem*backgrd_concnt**2
!
!      ENDDO
!      ENDDO
!      !$OMP END PARALLEL DO




       !====================================================================
       ! calculate the wind shear along plume corss-section
       ! clock-wise 90 degree from the plume length direction (box_alpha)
       ! calculate the diffusivity in horizontal and vertical direction
       !====================================================================

       ! calculate the wind_s shear along pressure direction
       wind_s_shear = Wind_shear_s(u, v, P_BXHEIGHT, box_alpha, i_lon, &
                        i_lat, i_lev,curr_lon, curr_lat, curr_pressure)
       ! *** attention *** ???
!       wind_s_shear = 0.004 ! 0.004


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
!       eddy_v = 0.1 !1.0 ! 0.1 ???

!       WRITE(6,*)'wind_s_shear, eddy_v: ', wind_s_shear, eddy_v

       ! Calculate horizontal eddy diffusivity:
!       Ch = 0.1
!       U_shear = Vertical_shear(u, P_BXHEIGHT, i_lon, i_lat, i_lev, curr_lon,
!       curr_lat, curr_pressure)
!       V_shear = Vertical_shear(v, P_BXHEIGHT, i_lon, i_lat, i_lev, curr_lon,
!       curr_lat, curr_pressure)
!       UV_shear = SQRT( U_shear**2 + V_shear**2 )
!       eddy_h(i_slab) = Ch*UV_shear*(Init_radius+(i_slab-1)*D_radius)**2
       eddy_h = 10.0  ! ???

! test for horizontal diffusivity
!       L_b = 2.0*PI*SQRT(curr_u**2+curr_v**2)/N_BV
!       Ee = 0.5 * Omega_N**2 * wind_s_shear**2
!       L_O = 2.0*PI*SQRT( Ee/(N_BV**3) )


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


           ! uodate box_extra(i_box)
!           V_grid_2D       = Pdx*Pdy*box_length*1.0e+6_fp
!           mass_plume      = V_grid_2D/9 *SUM(Pc(:,:))
!
!           mass_plume_new = V_grid_2D &
!                *SUM( box_concnt_2D(:,:))
!
!           D_mass_plume = mass_plume_new - mass_plume

!           IF(D_mass_plume>0)THEN
!             WRITE(6,*)'--- 9 to 1 grid cell:', i_box, mass_plume, mass_plume_new
!             box_extra = box_extra + D_mass_plume
!           ELSE
!             WRITE(6,*)'--- ERROR in 9 to 1 grid cell in 2D grid ---'
!           ENDIF



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

!       IF(box_label==1) WRITE(6,*)'Dt in 2D is:', Pdt, Pdx/Pu(1,2)
!       IF(box_label==1) WRITE(6,*)   Pdx/Pu(1,1), Pdy**2/(2*eddy_v)

       Concnt2D_bdy(:,:) = 0.0
       Concnt2D_bdy(2:n_x_max2-1,2:n_y_max2-1) = box_concnt_2D(:,:,i_species)

        
       IF( abs(Dt/Pdt-Nt) > 0.00001 ) THEN
         WRITE(6,*) "*** ERROR: Check Pdt ***"
         WRITE(6,*) Pdt, Nt, Dt
         WRITE(6,*) Pdx/Pu(1,1), Pdy**2/(2*eddy_v), Pdx**2/(2*eddy_h)
       ENDIF
       



       DO t1s = 1, NINT(Dt/Pdt)

         ! advection ----------------------------------------------

!         Pc_bdy(:,:) = 0.0
!         Pc_bdy(2:n_x_max2-1,2:n_y_max2-1) = Concnt2D_bdy(2:n_x_max2-1,2:n_y_max2-1)
!
!
!         ! Only calculate the vertical half 2D domain         
!
!         !$OMP PARALLEL DO           &
!         !$OMP DEFAULT( SHARED     ) &
!         !$OMP PRIVATE(i_y,i_x,CFL,Pc_middle,Pc_top,Pc_bottom,Pc_right,Pc_left)
!         DO i_y = 1, n_y_mid2, 1
!         DO i_x = 2, n_x_max2-1, 1
!
!           Pc_middle = Pc_bdy( i_x,   i_y  )
!           Pc_right  = Pc_bdy( i_x+1, i_y  )
!           Pc_left   = Pc_bdy( i_x-1, i_y  )
!
!           CFL       = Pdt*Pu(i_x,i_y)/Pdx
!
!           Concnt2D_bdy(i_x, i_y) = Pc_middle           &
!              - 0.5 * CFL    * ( Pc_right - Pc_left )   &
!              + 0.5 * CFL**2 * ( Pc_right - 2*Pc_middle + Pc_left ) 
!
!         ENDDO
!         ENDDO
!         !$OMP END PARALLEL DO
!
!         ! update the other half based on vertical symmetry 
!         Concnt2D_bdy(2:n_x_max2-1:1, n_y_mid2+1) = Concnt2D_bdy(n_x_max2-1:2:-1, n_y_mid2-1)


         ! diffusion ----------------------------------------------------

         Pc_bdy(:,:) = 0.0
         Pc_bdy(2:n_x_max2-1,2:n_y_max2-1) = Concnt2D_bdy(2:n_x_max2-1,2:n_y_max2-1)


         ! Only calculate the vertical half 2D domain         

         !$OMP PARALLEL DO           &
         !$OMP DEFAULT( SHARED     ) &
         !$OMP PRIVATE(i_y,i_x,CFL,Pc_middle,Pc_top,Pc_bottom,Pc_right,Pc_left)
         DO i_y = 2, n_y_mid2, 1
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


! !!! shw
!         box_concnt_2D(:,:,i_species) = Concnt2D_bdy(2:n_x_max2-1,2:n_y_max2-1)
!
!         !================================================================
!         ! Update the concentration in the background grid cell
!         ! after the interaction with 2D plume
!         !================================================================
!
!         grid_volumn     = State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp ! [cm3]
!         V_grid_2D       = Pdx*Pdy*box_length*1.0e+6_fp
!
!
!         ! the boundary always represents the background concentration
!         ! [molec]
!         mass_plume      = V_grid_2D*SUM(C2d_prev(:,:))
!
!         mass_plume_new = V_grid_2D * SUM(box_concnt_2D(:,:,i_species))
!
!
!         D_mass_plume = mass_plume_new - mass_plume
!
!
!         backgrd_concnt = State_Chm%Species(i_lon,i_lat,i_lev,id_GC(i_species))
!
!         backgrd_concnt = ( backgrd_concnt*grid_volumn &
!                                           - D_mass_plume) /grid_volumn
!
!         State_Chm%Species(i_lon,i_lat,i_lev,id_GC(i_species)) = backgrd_concnt
!
!
       ENDDO ! DO i_species=1,n_species,1



         !====================================================================
         ! Change from 2D to 1D, 
         ! once the tilting degree is bigger than 88 deg (88/180*3.14)
         !====================================================================

!         IF(Plume2d%LIFE>4.0*3600.0)THEN
         IF(Plume2d%LIFE>10.0*3600.0)THEN ! !!! shw

         frac_mass = 0.95


         Pc2 = box_concnt_2D(:,:,1)
         Xscale = Get_XYscale(Pc2, Pdx, Pdy, frac_mass, 2)
         Yscale = Get_XYscale(Pc2, Pdx, Pdy, frac_mass, 1)

         box_theta = ATAN( Xscale/Yscale )



!         IF( Plume2d%LIFE>18.0*3600.0 ) THEN
!           WRITE(6,*) "*** ERROR: 2D plume live too long! ***"
!           WRITE(6,*) box_label, Xscale, Yscale, Pdx, Pdy 
!           WRITE(6,*) wind_s_shear, eddy_h, eddy_v
!           WRITE(6,*) Pc2(n_x_mid-1,:)
!           WRITE(6,*) Pc2(:,n_y_mid)
!         ENDIF


!         IF( Xscale/Yscale>25.0 &
!                        .or. Plume2d%LIFE>12.0*3600.0 ) THEN !!! shw ???
!
         IF( Xscale/Yscale>25.0 &
                        .or. Plume2d%LIFE>12.0*3600.0 ) THEN !!! shw ???

           DO i_species=1,n_species,1

           Pc = box_concnt_2D(:,:,i_species)


           IF(box_theta>1.569) WRITE(6,*)&
            '*** ERROR *** i_box, box_theta:', i_box, box_theta,&
                                                     Plume2d%LIFE

           ! assign length/width and concentration to 1D slab model 
           ! return initial box_Ra(i_box), box_Rb(i_box), box_concnt_1D(i_box)

           Cslab = box_concnt_1D(:,i_species)
           CALL Slab_init(Pdx, Pdy, box_theta, Pc, Yscale, &
                               box_Ra, box_Rb, Cslab)
           box_concnt_1D(:,i_species) = Cslab


           V_grid_1D = box_Ra*box_Rb*box_length*1.0e+6_fp

           mass_plume_new = V_grid_1D &
                            * SUM(box_concnt_1D(1:n_slab_max,i_species))

           mass_plume = V_grid_2D * SUM(box_concnt_2D(:,:,i_species))



           !!! update the background concentration:
           D_mass_plume = mass_plume_new - mass_plume


!           IF(ABS(D_mass_plume/mass_plume)>0.01)THEN
!             WRITE(6,*) '*********************************************'
!             WRITE(6,*) '    more than 1% mass lost from 2-D to 1-D'
!             WRITE(6,*) ABS(D_mass_plume/mass_plume)
!             WRITE(6,*) '*********************************************'
!           ENDIF


!           i_advect = id_PASV_LA +i_species -1

           backgrd_concnt = State_Chm%Species(i_lon, i_lat, i_lev, id_GC(i_species))

           backgrd_concnt = ( backgrd_concnt*grid_volumn &
                                            - D_mass_plume ) /grid_volumn

           State_Chm%Species(i_lon, i_lat, i_lev, id_GC(i_species)) = backgrd_concnt

           ENDDO ! DO i_species=1,n_species,1


           ! ----------------------------------
           ! delete the node for 2D plume:
           ! ----------------------------------

           IF(.NOT.ASSOCIATED(Plume2d%next) .and. &
                          .NOT.ASSOCIATED(Plume2d_prev))THEN

             WRITE(6,*)'                '
             WRITE(6,*)'*** Delete the last Plume2d', Num_Plume2d
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

!            Plume2d => Plume2d_prev%next !!! ???
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


             ALLOCATE(Plume1d_tail%rWet(n_bins))
             Plume1d_tail%rWet  = 0.0e+0_fp

             ALLOCATE(Plume1d_tail%aWP(n_slab_max, n_bins))
             Plume1d_tail%aWP  = 0.0e+0_fp

    	     ALLOCATE(Plume1d_tail%aDen(n_slab_max, n_bins))
   	     Plume1d_tail%aDen = 0.0e+0_fp

    	     ALLOCATE(Plume1d_tail%Sfc_Ten(n_slab_max, n_bins))
    	     Plume1d_tail%Sfc_Ten = 0.0e+0_fp


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


             ALLOCATE(Plume1d_tail%rWet(n_bins))
             Plume1d_tail%rWet  = 0.0e+0_fp

             ALLOCATE(Plume1d_tail%aWP(n_slab_max, n_bins))
             Plume1d_tail%aWP  = 0.0e+0_fp

             ALLOCATE(Plume1d_tail%aDen(n_slab_max, n_bins))
             Plume1d_tail%aDen = 0.0e+0_fp

             ALLOCATE(Plume1d_tail%Sfc_Ten(n_slab_max, n_bins))
             Plume1d_tail%Sfc_Ten = 0.0e+0_fp


             ALLOCATE(Plume1d_new%CONCNT1d(n_slab_max,n_species))
             Plume1d_new%CONCNT1d = box_concnt_1D

             NULLIFY(Plume1d_new%next)
             Plume1d_tail%next => Plume1d_new
             Plume1d_tail => Plume1d_tail%next

!             WRITE(6,*) '*** CREATE other Plume1d nodes ***', i_box

           ENDIF

           Num_Plume1d = Num_Plume1d + 1

           IF(Stop_loop==1)THEN
             GOTO 900 ! skip the whole loop, go to Plume1d
           ELSE
             GOTO 800 ! skip this box, go to next box
           ENDIF

         ENDIF ! IF(box_theta(i_box)>98/180*3.14)
         ENDIF ! IF(Lifetime(i_box)>3*3600.0)THEN


!         !===================================================================
!         ! Third judge:
!         ! If this is the final time step, the model is going to end
!         ! dissolve all the plume into GCM in the last time step
!         !===================================================================
!         IF(MONTH==1 .AND. DAY==31)THEN
!         IF(HOUR==23 .AND. MINUTE>=50)THEN
!
!           Stop_inject = 1
!
!           i_box = i_box-1
!           Num_Plume2d  = Num_Plume2d - 1
!           Num_dissolve = Num_dissolve + 1
!
!
!           DO i_species=1,n_species
!
!           i_advect = id_PASV_LA +i_species -1
!
!           backgrd_concnt = State_Chm%Species(i_lon,i_lat,i_lev,i_advect)
!
!           backgrd_concnt = &
!                       ( SUM( box_concnt_2D(:,:,i_species) )*V_grid_2D  &
!                           + backgrd_concnt*grid_volumn ) / grid_volumn
!
!           State_Chm%Species(i_lon,i_lat,i_lev,i_advect) = backgrd_concnt
!
!           ENDDO
!
!           ! ----------------------------------
!           ! delete the node for 2D plume:
!           ! ----------------------------------
!
!           IF(.NOT.ASSOCIATED(Plume2d%next) .and. &
!                          .NOT.ASSOCIATED(Plume2d_prev))THEN
!
!              DEALLOCATE(Plume2d%CONCNT2d)
!              DEALLOCATE(Plume2d)
!              GOTO 900
!
!           ENDIF
!
!
!           IF(.NOT.ASSOCIATED(Plume2d%next))THEN ! delete the tail node
!              ! the tail node is also the head node, skip the loop
!              IF(.NOT.ASSOCIATED(Plume1d_prev)) GOTO 900
!
!              Plume2d_tail%next => Plume2d_prev
!              Plume2d_tail => Plume2d_tail%next
!!            Plume2d => Plume2d_prev%next !!! ???
!              IF(ASSOCIATED(Plume2d_prev%next)) NULLIFY(Plume2d_prev%next)
!              DEALLOCATE(Plume2d%CONCNT2d)
!              DEALLOCATE(Plume2d)
!              GOTO 900
!           ELSEIF(ASSOCIATED(Plume2d_prev))THEN ! delete node, not head/tail
!              Plume2d => Plume2d_prev%next
!              Plume2d_prev%next => Plume2d%next
!              DEALLOCATE(Plume2d%CONCNT2d)
!              DEALLOCATE(Plume2d)
!              Plume2d => Plume2d_prev%next
!           ELSE ! delete head node
!              ! the head node is also the tail node, skip the loop
!              IF(.NOT.ASSOCIATED(Plume2d%next)) GOTO 900
!              Plume2d => Plume2d_head
!              Plume2d_head => Plume2d%next
!              DEALLOCATE(Plume2d%CONCNT2d)
!              DEALLOCATE(Plume2d)
!              Plume2d => Plume2d_head
!           ENDIF
!
!
!
!
!           GOTO 800 ! skip this box, go to next box
!
!         ENDIF ! IF(MONTH==2 .AND. DAY==31)THEN
!         ENDIF ! IF(HOUR==23 .AND. MINUTE==50)THEN

         !===================================================================
         ! Fourth judge:
         ! If the plume touch the tropopause, dissolve the plume
         !===================================================================
         IF ( box_lev>State_Met%TROPP(i_lon,i_lat) ) THEN

           Num_dissolve = Num_dissolve+1


           DO i_species=1,n_species

!           i_advect = id_PASV_LA +i_species -1

           backgrd_concnt = State_Chm%Species(i_lon,i_lat,i_lev,id_GC(i_species))

           backgrd_concnt = ( SUM( box_concnt_2D(:,:,i_species) )*V_grid_2D &
                             + backgrd_concnt*grid_volumn ) / grid_volumn


           State_Chm%Species(i_lon,i_lat,i_lev,id_GC(i_species)) = backgrd_concnt

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




!       IF(mod(tt*NINT(Dt),24*60*60)==0)THEN   ! output once every day (24 hours)
!
!         OPEN( 485,      FILE=TRIM( FileConcnt   ), STATUS='OLD',  &
!           POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         OPEN( 486,      FILE=TRIM( FileXYZ   ), STATUS='OLD',  &
!           POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         Max_Concnt_Grid = MAXVAL(State_Chm%Species(:,:,:,id_PASV_LA))
!
!         DO i_y = 1,n_y_max,1
!         DO i_x = 1,n_x_max,1
!
!           Entropy_Concnt = box_concnt_2D(i_x,i_y,i_entropy) + &
!                        State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA)
!
!           tracer_mol = Entropy_Concnt * V_grid_2D / AVO
!           air_mol    = V_grid_2D/1e+6_fp*State_Met%AIRDEN(i_lon,i_lat,i_lev)*1000/AIRMW
!           mix_ratio  = tracer_mol/air_mol
!
!           IF(mix_ratio>0) Entropy = Entropy - BOLTZ*tracer_mol*log(mix_ratio)
!
!           i_cell = (i_lev-1)*IIPAR*JJPAR+(i_lat-1)*IIPAR+i_lon
!           Entropy_V(i_cell) = Entropy_V(i_cell) - V_grid_2D
!
!
!
!
!           ! Only keep the plume grid cell that has a larger concentration 
!           ! than the background concentration
!           IF( box_concnt_2D(i_x,i_y,i_tracer) > 100.0*Max_Concnt_Grid )THEN
!
!              WRITE(485,*) box_concnt_2D(i_x,i_y,i_tracer), &
!                              V_grid_2D, State_Met%AIRDEN(i_lon,i_lat,i_lev)
!
!              WRITE(486,*) i_lon,i_lat,i_lev
!
!           ELSE
!
!              backgrd_concnt = State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA4)
!
!              backgrd_concnt = &
!                   ( box_concnt_2D(i_x,i_y,i_tracer)*V_grid_2D  &
!                        + backgrd_concnt*grid_volumn ) / grid_volumn
!
!              State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA4) = backgrd_concnt
!
!           ENDIF
!
!         ENDDO
!         ENDDO
!
!
!
!         CLOSE(485)
!         CLOSE(486)
!
!       ENDIF



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
!       backgrd_concnt = State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA3)
!
!       backgrd_concnt = &
!                   ( SUM( box_concnt_2D(:,:,1) )*V_grid_2D  &
!                        + backgrd_concnt*grid_volumn ) / grid_volumn
!
!       State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA3) = backgrd_concnt
        

800     CONTINUE

    ENDDO ! DO WHILE(ASSOCIATED(Plume2d))

900     CONTINUE

!    WRITE(6,*)'=== loop for 2d plume in plume_run: ', i_box, Num_Plume2d



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
    


!     ! Compare SumV_Plume with Eulerian grid cell volume
!     ! determine the Eulerian grid cell containing too much plume segment
!     DO i_lev = 1, LLPAR
!     DO i_lat = 1, JJPAR
!     DO i_lon = 1, IIPAR
!
!        i_cell = (i_lev-1)*IIPAR*JJPAR+(i_lat-1)*IIPAR+i_lon
!
!        ! [cm3]
!        SumV_Plume(i_cell) = SumV_Plume(i_cell) - &
!                                   State_Met%AIRVOL(i_lon,i_lat,i_lev)*1e+6_fp
!
!        ! Now SumV_Plume>0 means too many plume segments in one Eulerian 
!        ! grid cell: accumalative plume volume exceed the volume of the host 
!        ! Eulerian grid cell.
!        ! We need to transfer the plume segment to the host Eulerian grid
!        ! cell to ensure SumV_Plume<0
!
!        ! First, select all the 1-D plume segment located in this host Eulerian
!        ! grid cell, label selected plume segment as Is_transfer=1
!        
!
!     ENDDO
!     ENDDO
!     ENDDO


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



!         call cpu_time(start)
!
!         !$OMP PARALLEL DO           &
!         !$OMP DEFAULT( SHARED     ) &
!         !$OMP PRIVATE( i_slab )
!         DO i_slab = 1,n_slab_max,1
!
!         ! run the fake 2nd order chemical reaction
!         box_concnt_1D(i_slab,2) = box_concnt_1D(i_slab,2) &
!                              + Dt* Kchem*box_concnt_1D(i_slab,1)**2
!
!         ENDDO
!         !$OMP END PARALLEL DO
!
!         call cpu_time(finish)
!         WRITE(6,*)'Time2 (finish-start) for 2D:', i_box, finish-start




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






!         ! run the fake 2nd order chemical reaction
!         backgrd_concnt = State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA)
!
!         box_concnt_1D(:,2) = box_concnt_1D(:,2) &
!             + Dt* Kchem* (box_concnt_1D(:,1)+backgrd_concnt)**2 &
!             - Dt* Kchem* backgrd_concnt**2



        !====================================================================
        ! calculate the wind shear along plume corss-section
        ! calculate the diffusivity in horizontal and vertical direction
        !====================================================================
        ! calculate the wind_s shear along pressure direction
        wind_s_shear = Wind_shear_s(u, v, P_BXHEIGHT, box_alpha, i_lon, &
                         i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
        ! *** attention *** ???
!        wind_s_shear = 0.002


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
!        eddy_v = 1.0 ! ???

       ! Calculate horizontal eddy diffusivity:
!       Ch = 0.1
!       U_shear = Vertical_shear(u, P_BXHEIGHT, i_lon, i_lat, i_lev, curr_lon,
!       curr_lat, curr_pressure)
!       V_shear = Vertical_shear(v, P_BXHEIGHT, i_lon, i_lat, i_lev, curr_lon,
!       curr_lat, curr_pressure)
!       UV_shear = SQRT( U_shear**2 + V_shear**2 )
!       eddy_h(i_slab) = Ch*UV_shear*(Init_radius+(i_slab-1)*D_radius)**2
        eddy_h = 10.0  ! ???

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





!       V_grid_1D        = box_Ra*box_Rb*box_length*1.0e+6_fp
!       IF(box_label==1)  WRITE(6,*) "Max mass cell 0:", &
!                MAXVAL(box_concnt_1D(:,i_tracer))*V_grid_1D/AVO





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



!       IF(box_label==1) WRITE(6,*) "Max mass cell 1:", &
!                MAXVAL(box_concnt_1D(:,i_tracer))*V_grid_1D/AVO


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

!       i_advect = id_PASV_LA +i_species -1

       backgrd_concnt = State_Chm%Species( i_lon, i_lat, i_lev, id_GC(i_species) )

       backgrd_concnt = ( backgrd_concnt*grid_volumn &
                                              - D_mass_plume) / grid_volumn

       State_Chm%Species(i_lon,i_lat,i_lev,id_GC(i_species)) = backgrd_concnt

       ENDDO ! DO i_species=1,n_species,1





!       IF(box_label==1) WRITE(6,*)  "Max mass cell 2:", &
!                MAXVAL(box_concnt_1D(:,i_tracer))*V_grid_1D/AVO






       !===================================================================
       ! Second judge:
       ! According to a 2-order process rate, if there is no big difference
       ! between PIG and Eulerian model, then dissolve the PIG into
       ! background grid cell.                                          shw
       ! test this once big time step
       !===================================================================

!---------------------------
! (1)
!---------------------------
!       ! Without consider the interaction between Plume model and Eulerian model
!       backgrd_concnt = State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA)
!
!       Prod_before = V_grid_1D * SUM( box_concnt_1D(1:n_slab_max,1)**2 ) &
!                + backgrd_concnt**2 * grid_volumn
!
!       Eul_concnt = ( SUM(box_concnt_1D(1:n_slab_max,1)) *V_grid_1D &
!                    + backgrd_concnt*grid_volumn ) / grid_volumn
!       Prod_after = Eul_concnt**2*grid_volumn
!
!
!       IF(box_label==1)THEN
!        WRITE(6,*)" "
!        WRITE(6,*)'Product criterion :'
!        WRITE(6,*)Prod_before, Prod_after, ABS(Prod_before-Prod_after)/Prod_after
!       ENDIF


!---------------------------
! (2)
!---------------------------
!       ! (1) Only consider the product producted by the tracer in Plume model
!       ! (2) Consider the interaction between Plume model and Eulerian model
!       ! Reference: Karamchandani et al., (2000)
!       !            Korsakissok and Mallet, (2010)

       backgrd_concnt = State_Chm%Species(i_lon,i_lat,i_lev,id_GC(i_inject)) ! shw, need change

       Prod_before = SUM( ( box_concnt_1D(1:n_slab_max,i_inject)+backgrd_concnt )**2 &
                                               - backgrd_concnt**2 ) *V_grid_1D

       Eul_concnt = ( SUM(box_concnt_1D(1:n_slab_max,i_inject)) *V_grid_1D &
                    + backgrd_concnt*grid_volumn ) / grid_volumn

       Prod_after = (Eul_concnt**2 - backgrd_concnt**2) *grid_volumn

!       IF(box_label==1) WRITE(6,*)Prod_before, Prod_after, ABS(Prod_before-Prod_after)/Prod_after
!       WRITE(6,*)' '
!
!       WRITE(6,*)'1d test:', i_box, Num_Plume1d


!---------------------------
! (3)
!---------------------------
!       ! Consider the background concentration;
!       ! plume segment is easier to be dissolved in the higher background
!       ! concentration;
!
!       backgrd_concnt = State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA)
!
!       Prod_before = SUM( ( box_concnt_1D(1:n_slab_max,1)+backgrd_concnt )**2 &
!                           - backgrd_concnt**2 ) *V_grid_1D &
!                     + backgrd_concnt**2 *grid_volumn
!
!       Eul_concnt = ( SUM(box_concnt_1D(1:n_slab_max,1)) *V_grid_1D &
!                    + backgrd_concnt*grid_volumn ) / grid_volumn
!
!       Prod_after = Eul_concnt**2 *grid_volumn



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


       ! !!! shw : Once plume1d is created, dissolve it immediately
       Plume1d%Is_transfer = 1


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

!         WRITE(6,*)'DISSOLVE1:', i_box, box_label, box_life/3600/24,&
!                              SUM(box_concnt_1D(1:n_slab_max)) * V_grid_1D

         mass_la2 = mass_la2 + SUM(box_concnt_1D(1:n_slab_max,1)) * V_grid_1D

         DO i_species=1,n_species,1
!         i_advect = id_PASV_LA +i_species -1

         backgrd_concnt = &
                  State_Chm%Species(i_lon,i_lat,i_lev,id_GC(i_species))

         backgrd_concnt = &
           ( SUM(box_concnt_1D(1:n_slab_max,i_species)) * V_grid_1D   &
                + backgrd_concnt*grid_volumn ) / grid_volumn

         State_Chm%Species(i_lon,i_lat,i_lev,id_GC(i_species)) = backgrd_concnt
         ENDDO

         ! ----------------------------------
         ! delete the node for 1D plume:
         ! ----------------------------------

         IF(.NOT.ASSOCIATED(Plume1d%next) .and. &
                        .NOT.ASSOCIATED(Plume1d_prev))THEN

            WRITE(6,*)'                 '
            WRITE(6,*)'*** Delete the last Plume1d', i_box, box_label, Num_Plume1d
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


       !===================================================================
       ! Third judge:
       ! If this is the final time step, the model is going to end
       ! dissolve all the plume into GCM in the last time step
       !===================================================================
!!       IF ( ITS_TIME_FOR_EXIT() ) THEN
!       IF(MONTH==1 .AND. DAY==31)THEN
!       IF(HOUR==23 .AND. MINUTE>=50)THEN
!
!         Stop_inject = 1
!
!         i_box = i_box-1
!         Num_Plume1d  = Num_Plume1d - 1
!         Num_dissolve = Num_dissolve+1
!!         WRITE(6,*) 'test:', box_label, eddy_v, eddy_h, box_theta
!
!
!         DO i_species=1,n_species,1
!         i_advect = id_PASV_LA +i_species -1
!
!         backgrd_concnt = &
!                  State_Chm%Species(i_lon,i_lat,i_lev,i_advect)
!
!         backgrd_concnt = &
!           ( SUM(box_concnt_1D(1:n_slab_max,i_species)) * V_grid_1D   &
!                + backgrd_concnt*grid_volumn ) / grid_volumn
!
!         State_Chm%Species(i_lon,i_lat,i_lev,i_advect) = backgrd_concnt
!         ENDDO
!
!
!
!         IF(.NOT.ASSOCIATED(Plume1d%next) .and. &
!                                .NOT.ASSOCIATED(Plume1d_prev))THEN
!            DEALLOCATE(Plume1d%CONCNT1d)
!            DEALLOCATE(Plume1d)
!            GOTO 500 ! Only 1 node left, skip the loop
!         ENDIF
!
!
!         IF(.NOT.ASSOCIATED(Plume1d%next))THEN ! delete the tail node
!            Plume1d_tail%next => Plume1d_prev
!            Plume1d_tail => Plume1d_tail%next
!
!!            Plume1d => Plume1d_prev%next !!! ???
!            IF(ASSOCIATED(Plume1d_prev%next)) NULLIFY(Plume1d_prev%next)
!            DEALLOCATE(Plume1d%CONCNT1d)
!            DEALLOCATE(Plume1d)
!            GOTO 500 ! delete the tail node, skip the loop
!         ELSEIF(ASSOCIATED(Plume1d_prev))THEN ! delete node, not head/tail
!            Plume1d => Plume1d_prev%next
!            Plume1d_prev%next => Plume1d%next
!            DEALLOCATE(Plume1d%CONCNT1d)
!            DEALLOCATE(Plume1d)
!            Plume1d => Plume1d_prev%next
!         ELSE ! delete head node
!            Plume1d => Plume1d_head
!            Plume1d_head => Plume1d%next
!            DEALLOCATE(Plume1d%CONCNT1d)
!            DEALLOCATE(Plume1d)
!            Plume1d => Plume1d_head
!         ENDIF
!
!           
!
!         WRITE(6,*)'333'
!
!
!         GOTO 200 ! skip this box, go to next box
!
!       ENDIF ! IF(MONTH==2 .AND. DAY==31)THEN
!       ENDIF ! IF(HOUR==23 .AND. MINUTE==50)THEN

       !===================================================================
       ! Fourth judge:
       ! If the plume touch the tropopause, dissolve the plume
       !===================================================================
!       IF ( box_lev>State_Met%TROPP(i_lon,i_lat) ) THEN
!
!
!         Num_dissolve = Num_dissolve+1
!
!         WRITE(6,*)'DISSOLVE 4:', i_box
!
!
!         DO i_species=1,n_species,1
!         i_advect = id_PASV_LA +i_species -1
!
!         backgrd_concnt = &
!                  State_Chm%Species(i_lon,i_lat,i_lev,i_advect)
!
!         backgrd_concnt = &
!           ( SUM(box_concnt_1D(1:n_slab_max,i_species)) * V_grid_1D   &
!                + backgrd_concnt*grid_volumn ) / grid_volumn
!
!         State_Chm%Species(i_lon,i_lat,i_lev,i_advect) = backgrd_concnt
!         ENDDO
!
!
!         WRITE(6,*) '*** Dissolved in 4th Judge ***', i_box
!
!
!
!         IF(.NOT.ASSOCIATED(Plume1d%next) .and. &
!                                .NOT.ASSOCIATED(Plume1d_prev))THEN
!
!            WRITE(6,*)'                 '
!            WRITE(6,*)'*** Only left the last Plume1d', Num_Plume1d
!            WRITE(6,*)'                 '
!            DEALLOCATE(Plume1d_head)
!
!            i_box = i_box-1
!            Num_Plume1d = Num_Plume1d - 1
!
!            GOTO 500 ! Only 1 plume1d left, skip loop
!
!         ENDIF
!
!
!         IF(.NOT.ASSOCIATED(Plume1d%next))THEN ! delete last node
!            Plume1d_tail%next => Plume1d_prev
!            Plume1d_tail => Plume1d_tail%next
!
!!            Plume1d => Plume1d_prev%next !!! ???
!            IF(ASSOCIATED(Plume1d_prev%next)) NULLIFY(Plume1d_prev%next)
!            DEALLOCATE(Plume1d%CONCNT1d)
!            DEALLOCATE(Plume1d)
!            i_box = i_box-1
!            Num_Plume1d = Num_Plume1d - 1
!            GOTO 500
!         ELSEIF(ASSOCIATED(Plume1d_prev))THEN ! delete node, not head/tail
!            Plume1d => Plume1d_prev%next
!            Plume1d_prev%next => Plume1d%next
!            DEALLOCATE(Plume1d%CONCNT1d)
!            DEALLOCATE(Plume1d)
!            Plume1d => Plume1d_prev%next
!         ELSE ! delete head node
!            Plume1d => Plume1d_head
!            Plume1d_head => Plume1d%next
!            DEALLOCATE(Plume1d%CONCNT1d)
!            DEALLOCATE(Plume1d)
!            Plume1d => Plume1d_head
!         ENDIF
!
!         i_box = i_box-1
!         Num_Plume1d = Num_Plume1d - 1
!
!
!
!         GOTO 200 ! skip this box, go to next box
!
!       ENDIF ! IF(ABS(Prod_before-Prod_after)/Prod_after<0.01)THEN


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



!       IF(mod(tt*NINT(Dt),24*60*60)==0)THEN   ! output once every day (24 hours)
!
!         OPEN( 485,      FILE=TRIM( FileConcnt   ), STATUS='OLD',  &
!           POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!         OPEN( 486,      FILE=TRIM( FileXYZ   ), STATUS='OLD',  &
!           POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!
!         Max_Concnt_Grid = MAXVAL(State_Chm%Species(:,:,:,id_PASV_LA))
!
!         DO i_slab = 1, n_slab_max, 1
!
!           Entropy_Concnt = box_concnt_1D(i_slab,i_tracer) + &
!                        State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA)
!
!           tracer_mol = Entropy_Concnt * V_grid_1D / AVO
!           air_mol    = V_grid_1D/1e+6_fp*State_Met%AIRDEN(i_lon,i_lat,i_lev)*1000/AIRMW
!           mix_ratio  = tracer_mol/air_mol
!              
!           IF(mix_ratio>0) Entropy = Entropy - BOLTZ*tracer_mol*log(mix_ratio)
!
!           i_cell = (i_lev-1)*IIPAR*JJPAR+(i_lat-1)*IIPAR+i_lon
!           Entropy_V(i_cell) = Entropy_V(i_cell) - V_grid_1D
!
!
!
!           ! Only keep the plume grid cell that has a larger concentration 
!           ! than the background concentration
!           IF( box_concnt_1D(i_slab,i_tracer) > 100.0*Max_Concnt_Grid )THEN
!
!              WRITE(485,*) box_concnt_1D(i_slab,i_tracer), &
!                              V_grid_1D, State_Met%AIRDEN(i_lon,i_lat,i_lev)
!
!              WRITE(486,*) i_lon,i_lat,i_lev
!
!           ELSE
!
!              backgrd_concnt = State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA4)
!
!              backgrd_concnt = &
!                   ( box_concnt_1D(i_slab,i_tracer)*V_grid_1D  &
!                        + backgrd_concnt*grid_volumn ) / grid_volumn
!
!              State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA4) = backgrd_concnt
!
!           ENDIF
!
!         ENDDO
!
!
!         CLOSE(485)
!         CLOSE(486)
!
!       ENDIF



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
!       backgrd_concnt = State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA3)
!
!       backgrd_concnt = &
!                     ( SUM(box_concnt_1D(1:n_slab_max,1)) * V_grid_1D  &
!                            + backgrd_concnt*grid_volumn ) / grid_volumn
!
!       State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA3) = backgrd_concnt


200     CONTINUE


     ENDDO ! DO WHILE(ASSOCIATED(Plume1d))


500     CONTINUE

!    WRITE(6,*)'=== loop for 1d plume in plume_run: ', i_box, Num_Plume1d




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

        Entropy_Concnt = box_concnt_2D(i_x,i_y,i_entropy) + &
                        State_Chm%Species(i_lon,i_lat,i_lev,id_GC(i_entropy))

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


!    WRITE(6,*)'Entropy1:', Entropy


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

        Entropy_Concnt = box_concnt_1D(i_slab,i_entropy) + &
                        State_Chm%Species(i_lon,i_lat,i_lev,id_GC(i_entropy))

        tracer_mol = Entropy_Concnt * V_grid_1D / AVO
        air_mol    = V_grid_1D/1e+6_fp*State_Met%AIRDEN(i_lon,i_lat,i_lev)*1000/AIRMW
        mix_ratio  = tracer_mol/air_mol

        IF(mix_ratio>0) Entropy = Entropy - BOLTZ*tracer_mol*log(mix_ratio)


!        IF(box_label==1)THEN
!          WRITE(6,*)"For box_label==1:"
!          WRITE(6,*) i_slab, box_concnt_1D(i_slab,i_tracer), tracer_mol, mix_ratio
!        ENDIF

!        IF(-BOLTZ*tracer_mol*log(mix_ratio)>1.0e+16)THEN
!          WRITE(6,*)"WARNING: super large Entropy happens!!!"
!          WRITE(6,*)air_mol, tracer_mol, mix_ratio
!        ENDIF


        i_cell = (i_lev-1)*IIPAR*JJPAR+(i_lat-1)*IIPAR+i_lon
        Entropy_V(i_cell) = Entropy_V(i_cell) - V_grid_1D



        ! Only output Day 1:
        IF(tt*NINT(Dt)==30*24*60*60 .AND. Calc_entropy==1)THEN
          WRITE(996,*)  1, tracer_mol, mix_ratio
        ENDIF




      ENDDO


!      IF(box_label==1)THEN
!         WRITE(6,*)"Total mass (molec) for box_label==1:", MAXVAL(box_concnt_1D(:,i_tracer))*V_grid_1D/AVO
!         WRITE(6,*) V_grid_1D, SUM(box_concnt_1D(:,i_tracer) )*V_grid_1D, &
!             State_Chm%Species(i_lon,i_lat,i_lev,id_PASV_LA)*V_grid_1D*n_slab_max
!      ENDIF


      Plume1d_prev => Plume1d
      Plume1d      => Plume1d%next


    ENDDO ! DO WHILE(ASSOCIATED(Plume1d))

    ENDIF ! IF(ASSOCIATED(Plume1d_head))THEN

!    WRITE(6,*)'Entropy2:', Entropy



    ! ------------------------------------------------------------------------------
    ! Third, for modified GEOS-Chem grid
    DO i_lon = 1, IIPAR
    DO i_lat = 1, JJPAR
    DO i_lev = 1, LLPAR

       i_cell = (i_lev-1)*IIPAR*JJPAR +(i_lat-1)*IIPAR +i_lon

       Entropy_Concnt = State_Chm%Species(i_lon,i_lat,i_lev,id_GC(i_entropy))

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
    IF(Num_Plume1d+Num_Plume2d==0) use_lagrange=88


    ! Everything is done, clean up pointers


    IF(ALLOCATED(box_lon_new))  deallocate(box_lon_new)
    IF(ALLOCATED(box_lat_new))  deallocate(box_lat_new)

    IF(ASSOCIATED(u))           nullify(u)
    IF(ASSOCIATED(v))           nullify(v)
    IF(ASSOCIATED(omeg))        nullify(omeg)

    IF(ASSOCIATED(Ptemp))       nullify(Ptemp)
    IF(ASSOCIATED(P_BXHEIGHT))  nullify(P_BXHEIGHT)

    IF(ASSOCIATED(X_edge))      nullify(X_edge)
    IF(ASSOCIATED(Y_edge))      nullify(Y_edge)


    IF(ASSOCIATED(Plume2d_new))  nullify(Plume2d_new)
    IF(ASSOCIATED(Plume2d))      nullify(Plume2d)
    IF(ASSOCIATED(Plume2d_prev)) nullify(Plume2d_prev)
    IF(ASSOCIATED(Plume1d_new))  nullify(Plume1d_new)
    IF(ASSOCIATED(Plume1d))      nullify(Plume1d)
    IF(ASSOCIATED(Plume1d_prev)) nullify(Plume1d_prev)

  END SUBROUTINE plume_run



  ! ==================================================================
  ! Calculate the aerosol process for the plume segemnt
  ! Based on the Debra's AER 2-D aerosol model
  ! ==================================================================
  SUBROUTINE Plume_Aerosol(am_I_Root, State_Chm, State_Grid, &
                                        State_Met, Input_Opt, RC)

    USE Plume_list

    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Met_Mod,   ONLY : MetState
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Grid_Mod,  ONLY : GrdState

    USE TIME_MOD,        ONLY : GET_TS_CHEM
    USE TIME_MOD,        ONLY : GET_TS_DYN
    USE UnitConv_Mod,    ONLY : Convert_Spc_Units

    USE State_Chm_Mod,   ONLY : Ind_

    Use Sect_Aer_Mod,    Only : Box_Sect_Aer

    logical, intent(in)           :: am_I_Root
    TYPE(MetState), intent(in)    :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(GrdState), INTENT(in)    :: State_Grid  ! Grid State objectgg
    TYPE(OptInput), intent(in)    :: Input_Opt

    INTEGER,        INTENT(out)   :: RC         ! Success or failure

    INTEGER     :: id_H2O
    INTEGER     :: i_box
    integer     :: i_x, i_y
    integer     :: i_lon, i_lat, i_lev

    integer	:: i_slab
    integer     :: i_species
    integer     :: ts_sec, ts_coag

    integer	:: axis

    real(fp)    :: box_lon, box_lat, box_lev
    real(fp)    :: box_length, box_alpha, box_theta
    real(fp)    :: Pdx, Pdy
    real(fp)    :: box_Ra, box_Rb
    real(fp)    :: box_extra, box_life, box_label

    real(fp), dimension(:,:,:), allocatable :: box_concnt_2D
    real(fp), dimension(:,:),   allocatable :: box_concnt_1D

    real(fp), pointer           :: X_edge(:), Y_edge(:)
    real(fp)                    :: X_edge2, Y_edge2

    ! Aerosol properties taken from input
    real(fp)                    :: T_K_Vec, p_hPa_Vec, ndens_Vec, vvH2O_Vec
    real(fp)                    :: MW_kg, MW_g

    real(fp)                    :: box2d_aWP     (n_x_max, n_y_max, n_bins)
    real(fp)                    :: box2d_aDen    (n_x_max, n_y_max, n_bins)
    real(fp)                    :: box2d_Sfc_Ten (n_x_max, n_y_max, n_bins)

    real(fp)                    :: box1d_aWP     (n_slab_max, n_bins)
    real(fp)                    :: box1d_aDen    (n_slab_max, n_bins)
    real(fp)                    :: box1d_Sfc_Ten (n_slab_max, n_bins)

    real(fp)                    :: aWP_arr     (n_bins)
    real(fp)                    :: aDen_arr    (n_bins)
    real(fp)                    :: Sfc_Ten_arr (n_bins)
    real(fp)                    :: box_rWet    (n_bins)
    real(fp)                    :: vvSul_Arr   (n_bins)

    real(fp)                    :: vvH2SO4_Vec

    real(fp)			:: Concnt_total(n_species)


    TYPE(Plume2d_list), POINTER :: Plume2d
    TYPE(Plume1d_list), POINTER :: Plume1d


    CHARACTER(LEN=255)          :: ErrMsg
    Character(Len=255)          :: Src_Unit
    CHARACTER(LEN=255)          :: File580

    Logical                     :: LAER_nuc
    Logical                     :: LAER_grow
    Logical                     :: LAER_coag
    Logical                     :: LAER_coag_imp
    Logical                     :: LAER_sed


    allocate(box_concnt_2D(n_x_max, n_y_max, n_species))
    allocate(box_concnt_1D(n_slab_max, n_species))

    ! Assume success
    RC = 0
    ErrMsg      = '' ! for Convert_Spc_Units()


    P_edge => State_Met%PEDGE(1,1,:)  ! Wet air press @ level edges [hPa]

    X_edge => State_Grid%XEdge(:,1) !XEDGE(:,1,1)   ! IIPAR+1 ! new
    Y_edge => State_Grid%YEdge(1,:) !YEDGE(1,:,1)  
    ! Use second YEDGE, because sometimes YMID(2)-YMID(1) is not DLAT

    X_edge2       = X_edge(2)
    Y_edge2       = Y_edge(2)


    ! Debug options
    LAER_nuc  = Input_Opt%LAER_nuc
    LAER_grow = Input_Opt%LAER_grow
    LAER_coag = Input_Opt%LAER_coag
    LAER_Coag_Imp =  .False.


    ts_sec  = Real(GET_TS_CHEM())
    ts_coag = ts_sec


    ! Convert Eulerian data to VMR
    Call Convert_Spc_Units(Input_Opt,  State_Chm, &
                           State_Grid, State_Met, 'v/v dry', &
                           RC, Src_Unit)


    !----------------------------------------------------------------
    ! Loop for 2-D plume segment
    !----------------------------------------------------------------

    ! If no 2-D plume segment left, skip 2-D plume segment loop
    IF(.NOT.ASSOCIATED(Plume2d_head)) GOTO 610

    Plume2d   => Plume2d_head
    i_box     = 0


    DO WHILE(ASSOCIATED(Plume2d))

      i_box         = i_box+1

      WRITE(6,*)'shw: begin Plume2d aerosol:', i_box

      box_lon       = Plume2d%LON
      box_lat       = Plume2d%LAT
      box_lev       = Plume2d%LEV
      box_length    = Plume2d%LENGTH
      box_alpha     = Plume2d%ALPHA

      box_label     = Plume2d%label
      box_life      = Plume2d%LIFE

      Pdx           = Plume2d%PDX
      Pdy           = Plume2d%PDY

      box_rWet      = Plume2d%rWet
      box2d_aWP       = Plume2d%aWP
      box2d_aDen      = Plume2d%aDen
      box2d_Sfc_Ten   = Plume2d%Sfc_Ten

      box_concnt_2D = Plume2d%CONCNT2d


      !------------------------------------------------------------
      ! For the center of the plume
      !------------------------------------------------------------

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


      ! Find the location (GEOS-Chem grid cell) of plume segment
      i_lon = Find_iLonLat(box_lon, DX, X_edge2)
      if(i_lon>IIPAR) i_lon=i_lon-IIPAR
      if(i_lon<1) i_lon=i_lon+IIPAR

      i_lat = Find_iLonLat(box_lat, DY, Y_edge2)
      if(i_lat>JJPAR) i_lat=JJPAR
      if(i_lat<1) i_lat=1

      i_lev = Find_iPLev(box_lev,P_edge)
      if(i_lev>LLPAR) i_lev=LLPAR


      ! Get properties of location
      T_K_Vec      = State_Met%T(i_lon,i_lat,i_lev)
      p_hPa_Vec    = State_Met%PMID(i_lon,i_lat,i_lev)
      ndens_Vec    = State_Met%AIRNUMDEN(i_lon,i_lat,i_lev)

!      Spc_vv      = State_Chm%Species(i_lon,i_lat,i_lev,:)
      id_H2O       = Ind_('H2O')
      vvH2O_Vec    = State_Chm%Species(i_lon,i_lat,i_lev,id_H2O)



      WRITE(6,*)'shw, PiG, T, P, ndens, H2O:'
      WRITE(6,*)T_K_Vec, p_hPa_Vec, ndens_Vec, vvH2O_Vec




      ! output once every day (24 hours)
      IF(mod(tt*INT(ts_sec),60*60)==0)THEN

        File580   = 'Plume_total_concentration.txt'
        OPEN( 580,      FILE=TRIM( File580 ), STATUS='OLD',  &
             POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

        DO i_species=1,n_species
          Concnt_total(i_species) = SUM(box_concnt_2D(:,:,i_species))
        END DO

!        WRITE(580,*) Concnt_total(3:42)
        WRITE(580,'(40E15.3E4)') box_concnt_2D(n_x_mid, n_y_mid, 3:42)
        WRITE(6,*)'shw: Concnt:', box_concnt_2D(n_x_mid-1, n_y_mid-1, n_species-15)

      ENDIF




      ! change unit from [molec/cm3] to [v/v] for species in plume segments
      DO i_species=1,n_species

        MW_kg = State_Chm%SpcData(id_GC(i_species))%Info%emMW_g * 1.e-3_fp
        MW_g  = State_Chm%SpcData(id_GC(i_species))%Info%emMW_g

        box_concnt_2D(:,:,i_species) = box_concnt_2D(:,:,i_species)     &
                                     /  AVO * MW_kg                     &
                                     / (State_Met%AIRDEN(i_lon,i_lat,i_lev)/1e6) &
                                     * ( AIRMW / MW_g)
      END DO


      DO i_y = 1, n_y_max, 1
      DO i_x = 1, n_x_max, 1

        aWP_arr         = box2d_aWP(i_x,i_y,:)
        aDen_arr        = box2d_aDen(i_x,i_y,:)
        Sfc_Ten_arr     = box2d_Sfc_Ten(i_x,i_y,:)

        vvSul_Arr       = box_concnt_2D(i_x,i_y,3:42)! 40 bins
        vvH2SO4_Vec     = box_concnt_2D(i_x,i_y,2)! 1


 !!! shw
	IF(MAXVAL(vvSul_Arr)>=0.0)THEN
          CALL Box_Sect_Aer(aWP_arr,aDen_arr,&
                            vvSul_Arr,Sfc_Ten_arr,vvH2O_Vec,&
                            vvH2SO4_Vec,box_rWet,T_K_Vec,p_hPa_Vec,&
                            ndens_Vec,ts_sec,ts_coag,&
                            LAER_Nuc,LAER_Grow,LAER_Coag,&
                            LAER_Coag_Imp, i_y, i_x, RC)
	ENDIF


        box2d_aWP(i_x,i_y,:)         = aWP_arr
        box2d_aDen(i_x,i_y,:)        = aDen_arr
        box2d_Sfc_Ten(i_x,i_y,:)     = Sfc_Ten_arr

        box_concnt_2D(i_x,i_y,3:42)   = vvSul_Arr     ! 40 bins
        box_concnt_2D(i_x,i_y,2)      = vvH2SO4_Vec   ! 1

      ENDDO
      ENDDO



      ! change unit from [v/v] to [molec/cm3] for species in plume segments
      DO i_species=1,n_species

        MW_kg = State_Chm%SpcData(id_GC(i_species))%Info%emMW_g * 1.e-3_fp
        MW_g  = State_Chm%SpcData(id_GC(i_species))%Info%emMW_g

        box_concnt_2D(:,:,i_species) = box_concnt_2D(:,:,i_species)     &
                                     *  AVO / MW_kg                     &
                                     * (State_Met%AIRDEN(i_lon,i_lat,i_lev)/1e6) &
                                     / ( AIRMW / MW_g )
      ENDDO



      ! -----------------------------------------------------------
      ! update     
      ! -----------------------------------------------------------
      Plume2d%LON       = box_lon
      Plume2d%LAT       = box_lat
      Plume2d%LEV       = box_lev
      Plume2d%LENGTH    = box_length
      Plume2d%ALPHA     = box_alpha
      Plume2d%label     = box_label
      Plume2d%LIFE      = box_life


      Plume2d%PDX       = Pdx
      Plume2d%PDY       = Pdy

      Plume2d%rWet      = box_rWet
      Plume2d%aWP       = box2d_aWP
      Plume2d%aDen      = box2d_aDen
      Plume2d%Sfc_Ten   = box2d_Sfc_Ten

      Plume2d%CONCNT2d  = box_concnt_2D

      Plume2d           => Plume2d%next

    ENDDO  ! DO WHILE(ASSOCIATED(Plume2d))

    610 CONTINUE




    !----------------------------------------------------------------
    ! Loop for 1-D plume segment
    !----------------------------------------------------------------

    IF(.NOT.ASSOCIATED(Plume1d_head)) GOTO 620

    Plume1d => Plume1d_head
    i_box = 0

    DO WHILE(ASSOCIATED(Plume1d))

      i_box = i_box+1

      WRITE(6,*)'shw: begin Plume1d aerosol:', i_box

      box_lon       = Plume1d%LON
      box_lat       = Plume1d%LAT
      box_lev       = Plume1d%LEV
      box_length    = Plume1d%LENGTH
      box_alpha     = Plume1d%ALPHA
      box_label     = Plume1d%label
      box_life      = Plume1d%LIFE

      box_Ra        = Plume1d%RA
      box_Rb        = Plume1d%RB

      box_rWet        = Plume1d%rWet
      box1d_aWP       = Plume1d%aWP
      box1d_aDen      = Plume1d%aDen
      box1d_Sfc_Ten   = Plume1d%Sfc_Ten

      box_concnt_1D = Plume1d%CONCNT1d



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


      ! Find the location (GEOS-Chem grid cell) of plume segment
      i_lon = Find_iLonLat(box_lon, DX, X_edge2)
      if(i_lon>IIPAR) i_lon=i_lon-IIPAR
      if(i_lon<1) i_lon=i_lon+IIPAR

      i_lat = Find_iLonLat(box_lat, DY, Y_edge2)
      if(i_lat>JJPAR) i_lat=JJPAR
      if(i_lat<1) i_lat=1

      i_lev = Find_iPLev(box_lev,P_edge)
      if(i_lev>LLPAR) i_lev=LLPAR


      ! Get properties of location
      T_K_Vec      = State_Met%T(i_lon,i_lat,i_lev)
      p_hPa_Vec    = State_Met%PMID(i_lon,i_lat,i_lev)
      ndens_Vec    = State_Met%AIRNUMDEN(i_lon,i_lat,i_lev)

!      Spc_vv      = State_Chm%Species(i_lon,i_lat,i_lev,:)
      id_H2O       = Ind_('H2O')
      vvH2O_Vec    = State_Chm%Species(i_lon,i_lat,i_lev,id_H2O)


      ! change unit from [molec/cm3] to [v/v] for species in plume segments
      DO i_species=1,n_species

        MW_kg = State_Chm%SpcData(id_GC(i_species))%Info%emMW_g * 1.e-3_fp
        MW_g  = State_Chm%SpcData(id_GC(i_species))%Info%emMW_g

        box_concnt_1D(:,i_species) = box_concnt_1D(:,i_species)     &
                                     /  AVO * MW_kg                     &
                                     / (State_Met%AIRDEN(i_lon,i_lat,i_lev)/1e6) &
                                     * ( AIRMW / MW_g)
      END DO


      DO i_slab = 1, n_slab_max, 1

        aWP_arr         = box1d_aWP(i_slab,:)
        aDen_arr        = box1d_aDen(i_slab,:)
        Sfc_Ten_arr     = box1d_Sfc_Ten(i_slab,:)

        vvSul_Arr       = box_concnt_1D(i_slab,3:42)! 40 bins
        vvH2SO4_Vec     = box_concnt_1D(i_slab,2)! 1

! !!! shw
!        CALL Box_Sect_Aer(aWP_arr,aDen_arr,&
!                          vvSul_Arr,Sfc_Ten_arr,vvH2O_Vec,&
!                          vvH2SO4_Vec,box_rWet,T_K_Vec,p_hPa_Vec,&
!                          ndens_Vec,ts_sec,ts_coag,&
!                          LAER_Nuc,LAER_Grow,LAER_Coag,&
!                          LAER_Coag_Imp, RC)


        box1d_aWP(i_slab,:)         = aWP_arr
        box1d_aDen(i_slab,:)        = aDen_arr
        box1d_Sfc_Ten(i_slab,:)     = Sfc_Ten_arr

        box_concnt_1D(i_slab,3:42)   = vvSul_Arr     ! 40 bins
        box_concnt_1D(i_slab,2)      = vvH2SO4_Vec   ! 1

      ENDDO




      Plume1d%LON       = box_lon
      Plume1d%LAT       = box_lat
      Plume1d%LEV       = box_lev
      Plume1d%LENGTH    = box_length
      Plume1d%ALPHA     = box_alpha
      Plume1d%label     = box_label
      Plume1d%LIFE      = box_life
      Plume1d%RA        = box_Ra
      Plume1d%RB        = box_Rb

      Plume1d%rWet      = box_rWet
      Plume1d%aWP       = box1d_aWP
      Plume1d%aDen      = box1d_aDen
      Plume1d%Sfc_Ten   = box1d_Sfc_Ten

      Plume1d%CONCNT1d  = box_concnt_1D


      Plume1d   => Plume1d%next

    ENDDO  ! DO WHILE(ASSOCIATED(Plume1d))


  620 CONTINUE


    ! Convert back to whatever we had before
    Call Convert_Spc_Units(Input_Opt,  State_Chm, &
                           State_Grid, State_Met, Src_Unit,  &
                           RC)




  END SUBROUTINE Plume_Aerosol


  !------------------------------------------------------------------
  ! Aerosol box model for one plume segment
  !------------------------------------------------------------------
!  SUBROUTINE Cal_Sect_Aer(aWP_Arr,aDen_Arr,&
!                         vvSO4_Arr,Sfc_Ten_Arr,vvH2O_Vec,&
!                         vvH2SO4_Vec,rWet_Arr,T_K_Vec,p_hPa_Vec,&
!                         ndens_Vec,ts_sec,ts_coag,&
!                         LAER_Nuc,LAER_Grow,LAER_Coag,&
!                         LAER_Coag_Imp,RC)
!
!  !
!  ! !USES:
!  !
!    ! For GEOS-Chem Precision (fp)
!    USE Precision_Mod,        ONLY : fp, &
!                                     dp=>f8
!
!    ! Data which will likely go to State_Chm
!    USE Sect_Aer_Data_Mod,    ONLY : n_aer_bin, &
!                                     aer_mass,  & ! Mass per particle in each bin
!                                     aer_molec, & ! Molec per particle in each bin
!
!    ! Arrays which should be put into State_Chm
!    USE Sect_Aer_Data_Mod,    ONLY : bvp,       & ! Eq. vapor pressure of H2SO4
!                                     st,        & ! Surface tension
!                                     airvdold,  & ! ???
!                                     ck           ! Coagulation kernel
!
!    ! Constants
!    USE PhysConstants,        ONLY : mw=>h2omw, &
!                                     av=>avo,   &
!                                     pi
!
!    USE Sect_Aer_Data_Mod,    ONLY : Ma=>aer_mwh2so4, &
!                                     den_h2so4,       &
!                                     Rgas=>aer_Rgas,  &
!                                     akb=>aer_akb,    &
!                                     aer_dry_rad,     &
!                                     aer_Vrat
!
!
!    USE sect_aer_mod,         ONLY : AER_wpden
!    USE sect_aer_mod,         ONLY : Cal_Coag_Kernel
!    USE sect_aer_mod,         ONLY : Cal_Coag_Kernel_Implicit
!    USE sect_aer_mod,         ONLY : AER_grow
!    USE sect_aer_mod,         ONLY : AER_nucleation
!    USE sect_aer_mod,         ONLY : AER_coagulation
!    USE sect_aer_mod,         ONLY : AER_coagulation_Implicit
!    USE sect_aer_mod,         ONLY : AER_mass_check
!
!
!    real(fp), intent(inout) :: aWP_Arr     (n_aer_bin)
!    real(fp), intent(inout) :: aDen_Arr    (n_aer_bin)
!    real(fp), intent(inout) :: vvSO4_Arr   (n_aer_bin)
!    real(fp), intent(inout) :: rWet_Arr    (n_aer_bin)
!    real(fp), intent(inout) :: Sfc_Ten_Arr (n_aer_bin)
!    real(fp), intent(inout) :: vvH2O_Vec   
!    real(fp), intent(inout) :: vvH2SO4_Vec 
!    real(fp), intent(in)    :: T_K_Vec    
!    real(fp), intent(in)    :: p_hPa_Vec  
!    real(fp), intent(in)    :: ndens_Vec 
!    integer, intent(in)     :: ts_sec
!    integer, intent(in)     :: ts_coag
!    logical, intent(in)     :: LAer_Coag_Imp
!    integer, intent(out)    :: RC
!
!
!    ! Atmospheric conditions
!    real(fp)                :: T_K, p_hPa, air_dens
!
!    ! Aerosol properties taken from input
!    real(fp)                :: aWP_Box      (n_aer_bin)
!    real(fp)                :: aDen_Box     (n_aer_bin)
!    real(fp)                :: vvSO4_Box    (n_aer_bin)
!    real(fp)                :: ST_Box       (n_aer_bin)
!    real(fp)                :: vvSO4_Box_0  (n_aer_bin)
!
!    real(fp)                :: vvH2O
!    real(fp)                :: vvH2SO4
!    real(fp)                :: vvH2SO4_0
!
!    ! Aerosol properties calculated during operation
!    real(fp)                :: BVP_Box (n_aer_bin)
!    real(fp)                :: rWet    (n_aer_bin)   
!
!    ! Coagulation kernel
!    real(fp), pointer       :: CK_Box(:,:)
!    real(fp)                :: fijk_Box(n_aer_bin,n_aer_bin,n_aer_bin)
!
!    ! Loop indices 
!    integer                 :: K,N
!    integer                 :: I_Bin
!
!    ! For time step control
!    Integer                 :: NAStep
!    Real(fp)                :: dt_substep
!
!    ! Diagnostics
!    Real(fp)                :: box_grow_d, box_nrate_d, box_wp_d
!    Real(fp)                :: S_Start, S_End
!
!    ! Logicals
!    logical, intent(in)     :: LAer_Nuc
!    logical, intent(in)     :: LAer_Grow
!    logical, intent(in)     :: LAer_Coag
!    logical                 :: Run_Micro
!
!
!
!    Run_Micro = (LAer_Nuc.or.LAer_Grow.or.LAer_Coag)
!
!    ! How many substeps?
!    NAStep = NINT(dble(ts_sec)/dble(ts_coag))
!    dt_substep = dble(ts_sec)/dble(NAStep)
!
!
!    ! Copy properties into temporary variables
!    vvH2O        = vvH2O_Vec 
!    vvH2SO4      = vvH2SO4_Vec
!    ST_Box(:)    = Sfc_Ten_Arr(:)
!    aWP_Box(:)   = aWP_Arr    (:)
!    aDen_Box(:)  = aDen_Arr   (:)
!    vvSO4_Box(:) = vvSO4_Arr  (:)
!
!    ! Store for later mass check
!    vvSO4_Box_0(:) = vvSO4_Box(:)
!    vvH2SO4_0      = vvH2SO4
!
!    ! Get properties of location
!    T_K          = T_K_Vec
!    P_hPa        = p_hPa_Vec
!
!    ! For safety's sake (these will be recalculated)
!    BVP_Box(:)   = 0.0e+0_fp
!    rWet(:)      = 0.0e+0_fp
!
!    ! Calculate updated weight % H2SO4 and density
!    Call AER_wpden(aWP_Box, aDen_Box, ST_Box, BVP_Box,&
!                   T_K,     P_hPa,    vvH2O,  rWet)
!
!    ! If we aren't running a microphysics step, stop here
!    If (Run_Micro) Then
!       ! Get more met properties
!       air_dens = ndens_vec
!
!       ! Calculate coagulation kernel
!       CK_Box => ck(1,1,:,:)
!       If (LAER_coag) &
!       Call Cal_Coag_Kernel(CK_Box,aWP_Box,aDen_Box,T_K,P_hPa)
!
!       ! If using implicit coagulation..
!       fijk_Box = 0.0e+0_fp
!       If (LAER_coag.and.LAER_Coag_Imp) &
!       Call Cal_Coag_Kernel_Implicit(fijk_Box,aWP_Box,aDen_Box)
!
!          ! Run microphysical processes
!          box_grow_d  = 0.0e+0_fp
!          box_nrate_d = 0.0e+0_fp
!
!          Do K=1,NAStep
!             If (LAER_grow) &
!             Call AER_grow(T_K,       air_dens,   dt_substep, aWP_Box, &
!                           aDen_box,  ST_Box,     BVP_Box,    vvH2So4, &
!                           vvSO4_Box, box_grow_d)
!             If (LAER_nuc) &
!             Call AER_nucleation(T_K,      P_hPa,       air_dens, dt_substep, &
!                                 aWP_Box,  aDen_Box,    vvH2O,    vvH2SO4,    &
!                                 vvSO4_Box,box_nrate_d, box_wp_d              )
!             If (LAER_coag.and.(.not.LAER_Coag_Imp)) &
!             Call AER_coagulation(vvSO4_Box,CK_box,air_dens,dt_substep)
!             If (LAER_coag.and.LAER_Coag_Imp) &
!             Call AER_Coagulation_Implicit(vvSO4_Box,fijk_Box,CK_box,aWP_Box, &
!                                           aDen_Box,air_dens,dt_substep       )
!          End Do
!          !aero_grow_d(I,J,L)  = box_grow_d/real(NAStep)
!          !aero_nrate_d(I,J,L) = box_nrate_d/real(NAStep)
!          CK_Box => Null()
!
!          ! Make sure that the mass hasn't drifted too far 
!          Call AER_mass_check(vvH2SO4, vvSO4_Box,  vvH2SO4_0, vvSO4_Box_0, &
!                              air_dens, s_start, s_end,        RC          )
!          If (RC.ne.0) Then
!             Call Error_Stop('Failure in strat aerosol mass balance', &
!                             'sect_aer_mod.F90')
!          End If
!          !Call AER_Washout(?,?,?,?)
!       End If
!
!       ! Copy results back out
!       vvH2O_Vec      = vvH2O
!       vvH2SO4_Vec    = vvH2SO4
!       Sfc_Ten_Arr(:) = ST_Box(:)
!       aWP_Arr    (:) = aWP_Box(:)
!       aDen_Arr   (:) = aDen_Box(:)
!       vvSO4_Arr  (:) = vvSO4_Box(:)
!       rWet_Arr   (:) = rWet(:)
!
!
!  END SUBROUTINE Cal_Sect_Aer


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

    ! identify the grid point located in the southwest of the particle or under
    ! the particle ??? 
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



!  REAL(fp) FUNCTION Find_theta(theta1, concnt1_2D, Height, Ibox)
!    ! n_x_max, n_y_max, Pdx, Pdy are global variables
!
!    IMPLICIT NONE
!
!    REAL(fp)    :: theta1, concnt1_2D(n_x_max, n_y_max), Height    
!    REAL(fp)    :: R
!    REAL(fp)    :: test(100), X(100), Y(100), Conc(100)
!
!    REAL(fp)    :: angle(1)
!
!    INTEGER     :: i, Ibox
!
!      R = Height/COS(theta1) ! 2.5 is a random number
!
!
!      DO i= 1, 100
!        test(i) = PI/2 - (PI/2-theta1)*2/100*i
!        Y(i)    = R *COS( test(i) )
!        X(i)    = R *SIN( test(i) )
!        Conc(i) = Interplt_2D(X(i), Y(i), concnt1_2D, Ibox, 1)
!      ENDDO
!
!
!      angle      = test( MAXLOC(Conc(:)) )
!      Find_theta = angle(1)
!
!    return 
!
!  END FUNCTION



  SUBROUTINE Slab_init(Pdx, Pdy, theta1, Pc_2D, Height1, box_Ra, box_Rb, Cslab)
 
    IMPLICIT NONE

    REAL(fp)    :: Pdx, Pdy, theta1, Height1
    REAL(fp), INTENT(INOUT)  :: box_Ra, box_Rb
    REAL(fp), INTENT(INOUT)  :: Cslab(n_slab_max)
    REAL(fp)    :: Pc_2D(n_x_max,n_y_max) !, Ec_2D(n_x_max,n_y_max)


    INTEGER, parameter     :: Nb = n_slab_max
    INTEGER, parameter     :: Na = 128

    INTEGER     :: Nb_mid, Na_mid

    REAL(fp)    :: X2d(Na,Nb), Y2d(Na,Nb), C2d(Na,Nb) !, Extra_C2d(Na,Nb)

    REAL(fp)    :: LenB, LenA
    REAL(fp)    :: Adx, Ady, Bdx, Bdy
    REAL(fp)    :: Prod, M, Lb, La

    REAL(fp)    :: C_slab(Nb) !, Extra_slab(Nb)

    real(fp)  :: start, finish

    INTEGER     :: i, j

      Nb_mid = INT(Nb/2)
      Na_mid = INT(Na/2)


      LenB = Pdy ! 8.0 *Height1 / Nb /SIN(theta1)
      LenA = 1.3 *Height1 *TAN(theta1) /Na /SIN(theta1)


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


!      call cpu_time(start)

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( i, j )
      DO i=1,Na,1
      DO j=1,Nb,1
        C2d(i,j)       = Interplt_2D(Pdx, Pdy, X2d(i,j), Y2d(i,j), Pc_2D, 2)
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

!      call cpu_time(finish)
!      WRITE(6,*)'Time2 (finish-start):', finish-start


      ! define the length/width of slab based on a seconde order chamical reaction

!      Prod = SUM( C2d(Na_mid,:)**2 ) *(LenA*LenB)
!      M    = SUM( C2d(Na_mid,:) ) *(LenA*LenB)
!      Lb   = LenB
!      La   = M**2 /Prod / Lb

!      DO i=1,Nb,1
!        C_slab(i)     = SUM(C2d(:,i)) *(LenA*LenB)/ (La*Lb)
!      ENDDO

        ! define the length/width of slab based on 95% total mass
        Lb = LenB
        La = Height1 *TAN(theta1)

        DO i=1,Nb,1
          C_slab(i)     = SUM(C2d(:,i)) *(LenA*LenB)/ (La*Lb)
        ENDDO


      ! Return:
      box_Ra = La
      box_Rb = Lb
      Cslab(1:n_slab_max) = C_slab(1:n_slab_max)

  END SUBROUTINE



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
!        IF(ids==1)WRITE(6,*)"*** ERROR: Find_theta() out bounds ***"

      ELSE IF((Iy0<1 .or. Iy0>=n_y_max))THEN
        Interplt_2D = 0.0
!        IF(ids==1)WRITE(6,*)"*** ERROR: Find_theta() out bounds***"

      ELSE
        C1 = Interplt_linear(x0, X1d(Ix0), X1d(Ix0+1), &
                C_2D(Ix0,Iy0),C_2D(Ix0+1,Iy0))

        C2 = Interplt_linear(x0, X1d(Ix0), X1d(Ix0+1), &
                C_2D(Ix0,Iy0+1), C_2D(Ix0+1,Iy0+1))

        Interplt_2D = Interplt_linear(y0, Y1d(Iy0), &
                                          Y1d(Iy0+1), C1, C2)

      ENDIF

    return

  END FUNCTION


  REAL(fp) FUNCTION Interplt_linear(xx, x1, x2, Cx1, Cx2)

    IMPLICIT NONE

    REAL(fp)    :: xx, x1, x2, Cx1, Cx2


    Interplt_linear = ( Cx1*(x2-xx) + Cx2*(xx-x1) ) /(x2-x1)

    return

  END FUNCTION



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


    TYPE(Plume2d_list), POINTER :: Plume2d_new, Plume2d, Plume2d_prev
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


    !-------------------------------------------------------------------
    ! Output plume location
    !-------------------------------------------------------------------
!    IF(mod(tt,6)==0)THEN     ! output once every hour
!    IF(mod(tt,14400000)==0)THEN   ! output once every day (24 hours)
!
!      FILENAME2   = 'Plume_concentration_molec_' // TRIM(ADJUSTL(YEAR_C)) // &
!         '-' //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) //      &
!         '-' //TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) //    &
!         ':' //TRIM(ADJUSTL(SECOND_C)) // '.txt'
!
!      OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='REPLACE', &
!            FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!      WRITE(262,*)'total plume number:', n_box
!
!
!      DO i_box = 1,n_box,1
!!      IF(i_box<3)THEN !!! shw
!
!        IF(Judge_plume(i_box)==0) THEN
!          WRITE(262,*) i_box, ': dissolved'
!        ENDIF
!
!        IF(Judge_plume(i_box)==2) THEN
!          WRITE(262,*) i_box, ': 2D model, total mass [molec]:'
!          WRITE(262,*) Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp &
!                * SUM(box_concnt_2D(:,:,N_species,i_box))
!        ENDIF
!
!        IF(Judge_plume(i_box)==1) THEN
!          WRITE(262,*) i_box, ': 1D model'
!          WRITE(262,*) SUM(box_concnt_1D(1:n_slab_max,N_species,i_box)) &
!                   *box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp
!        ENDIF
!
!        WRITE(262,*)' '
!
!!     ENDIF ! IF(i_box<3)THEN
!     ENDDO ! DO i_box = 1,n_box,1
!
!
!!       DO i_box = 1, n_box
!!          WRITE(262,*) i_box, SUM(box_concnt(i_box,:,N_species)*V_slab(i_box,:))
!!       ENDDO
!
!!        WRITE(262,*) box_concnt(1,:,N_species)
!
!    ENDIF ! IF(mod(tt,144)==0)THEN
!!    ENDIF
!
!

    CLOSE(261)

    IF(ASSOCIATED(Plume2d_new)) nullify(Plume2d_new)
    IF(ASSOCIATED(Plume2d)) nullify(Plume2d)
    IF(ASSOCIATED(Plume2d_prev)) nullify(Plume2d_prev)
    IF(ASSOCIATED(Plume1d_new)) nullify(Plume1d_new)
    IF(ASSOCIATED(Plume1d)) nullify(Plume1d)
    IF(ASSOCIATED(Plume1d_prev)) nullify(Plume1d_prev)


  END SUBROUTINE lagrange_write_std

!-------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------

  subroutine lagrange_cleanup()


!!! dissolve all the alive plume here ???


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

END FUNCTION




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

END FUNCTION
