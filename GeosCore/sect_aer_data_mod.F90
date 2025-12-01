module sect_aer_data_mod

  ! TEMPORARY - data module. Will eventually get moved to state_chm

  use precision_mod,    only: dp=>f8
  use precision_mod,    only: fp
  !use mo_socol_tracers, only: n_aer_bin, aer_bin
  use physconstants,    only: aer_pi=>pi
  use physconstants,    only: aer_av=>avo
  use physconstants,    only: aer_mwh2o=>h2omw

#if !defined( MDL_BOX )
  USE Input_Opt_Mod,    ONLY : OptInput
  USE State_Chm_Mod,    ONLY : ChmState
  USE State_Grid_Mod,   ONLY : GrdState
  USE Species_Mod                        ! For species database object
#endif

  implicit none
  public 

  public :: aer_allocate_ini
  public :: aer_cleanup

  real(dp), allocatable, target :: ck(:,:,:,:,:)    ! coagulation kernel
  real(dp), allocatable, target :: fijk(:,:,:)           ! volume fractions for coagulation scheme
  !!real(dp), allocatable, save :: aDen(:,:,:,:), aWP(:,:,:,:)   ! aerosol density, weight percent
  real(dp), parameter :: den_h2so4=1.8E-12_dp ! pure h2so4 density in g/um^3, used for calculating h2so4 mass/particle !eth_af_dryS
  real(dp), parameter :: aer_R0=3.9376E-4_dp !smallest bin's dry sulfate radius in um
  !real(dp), parameter ::  aer_Vrat=2.0_dp !multiplication factor between neighbouring bins
  real(dp) :: aer_Vrat !multiplication factor between neighbouring bins
  real(dp), allocatable :: aer_mass(:),aer_molec(:),aer_dry_rad(:)
  real(dp), allocatable, target :: st(:,:,:,:), bvp(:,:,:,:)
  real(dp), allocatable, target :: airvdold(:,:,:,:)

  ! contants used in AER module
  real(dp), parameter :: aer_akb     = 1.380662E-16_dp  ! [g cm2 s-2 K-1 molecule-1]
  real(dp), parameter :: aer_Rgas    = 8.3145E7_dp      ! [g cm2 s-2 mole-1 K-1]
  real(dp), parameter :: aer_mwh2so4 = 98.076_dp
  integer :: i_h2o, i_h2so4
  
  ! Variables for band_echam from BP Luo
  integer, parameter :: NT=8,Nwt=31, NR=40
  real(dp), save :: g3(nt,nwt,Nr,22), o3(nt,nwt,Nr,22), ex3(nt,nwt,Nr,22)

  integer :: n_aer_bin

  ! Number of microphysics substeps
  integer, parameter :: nastep = 1

contains

#if defined( MDL_BOX )
  subroutine AER_allocate_ini(n_boxes, RC)
    integer, intent(in)  :: n_boxes
    integer, intent(out) :: RC
#else
  subroutine AER_allocate_ini(Input_Opt,State_Chm,State_Grid,RC)

    TYPE(OptInput), INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm  ! Chemistry State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid ! Grid description
    integer, intent(out)        :: RC        

    TYPE(Species), POINTER      :: ThisSpc
#endif

    integer                     :: n_bins
    integer                     :: k, n
    integer                     :: it,iwt,ir,ioerror
    real                        :: tt, wtss

    integer                     :: nx, ny, nz

    integer                     :: as

    ! Assume success
    rc = 0

#if defined( MDL_BOX )
    nx = n_boxes
    ny = 1
    nz = 1
#else
    nx = State_Grid%NX
    ny = State_Grid%NY
    nz = State_Grid%NZ
#endif

    ! Needed for coagulation kernel
    allocate(ck(nx,ny,nz,n_aer_bin,n_aer_bin),stat=as)
    if (as.ne.0) then
       rc = -1
       return
    end if
    ! Constant - mass of SO4 per particle in each bin
    allocate(aer_mass(n_aer_bin),stat=as)
    if (as.ne.0) then
       rc = -1
       return
    end if
    ! Constant - # of SO4 molecules per particle in each bin
    allocate(aer_molec(n_aer_bin),stat=as)
    if (as.ne.0) then
       rc = -1
       return
    end if
    ! Constant - dry radius of each bin in um
    allocate(aer_dry_rad(n_aer_bin),stat=as)
    if (as.ne.0) then
       rc = -1
       return
    end if
    ! Needed for coagulation?
    allocate(fijk(n_aer_bin,n_aer_bin,n_aer_bin),stat=as)
    if (as.ne.0) then
       rc = -1
       return
    end if
    ! Surface tension
    allocate(st  (nx,ny,nz,n_aer_bin),stat=as)
    if (as.ne.0) then
       rc = -1
       return
    end if
    ! Equilibrium water vapor pressure
    allocate(bvp (nx,ny,nz,n_aer_bin),stat=as)
    if (as.ne.0) then
       rc = -1
       return
    end if
    allocate(airvdold(nx,ny,nz,1),stat=as)
    if (as.ne.0) then
       rc = -1
       return
    end if

    !!aden=0.
    !!awp=0.
    aer_mass=0.
    aer_molec=0.
    st=0.
    bvp=0.
    airvdold=0.

    ck = 0._dp
    fijk = 0._dp
    !calculate aerosol dry sulfate mass !eth_af_dryS

    if (n_aer_bin.eq.40) then
       ! Volume doubling
       aer_Vrat = 2.0e+0_dp
    !elseif (n_aer_bin.eq.150) then
    !   ! For consistency with AER code
    !   aer_Vrat = 1.2e+0_dp
    else
       ! Try to get the same range covered as the 40-bin model
       aer_Vrat = 2.0**(real(39,dp)/real(n_aer_bin-1,dp))
    end if
    
    write(6,*)'shw40: 40-bin radius [um]:'
    do k=1,n_aer_bin
      if (k.eq.1) then
         aer_mass(k)=den_h2so4*4./3.*aer_pi*aer_R0**3 !mass H2SO4/particle in g
         aer_dry_rad(k) = aer_R0
      else
         aer_mass(k)=aer_mass(k-1)*aer_Vrat !mass H2SO4/particle in g
         aer_dry_rad(k) = (3.0*aer_mass(k)/(4.0*aer_pi*den_h2so4))**(1.0/3.0)
      endif
      write(6,*)k, aer_dry_rad(k) !!! shw400
    enddo
    aer_molec=aer_mass/aer_mwh2so4*aer_av !molec H2SO4/particle 
    
    !ioerror=0
    !IF (p_parallel_io) THEN
    !! read band_echam
    !  open(100, file='band_echam5')
    !  ioerror=1
    !  DO IT=1,NT
    !     DO Iwt=1,Nwt
    !        read(100,*) tt,wtss
    !        DO ir=1,NR
    !        read(100,*) (ex3(IT,Iwt,ir,k),k=1,22), &
    !                    (o3(IT,Iwt,ir,k),k=1,22),  &
    !                    (g3(IT,Iwt,ir,k),k=1,22)               
    !        enddo
    !     enddo
    !  enddo
    !  o3=max(0.,o3)
    !  close(100)
    !  write(*,*) '*** Read radiation lookup table: done!'
    !ENDIF
    !IF (p_parallel) CALL p_bcast (ioerror, p_io)
    !IF (p_parallel) CALL p_bcast (ex3, p_io)
    !IF (p_parallel) CALL p_bcast (o3, p_io)
    !IF (p_parallel) CALL p_bcast (g3, p_io)   
    !if (ioerror==0) stop

    ! Success
    RC = 0

  end subroutine AER_allocate_ini

  subroutine AER_cleanup()

    If (allocated(ck))          deallocate(ck)
    If (allocated(aer_mass))    deallocate(aer_mass)
    If (allocated(aer_molec))   deallocate(aer_molec)
    If (allocated(aer_dry_rad)) deallocate(aer_dry_rad)
    If (allocated(fijk))        deallocate(fijk)
    If (allocated(st))          deallocate(st)
    If (allocated(bvp))         deallocate(bvp)
    If (allocated(airvdold))    deallocate(airvdold)
  end subroutine AER_cleanup
 

end module sect_aer_data_mod
