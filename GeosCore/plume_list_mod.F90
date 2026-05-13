MODULE Plume_list_mod

  USE PRECISION_MOD
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Plume2d_list, Plume1d_list

  TYPE :: Plume2d_list

    INTEGER :: IsNew     = MISSING_INT ! 1: the plume is new injected
    INTEGER :: label     = MISSING_INT! injected rank
    INTEGER :: lon_ind = MISSING_INT
    INTEGER :: lat_ind = MISSING_INT
    INTEGER :: lev_ind = MISSING_INT

    REAL(fp) :: LON = MISSING, LAT = MISSING, LEV = MISSING
    REAL(fp) :: LENGTH = MISSING, ALPHA = MISSING
    REAL(fp) :: LIFE = MISSING
    REAL(fp) :: PDX = MISSING, PDY = MISSING
    REAL(fp), DIMENSION(:,:,:), ALLOCATABLE :: CONCNT2d  ! [n_x_max, n_y_max, n_species]
    REAL(fp), DIMENSION(:), ALLOCATABLE :: MassRef2d  ! [n_species]

    TYPE(Plume2d_list), POINTER :: next => NULL()
  END TYPE


  TYPE :: Plume1d_list
    INTEGER :: label         = MISSING_INT ! injected rank
    INTEGER :: Is_transfer   = MISSING_INT ! if =1, transfer to the host Euletian model
    INTEGER :: lon_ind = MISSING_INT
    INTEGER :: lat_ind = MISSING_INT
    INTEGER :: lev_ind = MISSING_INT

    REAL(fp) :: LON = MISSING, LAT = MISSING, LEV = MISSING
    REAL(fp) :: LENGTH = MISSING, ALPHA = MISSING
    REAL(fp) :: LIFE = MISSING
    REAL(fp) :: RA = MISSING, RB = MISSING, THETA = MISSING
    REAL(fp), DIMENSION(:,:), ALLOCATABLE :: CONCNT1d ! [n_slab_max,n_species]
    REAL(fp), DIMENSION(:), ALLOCATABLE :: MassRef1d  ! [n_species]

    TYPE(Plume1d_list), POINTER :: next => NULL()
  END TYPE

END MODULE Plume_list_mod