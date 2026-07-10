MODULE Plume_list_mod

  USE PRECISION_MOD
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Plume2d_list, Plume1d_list

  TYPE :: Plume2d_list

    LOGICAL                                     :: IsTransfer     = .False.        ! 1: Transfer to 1-D
    LOGICAL                                     :: IsDissolve     = .False.        ! 1: Dissolve to the host Eulerian model

    INTEGER                                     :: label          = MISSING_INT    
    INTEGER                                     :: lon_ind        = MISSING_INT
    INTEGER                                     :: lat_ind        = MISSING_INT
    INTEGER                                     :: lev_ind        = MISSING_INT

    REAL(fp)                                    :: LON            = MISSING 
    REAL(fp)                                    :: LAT            = MISSING
    REAL(fp)                                    :: LEV            = MISSING
    REAL(fp)                                    :: LENGTH         = MISSING
    REAL(fp)                                    :: ALPHA          = MISSING
    REAL(fp)                                    :: LIFE           = MISSING
    REAL(fp)                                    :: PDX            = MISSING
    REAL(fp)                                    :: PDY            = MISSING
    REAL(fp), DIMENSION(:,:,:), ALLOCATABLE     :: CONCNT2d                         ! [n_x_max, n_y_max, n_species]
    REAL(fp), DIMENSION(:), ALLOCATABLE         :: MassRef2d                        ! [n_species]

    TYPE(Plume2d_list), POINTER                 :: next => NULL()
    
  END TYPE


  TYPE :: Plume1d_list

    LOGICAL                                      :: IsDissolve    = .False.          ! 1: Dissolve to the host Eulerian model

    INTEGER                                      :: label         = MISSING_INT 
    INTEGER                                      :: lon_ind       = MISSING_INT
    INTEGER                                      :: lat_ind       = MISSING_INT
    INTEGER                                      :: lev_ind       = MISSING_INT

    REAL(fp)                                     :: LON           = MISSING
    REAL(fp)                                     :: LAT           = MISSING
    REAL(fp)                                     :: LEV           = MISSING
    REAL(fp)                                     :: LENGTH        = MISSING
    REAL(fp)                                     :: ALPHA         = MISSING
    REAL(fp)                                     :: LIFE          = MISSING
    REAL(fp)                                     :: RA          = MISSING
    REAL(fp)                                     :: RB          = MISSING
    REAL(fp)                                     :: THETA          = MISSING
    !REAL(fp), DIMENSION(:), ALLOCATABLE          :: RB                                 ! [n_species]
    !REAL(fp), DIMENSION(:), ALLOCATABLE          :: RA                                 ! [n_species]
    !REAL(fp), DIMENSION(:), ALLOCATABLE          :: THETA                              ! [n_species]
    REAL(fp), DIMENSION(:,:), ALLOCATABLE        :: CONCNT1d                           ! [n_slab_max,n_species]
    REAL(fp), DIMENSION(:), ALLOCATABLE          :: MassRef1d                          ! [n_species]

    TYPE(Plume1d_list), POINTER                  :: next          => NULL()

  END TYPE

END MODULE Plume_list_mod