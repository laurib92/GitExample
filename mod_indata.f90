!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  sp3.f90  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

MODULE mod_indata
  USE mod_staticdata
  IMPLICIT  NONE
  SAVE

! data in the psi file
  INTEGER                 :: numx, numy, numz    ! grid points
  REAL*8                  :: deltax, deltay, deltaz
  COMPLEX*16, ALLOCATABLE :: psi(:,:,:)        ! sp eigenfunctions
! deriv const
  INTEGER                 :: numxyz              ! numx*numy*numz
  REAL*8                  :: deltakx, deltaky, deltakz

! indata_inoutput
  CHARACTER(80) :: bindirin      ! directory where to get bin files
  CHARACTER(80) :: psifilein     ! name of the psi file
  CHARACTER(80) :: psifileformat ! format of the psi file ("states" or "sp2")
  REAL*8 :: scalebox_x, scalebox_y, scalebox_z 
  CHARACTER(80) :: bindirout     ! dir in which to put bin files
  CHARACTER(80) :: ascdirout     ! dir in which to put bin files
  INTEGER :: loglevel            ! with 0 everything is logged
  CHARACTER(80) :: statusfile    ! name of status file
! reserved for future parallelization
! indata_parallel
  INTEGER :: npeoutput= 0, npeinput= 0, npelogga= 0  ! process for i/o
  INTEGER :: mypenum= 0

  NAMELIST /indata_inoutput/                   &
       &  bindirin,                            &
       &  psifilein,                           &
       &  psifileformat,                       &
       &  scalebox_x, scalebox_y, scalebox_z,  &
       &  bindirout,                           &
       &  ascdirout,                           &
       &  loglevel,                            &
       &  statusfile

CONTAINS

!===================================================================
SUBROUTINE INDATA_GET(nmlfile)
  IMPLICIT NONE
!  INTEGER*4           :: numx_read, numy_read, numz_read
!  REAL*4, ALLOCATABLE :: deltax_read(:), deltay_read(:), deltaz_read(:)
!  INTEGER*4           :: nconv_read
!  REAL*4, ALLOCATABLE :: energy_read(:)
!  COMPLEX*16, ALLOCATABLE :: psi_read(:,:)  ! tmp array for reading psi
  CHARACTER(*), INTENT(IN) :: nmlfile 
  REAL*8 :: cnormpsi
  INTEGER :: nn

  OPEN(11, FILE=TRIM(nmlfile), FORM="FORMATTED", ACTION="READ",  &
       &   STATUS="OLD")
  READ(11,NML=indata_inoutput)
  CLOSE(11)

  IF (psifileformat == "dbl") THEN
    CALL INDATA_PSIREAD_DBL()
  ELSE
    STOP "INDATA_GET:  Unknown psifileformat !"
  END IF

!!$  ! psi normalization 
!!$  DO nn= 1, nconv
!!$    cnormpsi= REAL( SUM( psi(:,:,:,nn)*CONJG(psi(:,:,:,nn)) ) )     &
!!$         &    * deltax * deltay * deltaz
!!$    PRINT*, nn, cnormpsi
!!$  END DO

  PRINT*, "numx,numy,numz: ", numx,numy,numz
  PRINT*, "deltax,deltay,deltaz: ", deltax, deltay, deltaz
  PRINT*
  PRINT*, "  now computing (see statusfile) ... "

END SUBROUTINE INDATA_GET

!===================================================================

SUBROUTINE INDATA_PSIREAD_DBL()
  IMPLICIT NONE
!...........................................psi data file read out
  OPEN(22, FILE=TRIM(bindirin)//"/"//TRIM(psifilein),           &
       &   STATUS="OLD", ACTION="READ", FORM="UNFORMATTED")

  READ(22)  numx, numy, numz, deltax, deltay, deltaz
  ALLOCATE(psi(numx, numy, numz))
  READ(22) psi
  CLOSE(22)
  !................................definition of derived global const
  numxyz= numx*numy*numz
  !psi(:,:,:)= psi(:,:,:)/SQRT(deltax*deltay*deltaz)
  deltakx= 2*PIG/(numx*deltax)
  deltaky= 2*PIG/(numy*deltay)
  deltakz= 2*PIG/(numz*deltaz)

END SUBROUTINE INDATA_PSIREAD_DBL


END MODULE mod_indata
