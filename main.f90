!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  sp3  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 
PROGRAM MAIN 
  USE mod_indata
  USE mod_logga
  USE mod_fftw3params
  IMPLICIT NONE
!*********************************************************************
!*  sp3    Computes the Fourier Transform of the Coulomb potential
!*  uses fftw 3 libraries - to be compiled with -r8
!*  from sp3 2003                               by Andrea Bertoni
!*  11 march 2020         PSXE_2019             by LB
!*********************************************************************

!............................................................Declarations
  
  INTEGER*8 :: planfw, planbk
  REAL*8 :: cnorm, fftnorm
  INTEGER :: nx, ny, nz, nn

!..........................................init vars definitions and readout
CALL INDATA_GET('fftw3.input')
CALL LOGGA(3, " == START ==")
cnorm= (deltax*deltay*deltaz) / ELCH   !  J -> eV
fftnorm= 1d0 / REAL(numx*numy*numz, 8)

OPEN(22, FILE=TRIM(bindirout)//"/psi_rspacenz2ini.dat", ACTION="WRITE", FORM="FORMATTED")
nz=2
  DO ny= 1, numy
    DO nx= 1, numx
      WRITE(22,*)  nx, ny, ABS(psi(nx, ny, nz))**2
    END DO
    WRITE(22,*) ''
  END DO
CLOSE(22)

!............................................................initializations

CALL LOGGA(2, "creating fftw3 forward plan...")
CALL dfftw_plan_dft_3d( planfw, numx, numy, numz, psi, psi,           &
     &                                   FFTW_FORWARD,  FFTW_MEASURE  )
CALL LOGGA(2, "...done")

CALL LOGGA(2, "computing fourier transform of the sp wavefunction...")
CALL dfftw_execute( planfw, psi, psi )
psi=psi*SQRT(fftnorm)

OPEN(21, FILE="psi_kspace.bin", ACTION="WRITE", FORM="UNFORMATTED")
WRITE(21)  numx, numy, numz
WRITE(21)  REAL(psi,4)
CLOSE(21)
OPEN(22, FILE=TRIM(bindirout)//"/psi_kspace_nz2.dat", ACTION="WRITE", FORM="FORMATTED")
nz= 2
  DO ny= 1, numy
    DO nx= 1, numx
      WRITE(22,*)  nx, ny, ABS(psi(nx, ny, nz))**2
    END DO
    WRITE(22,*) ''
  END DO
CLOSE(22)

CALL LOGGA(2, "deallocating plans and arrays...")
CALL dfftw_destroy_plan( planfw )
CALL LOGGA(2, "creating fftw3 backward plan...")
CALL dfftw_plan_dft_3d( planbk, numx, numy, numz, psi, psi,         &
     &                                   FFTW_BACKWARD, FFTW_MEASURE  )
CALL LOGGA(2, "...done")

CALL LOGGA(2, "computing anti fourier transform of the Coulomb potential...")
CALL dfftw_execute( planbk )
psi=psi*SQRT(fftnorm)

OPEN(22, FILE=TRIM(bindirout)//"/psi_rspacenz2fin.dat", ACTION="WRITE", FORM="FORMATTED")
nz=2
  DO ny= 1, numy
    DO nx= 1, numx
      WRITE(22,*)  nx, ny, ABS(psi(nx, ny, nz))**2
    END DO
    WRITE(22,*) ''
  END DO
CLOSE(22)


CALL LOGGA(2, "...done")

CALL LOGGA(2, "deallocating plans and arrays...")
CALL dfftw_destroy_plan( planbk )
DEALLOCATE( psi )
CALL LOGGA(2, "...done")

CALL LOGGA(3, "  == STOP ==")


END PROGRAM MAIN
