      SUBROUTINE READMOS(MATRIX,OCCN,OPTION1,OPTION2,OPTION3,DENS)
C
C     *******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
C
      INCLUDE 'fileio.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      CHARACTER*(*) :: OPTION1,OPTION2,OPTION3
C
      INTEGER :: DENS
C
      REAL, DIMENSION(MXSTO) :: OCCN
      REAL, DIMENSION(NSTO,NSTO) :: MATRIX
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Load MOs from Gaussian formatted files ___
C
      IF (READMO == 1) THEN
        CALL GAUSSIAN(MATRIX,OCCN,OPTION1,OPTION2,OPTION3,DENS)
      ELSE IF (READMO == 2) THEN
        CALL DEMON(MATRIX,OCCN,OPTION1,OPTION2,DENS)
      ELSE IF (READMO == 3) THEN
        CALL TURBOMOLE(MATRIX,OCCN,OPTION1,OPTION2,OPTION3,DENS)
C     ELSE IF (READMO == 4) THEN
C       CALL DALTON(MATRIX,OCCN,OPTION1,OPTION2,OPTION3,DENS)
C     ELSE IF (READMO == 5) THEN
C       CALL GAMESS(MATRIX,OCCN,OPTION1,OPTION2,'CARTESIAN',DENS)
      ELSE
        WRITE(OUT,"(/,T2,'UNKNOWN MOs OUTPUT FILES')")
        STOP
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE READMOS
