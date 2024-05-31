      SUBROUTINE BLDOVER 
C
C     ******************************************************************
C     ***                 Build overlap matrix                       ***
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
      INCLUDE 'moleculeD.h'
      INCLUDE 'molecule.h'
C
      INTEGER :: ISTO,JSTO
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:,:), ALLOCATABLE :: S
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Allocate local fields ___
C
      ALLOCATE(S(NSTO,NSTO),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'BLDOVER: ALLOCATION FAILED'
        STOP
      END IF
C
C     ___ Calculate overlap matrix ___
C
      CALL OVERLAP(S,NSTO)
C
C     ___ Store overlap matrix on tape ___
C
      CALL IOSCR(S,NSTO**2,'WRITE','ZERO',4)
      OPEN(UNIT=171716, FILE='molsoc_overlap.dat', STATUS='REPLACE')
      WRITE(171716, "(A15)") "AO_overlap"
      WRITE(171716, "(f19.9)") S(:,:)
      CLOSE(171716)
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(S,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'BLDOVER: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE BLDOVER
