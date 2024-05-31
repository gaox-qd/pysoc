      SUBROUTINE ALTMOS(CA,CB,STRWRT,DMAT) 
C
C     ******************************************************************
C     ***                       Alter MOs                            ***
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
      INCLUDE 'molecule.h'
C
      CHARACTER :: LINE*80,STRWRT*(*),STRING*80
C
      INTEGER :: DMAT,NALMOA,NALMOB
C
      INTEGER :: IORB
C
      REAL, DIMENSION(DMAT,DMAT) :: CA,CB
C
      INTEGER :: ALLOCATION
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ORBA,ORBB
      REAL, ALLOCATABLE, DIMENSION(:) :: MATRIX
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      IF (.NOT.ALTER) RETURN
C
      READ(INP,*) NALMOA,NALMOB
      IF (NALMOA == 0) GOTO 100
C
C     ___ Allocate local fileds ___
C
      ALLOCATE(ORBA(NALMOA),ORBB(NALMOA),MATRIX(DMAT),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'ALTMOS: ALLOCATION FAILED'
        STOP
      END IF
C
      STRING = STRWRT//' STATE ALPHA MOs ALTERED:'
      WRITE(OUT,"(/,T2,A80)") STRING
C
      DO IORB=1,NALMOA
        READ(INP,*) ORBA(IORB),ORBB(IORB)
        WRITE(OUT,1000) ORBA(IORB),ORBB(IORB)
      END DO
C
      DO IORB=1,NALMOA
        MATRIX(:) = CA(:,ORBB(IORB))
        CA(:,ORBB(IORB)) = CA(:,ORBA(IORB))
        CA(:,ORBA(IORB)) = MATRIX(:)
      END DO
C
C     ___ Deallocate local fields 
C
      DEALLOCATE(ORBA,ORBB,MATRIX,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'ALTMOS: DEALLOCATION FAILED'
        STOP
      END IF
C
  100 CONTINUE
C
      READ(INP,"(A80)") LINE
C
      IF (NALMOB == 0) RETURN 
C
C     ___ Allocate local fileds ___
C
      ALLOCATE(ORBA(NALMOB),ORBB(NALMOB),MATRIX(DMAT),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'ALTMOS: ALLOCATION FAILED'
        STOP
      END IF
C
      STRING = STRWRT//' STATE BETA MOs ALTERED:'
      WRITE(OUT,"(/,T2,A80)") STRING
C
      DO IORB=1,NALMOB
        READ(INP,*) ORBA(IORB),ORBB(IORB)
        WRITE(OUT,1000) ORBA(IORB),ORBB(IORB)
      END DO
C
      DO IORB=1,NALMOB
        MATRIX(:) = CB(:,ORBB(IORB))
        CB(:,ORBB(IORB)) = CB(:,ORBA(IORB))
        CB(:,ORBA(IORB)) = MATRIX(:)
      END DO
C
      READ(INP,"(A80)") LINE
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(ORBA,ORBB,MATRIX,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'ALTMOS: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
 1000 FORMAT(T2,I5,2X,I5)
C
C     O________________________________________________________________O
C
      END SUBROUTINE ALTMOS
