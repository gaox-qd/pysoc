      SUBROUTINE DISCONE(CA,CB,NOAF,NOBS,DMAT)
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
      INCLUDE 'fileio.h'
      INCLUDE 'math.h'
C
      INTEGER :: DMAT,NOAF,NOBS
      INTEGER :: IMAT,JMAT
      REAL :: PROD,PRODA,PRODB,SIG,SOMO
C
      REAL, DIMENSION(DMAT) :: CA,CB
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:,:), ALLOCATABLE :: MATRIX,SMOA,SMOB
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      SIG = 1.0
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Allocate local fields ___
C
      ALLOCATE(MATRIX(DMAT,DMAT),SMOA(DMAT,DMAT),SMOB(DMAT,DMAT),
     $         STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'DISCONE: ALLOCATION FAILED'
        STOP
      END IF
C
C     ___ Read from tape two-state MO overlap matrices ___
C
      CALL IOSCR(SMOA,DMAT**2,'READ','ZERO',5)
      CALL IOSCR(SMOB,DMAT**2,'READ','ZERO',6)
C
      PRODA = 1.0
      DO IMAT=1,NOAF
        PRODA = PRODA*SMOA(IMAT,IMAT)
      END DO
C
      PRODB = 1.0
      IF (NOBS > 0) THEN
        DO IMAT=1,NOBS
          PRODB = PRODB*SMOB(IMAT,IMAT)
        END DO
      END IF
C
      PROD = PRODA*PRODB
C
C     ___ Read from tape MO SOC integrals ___
C
      CALL IOSCR(MATRIX,DMAT**2,'READ','ZERO',1)
C
      SOMO = 0.0
      DO IMAT=1,DMAT
        DO JMAT=1,DMAT
          SOMO = SOMO + CB(IMAT)*CA(JMAT)*MATRIX(IMAT,JMAT)
        END DO
      END DO
C
      ENERD(1) = 0.5*FINESTRC**2*SOMO*PROD
      ENER(1) = ENERD(1)*CMM1
C
      CALL IOSCR(MATRIX,DMAT**2,'READ','ZERO',2)
C
      SOMO = 0.0
      DO IMAT=1,DMAT
        DO JMAT=1,DMAT
          SOMO = SOMO + CB(IMAT)*CA(JMAT)*MATRIX(IMAT,JMAT)
        END DO
      END DO
C
      ENERD(2) = 0.5*FINESTRC**2*SOMO*PROD
      ENER(2) = ENERD(2)*CMM1
C
      CALL IOSCR(MATRIX,DMAT**2,'READ','ZERO',3)
C
      SOMO = 0.0
      DO IMAT=1,DMAT
        DO JMAT=1,DMAT
          SOMO = SOMO + CB(IMAT)*CA(JMAT)*MATRIX(IMAT,JMAT)
        END DO
      END DO
C
      ENERD(3) = 0.5*FINESTRC**2*SOMO*PROD
      ENER(3) = ENERD(3)*CMM1
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(MATRIX,SMOA,SMOB,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'DISCONE: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE DISCONE
