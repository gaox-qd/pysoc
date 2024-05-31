      SUBROUTINE DISCZERO(CA1,CA2,NOA,NOB,DMAT)
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
      REAL, PARAMETER :: TOL = 1.0E-10
C
      INTEGER :: DMAT,NOA,NOB,LOWLIM
      INTEGER :: IMAT,JMAT
      REAL :: PROD,TPROD,SIG,SOMO
C
      REAL, DIMENSION(DMAT) :: CA1,CA2
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
        WRITE(OUT,*) 'DISCZERO: ALLOCATION FAILED'
        STOP
      END IF
C
C     ___ Read from tape two-state MO overlap matrices ___
C
      CALL IOSCR(SMOA,DMAT**2,'READ','ZERO',5)
      CALL IOSCR(SMOB,DMAT**2,'READ','ZERO',6)
C
      IF (ABS(SMOA(NOA,NOA)) > TOL) WRITE(OUT,1000),NOA,
     $                              SMOA(NOA,NOA)
C
      TPROD = 1.0
      IF (NOB > 0) THEN
        DO IMAT=1,NOB
          TPROD = TPROD*SMOB(IMAT,IMAT)
        END DO
      END IF
C
      PROD = 1.0
      IF (NOA-1 <= 0) STOP 'DISCZERO: NOA-1 <= 0'
      DO JMAT=1,NOA-1
        PROD = PROD*SMOA(JMAT,JMAT)
      END DO
C
C     ___ Read from tape MO SOC integrals ___
C
C     ___ X component ___
C
      CALL IOSCR(MATRIX,DMAT**2,'READ','ZERO',1)
      SOMO = 0.0
      DO IMAT=1,DMAT
        DO JMAT=1,DMAT
          SOMO = SOMO + CA1(IMAT)*CA2(JMAT)*MATRIX(IMAT,JMAT)
        END DO
      END DO
C
      ENERD(1) = SIG*0.5*SOMO*PROD
      ENERD(1) = ENERD(1)*TPROD*FINESTRC**2
      ENER(1) = ENERD(1)*CMM1
C
C     ___ Y component ___
C
      CALL IOSCR(MATRIX,DMAT**2,'READ','ZERO',2)
      SOMO = 0.0
      DO IMAT=1,DMAT
        DO JMAT=1,DMAT
          SOMO = SOMO + CA1(IMAT)*CA2(JMAT)*MATRIX(IMAT,JMAT)
        END DO
      END DO
C
      ENERD(2) = SIG*0.5*SOMO*PROD
      ENERD(2) = ENERD(2)*TPROD*FINESTRC**2
      ENER(2) = ENERD(2)*CMM1
C
C     ___ Z component ___
C
      CALL IOSCR(MATRIX,DMAT**2,'READ','ZERO',3)
      SOMO = 0.0
      DO IMAT=1,DMAT
        DO JMAT=1,DMAT
          SOMO = SOMO + CA1(IMAT)*CA2(JMAT)*MATRIX(IMAT,JMAT)
        END DO
      END DO
C
      ENERD(3) = SIG*0.5*SOMO*PROD
      ENERD(3) = ENERD(3)*TPROD*FINESTRC**2
      ENER(3) = ENERD(3)*CMM1 
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(MATRIX,SMOA,SMOB,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'DISCZERO: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
 1000 FORMAT(/,T2,55('*'),
     $       /,T2,'*** WARNING IN DISCZERO: ',
     $            'ABS(SA(NOA,NOA)) > 1.0E-10 ***',
     $       /,T2,'***',9X,I5,10X,1PE15.8,10X,'***',
     $       /,T2,55('*'))
C
C     O________________________________________________________________O
C
      END SUBROUTINE DISCZERO
