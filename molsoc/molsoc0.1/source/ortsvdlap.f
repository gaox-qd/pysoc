      SUBROUTINE ORTSVDLAP(CA,CB,SOMO,HOMOA,HOMOB,DMAT)
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
C
      INTEGER :: DMAT,HOMOA,HOMOB,DIMS
      INTEGER :: IMOA,IMOB,IMAT,JMAT,KMAT
      INTEGER :: LWORK,INFO
C
      REAL :: DET
      REAL, DIMENSION(DMAT,DMAT) :: CA,CB,SOMO
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:), ALLOCATABLE :: S,WORK
      REAL,DIMENSION(:,:), ALLOCATABLE :: U,V,A
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      LWORK = 5*MAX(3*MIN(HOMOA,HOMOB) + 
     $          MAX(HOMOA,HOMOB),5*MIN(HOMOA,HOMOB))
C
      DIMS = MIN(HOMOA,HOMOB) + 1
C
C     ___ Allocate local fields ___
C
      ALLOCATE(U(DMAT,DMAT),V(DMAT,DMAT),A(DMAT,DMAT),
     $         S(DIMS),WORK(LWORK),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'ORTSVDLAP: ALLOCATION FAILED'
        STOP
      END IF
C
C     ___ SVD matrices initialization ___
C
      U = 0.0
      V = 0.0
      A = 0.0
      S = 0.0
      WORK = 0.0
C
      DO IMOA=1,DMAT  
        DO IMOB=1,DMAT 
          A(IMOA,IMOB) = SOMO(IMOA,IMOB)
        END DO
      END DO
C
C     ___ Do singular value decomposition (Lapack) ___
C
      CALL DGESVD('A','A',HOMOA,HOMOB,A,DMAT,S,U,DMAT,V,DMAT,
     $                   WORK,LWORK,INFO)
C
      IF (INFO.EQ.0) THEN 
        WRITE(OUT,"(/,T2,'ORTSVDLAP: DGESVD SUCCESSFULL EXIT:',I5)") 
     $        INFO
      ELSE IF (INFO.LT.0) THEN
        WRITE(OUT,"(/,T2,'ORTSVDLAP: IN DGESVD THE ',I5,
     $              '-th ARGUMENT HAS AN ILLEGAL VALUE')") INFO
      ELSE IF (INFO.GT.0) THEN
        WRITE(OUT,"(/,T2,'ORTSVDLAP:',I5,'DGESVD DID NOT CONVERGE')") 
     $        INFO
      END IF
C
C     ___ Update MO alpha coefficients ___
C
      A = 0.0
      DO IMAT=1,DMAT
        DO JMAT=1,HOMOA
          DO KMAT=1,HOMOA
            A(IMAT,JMAT) = A(IMAT,JMAT) + CA(IMAT,KMAT)*U(KMAT,JMAT)
          END DO
        END DO
      END DO
C
      DO IMAT=1,DMAT
        DO JMAT=1,HOMOA
          CA(IMAT,JMAT) = A(IMAT,JMAT)
        END DO 
      END DO
C
C     ___ Update MO beta coefficients ___
C
      A = 0.0
      DO IMAT=1,DMAT
        DO JMAT=1,HOMOB
          DO KMAT=1,HOMOB
            A(IMAT,JMAT) = A(IMAT,JMAT) + CB(IMAT,KMAT)*V(JMAT,KMAT)
          END DO
        END DO
      END DO
C
      DO IMAT=1,DMAT
        DO JMAT=1,HOMOB
          CB(IMAT,JMAT) = A(IMAT,JMAT)
        END DO
      END DO
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(U,V,A,S,WORK,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'ORTSVDLAP: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O_______________________________________________________________O
C
      END SUBROUTINE ORTSVDLAP
