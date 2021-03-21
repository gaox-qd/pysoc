      SUBROUTINE MOS(C1,C2,SMO,DMAT)
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
C
      INTEGER :: DMAT,N
C
      REAL, DIMENSION(DMAT,DMAT) :: C1,C2,SMO
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:,:), ALLOCATABLE :: CMAT
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Allocate local fields ___
C
      ALLOCATE(CMAT(DMAT,DMAT),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'MOS: ALLOCATION FAILED'
        STOP
      END IF
C
      N = DMAT
C
C     ___ Read from tape overlap integrals ___
C
      SMO = 0.0
      CALL IOSCR(SMO,DMAT**2,'READ','ZERO',4)
C
C     ___ Calculate the two state MO overlap matrix ___
C
      CMAT = 0.0
      CALL DGEMM('N','N',N,N,N,1.0,SMO,DMAT,C2,DMAT,0.0,CMAT,DMAT)
C
      SMO = 0.0
      CALL DGEMM('T','N',N,N,N,1.0,C1,DMAT,CMAT,DMAT,0.0,SMO,DMAT)
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(CMAT,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'MOS: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE MOS
