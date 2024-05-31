      SUBROUTINE OVERFS(C1,C2,HOMO1,HOMO2,DMAT,NUM)
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
C
      INCLUDE 'fileio.h'
      INCLUDE 'moleculeD.h'
      INCLUDE 'molecule.h'
C
      INTEGER :: DMAT,NUM,HOMO1,HOMO2
C
      REAL :: CPUTOLD
      REAL, DIMENSION(DMAT,DMAT) :: C1,C2
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:,:), ALLOCATABLE :: SMO
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Allocate local fields ___
C
      ALLOCATE(SMO(DMAT,DMAT),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'OVERFS: ALLOCATION FAILED'
        STOP
      END IF
C
      IF (NOBIO) GOTO 100
C
C     ___ Calculate the two state MO overlap matrix  ___
C
      CALL MOS(C1,C2,SMO,DMAT)
C
C     ___ Get CPU time ___
C
      CALL CPUTIME
      CPUTOLD = CPUT
C
C     ___ SVD ___
C
      IF (HOMO1.LT.HOMO2) THEN
        SMO = TRANSPOSE(SMO)
        CALL ORTSVDLAP(C2,C1,SMO,HOMO2,HOMO1,DMAT)
      ELSE
        CALL ORTSVDLAP(C1,C2,SMO,HOMO1,HOMO2,DMAT)
      END IF
C
C     ___ Get CPU time ___
C
      CALL CPUTIME
      WRITE(OUT,1000) CPUT - CPUTOLD
C
  100 CONTINUE
C
C     ___ Calculate MO overlap integrals ___
C
      CALL MOS(C1,C2,SMO,DMAT)
C
C     CALL CHKSGMO(SMO,HOMO1,HOMO2,DMAT)
C
C     ___ Write on tape MO overlap integrals ___
C
      CALL IOSCR(SMO,DMAT**2,'WRITE','ZERO',NUM)
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(SMO,STAT=ALLOCATION)
      IF (ALLOCATION.GT.0) THEN 
        WRITE(OUT,*) 'OVERFS: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Format statements ___
C
 1000 FORMAT(/,T2,'BIORTHOGONALIZATION PERFORMED IN',F12.6,' SEC.')
C
C     O________________________________________________________________O
C
      END SUBROUTINE OVERFS
