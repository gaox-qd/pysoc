      SUBROUTINE MOOCC(OCCN,NOMO,LFOMO,HFOMO)
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
C
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
C
      INTEGER :: HFOMO,LFOMO,NOMO
      INTEGER :: IMO
C
      REAL, DIMENSION(MXSTO) :: OCCN
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     *** Set MO occupation ***
C
      DO IMO=1,NOMO
        OCCN(IMO) = OCCNUM
      END DO
C
      DO IMO=NOMO+1,NSTO
        OCCN(IMO) = 0.0
      END DO
C
      LFOMO = NOMO + 1
      HFOMO = NOMO
C
C
C     O________________________________________________________________O
C
      END SUBROUTINE MOOCC
