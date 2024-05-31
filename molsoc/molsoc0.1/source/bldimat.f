      SUBROUTINE BLDIMAT(ISHELL,JSHELL,NORM,MATRIX,BLOCKS,DMAT)
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
      INTEGER :: DMAT
      INTEGER :: ICO,I,J,IRO,ISHELL,JSHELL
C
      REAL, DIMENSION(DMAT,DMAT) :: MATRIX
      REAL, DIMENSION(DSHL,DSHL) :: BLOCKS
C
      REAL, DIMENSION(*) :: NORM
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      DO IRO=LLSTO(ISHELL),ULSTO(ISHELL)
        I = IRO - LLSTO(ISHELL) + 1
        DO ICO=LLSTO(JSHELL),ULSTO(JSHELL)
          J = ICO - LLSTO(JSHELL) + 1
          MATRIX(IRO,ICO) = MATRIX(IRO,ICO) +
     $                      NORM(IRO)*NORM(ICO)*BLOCKS(I,J)
        END DO
      END DO
C
C     O________________________________________________________________O
C
      END SUBROUTINE BLDIMAT
