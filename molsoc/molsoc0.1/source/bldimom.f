      SUBROUTINE BLDIMOM(ISHELL,JSHELL,NORM,MATRIX,BLOCKS,DMAT,DMOM)
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
      INTEGER :: DMAT,DMOM, counter=1
      INTEGER :: ICO,I,J,IRO,ISHELL,JSHELL
C
      REAL, DIMENSION(DMAT,DMAT,2:DMOM) :: MATRIX
      REAL, DIMENSION(DSHL,DSHL,2:DMOM) :: BLOCKS
C
      REAL, DIMENSION(*) :: NORM
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      DO IRO=LLSTO(ISHELL),ULSTO(ISHELL)
        I = IRO - LLSTO(ISHELL) + 1
        DO ICO=LLSTO(JSHELL),ULSTO(JSHELL)
          J = ICO - LLSTO(JSHELL) + 1
          MATRIX(IRO,ICO,2:DMOM) = MATRIX(IRO,ICO,2:DMOM) +
     $                             NORM(IRO)*NORM(ICO)*
     $                             BLOCKS(I,J,2:DMOM)
C        print *, "counter2", counter
C        print *, "IRO, ICO", IRO, ICO, I, J
C        counter = counter + 1
        END DO
      END DO
C
C     O________________________________________________________________O
C
      END SUBROUTINE BLDIMOM
