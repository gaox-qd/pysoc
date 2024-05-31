      SUBROUTINE SHELLOP
C
C     ******************************************************************
C     ***              Generate Shell orbital pointer                ***
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
C
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      INTEGER :: IS,MUL,L
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      SLL = 0
      SUL = 0
      NSPH = 0
C
      DO IS=1,NSHL
        SLL(IS) = NSPH + 1
        L = PSHELL(IS,2)
        DO MUL=1,2*L+1
          NSPH = NSPH + 1
        END DO  
        SUL(IS) = NSPH
      END DO  
C
C     O________________________________________________________________O
C
      END SUBROUTINE SHELLOP
