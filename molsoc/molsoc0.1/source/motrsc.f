      SUBROUTINE MOTRSC(SPHER,CART,DMAT)
C
C     ******************************************************************
C     ***       MO transformation from sphericals to cartesians      ***
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
C
      INCLUDE 'math.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      INTEGER :: DMAT
      INTEGER :: I,IMO,J,ISHELL,IS,ISLAT,L
C
      REAL, DIMENSION(DMAT,DMAT) :: CART,SPHER
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      CART = 0.0
      DO IMO=1,NSTO
        DO ISHELL=1,NSHL
          L = PSHELL(ISHELL,2)
          DO ISLAT=LLSTO(ISHELL),ULSTO(ISHELL)
            I = ISLAT - LLSTO(ISHELL) + 1
            DO IS=SLL(ISHELL),SUL(ISHELL)
              J = IS - SLL(ISHELL) + 1
              CART(ISLAT,IMO) = CART(ISLAT,IMO) +
     $                          TRMSC(J,I,L)*SPHER(IS,IMO)
            END DO
          END DO
        END DO
      END DO
C
C     O________________________________________________________________O
C
      END SUBROUTINE MOTRSC
