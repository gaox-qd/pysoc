      SUBROUTINE NORGTF
C
C     *******************************************************************
C     ***        Normalization of the Gaussian-type functions         ***
C     *******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
      INCLUDE 'math.h'
C
      INTEGER :: L,I,J
      REAL :: N
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      IF (TDTB) THEN 
        GCC = GCC
      ELSE
        DO I=1,NSHL
          L = PSHELL(I,2)
          N = 2.0**L*(2.0/PI)**0.75
          DO J=GLL(I),GUL(I)
            GCC(J) = N*(ZET(J)**(2*L+3))**0.25*GCC(J)
          END DO  
        END DO  
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE NORGTF
