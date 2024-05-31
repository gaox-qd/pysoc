      SUBROUTINE MATHCON
C
C     ******************************************************************
C     ***   Initialization of mathematical constants and functions   ***
C     ******************************************************************
C
C     List of mathematical constants:
C
C     PI: Mathematical constant Pi.
C
C     List of mathematical functions:
C
C     DFAC: Double factorial function.
C     FAC : Factorial function.
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'math.h'
C
      INTEGER :: I
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Definition of Pi ___
C
      PI = 2.0*ACOS(0.0)
C
C     ___ Definition of factorial function ___
C
      FAC(0) = 1.0
C
      DO I=1,MAXFAC
        FAC(I) = REAL(I)*FAC(I-1)
      END DO 
C
C     ___ Definition of double factorial function ___
C
      DFAC(-1) = 1.0
      DFAC(0) = 1.0
      DFAC(1) = 1.0
      DFAC(2) = 2.0
C
      DO I=3,2*MAXFAC+1
        DFAC(I) = REAL(I)*DFAC(I-2)
      END DO
C
      UM = 0.0
      DO I=1,3
        UM(I,I) = 1.0
      END DO
C
C     O________________________________________________________________O
C
      END SUBROUTINE MATHCON
