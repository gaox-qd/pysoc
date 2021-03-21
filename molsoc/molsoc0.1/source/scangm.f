      SUBROUTINE SCANGM(AX,AY,AZ,BX,BY,BZ,COOR,S,DELK)
C
C     Purpose: Get scaled ANGular Momentum indeces for 
C              cartesian gaussian functions.
C
C     ******************************************************************
C
C     List of local variables:
C
C     AX,AY,AZ : Angular momentum index of orbital A.
C     BX,BY,BZ : Angular momentum index of orbital B.
C     COOR     : Pointer position of the angular momentum to be scaled.
C     DELK     : Like Kronecker's delta variable.
C     S        : Scaling number.
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      IMPLICIT NONE
C
      INTEGER :: AX,AY,AZ,BX,BY,BZ,COOR,I
      INTEGER, DIMENSION(2) :: DELK
      INTEGER, DIMENSION(3) :: A,B,S
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Initialization ___
C
      A(1) = AX; A(2) = AY; A(3) = AZ
      B(1) = BX; B(2) = BY; B(3) = BZ
      COOR = 0
C
C     ___ Check pointer of the angular momentum ___
C
      DO I=1,3
        S(I) = 0
        IF (B(I) /= 0) COOR = I
      END DO
      S(COOR) = 1
C
C     ___ Check Kronecker's delta for left recurrence relation ___
C
      DELK(:) = 0
      IF (A(COOR) /= 0) DELK(1) = A(COOR)
C
C     ___ Check Kronecker's delta for right recurrence relation ___
C
      IF (B(COOR) > 1) DELK(2) = B(COOR) - 1
C
C     O________________________________________________________________O
C
      END SUBROUTINE SCANGM
