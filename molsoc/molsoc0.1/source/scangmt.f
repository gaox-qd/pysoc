      SUBROUTINE SCANGMT(AX,AY,AZ,BX,BY,BZ,CX,CY,CZ,COOR,S,DELK)
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
C     SA/B     : Scaling number.
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      IMPLICIT NONE
C
      INTEGER :: AX,AY,AZ,BX,BY,BZ,CX,CY,CZ,COOR,I
      INTEGER, DIMENSION(3) :: A,B,C,DELK,S
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Initialization ___
C
      A(1) = AX; A(2) = AY; A(3) = AZ
      B(1) = BX; B(2) = BY; B(3) = BZ
      C(1) = CX; C(2) = CY; C(3) = CZ
      COOR = 0
C
C     ___ Check pointer of the angular momentum ___
C
      DO I=1,3
        S(I) = 0
        IF (C(I) /= 0) COOR = I
      END DO
      S(COOR) = 1
C
C     ___ Check Kronecker's delta for recurrence relation ___
C
      DELK(:) = 0
C
      IF (A(COOR) /= 0) DELK(3) = A(COOR)
      IF (B(COOR) /= 0) DELK(2) = B(COOR) 
      IF (C(COOR) > 1) DELK(1) = C(COOR) - 1
C
C     O________________________________________________________________O
C
      END SUBROUTINE SCANGMT
