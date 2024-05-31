      REAL FUNCTION TRAATB(A,B,DMAT,N)
C
C     *****************************************************************
C
      IMPLICIT NONE
C
      INTEGER :: DMAT
      INTEGER :: I,J,N
C
      REAL, DIMENSION(DMAT,DMAT) :: A,B
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      TRAATB = 0.0
C
      DO I=1,N
        DO J=I,N
          TRAATB = TRAATB + A(I,J)*B(I,J)
        END DO
      END DO
C
C     O________________________________________________________________O
C
      END FUNCTION TRAATB
