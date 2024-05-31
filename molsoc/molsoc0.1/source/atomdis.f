      SUBROUTINE ATOMDIS(ATOMA,ATOMB,R,RSQUARE)
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
      INTEGER :: ATOMA,ATOMB,I
      REAL :: RSQUARE
C
      REAL, DIMENSION(3) :: R
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      RSQUARE = 0.0
      IF (ATOMA /= ATOMB) RSQUARE =  SDISTANCE(ATOMA,ATOMB)
C
      R(:) = C(:,ATOMA) - C(:,ATOMB)
C
C     O________________________________________________________________O
C
      END SUBROUTINE ATOMDIS
