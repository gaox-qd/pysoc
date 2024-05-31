      SUBROUTINE VECFSO(S,PDEL,U)
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'math.h'
C
      INTEGER :: U
C
      REAL, DIMENSION(3) :: S,PDEL
      REAL, DIMENSION(3,3) :: DE
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      CALL VPROD(S,UM(:,1),DE(:,1))
      CALL VPROD(S,UM(:,2),DE(:,2))
      CALL VPROD(S,UM(:,3),DE(:,3))
C
      PDEL(1) = DE(U,1)
      PDEL(2) = DE(U,2)
      PDEL(3) = DE(U,3)
C
C     O________________________________________________________________O
C
      END SUBROUTINE VECFSO
