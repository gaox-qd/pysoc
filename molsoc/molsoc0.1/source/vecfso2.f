      SUBROUTINE VECFSO2(S,DE)
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'math.h'
C
      REAL, DIMENSION(3) ::  S
      REAL, DIMENSION(3,3) :: DE,PT
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      CALL VPROD(S,UM(:,1),DE(:,1))
      CALL VPROD(S,UM(:,2),DE(:,2))
      CALL VPROD(S,UM(:,3),DE(:,3)) 
C
C     O________________________________________________________________O
C
      END SUBROUTINE VECFSO2
