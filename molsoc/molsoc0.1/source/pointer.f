      SUBROUTINE POINTER
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
      INCLUDE 'moleculeD.h'
C
      INTEGER :: AX,AY,AZ,I,IC,ICC
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      ICC = 0
      DO I=0,MXL
        IC = 0
        DO AX=I,0,-1
          DO AY=I-AX,0,-1
            AZ = I - AX - AY
            IC = IC + 1
            ICC = ICC + 1
            COP(AX,AY,AZ) = ICC
            SOP(AX,AY,AZ) = IC
C            print *, "ICC,IC", ICC,IC, AX,AY,AZ
          END DO
        END DO
      END DO
C      COP(2,0,0) = 5
C      COP(0,2,0) = 6
C      COP(0,0,2) = 7
C      COP(1,1,0) = 8
C      COP(1,0,1) = 9
C      COP(0,1,1) = 10
C      SOP(2,0,0) = 1
C      SOP(0,2,0) = 2
C      SOP(0,0,2) = 3
C      SOP(1,1,0) = 4
C      SOP(1,0,1) = 5
C      SOP(0,1,1) = 6
C
C     O________________________________________________________________O
C
      END SUBROUTINE POINTER
