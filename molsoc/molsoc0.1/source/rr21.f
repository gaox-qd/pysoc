      SUBROUTINE RR21(BLOCKS,RRINT,RRERI,IL,LA,LB,DPA,DPB,DQA,DQB,DWP,
     $                NIPR,RHO,ZETIG,ZETJG,ZETNI,ZETPR,DCPA,DCPC,DSH,
     $                ULLR,LBCD)
C
C      Creation (12.08.03, SC)
C               (05.12.04, SC)
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'moleculeD.h'
C
      INTEGER :: DCPA,DCPC,DSH,ULLR,LBCD
      INTEGER :: LA,LB,COOR,IL,N,X,Y,Z
      INTEGER :: IC,I1S,I2S,IS
      INTEGER, DIMENSION(3) :: S
      INTEGER, DIMENSION(4) :: DELK
C
      REAL :: NIPR,RHO,RRERI2,RRERI3,ZETIG,ZETJG,ZETNI,ZETPR
      REAL, DIMENSION(3) :: DPA,DPB,DQA,DQB,DWP,PDEL1,PDEL2,RRINT2,
     $                      RRINT3,RS
C
      REAL, DIMENSION(DSH,DSH,DSH,DSH,3) :: BLOCKS
      REAL, DIMENSION(DCPA,DCPA,DCPC,DCPC,0:ULLR,3) :: RRINT
      REAL, DIMENSION(DCPA,DCPA,DCPC,DCPC,0:ULLR+1) :: RRERI
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      DO 30 X=IL,0,-1
        DO 20 Y=IL-X,0,-1
          Z = IL - X - Y
          CALL SCANGMF(0,0,0,0,0,0,0,0,0,X,Y,Z,COOR,S,DELK)
          RS(:) = REAL(S(:))
          CALL VPROD(RS,DQB,PDEL1)
          CALL VPROD(RS,DQA,PDEL2)
          IC = COP(X,Y,Z)
          I1S = COP(X-S(1),Y-S(2),Z-S(3)) 
          DO 10 N=0,ULLR-IL
            RRINT2 = 0.0; RRERI2 = 0.0
C
            IF (DELK(1) /= 0) THEN
              I2S = COP(X-2*S(1),Y-2*S(2),Z-2*S(3))
              RRERI2 = RRERI(I2S,1,1,1,N+1) -  RHO/ZETPR*
     $                 RRERI(I2S,1,1,1,N+2)
              RRINT2(:) = RRINT(I2S,1,1,1,N,:) - RHO/ZETPR*
     $                    RRINT(I2S,1,1,1,N+1,:)
            END IF
            RRERI(IC,1,1,1,N+1) = DPA(COOR)*RRERI(I1S,1,1,1,N+1) +
     $                            DWP(COOR)*RRERI(I1S,1,1,1,N+2) +
     $                            0.5*DELK(1)/ZETPR*RRERI2
            RRINT(IC,1,1,1,N,:) = DPA(COOR)*RRINT(I1S,1,1,1,N,:) +
     $                            DWP(COOR)*RRINT(I1S,1,1,1,N+1,:) +
     $                            0.5*DELK(1)/ZETPR*RRINT2(:) - 
     $                            2.0*ZETJG*NIPR/ZETNI*PDEL1(:)*
     $                            RRERI(I1S,1,1,1,N+1)
C
            IF ((IL == LA).AND.(LBCD == 0)) THEN
              IS = SOP(X,Y,Z)
              BLOCKS(IS,1,1,1,:) = BLOCKS(IS,1,1,1,:) +
     $                             RRINT(IC,1,1,1,N,:)
              IF (N /= 0) THEN
                STOP 'ERROR IN RECUR. A'
              END IF
            END IF
C
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C
C     O_______________________________________________________________O
C
      END SUBROUTINE RR21
