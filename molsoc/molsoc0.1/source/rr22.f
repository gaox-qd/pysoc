      SUBROUTINE RR22(BLOCKS,RRINT,RRERI,ILA,ILB,LA,LB,DPA,DPB,DQA,DQB,
     $                DWP,NIPR,RHO,ZETIG,ZETJG,ZETNI,ZETPR,DCPA,DCPC,
     $                DSH,ULLR,LCD)
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
      INTEGER :: DCPA,DCPC,DSH,ULLR,LCD
      INTEGER :: COOR,ILA,ILB,LA,LB,N,U,AX,AY,AZ,BX,BY,BZ
      INTEGER :: IA,IB,I1SA,I1SB,I2SB,ISA,ISB
      INTEGER, DIMENSION(3) :: S
      INTEGER, DIMENSION(4) :: DELK
C
      REAL :: NIPR,RHO,RRERI2B,RRERI3B,RRERI2A,RRERI3A,ZETIG,ZETJG,
     $        ZETNI,ZETPR
      REAL, DIMENSION(3) :: DPA,DPB,DWP,DQA,DQB,DQQ,PDEL1,PDEL2,RS,
     $                      RRINT2B,RRINT3B,RRVEA,RRVEB,RRINT2A,RRINT3A
      REAL, DIMENSION(3,3) :: PDELV
      REAL, DIMENSION(DSH,DSH,DSH,DSH,3) :: BLOCKS
      REAL, DIMENSION(DCPA,DCPA,DCPC,DCPC,0:ULLR,3) :: RRINT
      REAL, DIMENSION(DCPA,DCPA,DCPC,DCPC,0:ULLR+1) :: RRERI 
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      DO 50 AX=ILA,0,-1
        DO 40 AY=ILA-AX,0,-1
          AZ = ILA - AX - AY
          DO 30 BX=ILB,0,-1
            DO 20 BY=ILB-BX,0,-1
              BZ = ILB - BX - BY
              CALL SCANGMF(0,0,0,0,0,0,AX,AY,AZ,BX,BY,BZ,COOR,S,DELK)
              RS(:) = REAL(S(:))
              CALL VPROD(RS,DQA,PDEL1)
              CALL VECFSO2(RS,PDELV)
C
              IA = COP(AX,AY,AZ)
              IB = COP(BX,BY,BZ)
              I1SB = COP(BX-S(1),BY-S(2),BZ-S(3))
C
              DO 10 N=0,ULLR-ILB-LA
                RRINT2B = 0.0; RRINT3B = 0.0
                RRERI2B = 0.0; RRERI3B = 0.0
                IF (DELK(1) /= 0) THEN
                  I2SB = COP(BX-2*S(1),BY-2*S(2),BZ-2*S(3))
                  RRERI2B = RRERI(IA,I2SB,1,1,N+1) - RHO/ZETPR*
     $                      RRERI(IA,I2SB,1,1,N+2)
                  RRINT2B(:) = RRINT(IA,I2SB,1,1,N,:) - RHO/ZETPR*
     $                         RRINT(IA,I2SB,1,1,N+1,:)
                END IF
                IF (DELK(2) /= 0) THEN
                  I1SA = COP(AX-S(1),AY-S(2),AZ-S(3))
                  RRERI3B = RRERI(I1SA,I1SB,1,1,N+1) - RHO/ZETPR*
     $                      RRERI(I1SA,I1SB,1,1,N+2)
                  RRINT3B(:) = RRINT(I1SA,I1SB,1,1,N,:) - RHO/ZETPR*
     $                         RRINT(I1SA,I1SB,1,1,N+1,:)
                END IF
                RRERI(IA,IB,1,1,N+1) = DPB(COOR)*RRERI(IA,I1SB,1,1,N+1)+
     $                                 DWP(COOR)*RRERI(IA,I1SB,1,1,N+2)+
     $                                 0.5*DELK(1)/ZETPR*RRERI2B + 
     $                                 0.5*DELK(2)/ZETPR*RRERI3B
C
                DO U=1,3
                  RRVEB(U) = 0.0
                  IF ((AX /= 0).AND.(PDELV(U,1) /= 0.0))
     $              RRVEB(U) = AX*PDELV(U,1)*
     $                         RRERI(COP(AX-1,AY,AZ),I1SB,1,1,N+1)
                  IF ((AY /= 0).AND.(PDELV(U,2) /= 0.0))
     $              RRVEB(U) = RRVEB(U) + AY*PDELV(U,2)*
     $              RRERI(COP(AX,AY-1,AZ),I1SB,1,1,N+1)
                  IF ((AZ /= 0).AND.(PDELV(U,3) /= 0.0))
     $              RRVEB(U) = RRVEB(U) + AZ*PDELV(U,3)*
     $              RRERI(COP(AX,AY,AZ-1),I1SB,1,1,N+1)
                END DO
C
                RRINT(IA,IB,1,1,N,:)=DPB(COOR)*RRINT(IA,I1SB,1,1,N,:) +
     $                               DWP(COOR)*RRINT(IA,I1SB,1,1,N+1,:)+
     $                               0.5*DELK(1)/ZETPR*RRINT2B(:) +
     $                               0.5*DELK(2)/ZETPR*RRINT3B(:) + 
     $                               2.0*ZETIG*NIPR/ZETNI*PDEL1(:)*
     $                               RRERI(IA,I1SB,1,1,N+1) - 
     $                               RHO/ZETPR*RRVEB(:)
C
                IF ((ILA == LA).AND.(ILB == LB).AND.(LCD == 0)) THEN
                  ISA = SOP(AX,AY,AZ)
                  ISB = SOP(BX,BY,BZ)
                  BLOCKS(ISA,ISB,1,1,:) = BLOCKS(ISA,ISB,1,1,:) +
     $                                    RRINT(IA,IB,1,1,N,:)
                  IF (N /= 0) THEN
                    STOP 'ERROR IN RECUR. B'
                  END IF
                END IF
C
   10         CONTINUE
   20       CONTINUE
   30      CONTINUE
   40    CONTINUE
   50  CONTINUE
C
C     O_______________________________________________________________O
C
      END SUBROUTINE RR22
