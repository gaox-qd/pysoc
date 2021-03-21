      SUBROUTINE RR23(BLOCKS,RRINT,RRERI,ILA,ILB,ILC,LA,LB,LC,LD,
     $                DPA,DPB,DQC,DQD,DWQ,NIPR,RHO,ZETIG,ZETJG,ZETNI,
     $                ZETPR,DCPA,DCPC,DSH,ULLR)
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
      INTEGER :: DCPA,DCPC,DSH
      INTEGER :: AX,AY,AZ,COOR,BX,BY,BZ,CX,CY,CZ,ILA,ILB,ILC,LA,LB,
     $           LC,LD,N,U,ULLR
      INTEGER :: IA,IB,IC,ISA,ISB,ISC,I1SA,I1SB,I1SC,I2SC,IIA,IIB
      INTEGER, DIMENSION(3) :: S
      INTEGER, DIMENSION(4) :: DELK
C
      REAL :: NIPR,RHO,RRERI2,RRERI3,RRERI4,ZETIG,ZETJG,ZETNI,ZETPR
      REAL, DIMENSION(3) :: DPA,DPB,DQC,DQD,DWQ,PDEL,RS,RRINTC2,RRINTC3,
     $                      RRINTC4,RRINTD2,RRINTD3,RRINTD4,RRVEA,RRVEB,
     $                      RRVEC,RRVED
      REAL, DIMENSION(3,3) :: PDELV
      REAL, DIMENSION(DSH,DSH,DSH,DSH,3) :: BLOCKS
      REAL, DIMENSION(DCPA,DCPA,DCPC,DCPC,0:ULLR,3) :: RRINT
      REAL, DIMENSION(DCPA,DCPA,DCPC,DCPC,0:ULLR+1) :: RRERI
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ C vertical recurrence relation steps ___
C
      DO 70 AX=ILA,0,-1
        DO 60 AY=ILA-AX,0,-1
          AZ = ILA - AX - AY
          DO 50 BX=ILB,0,-1
            DO 40 BY=ILB-BX,0,-1
              BZ = ILB - BX - BY
              IA = COP(AX,AY,AZ)
              IB = COP(BX,BY,BZ)
              DO 30 CX=ILC,0,-1
                DO 20 CY=ILC-CX,0,-1
                  CZ = ILC - CX - CY
                  CALL SCANGMF(0,0,0,AX,AY,AZ,BX,BY,BZ,CX,CY,CZ,COOR,
     $                         S,DELK)
                  RS(:) = REAL(S(:))
C
                  CALL VPROD(RS,DPB,PDEL)
                  CALL VECFSO2(RS,PDELV)
C
                  IC = COP(CX,CY,CZ)
                  I1SC = COP(CX-S(1),CY-S(2),CZ-S(3))
                  DO 10 N=0,ULLR-ILC-LB-LA
                    RRERI2 = 0.0
                    RRERI3 = 0.0
                    RRERI4 = 0.0
C
                    RRINTC2 = 0.0; RRINTC3 = 0.0; RRINTC4 = 0.0
C
                    IF (DELK(1) /= 0) THEN
                      I2SC = COP(CX-2*S(1),CY-2*S(2),CZ-2*S(3))
                      RRERI2 = RRERI(IA,IB,I2SC,1,N+1) - RHO/NIPR*
     $                         RRERI(IA,IB,I2SC,1,N+2)
                      RRINTC2(:) = RRINT(IA,IB,I2SC,1,N,:) - RHO/NIPR*
     $                             RRINT(IA,IB,I2SC,1,N+1,:)
                    END IF
                    IF (DELK(2) /= 0) THEN
                      I1SA = COP(AX,AY,AZ)
                      I1SB = COP(BX-S(1),BY-S(2),BZ-S(3))
                      RRERI3 = RRERI(I1SA,I1SB,I1SC,1,N+2)
                      RRINTC3(:) = RRINT(I1SA,I1SB,I1SC,1,N+1,:)
                    END IF
                    IF (DELK(3) /= 0) THEN
                      I1SA = COP(AX-S(1),AY-S(2),AZ-S(3))
                      I1SB = COP(BX,BY,BZ)
                      RRERI4 = RRERI(I1SA,I1SB,I1SC,1,N+2)
                      RRINTC4(:) = RRINT(I1SA,I1SB,I1SC,1,N+1,:)
                    END IF
                    RRERI(IA,IB,IC,1,N+1) = DQC(COOR)*
     $                                      RRERI(IA,IB,I1SC,1,N+1) +
     $                                      DWQ(COOR)*
     $                                      RRERI(IA,IB,I1SC,1,N+2) +
     $                                      0.5*DELK(1)/NIPR*RRERI2 +
     $                                      0.5*DELK(2)/ZETNI*RRERI3 +
     $                                      0.5*DELK(3)/ZETNI*RRERI4
C
                    DO U=1,3
                      RRVEA(U) = 0.0
                      IF ((AX /= 0).AND.(PDELV(U,1) /= 0.0)) THEN
                        IIA = COP(AX-1,AY,AZ)
                        IIB = COP(BX,BY,BZ)
                        RRVEA(U) = AX*PDELV(U,1)*
     $                  RRERI(IIA,IIB,I1SC,1,N+1)
                      END IF
                      IF ((AY /= 0).AND.(PDELV(U,2) /= 0.0)) THEN
                        IIA = COP(AX,AY-1,AZ)
                        IIB = COP(BX,BY,BZ)
                        RRVEA(U) = RRVEA(U) + AY*PDELV(U,2)*
     $                  RRERI(IIA,IIB,I1SC,1,N+1)
                      END IF
                      IF ((AZ /= 0).AND.(PDELV(U,3) /= 0.0)) THEN
                        IIA = COP(AX,AY,AZ-1)
                        IIB = COP(BX,BY,BZ)
                        RRVEA(U) = RRVEA(U) + AZ*PDELV(U,3)*
     $                  RRERI(IIA,IIB,I1SC,1,N+1)
                      END IF
C
                      RRVEB(U) = 0.0
                      IF ((BX /= 0).AND.(PDELV(U,1) /= 0.0)) THEN
                        IIA = COP(AX,AY,AZ)
                        IIB = COP(BX-1,BY,BZ)
                        RRVEB(U) = BX*PDELV(U,1)*
     $                  RRERI(IIA,IIB,I1SC,1,N+1)
                      END IF
                      IF ((BY /= 0).AND.(PDELV(U,2) /= 0.0)) THEN
                        IIA = COP(AX,AY,AZ)
                        IIB = COP(BX,BY-1,BZ)
                        RRVEB(U) = RRVEB(U) + BY*PDELV(U,2)*
     $                  RRERI(IIA,IIB,I1SC,1,N+1) 
                      END IF
                      IF ((BZ /= 0).AND.(PDELV(U,3) /= 0.0)) THEN
                        IIA = COP(AX,AY,AZ)
                        IIB = COP(BX,BY,BZ-1)
                        RRVEB(U) = RRVEB(U) + BZ*PDELV(U,3)*
     $                  RRERI(IIA,IIB,I1SC,1,N+1)
                      END IF
                    END DO
C
                    RRINT(IA,IB,IC,1,N,:)=DQC(COOR)*
     $                                    RRINT(IA,IB,I1SC,1,N,:) +
     $                                    DWQ(COOR)*
     $                                    RRINT(IA,IB,I1SC,1,N+1,:) +
     $                                    0.5*DELK(1)/NIPR*RRINTC2(:)+
     $                                    0.5*DELK(2)/ZETNI*RRINTC3(:)+
     $                                    0.5*DELK(3)/ZETNI*RRINTC4(:)+
     $                                    2.0*ZETJG*ZETPR/ZETNI*PDEL(:)*
     $                                    RRERI(IA,IB,I1SC,1,N+1) + 
     $                                    ZETJG/ZETNI*RRVEA(:) - 
     $                                    ZETIG/ZETNI*RRVEB(:)
C
                    IF ((ILA == LA).AND.(ILB == LB)
     $                               .AND.(ILC == LC)
     $                               .AND.(LD == 0)) THEN
                      ISA = SOP(AX,AY,AZ)
                      ISB = SOP(BX,BY,BZ)
                      ISC = SOP(CX,CY,CZ)
                      BLOCKS(ISA,ISB,ISC,1,:)=BLOCKS(ISA,ISB,ISC,1,:) +
     $                                        RRINT(IA,IB,IC,1,N,:)
                      IF (N /= 0) THEN 
                        STOP 'ERROR IN RECUR. C'
                      END IF
                    END IF
C
   10             CONTINUE
   20           CONTINUE
   30         CONTINUE
   40       CONTINUE
   50     CONTINUE
   60   CONTINUE
   70 CONTINUE
C
C     O_______________________________________________________________O
C
      END SUBROUTINE RR23
