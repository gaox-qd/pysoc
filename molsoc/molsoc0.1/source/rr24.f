      SUBROUTINE RR24(BLOCKS,RRINT,RRERI,ILA,ILB,ILC,ILD,LA,LB,LC,LD,
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
      INTEGER :: AX,AY,AZ,BX,BY,BZ,CX,CY,COOR,CZ,ILA,ILB,ILC,ILD,LA,LB,
     $           DX,DY,DZ,LC,LD,N,U,ULLR
      INTEGER :: IA,IB,IC,ID,I1SA,I1SB,I1SC,I1SD,I2SD,ISA,ISB,ISC,ISD,
     $           IIA,IIB
      INTEGER, DIMENSION(3) :: S
      INTEGER, DIMENSION(4) :: DELK
C
      REAL :: NIPR,RHO,RRERI2,RRERI3,RRERI4,RRERI5,ZETIG,ZETJG,ZETNI,
     $        ZETPR
      REAL, DIMENSION(3) :: DPA,DPB,DQC,DQD,DWQ,PDEL,RRINTD2,
     $                      RRINTD3,RRINTD4,RRINTD5,RRINTC2,RRINTC3,
     $                      RRINTC4,RRINTC5,RRVEA,RRVEB,RRVEC,RRVED,RS
      REAL, DIMENSION(3,3) :: PDELV
      REAL, DIMENSION(DSH,DSH,DSH,DSH,3) :: BLOCKS
      REAL, DIMENSION(DCPA,DCPA,DCPC,DCPC,0:ULLR,3) :: RRINT
      REAL, DIMENSION(DCPA,DCPA,DCPC,DCPC,0:ULLR+1) :: RRERI
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ D horizontal recurrence relation steps ___
C
      DO 90 AX=ILA,0,-1
        DO 80 AY=ILA-AX,0,-1
          AZ = ILA - AX - AY
          DO 70 BX=ILB,0,-1
            DO 60 BY=ILB-BX,0,-1
              BZ = ILB - BX - BY
              IA = COP(AX,AY,AZ)
              IB = COP(BX,BY,BZ)
              DO 50 CX=ILC,0,-1
                DO 40 CY=ILC-CX,0,-1
                  CZ = ILC - CX - CY
                  DO 30 DX=ILD,0,-1
                    DO 20 DY=ILD-DX,0,-1
                      DZ = ILD - DX - DY
                      CALL SCANGMF(AX,AY,AZ,BX,BY,BZ,CX,CY,CZ,DX,DY,DZ,
     $                             COOR,S,DELK)
                      RS(:) = REAL(S(:))
                      CALL VPROD(RS,DPB,PDEL)
                      CALL VECFSO2(RS,PDELV)
C
                      IC = COP(CX,CY,CZ)
                      ID = COP(DX,DY,DZ)
                      I1SD = COP(DX-S(1),DY-S(2),DZ-S(3))
                      DO 10 N=0,ULLR-ILD-LC-LB-LA
                        RRERI2 = 0.0; RRERI3 = 0.0
                        RRERI4 = 0.0; RRERI5 = 0.0
C
                        RRINTD2 = 0.0; RRINTD3 = 0.0
                        RRINTD4 = 0.0; RRINTD5 = 0.0
C
                        IF (DELK(1) /= 0) THEN
                          I2SD = COP(DX-2*S(1),DY-2*S(2),DZ-2*S(3))
                          RRERI2 = RRERI(IA,IB,IC,I2SD,N+1) - RHO/NIPR*
     $                             RRERI(IA,IB,IC,I2SD,N+2) 
                          RRINTD2(:) =RRINT(IA,IB,IC,I2SD,N,:)-RHO/NIPR*
     $                                RRINT(IA,IB,IC,I2SD,N+1,:)
                        END IF
                        IF (DELK(2) /= 0) THEN
                          I1SC = COP(CX-S(1),CY-S(2),CZ-S(3))
                          RRERI3 = RRERI(IA,IB,I1SC,I1SD,N+1) -RHO/NIPR*
     $                             RRERI(IA,IB,I1SC,I1SD,N+2)
                          RRINTD3(:) = RRINT(IA,IB,I1SC,I1SD,N,:)-
     $                                 RHO/NIPR*
     $                                 RRINT(IA,IB,I1SC,I1SD,N+1,:)
                        END IF
                        IF (DELK(3) /= 0) THEN
                          I1SA = COP(AX,AY,AZ)
                          I1SB = COP(BX-S(1),BY-S(2),BZ-S(3))
                          RRERI4 = RRERI(I1SA,I1SB,IC,I1SD,N+2)
                          RRINTD4(:) = RRINT(I1SA,I1SB,IC,I1SD,N+1,:)
                        END IF
                        IF (DELK(4) /= 0) THEN
                          I1SA = COP(AX-S(1),AY-S(2),AZ-S(3))
                          I1SB = COP(BX,BY,BZ)
                          RRERI5 = RRERI(I1SA,I1SB,IC,I1SD,N+2)
                          RRINTD5(:) = RRINT(I1SA,I1SB,IC,I1SD,N+1,:) 
                        END IF
                        RRERI(IA,IB,IC,ID,N+1)=DQD(COOR)*
     $                                         RRERI(IA,IB,IC,I1SD,N+1)+
     $                                         DWQ(COOR)*
     $                                         RRERI(IA,IB,IC,I1SD,N+2)+
     $                                         0.5*DELK(1)/NIPR*RRERI2 +
     $                                         0.5*DELK(2)/NIPR*RRERI3 +
     $                                         0.5*DELK(3)/ZETNI*RRERI4+
     $                                         0.5*DELK(4)/ZETNI*RRERI5
C
                        DO U=1,3
                          RRVEA(U) = 0.0
                          IF ((AX /= 0).AND.(PDELV(U,1) /= 0.0)) THEN
                            IIA = COP(AX-1,AY,AZ)
                            IIB = COP(BX,BY,BZ)
                            RRVEA(U) = AX*PDELV(U,1)*
     $                      RRERI(IIA,IIB,IC,I1SD,N+1)
                          END IF
                          IF ((AY /= 0).AND.(PDELV(U,2) /= 0.0)) THEN
                            IIA = COP(AX,AY-1,AZ)
                            IIB = COP(BX,BY,BZ)
                            RRVEA(U) = RRVEA(U) + AY*PDELV(U,2)*
     $                      RRERI(IIA,IIB,IC,I1SD,N+1)
                          END IF
                          IF ((AZ /= 0).AND.(PDELV(U,3) /= 0.0)) THEN
                            IIA = COP(AX,AY,AZ-1)
                            IIB = COP(BX,BY,BZ)
                            RRVEA(U) = RRVEA(U) + AZ*PDELV(U,3)*
     $                      RRERI(IIA,IIB,IC,I1SD,N+1)
                          END IF
                          RRVEB(U) = 0.0
                          IF ((BX /= 0).AND.(PDELV(U,1) /= 0.0)) THEN
                            IIA = COP(AX,AY,AZ)
                            IIB = COP(BX-1,BY,BZ)
                            RRVEB(U) = BX*PDELV(U,1)*
     $                      RRERI(IIA,IIB,IC,I1SD,N+1)
                          END IF
                          IF ((BY /= 0).AND.(PDELV(U,2) /= 0.0)) THEN
                            IIA = COP(AX,AY,AZ)
                            IIB = COP(BX,BY-1,BZ)
                            RRVEB(U) = RRVEB(U) + BY*PDELV(U,2)*
     $                      RRERI(IIA,IIB,IC,I1SD,N+1)
                          END IF
                          IF ((BZ /= 0).AND.(PDELV(U,3) /= 0.0)) THEN
                            IIA = COP(AX,AY,AZ)
                            IIB = COP(BX,BY,BZ-1)
                            RRVEB(U) = RRVEB(U) + BZ*PDELV(U,3)*
     $                      RRERI(IIA,IIB,IC,I1SD,N+1)
                          END IF
                        END DO
C
                        RRINT(IA,IB,IC,ID,N,:)=DQD(COOR)*
     $                                     RRINT(IA,IB,IC,I1SD,N,:)+
     $                                     DWQ(COOR)*
     $                                     RRINT(IA,IB,IC,I1SD,N+1,:)+
     $                                     0.5*DELK(1)/NIPR*RRINTD2(:)+
     $                                     0.5*DELK(2)/NIPR*RRINTD3(:)+
     $                                     0.5*DELK(3)/ZETNI*RRINTD4(:)+
     $                                     0.5*DELK(4)/ZETNI*RRINTD5(:)+
     $                                     2.0*ZETJG*ZETPR/ZETNI*
     $                                     PDEL(:)*
     R                                     RRERI(IA,IB,IC,I1SD,N+1) + 
     $                                     ZETJG/ZETNI*RRVEA(:) -
     $                                     ZETIG/ZETNI*RRVEB(:)
C
                        IF ((ILA == LA).AND.(ILB == LB)
     $                                   .AND.(ILC == LC)
     $                                   .AND.(ILD == LD)) THEN
                          ISA = SOP(AX,AY,AZ)
                          ISB = SOP(BX,BY,BZ)
                          ISC = SOP(CX,CY,CZ)
                          ISD = SOP(DX,DY,DZ)
                          BLOCKS(ISA,ISB,ISC,ISD,:) = 
     $                    BLOCKS(ISA,ISB,ISC,ISD,:) +
     $                    RRINT(IA,IB,IC,ID,N,:)
                          IF (N /= 0) THEN 
                            STOP 'ERROR IN RECUR. D'
                          END IF
                        END IF
C
   10                 CONTINUE
   20               CONTINUE
   30             CONTINUE
   40           CONTINUE
   50         CONTINUE
   60       CONTINUE
   70     CONTINUE
   80   CONTINUE
   90 CONTINUE
C
C     O_______________________________________________________________O
C
      END SUBROUTINE RR24
