      SUBROUTINE RRMOMINT(LA,LB,LU,SHELLA,SHELLB,BLOCKS,DMOM)
C
C     ******************************************************************
C     ***  Calculation of RRint field elements for MOMent INtegrals  ***
C     ******************************************************************
C
C     Creation (03.04.07, SC)
C
C     List of local variables:
C
C     AX,AY,AZ  : Angular momentum index of orbital A.
C     BX,BY,BZ  : Angular momentum index of orbital B.
C     COOR{A/B} : Pointer position of the angular momentum to be scaled.
C     DELK      : Like Kronecker's delta variable.
C     GLL{A/B}  : Lower GTF limit of shells.
C     GUL{A/B}  : Upper GTF limit of shells.
C     PCOV      : Primitive Cartesian overlap integrals.
C     SHELL{A/B}: Orbital shells of A and B.
C     BLOCKS    : Diatomic integral shell block.
C     S         : Scaling number.
C     ULAB      : Upper orbital L quantum number.
c     ULLR      : Upper limit LA + LB.
C     ZETPR     : Sum of Gaussian exponent. 
C
C     List of local dynamical fields:
C
C     RRINT : Work field for recurrence relations.
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
C
      INCLUDE 'fileio.h'
      INCLUDE 'math.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      INTEGER :: DMOM,DCLA,DCLB
      INTEGER :: IG,JG,ILA,ILU,ILB,LA,LB,LU,LLAB,GLLA,GLLB,N,SHELLA,
     $           SHELLB,ULAB,GULA,GULB,ULLR,AX,AY,AZ,BX,BY,BZ,UX,UY,UZ,
     $           COOR
      INTEGER, DIMENSION(3) :: S,DELK
C
      REAL :: PCSBOUND,PCSQ,PCOV,SPG,XI,T,ZETPR,RRINT2,RRINT3,RRINTU,
     $        ADDEXP,ADDFUN
      REAL, DIMENSION(3) :: DAB,DAC,DBC,DPA,DPB,DPC
      REAL, DIMENSION(DSHL,DSHL,2:DMOM) :: BLOCKS
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: RRINT
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Check lower and upper limits and distance vectors ___
C
      IF (LA >= LB) THEN
        GLLA = GLL(SHELLA)
        GULA = GUL(SHELLA)
        GLLB = GLL(SHELLB)
        GULB = GUL(SHELLB)
        LLAB   = LB
        ULAB   = LA
        DAB(:) = RAB(:)
        DAC(:) = RAC(:)
        DBC(:) = RBC(:)
        DCLA = (LA**3 - LA)/6 + LA**2 + 2*LA + 1
        DCLB = (LB**3 - LB)/6 + LB**2 + 2*LB + 1
      ELSE
        GLLA = GLL(SHELLB)
        GULA = GUL(SHELLB)
        GLLB = GLL(SHELLA)
        GULB = GUL(SHELLA)
        LLAB   = LA
        ULAB   = LB
        DAB(:) = -RAB(:)
        DAC(:) = RBC(:)
        DBC(:) = RAC(:)
        DCLA = (LB**3 - LB)/6 + LB**2 + 2*LB + 1
        DCLB = (LA**3 - LA)/6 + LA**2 + 2*LA + 1
      END IF
C
      ULLR = LA + LB
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Initialization ___
C
      BLOCKS = 0.0
C
C     ___ Allocate local field ___
C
      ALLOCATE(RRINT(DCLA,DCLB,DMOM),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'RRMOMINT: ALLOCATION FAILED'
        STOP
      END IF
C
C     ___ GTF loop ___
C
      DO 150 IG=GLLA,GULA
        DO 140 JG=GLLB,GULB
          ZETPR = ZET(IG) + ZET(JG)
          XI = ZET(IG)*ZET(JG)/ZETPR
          PCOV = 2.0*PI*EXP(-XI*RSQAB)/ZETPR
          PCSBOUND = GCC(IG)*GCC(JG)*PCOV 
          IF (ABS(PCSBOUND) < INTHRESH) GOTO 140
C
C     ___ Calculate the P - A and P - B factors ___
C
          DPA(:) = -ZET(JG)*DAB(:)/ZETPR 
          DPB(:) = ZET(IG)*DAB(:)/ZETPR 
          DPC(:) = (ZET(IG)*DAC(:) + ZET(JG)*DBC(:))/ZETPR
C
          RRINT(1,1,1) = 0.5*GCC(IG)*GCC(JG)*
     $                   PCOV*(PI/ZETPR)**(1.0/2.0)
C
C     ___ Calculate moment integrals over s functions ___
C
          DO ILU=1,LU 
            DO UX=ILU,0,-1
              DO UY=ILU-UX,0,-1
                UZ = ILU - UX - UY
                CALL SCANGMT(0,0,0,0,0,0,UX,UY,UZ,COOR,S,DELK)
                RRINTU = 0.0
                IF (DELK(1) /= 0) THEN
                  RRINTU = RRINT(1,1,COP(UX-2*S(1),
     $            UY-2*S(2),UZ-2*S(3)))
                END IF
                RRINT(1,1,COP(UX,UY,UZ)) = DPC(COOR)*
     $          RRINT(1,1,COP(UX-S(1),UY-S(2),UZ-S(3))) +
     $          0.5*DELK(1)/ZETPR*RRINTU 
                IF (ULLR == 0) THEN
                  BLOCKS(1,1,COP(UX,UY,UZ)) = 
     $            BLOCKS(1,1,COP(UX,UY,UZ)) + 
     $            RRINT(1,1,COP(UX,UY,UZ))
                END IF 
              END DO
            END DO
          END DO
C
          IF (ULLR == 0) GOTO 140
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Build basic cartesian integrals ___
C
C     ___ Left recurrence relation steps ___
C
          DO 60 ILA=1,ULAB
            DO 55 ILU=0,LU
              DO 50 AX=ILA,0,-1
                DO 40 AY=ILA-AX,0,-1
                  AZ = ILA - AX - AY
                  DO 35 UX=ILU,0,-1
                    DO 30 UY=ILU-UX,0,-1
                      UZ = ILU - UX - UY
                      CALL SCANGMT(0,0,0,UX,UY,UZ,AX,AY,AZ,COOR,
     $                            S,DELK)
                      RRINT2 = 0.0; RRINTU = 0.0
                      IF (DELK(1) /= 0) THEN
                        RRINT2 = RRINT(COP(AX-2*S(1),AY-2*S(2),
     $                  AZ-2*S(3)),1,COP(UX,UY,UZ)) 
                      END IF
                      IF (DELK(2) /= 0) THEN
                        RRINTU = RRINT(COP(AX-S(1),AY-S(2),AZ-S(3)),
     $                  1,COP(UX-S(1),UY-S(2),UZ-S(3)))
                      END IF
                      RRINT(COP(AX,AY,AZ),1,COP(UX,UY,UZ)) = 
     $                DPA(COOR)*RRINT(COP(AX-S(1),AY-S(2),
     $                AZ-S(3)),1,COP(UX,UY,UZ)) + 
     $                0.5*DELK(1)/ZETPR*RRINT2 +
     $                0.5*DELK(2)/ZETPR*RRINTU
C
C      ___ Build shell integral blocks ___
C
                      IF ((ILA == ULAB).AND.
     $                    (LLAB == 0).AND.
     $                    (ILU > 0)) THEN
                        BLOCKS(SOP(AX,AY,AZ),1,COP(UX,UY,UZ)) = 
     $                  BLOCKS(SOP(AX,AY,AZ),1,COP(UX,UY,UZ)) +
     $                  RRINT(COP(AX,AY,AZ),1,COP(UX,UY,UZ))
                      END IF
   30               CONTINUE
   35             CONTINUE
   40           CONTINUE
   50         CONTINUE
   55       CONTINUE
   60     CONTINUE
C
C     ___ Right recurrence relation steps ___
C
          DO 130 ILA=1,ULAB
C
            DO 120 ILB=1,ILA
              IF (ILB > LLAB) GOTO 130
C
              DO 115 ILU=0,LU
                DO 110 AX=ILA,0,-1
                  DO 100 AY=ILA-AX,0,-1
                    AZ = ILA - AX - AY
                    DO 90 BX=ILB,0,-1
                      DO 80 BY=ILB-BX,0,-1
                        BZ = ILB - BX - BY
                        DO 70 UX=ILU,0,-1
                          DO 65 UY=ILU-UX,0,-1
                            UZ = ILU - UX - UY
                            CALL SCANGMT(UX,UY,UZ,AX,AY,AZ,
     $                                   BX,BY,BZ,COOR,S,DELK)
                            RRINT2 = 0.0; RRINT3 = 0.0 
                            RRINTU = 0.0
                            IF (DELK(1) /= 0) THEN
                              RRINT2 = RRINT(COP(AX,AY,AZ),COP(BX-
     $                        2*S(1),BY-2*S(2),BZ-2*S(3)),
     $                        COP(UX,UY,UZ))
                            END IF
                            IF (DELK(2) /= 0) THEN
                              RRINT3 = RRINT(COP(AX-S(1),AY-S(2),AZ-
     $                        S(3)),COP(BX-S(1),BY-S(2),BZ-S(3)),
     $                        COP(UX,UY,UZ)) 
                            END IF
                            IF (DELK(3) /= 0) THEN
                              RRINTU = RRINT(COP(AX,AY,AZ),
     $                        COP(BX-S(1),BY-S(2),BZ-S(3)),
     $                        COP(UX-S(1),UY-S(2),UZ-S(3)))
                            END IF
                            RRINT(COP(AX,AY,AZ),COP(BX,BY,BZ),
     $                      COP(UX,UY,UZ)) = DPB(COOR)*
     $                      RRINT(COP(AX,AY,AZ),COP(BX-S(1),BY-S(2),
     $                      BZ-S(3)),COP(UX,UY,UZ)) + 
     $                      0.5*DELK(1)/ZETPR*RRINT2 +
     $                      0.5*DELK(2)/ZETPR*RRINT3 +
     $                      0.5*DELK(3)/ZETPR*RRINTU
C
C      ___ Build shell integral blocks ___
C
                            IF (ILA == ULAB.AND.ILB == LLAB.AND.
     $                          ILU > 0) THEN
                              BLOCKS(SOP(AX,AY,AZ),SOP(BX,BY,BZ),
     $                        COP(UX,UY,UZ)) = 
     $                        BLOCKS(SOP(AX,AY,AZ),SOP(BX,BY,BZ),
     $                        COP(UX,UY,UZ)) +
     $                        RRINT(COP(AX,AY,AZ),COP(BX,BY,BZ),
     $                        COP(UX,UY,UZ))
                            END IF
   65                     CONTINUE
   70                   CONTINUE
   80                 CONTINUE
   90               CONTINUE
  100             CONTINUE
  110           CONTINUE
  115         CONTINUE
  120       CONTINUE
  130     CONTINUE
C
C     ___ End of recurrence relation steps ___
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
  140   CONTINUE
  150 CONTINUE
C
C     ___ Deallocate local field ___
C
      DEALLOCATE(RRINT,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'RRMOMINT: DEALLOCATION FAILED'
        STOP
      END IF
C
C     ___ Transpose BLOCKS, if LA < LB ___
C
      IF (LA < LB) THEN
        DO ILU=2,DMOM
          BLOCKS(:,:,ILU) = TRANSPOSE(BLOCKS(:,:,ILU))
        END DO
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE RRMOMINT 
