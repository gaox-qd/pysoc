      SUBROUTINE RR1SOCIN(LA,LB,SHELLA,SHELLB,BLOCKS)
C
C     ******************************************************************
C     ***      Calculation of RR field elements for 1 electron       ***
C     ***      Spin Orbit Coupling INtegrals.                        ***
C     ******************************************************************
C
C     List of local variables:
C
C     AX,AY,AZ  : Angular momentum index of orbital A.
C     BX,BY,BZ  : Angular momentum index of orbital B.
C     COOR{A/B} : Pointer position of the angular momentum to be scaled.
C     DELK      : Like Kronecker's delta variable.
C     F(N)      : Incomplete Gamma function values.
C     GLL{A/B}  : Lower GTF limit of shells.
C     GUL{A/B}  : Upper GTF limit of shells.
C     PCOV      : Primitive Cartesian overlap integrals.
C     RRINT     : Work field for recurrence relations.
C     SHELL{A/B}: Orbital shells of A and B.
C     BLOCKS    : Diatomic integral shell block.
C     S         : Scaling number.
C     T         : Argument of incomplete Gamma function.
C     ULAB      : Upper orbital L quantum number.
c     ULLR      : Upper limit LA + LB.
C     ZETPRE    : Product of Gaussian exponent. 
C     ZETPR     : Sum of Gaussian exponent.
C
C     List of local dynamical fields:
C
C     RRINT   : Work field for SO recurrence relations.
C     RRINTNA : Work field for Nuclear attraction recurrence relations.
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
      INTEGER :: DCLA,DCLB
      INTEGER :: IG,JG,ILA,ILB,LA,LB,LLAB,GLLA,GLLB,N,U
      INTEGER :: SHELLA,SHELLB,ULAB,GULA,GULB,ULLR
      INTEGER :: AX,AY,AZ,BX,BY,BZ,COOR
      INTEGER, DIMENSION(3) :: S,DELK

      REAL :: PCSBOUND,PCSQ,PCOV,SPG,XI,T,ZETPR,ZETPRE,RRVEC,
     $        RRINTNA2,RRINTNA3
      REAL, DIMENSION(3) :: CX,DAB,DAC,DBC,DPA,DPB,DPC,PDEL,
     $                      PDEL2,DACR,DBCR,RS,RRINT2,RRINT3
      REAL, DIMENSION(DSHL,DSHL,3) :: BLOCKS
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:), ALLOCATABLE :: F
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: RRINTNA
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: RRINT
C
      LOGICAL :: LEFT,RIGHT
C
C     o~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~o
C
C     ___ Check lower and upper limits and distance vectors ___
C
      GLLA = GLL(SHELLA)
      GULA = GUL(SHELLA)
      GLLB = GLL(SHELLB)
      GULB = GUL(SHELLB)
      DAB(:) = RAB(:)
      DAC(:) = RAC(:)
      DBC(:) = RBC(:)
C
      LEFT = LA.GE.LB
      RIGHT = LA.LT.LB
C
      IF (LEFT) THEN
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
      ELSE IF (RIGHT) THEN
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
C     o~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~o
C
C     ___ Initialization ___
C
      BLOCKS = 0.0
C
C     ___ Allocate local fields ___
C
      ALLOCATE(RRINT(DCLA,DCLB,0:ULLR,3),
     $         RRINTNA(DCLA,DCLB,0:ULLR+1),F(0:ULLR+1),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'RR1SOCIN, ALLOCATION FAILED'
        STOP 
      END IF
C
C     ___ Fields initialization ___
C
      RRINT = 0.0; RRINTNA = 0.0
C
C     ___ GTF product loop ___
C
      DO 150 IG=GLLA,GULA
        DO 140 JG=GLLB,GULB
          ZETPRE = ZET(IG)*ZET(JG)
          ZETPR = ZET(IG) + ZET(JG)
          XI = ZET(IG)*ZET(JG)/ZETPR
          PCOV = 2.0*PI*EXP(-XI*RSQAB)/ZETPR
          PCSBOUND = GCC(IG)*GCC(JG)*PCOV 
          IF (ABS(PCSBOUND).LT.INTHRESH) GOTO 140
C
C     ___ Calculate the P - A and P - B factors ___
C
          DPA(:) = -ZET(JG)*DAB(:)/ZETPR 
          DPB(:) = ZET(IG)*DAB(:)/ZETPR 
C
C     ___ Calculate the P - C factor ___
C
          DPC(:) = (ZET(IG)/ZETPR*DAC(:)) +
     $             (ZET(JG)/ZETPR*DBC(:))
          PCSQ = SUM(DPC**2,1)
C
C     ___ Calculate the incomplete Gamma function ___
C
          T = ZETPR*PCSQ
          CALL GAMMAF(F,T,ULLR+1)
C
C     ___ Build primitive integrals over s functions ___
C
          RRINTNA(1,1,0) = GCC(IG)*GCC(JG)*F(0)*PCOV
          DO N=0,ULLR
            RRINTNA(1,1,N+1) = GCC(IG)*GCC(JG)*F(N+1)*PCOV
            CALL VPROD(DAC,DBC,CX)
            RRINT(1,1,N,:) = 4.0*ZETPRE*CX(:)*RRINTNA(1,1,N+1)
          END DO
C
          IF (ULLR == 0) THEN
            DO U=1,3
              BLOCKS(1,1,U) = BLOCKS(1,1,U) + RRINT(1,1,0,U)
            END DO
            GOTO 140
          END IF
C
C     o~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~o
C
C     ___ Build basic cartesian integrals ___
C
C     ___ Left or right recurrence relation steps ___
C
          DO 60 ILA=1,ULAB
            DO 50 AX=ILA,0,-1
              DO 40 AY=ILA-AX,0,-1
                AZ = ILA - AX - AY
                CALL SCANGM(0,0,0,AX,AY,AZ,COOR,S,DELK)
                RS(:) = REAL(S(:))
                CALL VPROD(RS,DBC,PDEL)
                DO 30 N=0,ULLR-ILA
                  RRINT2 = 0.0
                  RRINTNA2 = 0.0
                  IF (DELK(2) /= 0) THEN
                    RRINTNA2 = RRINTNA(COP(AX-2*S(1),AY-2*S(2),
     $              AZ-2*S(3)),1,N+1)-RRINTNA(COP(AX-2*S(1),AY-
     $              2*S(2),AZ-2*S(3)),1,N+2)
                  END IF
                  RRINTNA(COP(AX,AY,AZ),1,N+1) = DPA(COOR)*
     $            RRINTNA(COP(AX-S(1),AY-S(2),AZ-S(3)),1,N+1) -
     $            DPC(COOR)*RRINTNA(COP(AX-S(1),AY-S(2),AZ-S(3)),
     $            1,N+2) + 0.5*REAL(DELK(2))/ZETPR*RRINTNA2
                  DO U=1,3
                    IF (DELK(2) /= 0) THEN
                      RRINT2(U) = RRINT(COP(AX-2*S(1),AY-2*S(2),
     $                AZ-2*S(3)),1,N,U)-RRINT(COP(AX-2*S(1),AY-
     $                2*S(2),AZ-2*S(3)),1,N+1,U)
                    END IF
                    RRINT(COP(AX,AY,AZ),1,N,U) = DPA(COOR)*
     $              RRINT(COP(AX-S(1),AY-S(2),AZ-S(3)),1,N,U) -
     $              DPC(COOR)*RRINT(COP(AX-S(1),AY-S(2),AZ-S(3)),
     $              1,N+1,U) + 0.5*REAL(DELK(2))/ZETPR*RRINT2(U) +
     $              2.0*ZET(JG)*PDEL(U)*RRINTNA(COP(AX-S(1),
     $              AY-S(2),AZ-S(3)),1,N+1)
                  END DO
C
C      ___ Build shell integral blocks ___
C
                  IF (ILA == ULAB.AND.LLAB == 0) THEN
                    DO U=1,3
                      BLOCKS(SOP(AX,AY,AZ),1,U) = 
     $                BLOCKS(SOP(AX,AY,AZ),1,U) +
     $                RRINT(COP(AX,AY,AZ),1,N,U)
                    END DO
                  END IF
   30           CONTINUE
   40         CONTINUE
   50       CONTINUE
   60     CONTINUE
C
C     ___ Right or left recurrence relation steps ___
C
          DO 130 ILA=1,ULAB
            DO 120 ILB=1,ILA
              IF (ILB > LLAB) GOTO 130
              DO 110 AX=ILA,0,-1
                DO 100 AY=ILA-AX,0,-1
                  AZ = ILA - AX - AY
                  DO 90 BX=ILB,0,-1
                    DO 80 BY=ILB-BX,0,-1
                      BZ = ILB - BX - BY
                      CALL SCANGM(AX,AY,AZ,BX,BY,BZ,COOR,S,DELK)
                      RS(:) = REAL(S(:))
                      CALL VPROD(RS,DAC,PDEL)
                      DO 71 N=0,LLAB-ILB
                        RRINT2 = 0.0; RRINT3 = 0.0
                        RRINTNA2 = 0.0; RRINTNA3 = 0.0
                        IF (DELK(1).NE.0) THEN
                          RRINTNA2 = RRINTNA(COP(AX-S(1),AY-S(2),AZ-
     $                    S(3)),COP(BX-S(1),BY-S(2),BZ-S(3)),N+1)-
     $                    RRINTNA(COP(AX-S(1),AY-S(2),AZ-S(3)),
     $                    COP(BX-S(1),BY-S(2),BZ-S(3)),N+2)
                        END IF
                        IF (DELK(2) /= 0) THEN
                          RRINTNA3 = RRINTNA(COP(AX,AY,AZ),COP(BX-
     $                    2*S(1),BY-2*S(2),BZ-2*S(3)),N+1) -
     $                    RRINTNA(COP(AX,AY,AZ),COP(BX-2*S(1),BY-
     $                    2*S(2),BZ-2*S(3)),N+2)
                        END IF
                        RRINTNA(COP(AX,AY,AZ),COP(BX,BY,BZ),N+1) =
     $                  DPB(COOR)*RRINTNA(COP(AX,AY,AZ),COP(BX-S(1),
     $                  BY-S(2),BZ-S(3)),N+1) - DPC(COOR)*
     $                  RRINTNA(COP(AX,AY,AZ),COP(BX-S(1),
     $                  BY-S(2),BZ-S(3)),N+2) +
     $                  0.5*REAL(DELK(1))/ZETPR*RRINTNA2 +
     $                  0.5*REAL(DELK(2))/ZETPR*RRINTNA3
   71                 CONTINUE
                      DO 70 N=0,LLAB-ILB
                        DO U=1,3
                          IF (DELK(1) /= 0) THEN
                            RRINT2(U) = RRINT(COP(AX-S(1),AY-S(2),
     $                      AZ-S(3)),COP(BX-S(1),BY-S(2),BZ-S(3)),
     $                      N,U)-RRINT(COP(AX-S(1),AY-S(2),
     $                      AZ-S(3)),COP(BX-S(1),BY-S(2),BZ-S(3)),
     $                      N+1,U)
                          END IF
                          IF (DELK(2) /= 0) THEN
                            RRINT3(U) = RRINT(COP(AX,AY,AZ),COP(BX-
     $                      2*S(1),BY-2*S(2),BZ-2*S(3)),N,U)-
     $                      RRINT(COP(AX,AY,AZ),COP(BX-2*S(1),BY-
     $                      2*S(2),BZ-2*S(3)),N+1,U)
                          END IF
                          RS(:) = REAL(S(:))
                          CALL VECFSO(RS,PDEL2,U)
                          RRVEC = 0.0
                          IF ((AX /= 0).AND.(PDEL2(1) /= 0.0)) 
     $                    RRVEC = REAL(AX)*PDEL2(1)*
     $                    RRINTNA(COP(AX-1,AY,AZ),COP(BX-S(1),
     $                    BY-S(2),BZ-S(3)),N+1) 
                          IF ((AY /= 0).AND.(PDEL2(2) /= 0.0))
     $                    RRVEC = RRVEC + REAL(AY)*PDEL2(2)*
     $                    RRINTNA(COP(AX,AY-1,AZ),COP(BX-S(1),
     $                    BY-S(2),BZ-S(3)),N+1)
                          IF ((AZ /= 0).AND.(PDEL2(3) /= 0.0)) 
     $                    RRVEC = RRVEC + REAL(AZ)*PDEL2(3)*
     $                    RRINTNA(COP(AX,AY,AZ-1),COP(BX-S(1),
     $                    BY-S(2),BZ-S(3)),N+1)
C
                          RRINT(COP(AX,AY,AZ),COP(BX,BY,BZ),N,U) =
     $                    DPB(COOR)*RRINT(COP(AX,AY,AZ),COP(BX-S(1),
     $                    BY-S(2),BZ-S(3)),N,U) - DPC(COOR)*
     $                    RRINT(COP(AX,AY,AZ),COP(BX-S(1),BY-S(2),
     $                    BZ-S(3)),N+1,U) +
     $                    0.5*REAL(DELK(1))/ZETPR*RRINT2(U) +
     $                    0.5*REAL(DELK(2))/ZETPR*RRINT3(U) - 
     $                    2.0*ZET(IG)*PDEL(U)*RRINTNA(COP(AX,AY,AZ),
     $                    COP(BX-S(1),BY-S(2),BZ-S(3)),N+1) - 
     $                    RRVEC
                        END DO
C
C      ___ Build shell integral blocks ___
C
                        DO U=1,3
                          IF (ILA == ULAB.AND.ILB == LLAB) THEN
                            BLOCKS(SOP(AX,AY,AZ),SOP(BX,BY,BZ),U) = 
     $                      BLOCKS(SOP(AX,AY,AZ),SOP(BX,BY,BZ),U) + 
     $                      RRINT(COP(AX,AY,AZ),COP(BX,BY,BZ),N,U)
                          END IF
                        END DO
   70                 CONTINUE
   80               CONTINUE
   90             CONTINUE
  100           CONTINUE
  110         CONTINUE
  120       CONTINUE
  130     CONTINUE
C
C     ___ End of recurrence relation steps ___
C
C     o~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~o
C
  140   CONTINUE
  150 CONTINUE
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(RRINT,RRINTNA,F,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'RR1SOCIN, DEALLOCATION FAILED'
        STOP
      END IF
C 
C     ___ Transpose BLOCKS, if LA < LB ___
C
      IF (LA < LB) THEN
        DO U=1,3
          BLOCKS(:,:,U) = -TRANSPOSE(BLOCKS(:,:,U))
        END DO
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE RR1SOCIN
