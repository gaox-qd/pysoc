      SUBROUTINE RR2SOCIN(LA,LB,LC,LD,SHELLA,SHELLB,SHELLC,SHELLD,
     $                    BLOCKS)
C
C     ******************************************************************
C     ***      Calculation of RR field elements for 2 electron       ***
C     ***      Spin Orbit Coupling INtegrals.                        ***
C     ******************************************************************
C
C      Creation (12.08.03, SC)
C               (05.12.04, SC)
C
C     List of local variables:
C
C     DELK      : Like Kronecker's delta variable.
C     F(N)      : Incomplete Gamma function values.
C     GLL{A/B}  : Lower GTF limit of shells.
C     GUL{A/B}  : Upper GTF limit of shells.
C     NIPR      : Sum of Gaussian exponent C and D.
C     PCOV1/2   : Primitive Cartesian overlap integrals.
C     RRINT     : Work field for recurrence relations.
C     SHELL{A/B}: Orbital shells of A and B.
C     BLOCKS    : Diatomic integral shell block.
C     SA/B      : Scaling number.
C     T         : Argument of incomplete Gamma function.
c     ULLR      : Upper limit LA + LB + LC + LD.
C     ZETPRE    : Product of Gaussian exponent A and B. 
C     ZETPR     : Sum of Gaussian exponent A and B.
C
C     List of local dynamical fields:
C
C     RRINT    : Work field for SO recurrence relations.
C     RRERI    : Work field for Nuclear attraction recurrence relations.
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
      INTEGER :: DCPA,DCPB,DCPC,DCPD,THRESHC
      INTEGER :: IG,JG,KG,LG,ILA,ILB,LA,LB,LC,LD,GLLA,GLLB,GLLC,GLLD,N,
     $           SHELLA,SHELLB,SHELLC,SHELLD,GULA,GULB,GULC,GULD,ULLR,
     $           ILD,ILC,LBCD,LCD
      INTEGER, DIMENSION(3) :: SA,SB,SC,SD
      INTEGER, DIMENSION(4) :: DELK
C
      REAL :: PCSBOUNDA,PCSBOUNDB,PCSBOUND,PQSQ,PCOV1,PCOV2,XI,T,ZETPR,
     $        ZETPRE,NIPRE,NIPR,RHO,ZETNI
      REAL, DIMENSION(3) :: DAB,DCD,DPA,DPB,DWP,DWQ,DPQ,DQA,DQB,DQC,
     $                      DQD,TU1,TU2
      REAL, DIMENSION(DSHL,DSHL,DSHL,DSHL,3) :: BLOCKS
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:), ALLOCATABLE :: F
      REAL, DIMENSION(:,:,:,:,:), ALLOCATABLE :: RRERI
      REAL, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: RRINT
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Check lower and upper limits and distance vectors ___
C
      GLLA = GLL(SHELLA)
      GULA = GUL(SHELLA)
      GLLB = GLL(SHELLB)
      GULB = GUL(SHELLB)
      GLLC = GLL(SHELLC)
      GULC = GUL(SHELLC)
      GLLD = GLL(SHELLD)
      GULD = GUL(SHELLD)
      DAB(:) = RAB(:)
      DCD(:) = RCD(:)
C
      ULLR = LA + LB + LC + LD
C
      IF (LA >= LB) THEN
        DCPA = (LA**3 - LA)/6 + LA**2 + 2*LA + 1
      ELSE IF (LA < LB) THEN
        DCPA = (LB**3 - LB)/6 + LB**2 + 2*LB + 1
      END IF
C
      IF (LC >= LD) THEN
        DCPC = (LC**3 - LC)/6 + LC**2 + 2*LC + 1
      ELSE IF (LC < LD) THEN
        DCPC = (LD**3 - LD)/6 + LD**2 + 2*LD + 1
      END IF
C
      LCD = LC + LD
      LBCD = LB + LCD
      ULLR = LA + LB + LC + LD
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Initialization ___
C
      BLOCKS = 0.0
C
C     ___ Allocate local fields ___
C
      ALLOCATE(RRINT(DCPA,DCPA,DCPC,DCPC,0:ULLR,3),
     $         RRERI(DCPA,DCPA,DCPC,DCPC,0:ULLR+1),F(0:ULLR+1),
     $         STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'RR2SOCIN: ALLOCATION FAILED'
        STOP
      END IF
C
C     ___ Fields initialization ___
C
C     RRINT = 0.0; RRERI = 0.0
C
C     ___ GTF loop ___
C
      DO 150 IG=GLLA,GULA
        DO 140 JG=GLLB,GULB
          DO 132 KG=GLLC,GULC
            DO 131 LG=GLLD,GULD
              ZETPRE = ZET(IG)*ZET(JG)
              ZETPR = ZET(IG) + ZET(JG)
              XI = ZET(IG)*ZET(JG)/ZETPR
              PCOV1 = (PI/ZETPR)**(3.0/2.0)*EXP(-XI*RSQAB)
              PCSBOUNDA = GCC(IG)*GCC(JG)*PCOV1 
              NIPRE = ZET(KG)*ZET(LG)
              NIPR = ZET(KG) + ZET(LG)
              ZETNI = ZETPR + NIPR
              XI = ZET(KG)*ZET(LG)/NIPR
              PCOV2 = (PI/NIPR)**(3.0/2.0)*EXP(-XI*RSQCD)
              PCSBOUNDB = GCC(KG)*GCC(LG)*PCOV2
              PCSBOUND = PCSBOUNDA*PCSBOUNDB
              IF (ABS(PCSBOUND) < INTHRESH) GOTO 131
C
C     ___ Calculate the P - A, P - B, Q - C and Q - D factors ___
C
              DPA(:) = -ZET(JG)*DAB(:)/ZETPR 
              DPB(:) = ZET(IG)*DAB(:)/ZETPR 
              DQB(:) = (ZET(KG)*CAT(:) + ZET(LG)*DAT(:))/
     $                 NIPR - BAT(:)
              DQA(:) = (ZET(KG)*CAT(:) + ZET(LG)*DAT(:))/
     $                 NIPR - AAT(:)
              DQC(:) = -ZET(LG)*DCD(:)/NIPR
              DQD(:) = ZET(KG)*DCD(:)/NIPR
C
C     ___ Calculate the P - Q factor ___
C
              DPQ(:) = (ZET(IG)*AAT(:) + ZET(JG)*BAT(:))/ZETPR -
     $                 (ZET(KG)*CAT(:) + ZET(LG)*DAT(:))/NIPR
              DWQ(:) = ZETPR*DPQ(:)/ZETNI
              PQSQ = SUM(DPQ**2,1)
              RHO = ZETPR*NIPR/ZETNI
              DWP(:) = -NIPR/ZETNI*DPQ(:)
C
C     ___ Calculate the incomplete Gamma function ___
C
              T = RHO*PQSQ
              CALL GAMMAF(F,T,ULLR+1)
C
C     ___ Build primitive integrals over s functions ___
C
              RRERI(1,1,1,1,0) = GCC(IG)*GCC(JG)*PCOV1*
     $                           GCC(KG)*GCC(LG)*PCOV2*F(0)*
     $                         2.0*(RHO/PI)**(1.0/2.0)
              THRESHC = 0
              IF (ABS(RRERI(1,1,1,1,0)) <= INTHRESH2) THRESHC = 1
C
              DO N=0,ULLR
                RRERI(1,1,1,1,N+1) = GCC(IG)*GCC(JG)*PCOV1*
     $                               GCC(KG)*GCC(LG)*PCOV2*F(N+1)*
     $                             2.0*(RHO/PI)**(1.0/2.0)
                IF (ABS(RRERI(1,1,1,1,N+1)).LE.INTHRESH2)
     $          THRESHC = THRESHC + 1
C
                CALL VPROD(DPA,DPB,TU1)
                CALL VPROD(DAB,DWP,TU2)
                RRINT(1,1,1,1,N,:) = 4.0*ZETPRE*(TU1(:)*RRERI(1,1,1,1,N)
     $                               - TU2(:)*RRERI(1,1,1,1,N+1))
              END DO
C
              IF (THRESHC == ULLR+2) GOTO 131
C
              IF (ULLR == 0) THEN
                BLOCKS(1,1,1,1,:) = BLOCKS(1,1,1,1,:) + 
     $                              RRINT(1,1,1,1,0,:)
                GOTO 131
              END IF
C
C     ___ Build basic cartesian integrals ___
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ A vertical recurrence relation steps ___
C
              IF (LA == 0) GOTO 510 
C
              DO 60 ILA=1,LA
                CALL RR21(BLOCKS,RRINT,RRERI,ILA,LA,LB,DPA,DPB,DQA,DQB,
     $                    DWP,NIPR,RHO,ZET(IG),ZET(JG),ZETNI,ZETPR,DCPA,
     $                    DCPC,DSHL,ULLR,LBCD)
   60         CONTINUE
C
  510         CONTINUE
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ B vertical recurrence relation steps ___
C
              IF (LB == 0) GOTO 520
C
              DO 80 ILA=0,LA
                DO 70 ILB=1,LB
                  CALL RR22(BLOCKS,RRINT,RRERI,ILA,ILB,LA,LB,DPA,DPB,
     $                      DQA,DQB,DWP,NIPR,RHO,ZET(IG),ZET(JG),
     $                      ZETNI,ZETPR,DCPA,DCPC,DSHL,ULLR,LCD)
   70           CONTINUE
   80         CONTINUE
C
  520         CONTINUE
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ C vertical recurrence relation steps ___
C
              IF (LC == 0) GOTO 530
C
              DO 110 ILA=0,LA  
                DO 100 ILB=0,LB
                  DO 90 ILC=1,LC
                    CALL RR23(BLOCKS,RRINT,RRERI,ILA,ILB,ILC,LA,LB,LC,
     $                        LD,DPA,DPB,DQC,DQD,DWQ,NIPR,RHO,ZET(IG),
     $                        ZET(JG),ZETNI,ZETPR,DCPA,DCPC,DSHL,ULLR)
   90             CONTINUE
  100           CONTINUE
  110         CONTINUE
C
  530         CONTINUE
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ D horizontal recurrence relation steps ___
C
              IF (LD == 0) GOTO 540
C
              DO 123 ILA=0,LA
                DO 122 ILB=0,LB
                  DO 121 ILC=0,LC
                    DO 120 ILD=1,LD
                      CALL RR24(BLOCKS,RRINT,RRERI,ILA,ILB,ILC,ILD,
     $                          LA,LB,LC,LD,DPA,DPB,DQC,DQD,DWQ,NIPR,
     $                          RHO,ZET(IG),ZET(JG),ZETNI,ZETPR,DCPA,
     $                          DCPC,DSHL,ULLR)
  120               CONTINUE
  121             CONTINUE
  122           CONTINUE
  123         CONTINUE
C
  540         CONTINUE
C
C     ___ End of recurrence relation steps ___
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
  131       CONTINUE
  132     CONTINUE
  140   CONTINUE
  150 CONTINUE
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(RRINT,RRERI,F,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'RR2SOCIN: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O_______________________________________________________________O
C
      END SUBROUTINE RR2SOCIN
