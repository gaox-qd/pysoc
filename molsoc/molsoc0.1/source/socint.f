      SUBROUTINE SOCINT(MATRIX,CENTS)
C
C     ******************************************************************
C     ***        Calculate the one-electron SOC INTegrals            ***
C     ******************************************************************
C
C     Lit.: S.Obara, A.Saika, J.Chem.Phys. 84, 3963 (1986)
C
C     History: - Creation (14.08.03, SC)
C
C     List of local dimensions:
C
C     List of local variables:
C
C     CENTS  : Atomic center.
C     MATRIX : The <Chimu | SOC | Chini > integral matrix.
C     L{A/B} : Orbital L quantum number.
C
C     List of local dynamical fields:
C
C     BLOCKS: Diatomic integral shell block.
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
C
      INCLUDE 'fileio.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      INTEGER :: IATOM,ISHL,JATOM,JSHL,LA,LB,ISTO,JSTO,IVEC
C
      REAL, DIMENSION(3) :: CENTS
      REAL, DIMENSION(NSTO,NSTO,3) :: MATRIX
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: BLOCKS
C
C     o~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~o
C
C     ___ Allocate local fields ___
C
      ALLOCATE(BLOCKS(DSHL,DSHL,3),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'SOCINT: ALLOCATION FAILED'
        STOP
      END IF
C
C     ___ Initialize integral matrix ___
C
      BLOCKS = 0.0
      MATRIX = 0.0
C
C     ___ Calculate the <Chimu | SOC | Chini> integral matrix ___
C
      DO 40 IATOM=1,NATOM
        DO 30 JATOM=IATOM,NATOM
          CALL ATOMDIS(IATOM,JATOM,RAB,RSQAB)
          RAC(:) = C(:,IATOM) - CENTS(:)
          RBC(:) = C(:,JATOM) - CENTS(:)
          DO 20 ISHL=LLS(IATOM),ULS(IATOM)
            DO 10 JSHL=MAX(ISHL,LLS(JATOM)),ULS(JATOM)
              LA = PSHELL(ISHL,2)
              LB = PSHELL(JSHL,2)
              CALL RR1SOCIN(LA,LB,ISHL,JSHL,BLOCKS)
              CALL BLDIMSO(ISHL,JSHL,NCSTO,MATRIX,BLOCKS,NSTO)
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
   40 CONTINUE
C
C     ___ Antisymmetrize SOC matrix ___
C
      DO ISTO=1,NSTO-1
        DO JSTO=ISTO+1,NSTO
          MATRIX(JSTO,ISTO,1:3) = -MATRIX(ISTO,JSTO,1:3)
        END DO
      END DO
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(BLOCKS,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'SOCINT: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE SOCINT
