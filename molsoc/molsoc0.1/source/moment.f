      SUBROUTINE MOMENT
C
C     ******************************************************************
C     ***             Calculate the MOMENT integrals                 ***
C     ******************************************************************
C
C     Lit.: S.Obara, A.Saika, J.Chem.Phys. 84, 3963 (1986)
C
C     Creation (03.04.07, SC)
C
C     List of local variables:
C
C     CENTS  : Atomic center.
C     MATRIX : The moment integral matrix.
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
      INTEGER :: DMOM, counter
      INTEGER :: IATOM,ISHL,JATOM,JSHL,MATOM,LA,LB,U,ISTO,JSTO
C
      INTEGER :: ALLOCATION
      REAL,DIMENSION(:,:,:), ALLOCATABLE :: BLOCKS,MATRIX
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      IF (LM == 0) RETURN
C
C     ___ Allocate local fields ___
C
      DMOM = (LM**3 - LM)/6 + LM**2 + 2*LM + 1
C
      ALLOCATE(BLOCKS(DSHL,DSHL,2:DMOM),MATRIX(NSTO,NSTO,2:DMOM),
     $         STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'MOMENT: ALLOCATION FAILED'
        STOP
      END IF
C
C     ___ Calculate the moment integral matrix ___
C
      MATRIX = 0.0
      BLOCKS = 0.0
      counter = 1
C
      DO IATOM=1,NATOM
        DO JATOM=IATOM,NATOM
C          print *, "IATOM,JATOM", IATOM,JATOM
          CALL ATOMDIS(IATOM,JATOM,RAB,RSQAB)
          RAC(:) = C(:,IATOM)
          RBC(:) = C(:,JATOM)
          DO ISHL=LLS(IATOM),ULS(IATOM)
            DO JSHL=MAX(ISHL,LLS(JATOM)),ULS(JATOM)
              LA = PSHELL(ISHL,2)
              LB = PSHELL(JSHL,2)
C              print *, "counter", counter
C              print *, "LA,LB", LA,LB
              CALL RRMOMINT(LA,LB,LM,ISHL,JSHL,BLOCKS,DMOM)
              CALL BLDIMOM(ISHL,JSHL,NCSTO,MATRIX,BLOCKS,NSTO,DMOM)
              counter = counter + 1
            END DO  
          END DO  
        END DO  
      END DO  
C
C     ___ symmetrize moment matrix ___
C
      DO ISTO=1,NSTO-1
        DO JSTO=ISTO+1,NSTO
          MATRIX(JSTO,ISTO,2:DMOM) = MATRIX(ISTO,JSTO,2:DMOM)
        END DO
      END DO
C
C     ___ Write on tape moment integrals ___
C
      OPEN(UNIT=171717, FILE='molsoc_dipole.dat', STATUS='REPLACE') 
      DO U=2,DMOM
        CALL IOSCR(MATRIX(:,:,U),NSTO**2,'WRITE','MOME',U)
        WRITE(171717, "(A4,I1)") "DIM=", U-1
        WRITE(171717, "(f19.9)") MATRIX(:,:,U)
      END DO
      CLOSE(171717)
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(BLOCKS,MATRIX,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'MOMENT: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE MOMENT
