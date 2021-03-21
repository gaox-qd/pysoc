      SUBROUTINE OVERLAP(S,DSTO)
C
C     ******************************************************************
C     ***               Calculate the OVERLAP matrix                 ***
C     ******************************************************************
C
C     Lit.: S.Obara, A.Saika, J.Chem.Phys. 84, 3963 (1986)
C
C     Creation (16.01.05, SC)
C
C     List of local dimensions:
C
C     DSPP: Dimension of scaling parameters.
C
C     List of local variables:
C
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
      INTEGER :: DSTO
      INTEGER :: IATOM,ISHL,JATOM,JSHL,LA,LB,ISTO,JSTO
C
      REAL, DIMENSION(DSTO,DSTO) :: S
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:,:), ALLOCATABLE :: BLOCKS
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Allocate local fields ___
C
      ALLOCATE(BLOCKS(DSHL,DSHL),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'OVERLAP: ALLOCATION FAILED'
        STOP
      END IF
C
C     ___ Initialize overlap matrix ___
C
      S = 0.0
C
C     ___ Calculate the overlap matrix ___
C
      DO IATOM=1,NATOM
        DO JATOM=IATOM,NATOM
          CALL ATOMDIS(IATOM,JATOM,RAB,RSQAB)
          DO ISHL=LLS(IATOM),ULS(IATOM)
            DO JSHL=MAX(ISHL,LLS(JATOM)),ULS(JATOM)
              LA = PSHELL(ISHL,2)
              LB = PSHELL(JSHL,2)
C              print *, "LA,LB", LA, LB
              CALL RROVINT(LA,LB,ISHL,JSHL,BLOCKS)
              CALL BLDIMAT(ISHL,JSHL,NCSTO,S,BLOCKS,DSTO)
            END DO  
          END DO  
        END DO  
      END DO  
C
C     ___ Symmetrize overlap matrix ___
C
      DO ISTO=1,NSTO-1
        DO JSTO=ISTO+1,NSTO
          S(JSTO,ISTO) = S(ISTO,JSTO)
        END DO
      END DO
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(BLOCKS,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'OVERLAP: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE OVERLAP
