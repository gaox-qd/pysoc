      SUBROUTINE KSSOC
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      CHARACTER :: STRMO*9
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:,:), ALLOCATABLE :: COEFA1,COEFB1,COEFA2,COEFB2
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Allocate local fields ___
C
      ALLOCATE(COEFA1(NSTO,NSTO),COEFB1(NSTO,NSTO),
     $         COEFA2(NSTO,NSTO),COEFB2(NSTO,NSTO),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'KSSOC: ALLOCATION FAILED'
        STOP
      END IF
C
      IF (SPHEORB) THEN
        STRMO = 'SPHERICAL'
      ELSE 
        STRMO = 'CARTESIAN'
      END IF
C
C     ___ Load coefficient's matrix for I and II states ___
C
      IF (UKS) THEN
C
        IF (DIRECTA) THEN
          CALL READMOS(COEFA1,OCCNAF,'ALPHA','FIRST',STRMO,DEN1)
          CALL READMOS(COEFB1,OCCNBF,'BETA','FIRST',STRMO,DEN1)
        ELSE
          CALL READMOS(COEFB1,OCCNBF,'ALPHA','FIRST',STRMO,DEN1)
          CALL READMOS(COEFA1,OCCNAF,'BETA','FIRST',STRMO,DEN1)
        END IF
        IF (DIRECTB) THEN
          CALL READMOS(COEFA2,OCCNAS,'ALPHA','SECOND',STRMO,DEN2)
          CALL READMOS(COEFB2,OCCNBS,'BETA','SECOND',STRMO,DEN2)
        ELSE
          CALL READMOS(COEFB2,OCCNBS,'ALPHA','SECOND',STRMO,DEN2)
          CALL READMOS(COEFA2,OCCNAS,'BETA','SECOND',STRMO,DEN2)
        END IF
C
      ELSE IF (SROKS) THEN
C
        CALL READMOS(COEFA1,OCCNF,'ALPHA','FIRST',STRMO,DEN1)
        CALL READMOS(COEFA2,OCCNS,'ALPHA','SECOND',STRMO,DEN2)
        COEFB1 = COEFA1
        COEFB2 = COEFA2
C
      END IF
C
C     ___ Alter the MOs occupation ___
C
      CALL ALTMOS(COEFA1,COEFB1,'FIRST',NSTO)
      CALL ALTMOS(COEFA2,COEFB2,'SECOND',NSTO)
C
C     ___ Calculate density matrix before biorthogonalization ___
C
      CALL MULDENS(COEFA1,NSTO,HFOMOAF,OCCNAF,'FIRST')
      CALL MULDENS(COEFB1,NSTO,HFOMOBF,OCCNBF,'NONE')
      CALL MULDENS(COEFA2,NSTO,HFOMOAS,OCCNAS,'NONE')
      CALL MULDENS(COEFB2,NSTO,HFOMOBS,OCCNBS,'NONE')
C
C     ___ Calculate multipole moments ___
C
      CALL MULTIPOLE(COEFA1,COEFB1,COEFA2,COEFB2,HFOMOAF,
     $               HFOMOBF,HFOMOAS,HFOMOBS,NSTO)
C
C     ___ Calculate I-II state MO overlaps (biorthogonalization) ___
C
      CALL OVERFS(COEFA1,COEFA2,HFOMOAF,HFOMOAS,NSTO,5)
      CALL OVERFS(COEFB1,COEFB2,HFOMOBF,HFOMOBS,NSTO,6)
C
C     ___ Disable orthogonalization ___
C
      GRAMS = .FALSE. 
C
C     ___ Calculate density matrix after biorthogonalization ___
C
      IF (NOBIO) GOTO 100
C
      CALL MULDENS(COEFA1,NSTO,HFOMOAF,OCCNAF,'SECOND')
      CALL MULDENS(COEFB1,NSTO,HFOMOBF,OCCNBF,'NONE')
      CALL MULDENS(COEFA2,NSTO,HFOMOAS,OCCNAS,'NONE')
      CALL MULDENS(COEFB2,NSTO,HFOMOBS,OCCNBS,'NONE')
C
C     ___ Calculate multipole moments ___
C
      CALL MULTIPOLE(COEFA1,COEFB1,COEFA2,COEFB2,HFOMOAF,
     $               HFOMOBF,HFOMOAS,HFOMOBS,NSTO)
C
  100 CONTINUE
C
C     ___ Calculate MO I-state II-state SOC 1e-contributions ___
C
      IF (DISCIN == 2) THEN
        WRITE(OUT,2000) 1,0,'BETA ALPHA '
        IF (DIRECTA.AND.DIRECTB) THEN
          CALL DISCONE(COEFA2(:,HFOMOAS),COEFB1(:,HFOMOBF),
     $                 HFOMOAF,HFOMOBS,NSTO)
          CALL BLDPFS2(COEFA1,COEFB1,COEFA2,COEFB2,
     $                 HFOMOAF,HFOMOBF,HFOMOAS,HFOMOBS,NSTO)
        ELSE IF ((.NOT.DIRECTA).AND.(.NOT.DIRECTB)) THEN
          CALL DISCONE(COEFB2(:,HFOMOBS),COEFA1(:,HFOMOAF),
     $                 HFOMOBF,HFOMOAS,NSTO)
          CALL BLDPFS2(COEFB1,COEFA1,COEFB2,COEFA2,
     $                 HFOMOBF,HFOMOAF,HFOMOBS,HFOMOAS,NSTO)
        END IF
      ELSE IF (DISCIN == 0) THEN
        IF ((.NOT.DIRECTA).AND.(DIRECTB)) THEN
          WRITE(OUT,2000) 1,1,'BETA ALPHA '
          CALL DISCONE(COEFA2(:,HFOMOAS),COEFB1(:,HFOMOBF),
     $                 HFOMOAF,HFOMOBS,NSTO)
          CALL BLDPFS2(COEFA1,COEFB1,COEFA2,COEFB2,
     $                 HFOMOAF,HFOMOBF,HFOMOAS,HFOMOBS,NSTO)
        ELSE IF (DIRECTA.AND.(.NOT.DIRECTB)) THEN
          WRITE(OUT,2000) 1,1,'ALPHA BETA'
          CALL DISCONE(COEFB2(:,HFOMOBS),COEFA1(:,HFOMOAF),
     $                 HFOMOBF,HFOMOAS,NSTO)
          CALL BLDPFS2(COEFB1,COEFA1,COEFB2,COEFA2,
     $                 HFOMOBF,HFOMOAF,HFOMOBS,HFOMOAS,NSTO)
        ELSE
          IF ((.NOT.DIRECTA).AND.(.NOT.DIRECTB)) THEN
            WRITE(OUT,2000) 0,1,'BETA BETA  '
          ELSE
            WRITE(OUT,2000) 0,1,'ALPHA ALPHA'
          END IF
          CALL DISCZERO(COEFA1(:,HFOMOAF),COEFA2(:,HFOMOAS),
     $                  HFOMOAF,HFOMOBF,NSTO)
          CALL BLDPFS0(COEFA1,COEFB1,COEFA2,COEFB2,
     $                 HFOMOAF,HFOMOAS,HFOMOBS,NSTO)
        END IF
      END IF
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(COEFA1,COEFB1,COEFA2,COEFB2,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'KSSOC: DEALLOCATION FAILED'
        STOP
      END IF
C
C     ___ Calculate 2e-electron contributions ___
C
      IF (DISCIN == 2) THEN
        CALL SOCFST2
      ELSE IF (DISCIN == 0) THEN
        IF ((.NOT.DIRECTA).AND.(DIRECTB)) THEN
          CALL SOCFST2
        ELSE IF (DIRECTA.AND.(.NOT.DIRECTB)) THEN
          CALL SOCFST2
        ELSE
          CALL SOCFST0
        END IF
      END IF
C
      WRITE(OUT,1000) HFOMOAF,HFOMOBF,
     $                HFOMOAS,HFOMOBS
C
      CALL XYZCOMP
C
C     ___ Write I-II state MO overlaps ___
C
      CALL OVERABU(NSTO,HFOMOAF,HFOMOBF,HFOMOAS,HFOMOBS)
C
      CALL CPUTIME
      WRITE(OUT,3000) CPUT
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     *** Format statements ***
C
 1000 FORMAT(/,T2,'FIRST STATE:',
     $      //,T3,'ALPHA HOMO NO. =', I5,
     $       /,T3,'BETA HOMO NO.  =', I5,
     $      //,T2,'SECOND STATE:',
     $      //,T3,'ALPHA HOMO NO. =', I5,
     $       /,T3,'BETA HOMO NO.  =', I5)
C
 2000 FORMAT(/,T2,'SPIN DISCOINCIDENCE    =', I3,
     $       /,T2,'ORBITAL DISCOINCIDENCE =', I3,
     $       /,T2,'COUPLING DETERMINANTS  =',2X,A11)
 3000 FORMAT(/,T2,'THE TOTAL CPU TIME IS:', F12.6,X,'SEC.') 
C
C     O________________________________________________________________O
C
      END SUBROUTINE KSSOC
