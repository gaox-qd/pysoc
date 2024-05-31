      SUBROUTINE OVERABU(DMAT,NOAF,NOBF,NOAS,NOBS)
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'molecule.h'
      INCLUDE 'fileio.h'
      INCLUDE 'math.h'
C
      REAL, PARAMETER :: OVERTOL = 1.0E-10
C
      INTEGER :: DMAT,NOAF,NOBF,NOAS,NOBS,IMAT,JMAT,NOLIM
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:,:), ALLOCATABLE :: SMOA,SMOB
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Allocate local fields ___
C
      ALLOCATE(SMOA(DMAT,DMAT),SMOB(DMAT,DMAT),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'OVERABU: ALLOCATION FAILED'
        STOP
      END IF
C
C     ___ Read from tape two-state MO overlap integrals ___
C
      CALL IOSCR(SMOA,DMAT**2,'READ','ZERO',5)
      CALL IOSCR(SMOB,DMAT**2,'READ','ZERO',6)
C
C     ___ Write two-state alpha MO overlap integrals ___
C
      WRITE(IOAL,"(/,'TWO-STATE MO OVERLAP INTEGRALS',/)")
C
      DO IMAT=1,NOAF
        DO JMAT=1,NOAS
          WRITE(IOAL,1000) IMAT,JMAT,SMOA(IMAT,JMAT)
          IF (IMAT /= JMAT) THEN
            IF (SMOA(IMAT,JMAT) > OVERTOL) THEN
              WRITE(IOAL,"('OFF DIAGONAL ELEMENT >',F20.9)") OVERTOL
            END IF
          END IF
        END DO
      END DO
C
C     ___ Write alpha diagonal elements ___
C
      WRITE(IOAL,"(/,'DIAGONAL ELEMENTS',/)")
C
      NOLIM = NOAF
      IF (NOAS > NOAF) NOLIM = NOAS
      DO IMAT=1,NOLIM
        WRITE(IOAL,1000) IMAT,IMAT,SMOA(IMAT,IMAT)
      END DO
C
C     ___ Write two-state beta MO overlap integrals ___
C
      WRITE(IOBE,"(/,'TWO-STATE MO OVERLAP INTEGRALS',/)")
C
      DO IMAT=1,NOBF
        DO JMAT=1,NOBS
          WRITE(IOBE,1000) IMAT,JMAT,SMOB(IMAT,JMAT)
          IF (IMAT /= JMAT) THEN
            IF (SMOB(IMAT,JMAT) > OVERTOL) THEN
              WRITE(IOAL,"('OFF DIAGONAL ELEMENT >',F20.9)") OVERTOL
            END IF
          END IF
        END DO
      END DO
C
C     ___ Write beta diagonal elements ___
C
      WRITE(IOBE,"(/,'DIAGONAL ELEMENTS',/)")
C
      NOLIM = NOBF
      IF (NOBS.GT.NOBF) NOLIM = NOBS
      DO IMAT=1,NOLIM
        WRITE(IOBE,1000) IMAT,IMAT,SMOB(IMAT,IMAT)
      END DO
C
      DEALLOCATE(SMOA,SMOB,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'OVERABU: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Format statements ___
C
 1000 FORMAT(I7,4X,I7,10X,F20.9)
C
C     O_______________________________________________________________O
C
      END SUBROUTINE OVERABU
