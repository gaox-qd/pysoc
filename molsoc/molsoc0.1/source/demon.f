      SUBROUTINE DEMON(MATRIX,OCCN,OPTION1,OPTION2,DENS)
C
C     *******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
C
      INCLUDE 'fileio.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      CHARACTER*(*) :: OPTION1,OPTION2
      CHARACTER :: PRSTR*80,PRSTR1*29,FSTRING*21
      INTEGER :: START,I,J,LOOP,NTIM,ARR,DENS,DIM,JCOUN
      INTEGER, DIMENSION(5) :: ORBIT
C
      REAL, DIMENSION(NSTO,NSTO) :: MATRIX
      REAL, DIMENSION(MXSTO) :: OCCN
      REAL, DIMENSION(5) :: VECTOR
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:), ALLOCATABLE :: EIGVAL
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      PRSTR = OPTION1//' MO COEFFICIENTS OF THE '//OPTION2//' STATE'
      PRSTR1 = 'COEFFICIENTS OF CARTESIAN MOs'
C
C     ___ Allocate local fields ___
C
      ALLOCATE(EIGVAL(NSTO),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'ALLOCATION FAILED'
        STOP
      END IF
C
C     ___ Read deMon MO coefficients ___
C
      DIM = NSTO
      READ(DENS,END=10010,ERR=10000) ((MATRIX(J,I),J=1,DIM),I=1,DIM)
      READ(DENS,END=10010,ERR=10000) (EIGVAL(I),I=1,DIM)
C
      WRITE(COESP,*)  
      WRITE(COESP,*) PRSTR
      WRITE(COESP,*) PRSTR1
C
      NTIM = (DIM-1)/5 + 1
      START = 0
      ARR = 0
      DO LOOP=1,NTIM
        START = ARR + 1
        ARR = ARR + 5
        ARR = MIN(ARR,DIM)
        JCOUN = 0
        DO J=START,ARR
          JCOUN = JCOUN + 1
          ORBIT(JCOUN) = J
        END DO
        WRITE(COESP,"(4X,5I9)") (ORBIT(J),J=1,JCOUN) 
        WRITE(COESP,"(4X,5F10.5)") (OCCN(J),J=START,ARR)
        WRITE(COESP,"(4X,5F10.5)") (EIGVAL(J),J=START,ARR)
        WRITE(COESP,*) 
        DO I=1,DIM 
          WRITE(COESP,"(A4,5F10.5)") ORSYM(I),(MATRIX(I,J),J=START,ARR)
        END DO
      END DO
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(EIGVAL,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'DEMON: DEALLOCATION FAILED'
        STOP
      END IF
C
      RETURN
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Error handling ___
C
10000 CONTINUE
      WRITE(OUT,"(/,T2,'RSTIO: ERROR READING FROM RESTART FILE')")
      STOP
C
10010 CONTINUE
      WRITE(OUT,"(/,T2,'RSTIO: UNEXPECTED END OF RESTART FILE')")
      STOP
C
C     O________________________________________________________________O
C
      END SUBROUTINE DEMON
