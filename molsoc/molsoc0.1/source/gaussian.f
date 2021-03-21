      SUBROUTINE GAUSSIAN(MATRIX,OCCN,OPTION1,OPTION2,OPTION3,DENS)
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
      CHARACTER :: OPTION1*(*),OPTION2*(*),OPTION3*(*),PRSTR*80,
     $             PRSTR1*29,LINE*80,FSTRING*21
C
      INTEGER :: START,ISTO,JSTO,LOOP,NTIM,ARR,DENS,DIM,ILINE
      INTEGER, DIMENSION(5) :: ORBIT
C
      REAL, DIMENSION(NSTO,NSTO) :: MATRIX
      REAL, DIMENSION(5) :: VECTOR
      REAL, DIMENSION(MXSTO) :: OCCN
C
      INTEGER :: DIFF,COUNT,LIMIT,J,ISTOP,JSTOP,JCOUN
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:,:), ALLOCATABLE :: CARTESIAN
C
      LOGICAL :: ISTAT
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      PRSTR = OPTION1//' MO COEFFICIENTS OF THE '//OPTION2//' STATE'
      PRSTR1 = 'COEFFICIENTS OF '//OPTION3//' MOs'
C
      REWIND(DENS)
C
      IF (OPTION1 == 'ALPHA') THEN
        FSTRING = 'Alpha MO coefficients'
      ELSE IF (OPTION1 == 'BETA') THEN
        FSTRING = 'Beta MO coefficients '
      ELSE
        WRITE(OUT,*) 'RADMOS: UNKNOWN READING MODE'
        STOP
      END IF
C
      IF (OPTION3 == 'SPHERICAL') THEN
        DIM = NSPH
      ELSE IF (OPTION3 == 'CARTESIAN') THEN
        DIM = NSTO
      ELSE
        WRITE(OUT,*) 'GAUSSIAN: UNKNOWN MOs OPTION'
        STOP
      END IF
C
C     ___ Read Gaussian MO coefficients ___
C
      ILINE = 0
  200 CONTINUE
      ILINE = ILINE + 1
      READ(DENS,"(A80)") LINE
      IF (LINE(1:21) == FSTRING) GOTO 300
      GOTO 200
  300 CONTINUE
C
      MATRIX = 0.0
      NTIM = (DIM-1)/5 + 1
      START = 0
      ARR = 0
      COUNT = 0
      LIMIT = DIM*DIM 
      ISTAT = .FALSE.
      VECTOR = 0.0
C
C     ___ Read coefficients matrix ___
C
      ISTOP = 1
      JSTOP = 0
   50 CONTINUE
      COUNT = COUNT + 5
      START = 1
      ARR = 5
      IF (COUNT >= LIMIT) THEN 
        ARR = LIMIT - (COUNT - 5)
        ISTAT = .TRUE.
      END IF
      READ(DENS,*) (VECTOR(J),J=START,ARR)
      DO J=START,ARR
        JSTOP = JSTOP + 1
        IF (JSTOP > DIM) THEN
          ISTOP = ISTOP + 1
          JSTOP = 1
        ELSE IF (ISTOP > DIM) THEN
          WRITE(OUT,*) 'GAUSSIAN: ISTOP > DIM'
          STOP 
        END IF
        MATRIX(JSTOP,ISTOP) = VECTOR(J)
      END DO
      IF (ISTAT) GOTO 51
      GOTO 50
   51 CONTINUE
C
      CALL ORBORD(MATRIX,NSTO,OPTION3)
C
C     ___ Write coefficients matrix ___
C
      WRITE(COESP,*)  
      WRITE(COESP,"(T2,A80)") PRSTR
      WRITE(COESP,"(T2,A29)") PRSTR1
C
      NTIM = (DIM-1)/5 + 1
      START = 0
      ARR = 0
      DO LOOP=1,NTIM
        START = ARR + 1
        ARR = ARR + 5
        ARR = MIN(ARR,DIM)
        JCOUN = 0
        DO JSTO=START,ARR
          JCOUN = JCOUN + 1
          ORBIT(JCOUN) = JSTO
        END DO
        WRITE(COESP,"(6X,5I10)") (ORBIT(JSTO),JSTO=1,JCOUN) 
        WRITE(COESP,"(9X,5F10.5)") (OCCN(JSTO),JSTO=START,ARR)
        WRITE(COESP,*) 
        DO ISTO=1,DIM 
          IF (OPTION3.EQ.'SPHERICAL') THEN
            WRITE(COESP,"(I5,A4,5F10.5)") ISTO,ORSYMS(ISTO),
     $                                    (MATRIX(ISTO,JSTO),
     $                                     JSTO=START,ARR)
          ELSE IF (OPTION3.EQ.'CARTESIAN') THEN
            WRITE(COESP,"(I5,A4,5F10.5)") ISTO,ORSYM(ISTO),
     $                                    (MATRIX(ISTO,JSTO),
     $                                     JSTO=START,ARR)
          END IF
        END DO
      END DO
C
C     ___ Allocate local fields ___
C
      IF (OPTION3 /= 'SPHERICAL') RETURN
C
      ALLOCATE(CARTESIAN(NSTO,NSTO),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'GAUSSIAN: ALLOCATION FAILED'
        STOP
      END IF
C
C     ___ MO transformation sphericals to cartesians ___
C
      CARTESIAN = 0.0
      CALL MOTRSC(MATRIX,CARTESIAN,NSTO)
C
      MATRIX = 0.0
      MATRIX = CARTESIAN
C
C     ___ Write cartesian coefficients ___
C
      PRSTR = OPTION1//' MO COEFFICIENTS OF THE '//OPTION2//' STATE'
      PRSTR1 = 'COEFFICIENTS OF CARTESIAN MOs'
      WRITE(COESP,*)
      WRITE(COESP,"(T2,A80)") PRSTR
      WRITE(COESP,"(T2,A25)") PRSTR1
C
      NTIM = (NSTO-1)/5 + 1
      START = 0
      ARR = 0
      DO LOOP=1,NTIM
        START = ARR + 1
        ARR = ARR + 5
        ARR = MIN(ARR,NSTO)
        JCOUN = 0
        DO JSTO=START,ARR
          JCOUN = JCOUN + 1
          ORBIT(JCOUN) = JSTO
        END DO
        WRITE(COESP,"(6X,5I10)") (ORBIT(JSTO),JSTO=1,JCOUN)
        WRITE(COESP,"(9X,5F10.5)") (OCCN(JSTO),JSTO=START,ARR)
        WRITE(COESP,*)
        DO ISTO=1,NSTO
          WRITE(COESP,"(I5,A4,5F10.5)") ISTO,ORSYM(ISTO),
     $                                  (MATRIX(ISTO,JSTO),
     $                                   JSTO=START,ARR)
        END DO
      END DO
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(CARTESIAN,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'GAUSSIAN: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE GAUSSIAN
