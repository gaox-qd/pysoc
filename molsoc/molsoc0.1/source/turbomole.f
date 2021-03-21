      SUBROUTINE TURBOMOLE(MATRIX,OCCN,OPTION1,OPTION2,OPTION3,DENS)
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
     $             PRSTR1*29,LINE*80,FSTRING*21,PSTR*5
C
      INTEGER :: START,ISTO,JSTO,LOOP,TIMES,ARR,DENS,DIM,ILINE,DIFF,
     $           COUNT,LIMIT,J,ISTOP,JSTOP,JCOUN
      INTEGER, DIMENSION(5) :: ORBIT
C
      REAL, DIMENSION(5) :: VECTOR
      REAL, DIMENSION(NSTO,NSTO) :: MATRIX
      REAL, DIMENSION(MXSTO) :: OCCN
C
      INTEGER ALLOCATION
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
      IF (UKS) THEN
C
        IF (OPTION1 == 'ALPHA') THEN 
          PSTR = 'alpha'
        ELSE IF (OPTION1 == 'BETA') THEN
          PSTR = 'beta '
        ELSE
          WRITE(OUT,*) 'TURBOMOLE: UNKNOWN READING MODE'
          STOP
        END IF
C
        ILINE = 0
   20   CONTINUE
        ILINE = ILINE + 1
        READ(DENS,"(A80)") LINE
        IF (LINE(8:12) /= PSTR) GOTO 20
C
        DO ILINE=1,2
          READ(DENS,"(A80)") LINE
        END DO
C
      ELSE IF (SROKS) THEN
C
        DO ILINE=1,3
          READ(DENS,"(A80)") LINE
        END DO
C
      END IF
C 
      IF (OPTION3 == 'SPHERICAL') THEN
        DIM = NSPH
      ELSE IF (OPTION3 == 'CARTESIAN') THEN
        WRITE(OUT,*) 'CARTESIAN CONVERSION AOs IS NOT AVAILABLE'
        STOP
      ELSE
        WRITE(OUT,*) 'TURBOMOLE: UNKNOWN MOs OPTION'
        STOP
      END IF
C
      MATRIX = 0.0
      LIMIT = DIM*DIM 
      VECTOR = 0.0
C
C     ___ Read coefficients matrix ___
C
      DO JSTOP=1,DIM
        READ(DENS,"(A80)") LINE
        ISTAT = .FALSE.
        COUNT = 0
        START = 0
        ARR = 0
        ISTOP = 0
   50   CONTINUE
        COUNT = COUNT + 4
        START = 1
        ARR = 4
        IF (COUNT >= DIM) THEN
          ARR = DIM - (COUNT - 4)
          ISTAT = .TRUE.
        END IF
        READ(DENS,"(4D20.14)") (VECTOR(J),J=START,ARR)
        DO J=START,ARR
          ISTOP = ISTOP + 1
          MATRIX(ISTOP,JSTOP) = VECTOR(J)
        END DO
        IF (ISTAT) GOTO 51
        GOTO 50
   51   CONTINUE
      END DO
C
      CALL TURBORD(MATRIX,NSTO,OPTION3)
C
C     ___ Write coefficient's matrix ___
C
      WRITE(COESP,*)  
      WRITE(COESP,"(T2,A80)") PRSTR
      WRITE(COESP,"(T2,A29)") PRSTR1
C
      TIMES = (DIM-1)/5 + 1
      START = 0
      ARR = 0
      DO LOOP=1,TIMES
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
          IF (OPTION3 == 'SPHERICAL') THEN
            WRITE(COESP,"(I5,A4,5F10.5)") ISTO,ORSYMS(ISTO),
     $                                    (MATRIX(ISTO,JSTO),
     $                                     JSTO=START,ARR)
          ELSE IF (OPTION3 == 'CARTESIAN') THEN
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
        WRITE(OUT,*) 'TURBOMOLE: ALLOCATION FAILED'
        STOP
      END IF
C
C     ___ MO transformation from sphericals to cartesians ___
C
      CALL MOTRSC(MATRIX,CARTESIAN,NSTO)
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
      TIMES = (NSTO-1)/5 + 1
      START = 0
      ARR = 0
      DO LOOP=1,TIMES
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
     $                                  JSTO=START,ARR)
        END DO
      END DO
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(CARTESIAN,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'TURBOMOLE: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE TURBOMOLE
