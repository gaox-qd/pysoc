      SUBROUTINE FINSTR(LINE,SUCCESS,STRING) 
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
      CHARACTER*(*) :: LINE
C
      INTEGER :: I,ILINE,ILOW,COUNC
      INTEGER, DIMENSION(80) :: LENSTR,LLSTR
C
      LOGICAL :: SUCCESS
C
      CHARACTER :: UCC*7,STRING*(*)
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      ILINE = 1
      ILOW = 0
      LENSTR = 0
      LLSTR = 0
      SUCCESS = .FALSE.
C
      DO I=1,LEN(LINE)
        IF (INDEX(' ',LINE(I:I)) == 0) THEN
          ILOW = ILOW + 1
          LENSTR(ILINE) = LENSTR(ILINE) + 1
        ELSE
          ILINE = ILINE + 1
          ILOW = 0
        END IF
        IF (ILOW == 1) LLSTR(ILINE) = I
      END DO
C
      COUNC = 0
      DO I=1,ILINE
        IF (LENSTR(I) /= 0) THEN
          IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == STRING(1:3)) 
     $       SUCCESS = .TRUE.
          COUNC = COUNC + 1
        END IF
      END DO
C
      IF ((UCC(STRING(1:3)) == 'END').AND.(COUNC == 0)) SUCCESS = .TRUE.
C
C     O________________________________________________________________O
C
      END SUBROUTINE FINSTR
      
