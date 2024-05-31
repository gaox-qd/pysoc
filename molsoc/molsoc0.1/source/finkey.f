      SUBROUTINE FINKEY(LINE,SUCCESS) 
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
      INTEGER :: I,ILINE,ILOW
      INTEGER, DIMENSION(80) :: LENSTR,LLSTR
C
      LOGICAL :: SUCCESS
C
      CHARACTER :: UCC*7
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Find keywords decoding the line ___
C
      ILINE = 1
      ILOW = 0
      LENSTR = 0
      LLSTR = 0
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
      DO I=1,ILINE
        IF (LENSTR(I) /= 0) THEN
          IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'ANG') THEN
            TRANS = BOHR
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'BOH') THEN
            TRANS = 1.0
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'ORT') THEN
            GRAMS = .TRUE.
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'UKS') THEN
            UKS = .TRUE.
            SROKS = .FALSE.
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'RKS') THEN
            SROKS = .TRUE.
            UKS = .FALSE.
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'SPH') THEN
            SPHEORB = .TRUE.
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'CAR') THEN
            SPHEORB = .FALSE.
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'ONE') THEN
            ZTYPE = 'ATOM'
            SOCTWOM = .FALSE.
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'ZEF') THEN
            ZTYPE = 'ZEFF'
            SOCTWOM = .FALSE.
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'TWO') THEN
            ZTYPE = 'ATOM'
            SOCTWOM = .TRUE.
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'DIP') THEN
            LM = 1
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'QUA') THEN
            LM = 2
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'OCT') THEN
            LM = 3
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'GAU') THEN
            READMO = 1
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'DEM') THEN
            READMO = 2
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'TUR') THEN
            READMO = 3
            SUCCESS = .TRUE.
C         ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'DAL') THEN
C           READMO = 4
C           SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'ALT') THEN
            ALTER = .TRUE.
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'GOR') THEN
            MELEM = 1
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'MAL') THEN
            MELEM = 2
            SUCCESS = .TRUE.
C         ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'GAM') THEN
C           READMO = 5
C           SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'NOB') THEN
            NOBIO = .TRUE.
            SUCCESS = .TRUE.
          ELSE IF (UCC(LINE(LLSTR(I):LLSTR(I)+2)) == 'TDB') THEN
            TDTB = .TRUE.
            SUCCESS = .TRUE.
          END IF
        END IF
      END DO
C
C     O________________________________________________________________O
C
      END SUBROUTINE FINKEY
      
