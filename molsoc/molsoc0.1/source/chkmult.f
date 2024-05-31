      SUBROUTINE CHKMULT
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
      CHARACTER :: STRING*80
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Check the right multiplicity ___
C
      IF (MULTF < 1) THEN 
        STRING = 'MULTIPLICITY OF THE I-STATE MUST BE > 0'
        GOTO 100
      ELSE IF (MULTS < 1) THEN
        STRING = 'MULTIPLICITY OF THE II-STATE MUST BE > 0'
        GOTO 100
      ELSE IF ((MULTF - 1) > NELEC) THEN
        STRING = 'MULTIPLICITY OF THE I-STATE CANNOT BE OBTAINED '//
     $           'BY THE GIVEN NUMBER OF ELECTRONS'
        GOTO 100
      ELSE IF ((MULTS - 1) > NELEC) THEN
        STRING = 'MULTIPLICITY OF THE II-STATE CANNOT BE OBTAINED '//
     $           'BY THE GIVEN NUMBER OF ELECTRONS'
        GOTO 100
      ELSE IF ((MOD(MULTF,2) == 0).AND.(MOD(NELEC,2) == 0)) THEN
        STRING = 'I-STATE: AN EVEN MULTIPLICITY IMPLIES AN ODD '//
     $           'NUMBER OF ELECTRONS'
        GOTO 100
      ELSE IF ((MOD(MULTS,2) == 0).AND.(MOD(NELEC,2) == 0)) THEN
        STRING = 'II-STATE: AN EVEN MULTIPLICITY IMPLIES AN ODD '//
     $           'NUMBER OF ELECTRONS'
        GOTO 100
      ELSE IF ((MOD(MULTF,2) /= 0).AND.(MOD(NELEC,2) /= 0)) THEN  
        STRING = 'I-STATE: AN ODD MULTIPLICITY IMPLIES AN EVEN '//
     $           'NUMBER OF ELECTRONS'
        GOTO 100
      ELSE IF ((MOD(MULTS,2) /= 0).AND.(MOD(NELEC,2) /= 0)) THEN
        STRING = 'II-STATE: AN ODD MULTIPLICITY IMPLIES AN EVEN '//
     $           'NUMBER OF ELECTRONS'
        GOTO 100
      ELSE IF (NELEC > 2*NSTO) THEN
        STRING = 'NOT ENOUGH ORBITAL FUNCTIONS'
        GOTO 100
      END IF
C
      RETURN
C
  100 CONTINUE
C
      WRITE(OUT,"(/,A80)") STRING
      STOP
C
C     O________________________________________________________________O
C
      END SUBROUTINE CHKMULT
