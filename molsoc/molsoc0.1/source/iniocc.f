      SUBROUTINE INIOCC
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      INTEGER :: IMO,NOMOUA,NOMOLA,NOMOUB,NOMOLB
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Initialize orbital occupation ___
C
      OCCNAF = 0.0
      OCCNBF = 0.0
      OCCNAS = 0.0
      OCCNBS = 0.0
      OCCNF = 0.0
C
      OCCNUM = 1.0
      IF (DIRECTA) THEN
        NOMOAF = (NELEC + MULTF - 1)/2
        NOMOBF = (NELEC - MULTF + 1)/2
        NOMOUA = NOMOAF
        NOMOLA = NOMOBF
      ELSE
        NOMOAF = (NELEC - MULTF + 1)/2
        NOMOBF = (NELEC + MULTF - 1)/2
        NOMOLA = NOMOAF
        NOMOUA = NOMOBF
      END IF
      IF (DIRECTB) THEN
        NOMOAS = (NELEC + MULTS - 1)/2
        NOMOBS = (NELEC - MULTS + 1)/2
        NOMOUB = NOMOAS
        NOMOLB = NOMOBS
      ELSE
        NOMOAS = (NELEC - MULTS + 1)/2
        NOMOBS = (NELEC + MULTS - 1)/2
        NOMOLB = NOMOAS
        NOMOUB = NOMOBS
      END IF
C
      CALL MOOCC(OCCNAF,NOMOAF,LFOMOAF,HFOMOAF)
      CALL MOOCC(OCCNBF,NOMOBF,LFOMOBF,HFOMOBF)
      CALL MOOCC(OCCNAS,NOMOAS,LFOMOAS,HFOMOAS)
      CALL MOOCC(OCCNBS,NOMOBS,LFOMOBS,HFOMOBS)
C
      DO IMO=1,NOMOLA
        OCCNF(IMO) = 2.0
      END DO  
C
      IF (NOMOUA.GT.NOMOLA) THEN
        DO IMO=NOMOLA+1,NOMOUA
          OCCNF(IMO) = 1.0
        END DO
      END IF
C
      DO IMO=1,NOMOLB
        OCCNS(IMO) = 2.0 
      END DO 
C
      IF (NOMOUB.GT.NOMOLB) THEN
        DO IMO=NOMOLB+1,NOMOUB
          OCCNS(IMO) = 1.0
        END DO
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE INIOCC
