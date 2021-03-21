      REAL FUNCTION SOZEFF(ATOM)
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'fileio.h'
C
      INTEGER, PARAMETER :: MAXAN = 54
C
      INTEGER ATOM
      INTEGER, DIMENSION(MAXAN) :: NEVAL
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      DATA NEVAL/1,2, 
     $           1,2,3,4,5,6,7,8,
     $           1,2,3,4,5,6,7,8, 
     $           1,2,3,4,5,6,7,8,9,10,11,12,
     $                   3,4,5,6,7,8,
     $           1,2,3,4,5,6,7,8,9,10,11,12,
     $                   3,4,5,6,7,8/
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      IF (ATOM == 1) THEN
        SOZEFF = 1.0
      ELSE IF (ATOM == 2) THEN
        SOZEFF = 2.0
      ELSE IF ((ATOM > 2).AND.(ATOM <= 10)) THEN
        SOZEFF = (0.2517 + 0.0626*NEVAL(ATOM))*ATOM
      ELSE IF ((ATOM > 10).AND.(ATOM <= 18)) THEN
        SOZEFF = (0.7213 + 0.0144*NEVAL(ATOM))*ATOM
      ELSE IF (((ATOM > 18).AND.(ATOM <= 20)).OR.((ATOM > 30).
     $                      AND.(ATOM <= 36))) THEN
        SOZEFF = (0.8791 + 0.0039*NEVAL(ATOM))*ATOM
      ELSE IF (((ATOM > 36).AND.(ATOM <= 38)).OR.((ATOM > 48).
     $                      AND.(ATOM <= 54))) THEN
        SOZEFF = (0.9228 + 0.0017*NEVAL(ATOM))*ATOM
      ELSE IF (ATOM == 26) THEN
        SOZEFF = 0.583289*ATOM
      ELSE IF (ATOM == 77) THEN
        SOZEFF = 1234
      ELSE IF (ATOM == 30) THEN
        SOZEFF = 330
C	 Zeff taken from J. Phys. Chem. A 1998, 102, 50, 10430–10435:
C    First row transition metals.
      ELSE IF (ATOM >= 21 .AND. ATOM <= 30) THEN
        SOZEFF = ATOM * (0.385 + 0.025 * (NEVAL(ATOM) -2))
C    Second row transition metals.
      ELSE IF (ATOM >= 39 .AND. ATOM <= 48) THEN
        SOZEFF = ATOM * (4.680 + 0.060 * (NEVAL(ATOM) -2))
C    Third row transition metals.
      ELSE IF (ATOM == 72) THEN
        SOZEFF = 1025.28
      ELSE IF (ATOM == 73) THEN
        SOZEFF = 1049.74
      ELSE IF (ATOM == 74) THEN
        SOZEFF = 1074.48
      ELSE IF (ATOM == 75) THEN
        SOZEFF = 1099.5
      ELSE IF (ATOM == 76) THEN
        SOZEFF = 1124.8
      ELSE IF (ATOM == 77) THEN
        SOZEFF = 1150.38
      ELSE IF (ATOM == 78) THEN
        SOZEFF = 1176.24
      ELSE IF (ATOM == 79) THEN
        SOZEFF = 1202.38
      ELSE IF (ATOM == 80) THEN
        SOZEFF = 1228.8
      ELSE
        WRITE(OUT,1000) ATOM
        STOP
      END IF
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Format statements ___
C
 1000 FORMAT(/,T2,'SOZEFF: SOZEFF IS NOT AVAILABLE FOR ATOMIC NUMBER',
     $         I5)
C
C     O________________________________________________________________O
C
      END FUNCTION SOZEFF
