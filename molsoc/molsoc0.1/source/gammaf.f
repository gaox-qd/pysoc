      SUBROUTINE GAMMAF(F,T,NMAX)
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'math.h'
C
      REAL, PARAMETER :: TEPS = 1.E-13
C
      INTEGER :: NMAX
C
      INTEGER :: I,ITAB,K
C
      REAL :: TDELTA,EXPO,G,TMP,TTAB,T
      REAL, DIMENSION(0:NMAX) :: F
C
C     ------------------------------------------------------------------
C
C     ___ Calculate F(t) ___
C
      IF (T < TEPS) THEN
C
C     ___ Special cases: t = 0 ___
C
        DO I=0,NMAX
          F(I) = 1.0/(2.0*I + 1.0)
        END DO
C
      ELSE IF (T <= 12.0) THEN
C
C     ___ 0 < T < 12 --> Taylor expansion ___
C
C     ___ Use the pretabulation of the F(t) function ___
C     ___ for the Taylor series expansion            ___
C
        TDELTA = 0.1
        ITAB = NINT(T/TDELTA)
        TTAB = REAL(ITAB)*TDELTA
C
        F(NMAX) = FTABLE(ITAB,NMAX)
C
        TMP = 1.0
        DO K=1,6
          TMP = TMP*(TTAB - T)
          F(NMAX) = F(NMAX) + FTABLE(ITAB,NMAX+K)*TMP*FAC(K)
        END DO
C
        EXPO = EXP(-T)
C
C     ___ With the downward recursive relation ___
C     ___ are create the other F(t) values     ___
C
        DO I=NMAX-1,0,-1
          F(I) = (2.0*T*F(I+1) + EXPO)/(2.0*I + 1.0)
        END DO
C
      ELSE
C
C     ___ T > 12 ___
C
        IF (T <= 15.0) THEN
C
C     ___ 12 < T <= 15 --> Four term polynom expansion ___
C
          G = 0.4999489092 - 0.2473631686/T +
     $        0.321180909/T**2 - 0.3811559346/T**3
          F(0) = 0.5*SQRT(PI/T) - G*EXP(-T)/T
C
        ELSE IF (T <= 18.0) THEN
C
C     ___ 15 < T <= 18 --> three term polynom expansion ___
C
          G = 0.4998436875 - 0.24249438/T + 0.24642845/T**2
          F(0) = 0.5*SQRT(PI/T) - G*EXP(-T)/T
C
        ELSE IF (T <= 24.0) THEN
C
C     ___ 18 < t <= 24 --> Two term polynom expansion ___
C
          G = 0.499093162 - 0.2152832/T
          F(0) = 0.5*SQRT(PI/T) - G*EXP(-T)/T
C
        ELSE IF (T <= 30.0) THEN
C
C     ___ 24 < t <= 30 --> 1 term polynom expansion ___
C
          G = 0.49
          F(0) = 0.5*SQRT(PI/T) - G*EXP(-T)/T
C
        ELSE
C
C     ___ t > 30 -> Asymptotic formula ___
C
          F(0) = 0.5*SQRT(PI/T)
C
        END IF
C
        IF (T.GT.(2.0*NMAX + 36.0)) THEN
          EXPO = 0.0
        ELSE
          EXPO = EXP(-T)
        END IF
C
C     ___ With the upward recursion relation ___
C     ___ are created the other F(t) values   ___
C
        DO I=1,NMAX
          F(I) = 0.5*((2.0*I - 1)*F(I-1) - EXPO)/T
        END DO
C
      END IF
C
C     ------------------------------------------------------------------
C
      END SUBROUTINE GAMMAF
