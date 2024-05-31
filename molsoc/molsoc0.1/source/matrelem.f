      SUBROUTINE MATRELEM
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
      INTEGER :: ISPIN1,ISPIN2,J1,J2,M1,M2,J,M,I,II
      REAL :: MS1,MS2,S1,S2,RJ1,RJ2,RM1,RM2,RJ,RM,CLEB,
     $        SUMCONTX,SUMCONTY,SUMCONTZ,MSP1,MSP2,
     $        SUMCONTAUX,SUMCONTAUY,SUMCONTAUZ
      REAL, DIMENSION(3) :: NORMA
C
      COMPLEX :: MAT1E,MAT2E,MATT 
      COMPLEX, DIMENSION(3) :: MATRE1E(3),MATRE2E(3),MATRET(3)
C
      REAL :: CLEBSCH
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      S1 = (REAL(MULTF) - 1.0)/2.0
      S2 = (REAL(MULTS) - 1.0)/2.0
C
      NORMA = 0.0
      IF (DISCIN == 0) THEN
        MSP1 = S1
        MSP2 = S2
C
        J1 = INT(2.0*S2)
        J2 = 2
        M1 = INT(2.0*S2)
        M2 = 0
        J = INT(2.0*S1)
        M = INT(2.0*S1)
        NORMA(3) = CLEBSCH(J1,J2,M1,M2,J,M)
C
        MSP1 = S1-1.0
        MSP2 = S2
        J1 = INT(2.0*S2)
        J2 = 2
        M1 = INT(2.0*MSP2)
        M2 = -2
        J = INT(2.0*S1)
        M = INT(2.0*MSP1)
        NORMA(1) = -CLEBSCH(J1,J2,M1,M2,J,M)
        norma(1) = NORMA(1)*sqrt(S1*2.0)
        NORMA(2) = NORMA(1)
C
      ELSE IF (DISCIN == 2) THEN
C
        MSP1 = S1
        MSP2 = S2
        J1 = INT(2.0*S2)
        J2 = 2
        M1 = INT(2.0*S2)
        M2 = -2
        J = INT(2.0*S1)
        M = INT(2.0*S1)
        NORMA(1) = CLEBSCH(J1,J2,M1,M2,J,M)
        NORMA(2) = NORMA(1)
C
        J1 = INT(2.0*S2)
        J2 = 2
        M1 = INT(2.0*S1)
        M2 = 0 
        J = INT(2.0*S1)
        M = INT(2.0*S1)
        NORMA(3) = CLEBSCH(J1,J2,M1,M2,J,M)
        norma(3) = norma(3)*sqrt(S2*2.0)/2.0
      END IF
C
      MATRE1E(1) = CMPLX(ENER(1),ENER(2))
      MATRE1E(2) = CMPLX(-ENER(1),ENER(2))
      MATRE1E(3) = CMPLX(ENER(3),0.0)
C
      MATRE2E(1) = CMPLX(ENERTWO(1),ENERTWO(2))
      MATRE2E(2) = CMPLX(-ENERTWO(1),ENERTWO(2))
      MATRE2E(3) = CMPLX(ENERTWO(3),0.0)
C
      SUMCONTX = ENER(1) + ENERTWO(1)
      SUMCONTY = ENER(2) + ENERTWO(2)
      SUMCONTZ = ENER(3) + ENERTWO(3)
C
      SUMCONTAUX = ENERD(1) + ENERTWOAU(1)
      SUMCONTAUY = ENERD(2) + ENERTWOAU(2)
      SUMCONTAUZ = ENERD(3) + ENERTWOAU(3)
C
      MATRET(1) = CMPLX(SUMCONTX,SUMCONTY)
      MATRET(2) = CMPLX(-SUMCONTX,SUMCONTY)
      MATRET(3) = CMPLX(SUMCONTZ,0.0)
C
      WRITE(OUT,1000)
      WRITE(OUT,2000)
C
      DO II=1,3
        IF (NORMA(II) == 0.0) THEN
          WRITE(OUT,"(T2,'MATRELEM: NORMALIZATION = 0')")
          STOP 
        END IF
      END DO
C
      WRITE(OUT,*)
      CLEB = 0.0
      MS1 = S1
      DO ISPIN1=1,MULTF
        MS2 = S2
        DO ISPIN2=1,MULTS
          I = INT(MS2 - MS1)
          IF (ABS(I) == 0) THEN
            J1 = INT(2.0*S2)
            J2 = 2 
            M1 = INT(2.0*MS2)
            M2 = -2*I
            J = INT(2.0*S1)
            M = INT(2.0*MS1)
            CLEB = CLEBSCH(J1,J2,M1,M2,J,M)/NORMA(3)
            MAT1E = CLEB*MATRE1E(3)
            MAT2E = CLEB*MATRE2E(3)
            MATT = CLEB*MATRET(3)
            WRITE(OUT,3000) MS1,MS2,MAT1E,MAT2E,MATT
          END IF
          MS2 = MS2 - 1.0
        END DO
        MS1 = MS1 - 1.0
      END DO
C
      WRITE(OUT,*)
      CLEB = 0.0
      MS1 = S1
      DO ISPIN1=1,MULTF
        MS2 = S2
        DO ISPIN2=1,MULTS
          I = INT(MS2 - MS1)
          IF (ABS(I) == 1) THEN
C           IF (I == 1) II = 2
C           IF (I == -1) II = 1
            II = (3 + I)/2
            J1 = INT(2.0*S2)
            J2 = 2
            M1 = INT(2.0*MS2)
            M2 = -2*I
            J = INT(2.0*S1)
            M = INT(2.0*MS1)
            CLEB = (-1)**I*CLEBSCH(J1,J2,M1,M2,J,M)/NORMA(II)
            MAT1E = CLEB*MATRE1E(II)
            MAT2E = CLEB*MATRE2E(II)
            MATT = CLEB*MATRET(II)
            WRITE(OUT,3000) MS1,MS2,MAT1E,MAT2E,MATT
          END IF
          MS2 = MS2 - 1.0
        END DO
        MS1 = MS1 - 1.0
      END DO
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Format statements ___
C
 1000 FORMAT(//,T2,30('*'),
     $       /,T2,'*** MATRIX ELEMENTS (CM-1) ***',
     $       /,T2,30('*'))
 2000 FORMAT(/,T3,'(I-STATE)',X,'(II-STATE)',12X,'1-el',20X,'2-el',19X,
     $            'TOTAL',
     $       /,T7,'MS',8X,'MS',12X,'R',11X,'I',11X,'R',11X,'I',
     $                         11X,'R',11X,'I',
     $       /,T2,94('-'))
 3000 FORMAT(T5,F5.2,5X,F5.2,3X,F12.4,2F12.4,2F12.4,2F12.4)
C
C     O________________________________________________________________O
C     
      END SUBROUTINE MATRELEM
