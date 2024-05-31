      SUBROUTINE WRTHEAD
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
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      WRITE(OUT,1000) 
C
C     ___ Write input header ___
C
      IF (TRANS.EQ.BOHR) THEN
        WRITE(OUT,2000) 
      ELSE IF (TRANS.EQ.1.0) THEN
        WRITE(OUT,3000) 
      END IF
C
      IF (GRAMS) THEN
        WRITE(OUT,4000)
      END IF
C
      IF (SROKS) THEN
        WRITE(OUT,5000) 
      ELSE
        WRITE(OUT,6000)
      END IF
C
      IF (SPHEORB) THEN
        WRITE(OUT,7000) 'SPHERICAL'
      ELSE
        WRITE(OUT,7000) 'CARTESIAN'
      END IF
C
      IF ((ZTYPE == 'ATOM').AND.(.NOT.SOCTWOM)) THEN
        WRITE(OUT,8000) 
      ELSE IF ((ZTYPE == 'ZEFF').AND.(.NOT.SOCTWOM)) THEN
        WRITE(OUT,9000)
      ELSE IF ((ZTYPE == 'ATOM').AND.(SOCTWOM)) THEN
        WRITE(OUT,10000)
      END IF
C
      IF (LM == 1) THEN
        WRITE(OUT,11000) 
      ELSE IF (LM == 2) THEN
        WRITE(OUT,12000)
      ELSE IF (LM == 3) THEN
        WRITE(OUT,13000)
      END IF
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Format statements ___
C
 1000 FORMAT(/,T2,'*** SPIN-ORBIT CALCULATIONS ***')
 2000 FORMAT(/,T2,'GEOMETRY IN ANGSTROM')
 3000 FORMAT(/,T2,'GEOMETRY IN BOHR')
 4000 FORMAT(T2,'GRAM-SCHMIDT ORTHOGONALIZATION')
 5000 FORMAT(T2,'RESTRICTED MOs')
 6000 FORMAT(T2,'UNRESTRICTED MOs')
 7000 FORMAT(T2,'COEFFICIENTS OF ',A9,' MOs') 
 8000 FORMAT(T2,'ONE-ELECTRON SOC CALCULATIONS')
 9000 FORMAT(T2,'ONE-ELECTRON SOC CALCULATIONS WITH Z_eff')
10000 FORMAT(T2,'ONE AND TWO-ELECTRON SOC CALCULATIONS')
11000 FORMAT(T2,'DIPOLE MOMENT COMPONENTS')
12000 FORMAT(T2,'DIPOLE AND QUADRUPOLE MOMENT COMPONENTS')
13000 FORMAT(T2,'DIPOLE, QUADRUPOLE, AND OCTUPOLE MOMENT COMPONENTS') 
C
C     O________________________________________________________________O
C
      END SUBROUTINE WRTHEAD
      
