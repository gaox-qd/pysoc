      SUBROUTINE ATOMDATA(ATOM,ILINE)
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
      INCLUDE 'molecule.h'
C
      INTEGER :: I,ATOM,ILINE
      CHARACTER*2, DIMENSION(92) :: ELEM
C
      REAL :: SOZEFF
C
      LOGICAL :: FAILED
C
      CHARACTER :: UCC*7
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      DATA ELEM /'H ','He',
     $          'Li','Be','B ','C ','N ','O ','F ','Ne',
     $          'Na','Mg','Al','Si','P ','S ','Cl','Ar',
     $          'K ','Ca',
     $          'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     $                    'Ga','Ge','As','Se','Br','Kr',
     $          'Rb','Sr',
     $          'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
     $                    'In','Sn','Sb','Te','I ','Xe',
     $          'Cs','Ba',
     $          'La',
     $          'Ce','Pr','Nd','Pm','Sm','Eu','Gd',
     $                         'Tb','Dy','Ho','Er','Tm','Yb','Lu',
     $               'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     $                     'Tl','Pb','Bi','Po','At','Rn',
     $          'Fr','Ra',
     $          'Ac',
     $          'Th','Pa','U '/
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Load atomic number ___
C
      FAILED = .TRUE.
      DO I=1,92
        IF (UCC(ELEM(I)(1:2)) == UCC(ELEMENT(ATOM)(1:2))) THEN
          ZATOM(ATOM) = REAL(I)
          FAILED = .FALSE.
        END IF
      END DO
C
      IF (FAILED) THEN 
        WRITE(OUT,"(/,'ERROR READING INPUT LINE ',I5)") ILINE
        WRITE(OUT,"('ATOMDATA: NO ATOMIC NUMBER FOUND')")
        STOP
      END IF
C
C     ___ Calculate Zeff ___
C
      IF (ZTYPE == 'ZEFF') THEN
        ZEFF(ATOM) = SOZEFF(INT(ZATOM(ATOM)))
        IF (ZEFF(ATOM) == 0.0) THEN
          WRITE(OUT,"('ATOMDATA: CANNOT BE FOUND Z_eff')")
          STOP
        END IF
C
C     ___ ZATOM ___
C
      ELSE IF (ZTYPE == 'ATOM') THEN
        ZEFF(ATOM) = ZATOM(ATOM)*SCAL(ATOM)
      ELSE
        STOP 'ATOMDATA: Z_eff IS NOT AVAILABLE'
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE ATOMDATA
