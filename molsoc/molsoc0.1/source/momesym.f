      SUBROUTINE MOMESYM
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'molecule.h'
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      MOMSYM(2) = 'X  ' 
      MOMSYM(3) = 'Y  '
      MOMSYM(4) = 'Z  '
      MOMSYM(5) = 'XX '
      MOMSYM(6) = 'XY '
      MOMSYM(7) = 'XZ '
      MOMSYM(8) = 'YY '
      MOMSYM(9) = 'YZ '
      MOMSYM(10) = 'ZZ '
      MOMSYM(11) = 'XXX'
      MOMSYM(12) = 'XXY'
      MOMSYM(13) = 'XXZ'
      MOMSYM(14) = 'XYY'
      MOMSYM(15) = 'XYZ'
      MOMSYM(16) = 'XZZ'
      MOMSYM(17) = 'YYY'
      MOMSYM(18) = 'YYZ'
      MOMSYM(19) = 'YZZ'
      MOMSYM(20) = 'ZZZ'
C
C     O________________________________________________________________O
C
      END SUBROUTINE MOMESYM 
