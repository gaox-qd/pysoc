      SUBROUTINE INITIALIZE
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'math.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      INTEGER :: I
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      GCC = 0.0
      C = 0.0 
      MULTF = 0 
      MULTS = 0
      NAT = 0
      NGTF = 0
      NSHL = 0
      ZATOM = 0.0
      NATOM = 0
      TRMSC = 0.0
C
      PSHELL = 0
      NSTO = 0
      COP = 0
      SOP = 0
C
      MS = 1.0
C
      BOHR = 1.88972595820018
      FINESTRC = 7.297352568E-3
      CMM1 = 219474.6
      DEBYE = 2.54174780119995
C     CMM1 = 219474.625390476
      ZTYPE = 'ATOMS'
      TOLTIME = 1.0E-14
C
      ENER = 0.0
      ENERD = 0.0
      ENERTWOAU = 0.0
      ENERTWO = 0.0
C
C     ___ Initialize input keywords ___
C
      TRANS = 1.0
      GRAMS = .FALSE.
      UKS = .FALSE.
      SROKS = .TRUE.
      SPHEORB = .TRUE.
      ZTYPE = 'ATOM'
      SOCTWOM = .TRUE.
      LM = 0
      READMO = 1
      ALTER = .FALSE.
      MELEM = 1
      NOBIO = .FALSE.
      TDTB = .FALSE.
C
C     O________________________________________________________________O
C
      END SUBROUTINE INITIALIZE
