C     PARAMETER DEFINITION
C
C     *** MolSOC parameters ***
C
      INTEGER MXATM,MXSH,MXGTF,MXSTO
C
      PARAMETER (MXATM = 5000,
     $           MXSH  = 10000,
     $           MXGTF = 50000,
     $           MXSTO = 50000)
C
      INTEGER MXL,MXLT,MXSO,MXCO
C
      PARAMETER (MXL = 3,
     $           MXLT = 4*MXL+1,
     $           MXSO  = 2*MXL+1,
     $           MXCO  = (MXL+1)*(MXL+2)/2)
C
      INTEGER MAXFAC,MAXMOM,MAXITAB
C
      PARAMETER (MAXFAC = 30,
     $           MAXMOM = 20,
     $           MAXITAB = 120)
C
