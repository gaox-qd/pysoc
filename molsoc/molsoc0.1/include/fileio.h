C     FILEIO COMMON
C
C     Purpose: Common definition for FILE Input/Output routines.
C
C     ------------------------------------------------------------------
C
C     *** Input/Output parameters ***
C
      INTEGER MINPER,MAXPER,MINTEM,MAXTEM,INP,BASIS,OUT,DEN1,DEN2,SOCA,
     $        COESP,SOMAT,IOSO,IOAL,IOBE,IOMO,IOPM
C
      PARAMETER (MINPER = 1,
     $           INP    = 1,
     $           BASIS  = 2,
     $           OUT    = 3,
     $           DEN1   = 4,
     $           DEN2   = 5,
     $           SOCA   = 6,
     $           COESP  = 7,
     $           SOMAT  = 8,
     $           IOAL   = 9,
     $           IOBE   = 10,
     $           MAXPER = 10,
     $           MINTEM = 11,
     $           IOSO   = 11,
     $           IOMO   = 12,
     $           IOPM   = 13,
     $           MAXTEM = 13)
   
C
C     *** Common blocks for FILE Input/Output variables ***
C
      CHARACTER*10 FNAME(MAXTEM)
C
      COMMON /CFILEIO/ FNAME
C
C     ------------------------------------------------------------------
