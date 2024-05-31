C     MATH COMMON
C
C     ------------------------------------------------------------------
C
C     *** Common block for mathematical constants and functions ***
C
      REAL PI,UM(3,3)
C
      REAL DFAC(-1:2*MAXFAC+1),FAC(0:MAXFAC),FTABLE(MAXITAB,0:MXLT+6),
     $     TRMSC(MXSO,MXCO,0:MXL)
C
      COMMON /RMATH/ DFAC,FAC,FTABLE,PI,UM,TRMSC
C
C     ------------------------------------------------------------------
