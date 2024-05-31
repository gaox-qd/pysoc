      SUBROUTINE GETDATE
C
C     *****************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'fileio.h'
C
      CHARACTER :: DATE*8,TIME*10,ZONE*5
      CHARACTER*8, DIMENSION(12) :: MONTH
      INTEGER :: I
      INTEGER, DIMENSION(8) :: TIMEDATE
C
C     ------------------------------------------------------------------
C
      MONTH(1)  = 'Jan'
      MONTH(2)  = 'Feb'
      MONTH(3)  = 'Mar'
      MONTH(4)  = 'Apr'
      MONTH(5)  = 'May'
      MONTH(6)  = 'Jun'
      MONTH(7)  = 'Jul'
      MONTH(8)  = 'Aug'
      MONTH(9)  = 'Sep'
      MONTH(10) = 'Oct'
      MONTH(11) = 'Nov'
      MONTH(12) = 'Dec'
C
      CALL DATE_AND_TIME(DATE,TIME,ZONE,TIMEDATE)
C
      WRITE(OUT,1000) TIMEDATE(3),MONTH(TIMEDATE(2)),TIMEDATE(1),
     $                TIME(1:2),TIME(3:4),TIME(5:6)
C
C     ------------------------------------------------------------------
C
C     *** Format statements ***
C
 1000 FORMAT(/,T2,3("*"),I3,X,A3,I5,X,("-"),X,A2,':',A2,':',A2,X,3("*"))
C
C     ------------------------------------------------------------------
C
      END SUBROUTINE GETDATE 
