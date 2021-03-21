      SUBROUTINE BLDSOC 
C
C     ******************************************************************
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Calculate the occupation numbers ___
C
      CALL INIOCC
C
C     ___ Calculate UKS or SROKS SOC contributions ___
C
      CALL KSSOC
C
C     O_________________________________________________________________
C
      END SUBROUTINE BLDSOC
