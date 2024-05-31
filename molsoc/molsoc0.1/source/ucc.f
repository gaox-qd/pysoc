      CHARACTER*(*) FUNCTION UCC(STR)
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      CHARACTER :: STR*(*)
C
      INTEGER :: I
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      UCC = STR
      DO I=1,LEN(STR)
        IF ((ICHAR(STR(I:I)).GE.97).AND.(ICHAR(STR(I:I)).LE.122))
     $  UCC(I:I) = CHAR(ICHAR(STR(I:I)) - 32)
      END DO  
C
C     O________________________________________________________________O
C
      END FUNCTION UCC
