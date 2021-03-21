      SUBROUTINE BLDSOC0(ISHELL,JSHELL,KSHELL,LSHELL,TWOSOA,TWOSOB,
     $                   CA1,CA2,PABA,PABB,BLOCKS)
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
C
      INCLUDE 'fileio.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      REAL, DIMENSION(NSTO) :: CA1,CA2
      REAL, DIMENSION(NSTO,NSTO) :: PABA,PABB
      REAL, DIMENSION(2,3) :: SOCTWOB
      REAL, DIMENSION(4,3) :: MATR,SOCTWOA
      REAL, DIMENSION(DSHL,DSHL,DSHL,DSHL,3) :: BLOCKS
      REAL, DIMENSION(3) :: MATRIX,TWOSOA,TWOSOB
C
      INTEGER :: FUNAB,FUNCD
      INTEGER :: I,J,K,L,ISHELL,JSHELL,KSHELL,LSHELL,INA,INB,INC,IND
      INTEGER :: ILLSTO,IULSTO,JLLSTO,JULSTO,KLLSTO,KULSTO,LLLSTO,LULSTO
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      ILLSTO = LLSTO(ISHELL)
      IULSTO = ULSTO(ISHELL)
      JLLSTO = LLSTO(JSHELL)
      JULSTO = ULSTO(JSHELL)
      KLLSTO = LLSTO(KSHELL)
      KULSTO = ULSTO(KSHELL)
      LLLSTO = LLSTO(LSHELL)
      LULSTO = ULSTO(LSHELL)
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      DO INA=ILLSTO,IULSTO
        I = INA - LLSTO(ISHELL) + 1
        IF (ISHELL == JSHELL) JLLSTO = INA
        DO INB=JLLSTO,JULSTO
          J = INB - LLSTO(JSHELL) + 1
          DO INC=KLLSTO,KULSTO
            K = INC - LLSTO(KSHELL) + 1
            IF (KSHELL == LSHELL) LLLSTO = INC
            DO IND=LLLSTO,LULSTO
              L = IND - LLSTO(LSHELL) + 1
              IF ((INA.EQ.INB).AND.(INC.EQ.IND)) GOTO 100
              MATRIX(:) = NCSTO(INA)*NCSTO(INB)*NCSTO(INC)*NCSTO(IND)*
     $                    BLOCKS(I,J,K,L,:)
              MATR(1,:) = MATRIX(:)
              MATR(2,:) = -MATRIX(:)
              MATR(3,:) = MATRIX(:)
              MATR(4,:) = -MATRIX(:)
              IF (INA.EQ.INB) MATR(2,:) = 0.0
              IF (INC.EQ.IND) MATR(3,:) = 0.0
              IF (INC.EQ.IND) MATR(4,:) = 0.0
              SOCTWOA(1,:) = CA1(INC)*CA2(IND)*
     $                       PABA(INA,INB)*MATR(1,:) +
C
     $                       CA1(INC)*CA2(IND)*
     $                       PABA(INB,INA)*MATR(2,:) +
C
     $                       CA1(IND)*CA2(INC)*
     $                       PABA(INA,INB)*MATR(3,:) +
C
     $                       CA1(IND)*CA2(INC)*
     $                       PABA(INB,INA)*MATR(4,:)
C
              SOCTWOA(2,:) = CA1(INA)*CA2(INB)*
     $                       PABA(INC,IND)*MATR(1,:) +
C
     $                       CA1(INB)*CA2(INA)*
     $                       PABA(INC,IND)*MATR(2,:) +
C
     $                       CA1(INA)*CA2(INB)*
     $                       PABA(IND,INC)*MATR(3,:) +
C
     $                       CA1(INB)*CA2(INA)*
     $                       PABA(IND,INC)*MATR(4,:)
C
              SOCTWOA(3,:) = -CA1(INC)*CA2(INB)*
     $                        PABA(INA,IND)*MATR(1,:) -
C
     $                        CA1(INC)*CA2(INA)*
     $                        PABA(INB,IND)*MATR(2,:) -
C
     $                        CA1(IND)*CA2(INB)* 
     $                        PABA(INA,INC)*MATR(3,:) -
C
     $                        CA1(IND)*CA2(INA)*
     $                        PABA(INB,INC)*MATR(4,:)
C
              SOCTWOA(4,:) = -CA1(INA)*CA2(IND)*
     $                        PABA(INC,INB)*MATR(1,:) -
C
     $                        CA1(INB)*CA2(IND)*
     $                        PABA(INC,INA)*MATR(2,:) -
C
     $                        CA1(INA)*CA2(INC)*
     $                        PABA(IND,INB)*MATR(3,:) -
C
     $                        CA1(INB)*CA2(INC)*
     $                        PABA(IND,INA)*MATR(4,:)
C
              TWOSOA(:) = TWOSOA(:) + 3.0*(SOCTWOA(1,:) +
     $                                     SOCTWOA(2,:) +
     $                                     SOCTWOA(3,:) +
     $                                     SOCTWOA(4,:))
C
              SOCTWOB(1,:) = CA1(INC)*CA2(IND)*
     $                       PABB(INA,INB)*MATR(1,:) +
C
     $                       CA1(INC)*CA2(IND)*
     $                       PABB(INB,INA)*MATR(2,:) +
C
     $                       CA1(IND)*CA2(INC)*
     $                       PABB(INA,INB)*MATR(3,:) +
C
     $                       CA1(IND)*CA2(INC)*
     $                       PABB(INB,INA)*MATR(4,:)
C
              SOCTWOB(2,:) = -CA1(INA)*CA2(INB)*
     $                        PABB(INC,IND)*MATR(1,:) -
C
     $                        CA1(INB)*CA2(INA)*
     $                        PABB(INC,IND)*MATR(2,:) -
C
     $                        CA1(INA)*CA2(INB)*
     $                        PABB(IND,INC)*MATR(3,:) -
C
     $                        CA1(INB)*CA2(INA)*
     $                        PABB(IND,INC)*MATR(4,:)
C
              TWOSOB(:) = TWOSOB(:) + SOCTWOB(1,:) + SOCTWOB(2,:)
C
  100         CONTINUE
            END DO
          END DO
        END DO
      END DO
C
C     O________________________________________________________________O
C
      END SUBROUTINE BLDSOC0
