      SUBROUTINE IOSCR(VECT,DIM,OPTION,KEYWORD,NUM)
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INCLUDE 'fileio.h'
C
      CHARACTER*(*) :: OPTION,KEYWORD
      CHARACTER :: ACCE*10,FRM*11,UCC*20,FILNAME*10
C
      INTEGER :: DIM,IONUM,I,NUM
      REAL, DIMENSION(DIM) :: VECT
C
      LOGICAL :: EXIS,OPE
C
C     ------------------------------------------------------------------
C
      SELECT CASE (KEYWORD)
C
      CASE ('ZERO')
        IONUM = IOSO
      CASE ('MOME')
        IONUM = IOMO
      CASE ('PMAT')
        IONUM = IOPM
      CASE DEFAULT
        WRITE(OUT,*) 'IOSCR:','UNKNOWN KEYWORD '//KEYWORD
        STOP
      END SELECT
C
      FILNAME = FNAME(IONUM)
C
      INQUIRE(FILE=FILNAME,EXIST=EXIS,OPENED=OPE,ERR=99910)
C
      IF (OPE) THEN
C
        INQUIRE(FILE=FILNAME,ACCESS=ACCE,FORM=FRM,ERR=99910)
C
        IF (UCC(ACCE) /= 'DIRECT') THEN
          WRITE(OUT,"('THE FILE ',A,
     $                ' (UNIT',I5,') WAS NOT OPENED WITH ',
     $                'DIRECT ACCESS')") FILNAME,IONUM
          STOP
        END IF
C
        IF (UCC(FRM) /= 'UNFORMATTED') THEN
          WRITE(OUT,"('THE FILE ',A,
     $                ' (UNIT',I5,') WAS NOT OPENED AS ',
     $                'UNFORMATTED')") FILNAME,IONUM
          STOP
        END IF
C
      END IF
C
      IF (OPTION == 'WRITE') THEN
C
        IF (.NOT.OPE) THEN
          OPEN(UNIT=IONUM,FILE=FILNAME,STATUS='UNKNOWN',
     $         ACCESS='DIRECT',FORM='UNFORMATTED',RECL=8*DIM,ERR=99915)
        END IF
C
        WRITE(UNIT=IONUM,REC=NUM,ERR=99920) (VECT(I),I=1,DIM)


C
      ELSE IF (OPTION == 'READ') THEN
C
        IF (OPE) THEN
          READ(UNIT=IONUM,REC=NUM,ERR=99925) (VECT(I),I=1,DIM)
        ELSE IF (EXIS) THEN
          OPEN(UNIT=IONUM,FILE=FILNAME,STATUS='UNKNOWN',
     $          ACCESS='DIRECT',FORM='UNFORMATTED',RECL=8*DIM,ERR=99915)
          READ(UNIT=IONUM,REC=NUM,ERR=99925) (VECT(I),I=1,DIM)
        ELSE
          WRITE(OUT,"('THE FILE ',A,' (UNIT',I5,') CANNOT BE FOUND ',
     $                'TO BE READ')") FILNAME,IONUM
          STOP
        END IF
C
      ELSE 
        WRITE(OUT,"('UNKNOWN IOSCR OPTION')")
        STOP
C
      END IF
C
      RETURN
C
C     ------------------------------------------------------------------
C
99910 CONTINUE
      WRITE(OUT,"('ERROR TRYING TO INQUIRE THE FILE ',A,
     $            ' (Unit',I5,')')") FILNAME,IONUM
      STOP
C
99915 CONTINUE
      WRITE(OUT,"('ERROR TRYING TO OPEN THE FILE ',A,
     $            ' (Unit',I5,')')") FILNAME,IONUM
      STOP
C
99920 CONTINUE
      WRITE(OUT,"('ERROR TRYING TO WRITE IN THE FILE ',A,
     $            ' (Unit',I5,')')") FILNAME,IONUM
      STOP
C
99925 CONTINUE
      WRITE(OUT,"('ERROR TRYING TO READ FROM THE FILE ',A,
     $            ' (Unit',I5,')')") FILNAME,IONUM
      STOP
C
C     ------------------------------------------------------------------
C
C     *** End of SUBROUTINE IOSCR ***
C
      END
