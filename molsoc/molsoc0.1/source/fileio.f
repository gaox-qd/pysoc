      SUBROUTINE FILEIO(OPTION)
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
      INCLUDE 'molecule.h'
C
      CHARACTER :: OPTION*(*)
C
      CHARACTER :: NAME*10,STRING*80
      INTEGER :: ITA
      LOGICAL :: EXIS,NMD,OPE
C
C     ------------------------------------------------------------------
C
      SELECT CASE (OPTION)
C
C     *** Open permanent files ***
C
      CASE ('OPEN PERMANENT FILES')
C
        OPEN(UNIT=INP,FILE='molsoc.inp',STATUS='UNKNOWN',
     $       ACCESS='SEQUENTIAL',FORM='FORMATTED')
        OPEN(UNIT=BASIS,FILE='molsoc_basis',STATUS='UNKNOWN',
     $       POSITION='REWIND',ACCESS='SEQUENTIAL',FORM='FORMATTED')
        OPEN(OUT,FILE='molsoc.out',STATUS='UNKNOWN',POSITION='REWIND',
     $       ACCESS='SEQUENTIAL',FORM='FORMATTED')
        OPEN(SOCA,FILE='soint',STATUS='UNKNOWN',POSITION='REWIND',
     $       ACCESS='SEQUENTIAL',FORM='FORMATTED')
C        OPEN(COESP,FILE='coef',STATUS='UNKNOWN',POSITION='REWIND',
C     $       ACCESS='SEQUENTIAL',FORM='FORMATTED')
C        OPEN(IOAL,FILE='aover',STATUS='UNKNOWN',POSITION='REWIND',
C     $       ACCESS='SEQUENTIAL',FORM='FORMATTED')
C        OPEN(IOBE,FILE='bover',STATUS='UNKNOWN',POSITION='REWIND',
C     $       ACCESS='SEQUENTIAL',FORM='FORMATTED')
C
      CASE ('OPEN MOS PERMANENT FILES')
C
C        IF ((READMO == 1).OR.(READMO == 3).OR.(READMO == 4)
C     $      .OR.(READMO == 5)) THEN
C          OPEN(DEN1,FILE='mos1',STATUS='UNKNOWN',POSITION='REWIND',
C     $         ACCESS='SEQUENTIAL',FORM='FORMATTED')
C          OPEN(DEN2,FILE='mos2',STATUS='UNKNOWN',POSITION='REWIND',
C     $         ACCESS='SEQUENTIAL',FORM='FORMATTED')
C        ELSE IF (READMO == 2) THEN
C          OPEN(DEN1,FILE='mos1',STATUS='UNKNOWN',POSITION='REWIND',
C     $         ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
C          OPEN(DEN2,FILE='mos2',STATUS='UNKNOWN',POSITION='REWIND',
C     $         ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
C        END IF
C
C     *** Temporary file names ***
C
        FNAME(IOSO) = 'ioso.scr'
        FNAME(IOMO) = 'iomom.scr'
        FNAME(IOPM) = 'iopm.scr'
C
      CASE ('CLOSE ALL FILES')
C
        DO ITA=MINPER,MAXPER
          CLOSE (UNIT=ITA,STATUS='KEEP')
        END DO
C
        DO ITA=MINTEM,MAXTEM
          CLOSE (UNIT=ITA,STATUS='DELETE')
        END DO
C
C     *** Error handling ***
C
      CASE DEFAULT
C
        STRING = 'UNKNOWN OPTION '//OPTION
        WRITE(OUT,*) STRING
C
      END SELECT
C
C     ------------------------------------------------------------------
C
      END SUBROUTINE FILEIO
