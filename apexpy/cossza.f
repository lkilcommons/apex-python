      SUBROUTINE COSSZA (GLAT,GLON,SBSLLAT,SBSLLON,CSZA)
C          Compute cosine of solar zenith angle from geographic coordinates of
C          point in question and of subsolar point.
C
C          INPUTS:
C            GLAT    = geographic latitude of point (degrees)
C            GLON    = geographic longitude of point (degrees)
C            SBSLLAT = geographic latitude of subsolar point (degrees)
C            SBSLLON = geographic longitude of subsolar point (degrees)
C          OUTPUT:
C            CSZA    = cosine of solar zenith angle
C
C          HISTORY:
C          Jul 1994: Completed on the 30th.  A. D. Richmond, NCAR
      implicit none
C
C     Parameter variables
C
      REAL                DTOR
      PARAMETER           (DTOR = 0.01745329251994330)
C
C     Argument variables
C
      REAL                CSZA,        GLAT,        GLON,        SBSLLAT
      REAL                SBSLLON
C
C     Local variables
C
      REAL                ANG,         CS,          CT,          SS
      REAL                ST
C
      CT   = SIN (GLAT*DTOR)
      ST   = SQRT (1. - CT*CT)
      CS   = SIN (SBSLLAT*DTOR)
      SS   = SQRT (1. - CS*CS)
      ANG  = (GLON - SBSLLON)*DTOR
      CSZA = CT*CS + ST*SS*COS(ANG)
      RETURN
      END
