      SUBROUTINE MAGLOCTM (ALON,SBSLLAT,SBSLLON,CLATP,POLON,MLT)
C  Computes magnetic local time from magnetic longitude, subsolar coordinates,
C   and geomagnetic pole coordinates.
C  950302 A. D. Richmond, NCAR
C  Algorithm:  MLT is calculated from the difference of the apex longitude,
C   alon, and the geomagnetic dipole longitude of the subsolar point.
C
C   Inputs:
C    ALON    = apex magnetic longitude of the point (deg)
C    SBSLLAT = geographic latitude of subsolar point (degrees)
C    SBSLLON = geographic longitude of subsolar point (degrees)
C    CLATP   = Geocentric colatitude of geomagnetic dipole north pole (deg)
C    POLON   = East longitude of geomagnetic dipole north pole (deg)
C
C   Output:
C    mlt (real) = magnetic local time for the apex longitude alon (hours)
C
C To go from mlt to alon (see comments following Entry mlt2alon for definition 
C  of variables), use:
C
C     CALL MLT2ALON (MLT,SBSLLAT,SBSLLON,CLATP,POLON,ALON)
C
C  NOTE: If the calling routine uses subroutine magloctm in conjunction with 
C   file magfld.f (which is used by subroutine APEX), then clatp and polon can 
C   be found by invoking
C
C     CALL DYPOL (CLATP,POLON,VP)
C
C   where vp is an unneeded variable.  (Note that subroutine COFRM must have
C   been called before DYPOL, in order to set up the coefficients for the
C   desired epoch.)  Alternatively, if subroutine apxntrp is
C   used to get alon from previously computed arrays, then
C   clatp and polon can be obtained for use in magloctm by adding
C
C     COMMON /APXDIPL/ CLATP,POLON,DUM1,DUM2,DUM3
C
C   to the calling routine (where DUM1,DUM2,DUM3 are unneeded dummy variables).
C
      implicit none
C
C     Argument variables
C
      REAL                ALON,        ALONX,       CLATP,       MLT
      REAL                POLON,       SBSLLAT,     SBSLLON,     XMLT
C
C     Local variables
C
      REAL                SMLON
C
	CALL SOLGMLON (SBSLLAT,SBSLLON,CLATP,POLON,SMLON)
	MLT = (ALON - SMLON)/15.0 + 12.0
	IF (MLT .GE. 24.0) MLT = MLT - 24.0
	IF (MLT .LT.   0.) MLT = MLT + 24.0
	RETURN
C
      ENTRY MLT2ALON (XMLT,SBSLLAT,SBSLLON,CLATP,POLON,ALONX)
C
C   Inputs:
C    XMLT (real) = magnetic local time for the apex longitude alon (hours,
C                 0. to 24.)
C    SBSLLAT     = geographic latitude of subsolar point (degrees)
C    SBSLLON     = geographic longitude of subsolar point (degrees)
C    CLATP       = Geocentric colatitude of geomagnetic dipole north pole (deg)
C    POLON       = East longitude of geomagnetic dipole north pole (deg)
C
C   Output:
C    ALONX       = apex magnetic longitude of the point (deg, -180. to 180.)
C
	CALL SOLGMLON (SBSLLAT,SBSLLON,CLATP,POLON,SMLON)
	ALONX = 15.*(XMLT - 12.0) + SMLON
	IF (ALONX .GT.  180.) ALONX = ALONX - 360.0
	IF (ALONX .LE. -180.) ALONX = ALONX + 360.0
	RETURN
	END

      SUBROUTINE SOLGMLON (XLAT,XLON,COLAT,ELON,MLON)
C Computes geomagnetic longitude of the point with geocentric spherical
C  latitude and longitude of XLAT and XLON, respectively.
C 940719 A. D. Richmond, NCAR
C Inputs:
C   XLAT  = geocentric spherical latitude (deg)
C   XLON  = geocentric spherical longitude (deg)
C   COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
C   ELON  = East longitude of geomagnetic dipole north pole (deg)
C Output:
C   MLON  = Geomagnetic dipole longitude of the point (deg, -180. to 180.)
      implicit none
C
C     Parameter variables
C
      REAL                RTOD
      PARAMETER           (RTOD = 5.72957795130823E1)
      REAL                DTOR
      PARAMETER           (DTOR = 1.745329251994330E-2)
C
C     Argument variables
C
      REAL                COLAT,       ELON,        MLON,        XLAT
      REAL                XLON
C
C     Local variables
C
      REAL                ANG,         CANG,        CTE,         CTP
      REAL                SANG,        STE,         STFCPA,      STFSPA
      REAL                STP
C
C Algorithm:
C   Use spherical coordinates.
C   Let GP be geographic pole.
C   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
C   Let XLON be longitude of point P.
C   Let TE be colatitude of point P.
C   Let ANG be longitude angle from GM to P.
C   Let TP be colatitude of GM.
C   Let TF be arc length between GM and P.
C   Let PA = MLON be geomagnetic longitude, i.e., Pi minus angle measured
C     counterclockwise from arc GM-P to arc GM-GP.
C   Then, using notation C=cos, S=sin, spherical-trigonometry formulas
C     for the functions of the angles are as shown below.  Note: STFCPA,
C     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.

      CTP = COS(COLAT*DTOR)
      STP = SQRT(1. - CTP*CTP)
      ANG = (XLON-ELON)*DTOR
      CANG = COS(ANG)
      SANG = SIN(ANG)
      CTE = SIN(XLAT*DTOR)
      STE = SQRT(1.-CTE*CTE)
      STFCPA = STE*CTP*CANG - CTE*STP
      STFSPA = SANG*STE
      MLON = ATAN2(STFSPA,STFCPA)*RTOD
      RETURN
      END
