      SUBROUTINE SUBSOL (IYR,IDAY,IHR,IMN,SEC,SBSLLAT,SBSLLON)
C          Find subsolar geographic latitude and longitude given the
C          date and time (Universal Time).
C
C          This is based on formulas in Astronomical Almanac for the
C          year 1996, p.  C24. (U.S.  Government Printing Office,
C          1994).  According to the Almanac, results are good to at
C          least 0.01 degree latitude and 0.025 degree longitude
C          between years 1950 and 2050.  Accuracy for other years has
C          not been tested although the algorithm has been designed to
C          accept input dates from 1601 to 2100.  Every day is assumed
C          to have exactly 86400 seconds; thus leap seconds that
C          sometimes occur on June 30 and December 31 are ignored:
C          their effect is below the accuracy threshold of the algorithm.
C
C          961026 A. D. Richmond, NCAR
C
C          INPUTS:
C            IYR  = Year (e.g., 1994). IYR must be in the range: 1601 to 2100.
C            IDAY = Day number of year (e.g., IDAY = 32 for Feb 1)
C            IHR  = Hour of day    (e.g., 13 for 13:49)
C            IMN  = Minute of hour (e.g., 49 for 13:49)
C            SEC  = Second and fraction after the hour/minute.
C          Note:  While IYR is bounds tested, there is no constraint
C                 placed on values: IDAY,IHR,IMN,SEC; e.g., IHR=25 is
C                 properly interpreted.
C          RETURNS:
C            SBSLLAT = geographic latitude of subsolar point (degrees)
C            SBSLLON = geographic longitude of subsolar point (degrees,
C                      between -180 and +180)
C 
      implicit none
C
      INTEGER             NINT
C
      REAL                FLOAT
C
C     Parameter variables
C
      INTEGER             MSGUN
      PARAMETER           (MSGUN = 6)
C
      REAL                D2R
      PARAMETER           (D2R = 0.0174532925199432957692369076847)
      REAL                R2D
      PARAMETER           (R2D = 57.2957795130823208767981548147)
C
C     Argument variables
C
      INTEGER             IDAY,        IHR,         IMN,         IYR
C
      REAL                SBSLLAT,     SBSLLON,     SEC
C
C     Local variables
C
      INTEGER             NCENT,       NLEAP,       NROT,        YR
C
      REAL                ALPHA,       APTIME,      DELTA,       DF
      REAL                EPSILON,     EPSRAD,      ETDEG,       G
      REAL                G0,          GF,          GRAD,        L
      REAL                L0,          LAMBDA,      LAMRAD,      LF
      REAL                N,           SINLAM,      UT
C
C Number of years from 2000 to IYR (negative if IYR < 2000):
      YR = IYR - 2000
C
C NLEAP (final) = number of leap days from (2000 January 1) to (IYR January 1)
C                 (negative if IYR is before 1997)
      NLEAP = (IYR-1601)/4
      NLEAP = NLEAP - 99
      IF (IYR.LE.1900) THEN
	IF (IYR.LE.1600) THEN
	 WRITE(MSGUN,*) 'SUBSOLR INVALID BEFORE 1601: INPUT YEAR = ',IYR
	 STOP
	ENDIF
	NCENT = (IYR-1601)/100
	NCENT = 3 - NCENT 
	NLEAP = NLEAP + NCENT
      ENDIF
      IF (IYR.GE.2101) THEN
	WRITE(MSGUN,*) 'SUBSOLR INVALID AFTER 2100:  INPUT YEAR = ',IYR
	STOP
      ENDIF
C
C L0 = Mean longitude of Sun at 12 UT on January 1 of IYR:
C     L0 = 280.461 + .9856474*(365*(YR-NLEAP) + 366*NLEAP) 
C	   - (ARBITRARY INTEGER)*360.
C        = 280.461 + .9856474*(365*(YR-4*NLEAP) + (366+365*3)*NLEAP) 
C	   - (ARBITRARY INTEGER)*360.
C        = (280.461 - 360.) + (.9856474*365 - 360.)*(YR-4*NLEAP) 
C	   + (.9856474*(366+365*3) - 4*360.)*NLEAP,
C  where ARBITRARY INTEGER = YR+1.  This gives:
      L0 = -79.549 + (-.238699*(YR-4*NLEAP) + 3.08514E-2*NLEAP)
C
C G0 = Mean anomaly at 12 UT on January 1 of IYR:
C     G0 = 357.528 + .9856003*(365*(YR-NLEAP) + 366*NLEAP) 
C	   - (ARBITRARY INTEGER)*360.
C        = 357.528 + .9856003*(365*(YR-4*NLEAP) + (366+365*3)*NLEAP) 
C	   - (ARBITRARY INTEGER)*360.
C        = (357.528 - 360.) + (.9856003*365 - 360.)*(YR-4*NLEAP) 
C	   + (.9856003*(366+365*3) - 4*360.)*NLEAP,
C  where ARBITRARY INTEGER = YR+1.  This gives:
      G0 = -2.472 + (-.2558905*(YR-4*NLEAP) - 3.79617E-2*NLEAP)
C
C Universal time in seconds:
      UT = FLOAT(IHR*3600 + IMN*60) + SEC
C
C Days (including fraction) since 12 UT on January 1 of IYR:
      DF = (UT/86400. - 1.5) + IDAY
C
C Addition to Mean longitude of Sun since January 1 of IYR:
      LF = .9856474*DF
C
C Addition to Mean anomaly since January 1 of IYR:
      GF = .9856003*DF
C
C Mean longitude of Sun:
      L = L0 + LF
C
C Mean anomaly:
      G = G0 + GF
      GRAD = G*D2R
C
C Ecliptic longitude:
      LAMBDA = L + 1.915*SIN(GRAD) + .020*SIN(2.*GRAD)
      LAMRAD = LAMBDA*D2R
      SINLAM = SIN(LAMRAD)
C
C Days (including fraction) since 12 UT on January 1 of 2000:
      N = DF + FLOAT(365*YR + NLEAP)
C
C Obliquity of ecliptic:
      EPSILON = 23.439 - 4.E-7*N
      EPSRAD = EPSILON*D2R
C
C Right ascension:
      ALPHA = ATAN2(COS(EPSRAD)*SINLAM,COS(LAMRAD))*R2D
C
C Declination:
      DELTA = ASIN(SIN(EPSRAD)*SINLAM)*R2D
C
C Subsolar latitude:
      SBSLLAT = DELTA
C
C Equation of time (degrees):
      ETDEG = L - ALPHA
      NROT = NINT(ETDEG/360.)
      ETDEG = ETDEG - FLOAT(360*NROT)
C
C Apparent time (degrees):
      APTIME = UT/240. + ETDEG
C          Earth rotates one degree every 240 s.
C
C Subsolar longitude:
      SBSLLON = 180. - APTIME
      NROT = NINT(SBSLLON/360.)
      SBSLLON = SBSLLON - FLOAT(360*NROT)
C
      RETURN
      END
