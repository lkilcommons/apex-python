      SUBROUTINE APEX (DATE,DLAT,DLON,ALT,
     +                 A,ALAT,ALON,BMAG,XMAG,YMAG,ZMAG,V)
C          Calculate apex radius, latitude, longitude; and magnetic field and
C          scalar magnetic potential.
C
C          INPUTS:
C            DATE = Year and fraction (1990.0 = 1990 January 1, 0 UT)
C            DLAT = Geodetic latitude in degrees
C            DLON = Geodetic longitude in degrees
C            ALT = Altitude in km
C
C          RETURNS:
C            A    = (Apex height + REQ)/REQ, where REQ = equatorial Earth radius.
C                   A is analogous to the L value in invariant coordinates.
C            ALAT = Apex latitude in degrees (negative in S. magnetic hemisphere)
C            ALON = Apex longitude (geomagnetic longitude of apex) in degrees
C            BMAG = geomagnetic field magnitude (nT)
C            XMAG = geomagnetic field component (nT): north
C            YMAG = geomagnetic field component (nT): east
C            ZMAG = geomagnetic field component (nT): downward
C            V    = geomagnetic potential (T-m)
C
C          COMMON BLOCKS:
C            COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
C
C          DIPOLE has IGRF variables obtained from routines in magfld.f:
C            COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
C            ELON  = East longitude of geomagnetic dipole north pole (deg)
C            VP    = Magnitude (T-m) of dipole component of magnetic potential at
C                    geomagnetic pole and geocentric radius of 6371.2 km
C            CTP   = cosine of COLAT
C            STP   = sine   of COLAT
C
C------------------------------------------------------------------------------
C          HISTORY:
C          Aug 1994: First version completed on the 22nd by A.D. Richmond.
C          May 1999: Revise DS calculation in LINAPX to avoid divide by zero.
C          Apr 2004: - Change definition of earth's equatorial radius (REQ)
C                      from the IAU-1966 spheroid (6378.160 km) to the WGS-1984
C                      spheroid (6378.137 km); see description in
C                      '$APXROOT/docs/igrf.2004.spheroid.desc'.
C                    - Revise comments toward a consistent format so they are
C                      easy to read.
C                    - Replace computed GO TO in ITRACE with IF blocks.
C                    - Refine FNDAPX to insure |Bdown/Btot| < 1.E-6 at apex

      implicit none
C
C     Parameter variables
C
      REAL                RE
      PARAMETER           (RE = 6371.2)
      REAL                DTOR
      PARAMETER           (DTOR = .01745329251994330)
C
C     Argument variables
C
      REAL                A,           ALAT,        ALON,        ALT
      REAL                BMAG,        DATE,        DLAT,        DLON
      REAL                V,           XMAG,        YMAG,        ZMAG
C
C     Local variables
C
      REAL                BX,          BY,          BZ,          CLATP
      REAL                POLON,       VPOL,        X,           Y
      REAL                Z
C
C     Common variables
C
      REAL                COLAT,       CTP,         ELON,        STP
      REAL                VP
C
      COMMON / DIPOLE /   COLAT,       ELON,        VP,          CTP
      COMMON / DIPOLE /   STP

      CALL COFRM (DATE)
      CALL DYPOL (CLATP,POLON,VPOL)
      COLAT = CLATP
      CTP   = COS(CLATP*DTOR)
      STP   = SQRT(1. - CTP*CTP)
      ELON  = POLON
      VP    = VPOL
      CALL LINAPX (DLAT,DLON,ALT, A,ALAT,ALON,XMAG,YMAG,ZMAG,BMAG)
      XMAG = XMAG*1.E5                 ! convert from gauss to nT
      YMAG = YMAG*1.E5
      ZMAG = ZMAG*1.E5
      BMAG = BMAG*1.E5
      CALL GD2CART (DLAT,DLON,ALT,X,Y,Z)
      CALL FELDG (3, X/RE,Y/RE,Z/RE, BX,BY,BZ,V)
      RETURN
      END

      SUBROUTINE LINAPX (GDLAT,GLON,ALT, A,ALAT,ALON,XMAG,YMAG,ZMAG,F)

C          Transform geographic coordinates to Apex coordinates.
C
C          INPUTS:
C            GDLAT = Latitude  (degrees, positive northward)
C            GLON  = Longitude (degrees, positive eastward)
C            ALT   = Height of starting point (km above mean sea level)
C
C          OUTPUTS:
C            A     = (Apex height + REQ)/REQ, where REQ = equatorial Earth radius.
C                    A is analogous to the L value in invariant coordinates.
C            ALAT  = Apex Lat. (deg)
C            ALON  = Apex Lon. (deg)
C            XMAG  = Geomagnetic field component (gauss): north
C            YMAG  = Geomagnetic field component (gauss): east
C            ZMAG  = Geomagnetic field component (gauss): down
C            F     = Geomagnetic field magnitude (gauss)
C
C          Trace the geomagnetic field line from the given location to find the
C          apex of the field line.  Before starting iterations to trace along
C          the field line: (1) Establish a step size (DS, arc length in km)
C          based on the geomagnetic dipole latitude; (2) determine the step
C          direction from the sign of the vertical component of the geomagnetic
C          field; and (3) convert to geocentric cartesian coordinates.  Each
C          iteration increments a step count (NSTP) and calls ITRACE to move
C          along the the field line until reaching the iteration count limit
C          (MAXS) or passing the apex (IAPX=2) and then calling FNDAPX to
C          determine the apex location from the last three step locations
C          (YAPX); however, if reaching the iteration limit, apex coordinates
C          are calculated by DIPAPX which assumes a simplified dipole field.
C
C          COMMON BLOCKS:
C            COMMON /APXIN/   YAPX(3,3)
C            COMMON /DIPOLE/  COLAT,ELON,VP,CTP,STP
C            COMMON /FLDCOMD/ BX, BY, BZ, BB
C            COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
C
C          APXIN has step locations determined in ITRACE:
C            YAPX  = Matrix of cartesian coordinates (loaded columnwise) of the
C                    three points about the apex.  Set in subroutine ITRACE.
C                                                                               
C          DIPOLE has IGRF variables obtained from routines in magfld.f:
C            COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
C            ELON  = East longitude of geomagnetic dipole north pole (deg)
C            VP    = Magnitude (T-m) of dipole component of magnetic potential at
C                    geomagnetic pole and geocentric radius of 6371.2 km
C            CTP   = cosine of COLAT
C            STP   = sine   of COLAT
C                                                                               
C          FLDCOMD has geomagnetic field at current trace point:
C            BX    = X component (Gauss)
C            BY    = Y component (Gauss)
C            BZ    = Z component (Gauss)
C            BB    = Magnitude   (Gauss)
C
C          ITRA has field line tracing variables determined in LINAPX:
C            NSTP  = Step count.
C            Y     = Array containing current tracing point cartesian coordinates.
C            YOLD  = Array containing previous tracing point cartesian coordinates.
C            SGN   = Determines direction of trace.
C            DS    = Step size (arc length in km).
C                                                                               
C          REFERENCES:
C            Stassinopoulos E. G. , Mead Gilbert D., X-841-72-17 (1971) GSFC,
C            Greenbelt, Maryland
C                                                                               
C          EXTERNALS:
C            GD2CART = Convert geodetic to geocentric cartesian coordinates (in magfld.f)
C            CONVRT  = Convert geodetic to geocentric cylindrical or geocentric spherical
C                      and back (in magfld.f).
C            FELDG   = Obtain IGRF magnetic field components (in magfld.f).
C            ITRACE  = Follow a geomagnetic field line
C            DIPAPX  = Compute apex coordinates assuming a geomagnetic dipole field
C            FNDAPX  = Compute apex coordinates from the last three traced field line points
C
C------------------------------------------------------------------------------
C          HISTORY:
C          Oct 1973: Initial version completed on the 29th by Wally Clark, NOAA
C                    ERL Lab.
C          Feb 1988: Revised on the 1st by Harsh Anand Passi, NCAR.
C          Aug 1994: Revision by A. D. Richmond, NCAR.

      implicit none
C
      REAL                AMAX1
C
C     Parameter variables
C
      INTEGER             MAXS
      PARAMETER           (MAXS = 200)
C
      REAL                RTOD
      PARAMETER           (RTOD = 57.2957795130823)
      REAL                RE
      PARAMETER           (RE = 6371.2)
      REAL                DTOR
      PARAMETER           (DTOR = .01745329251994330)
      REAL                REQ
      PARAMETER           (REQ = 6378.137)
C
C     Argument variables
C
      REAL                A,           ALAT,        ALON,        ALT
      REAL                F,           GDLAT,       GLON,        XMAG
      REAL                YMAG,        ZMAG
C
C     Local variables
C
      INTEGER             I,           IAPX,        J
C
      REAL                BABS,        BDOWN,       BEAST,       BNRTH
      REAL                CGML2,       GCLAT,       HT,          R
      REAL                RHO,         SINGML,      XLAT,        XLON
C
C     Common variables
C
      REAL                YAPX(3,3)
C
      COMMON / APXIN  /   YAPX
C
C     Common variables
C
      REAL                COLAT,       CTP,         ELON,        STP
      REAL                VP
C
      COMMON / DIPOLE /   COLAT,       ELON,        VP,          CTP
      COMMON / DIPOLE /   STP
C
C     Common variables
C
      REAL                BB,          BX,          BY,          BZ
C
      COMMON / FLDCOMD/   BX,          BY,          BZ,          BB
C
C     Common variables
C
      INTEGER             NSTP
C
      REAL                DS,          SGN,         Y(3),        YP(3)
C
      COMMON / ITRA   /   NSTP,        Y,           YP,          SGN
      COMMON / ITRA   /   DS


C          Set step size based on the geomagnetic dipole latitude of the starting point
      CALL CONVRT (2,GDLAT,ALT,GCLAT,R)
      SINGML = CTP*SIN(GCLAT*DTOR) + STP*COS(GCLAT*DTOR)*
     +                                             COS((GLON-ELON)*DTOR)
C          May 1999: avoid possible divide by zero (when SINGML = 1.): the old version
C          limited DS to its value at 60 deg GM latitude with: DS = .06*R/(1.-SINGML*SINGML) - 370.
C                                                              IF (DS .GT. 1186.) DS = 1186.
      CGML2 = AMAX1 (0.25,1.-SINGML*SINGML)
      DS = .06*R/CGML2 - 370.

C          Initialize YAPX array
      DO 4 J=1,3
      DO 4 I=1,3
    4 YAPX(I,J) = 0.

C          Convert from geodetic to earth centered cartesian coordinates
      CALL GD2CART (GDLAT,GLON,ALT,Y(1),Y(2),Y(3))
      NSTP = 0

C          Get magnetic field components to determine the direction for
C          tracing the field line
      CALL FELDG (1,GDLAT,GLON,ALT,XMAG,YMAG,ZMAG,F)
      SGN = SIGN (1.,-ZMAG)

C          Use cartesian coordinates to get magnetic field components
C          (from which gradients steer the tracing)
   10 CALL FELDG (2, Y(1)/RE,Y(2)/RE,Y(3)/RE, BX,BY,BZ,BB)
      NSTP = NSTP + 1

      IF (NSTP .LT. MAXS) THEN
        ! trace along field line
	CALL ITRACE (IAPX)
	IF (IAPX .EQ. 1) GO TO 10
        ! (IAPX=2) => passed max radius; find its coordinates
	CALL FNDAPX (ALT,ZMAG,A,ALAT,ALON)
      ELSE
	RHO = SQRT (Y(1)*Y(1) + Y(2)*Y(2))               ! too many steps; get apex from dipole approximation
	CALL CONVRT (3,XLAT,HT,RHO,Y(3))
	XLON = RTOD*ATAN2 (Y(2),Y(1))
	CALL FELDG (1,XLAT,XLON,HT,BNRTH,BEAST,BDOWN,BABS)
	CALL DIPAPX  (XLAT,XLON,HT,BNRTH,BEAST,BDOWN,A,ALON)
	ALAT = -SGN*RTOD*ACOS (SQRT(1./A))
      ENDIF

      RETURN
      END

      SUBROUTINE ITRACE (IAPX)

C          Follow a geomagnetic field line until passing its apex
C
C          INPUTS:
C            (all are in common blocks)
C          OUTPUTS:
C            IAPX = 2 (when apex passed) or 1 (not)
C                                                                               
C          This uses the 4-point Adams formula after initialization.
C          First 7 iterations advance point by 3 steps.
C                                                                               
C          COMMON BLOCKS:
C            COMMON /APXIN/   YAPX(3,3)
C            COMMON /FLDCOMD/ BX, BY, BZ, BB
C            COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
C
C          APXIN has step locations determined in ITRACE:
C            YAPX  = Matrix of cartesian coordinates (loaded columnwise) of the
C                    three points about the apex.  Set in subroutine ITRACE.
C                                                                               
C          FLDCOMD has geomagnetic field at current trace point:
C            BX    = X component (Gauss)
C            BY    = Y component (Gauss)
C            BZ    = Z component (Gauss)
C            BB    = Magnitude   (Gauss)
C                                                                               
C          ITRA has field line tracing variables determined in LINAPX:
C            NSTP  = Step count.
C            Y     = Array containing current tracing point cartesian coordinates.
C            YOLD  = Array containing previous tracing point cartesian coordinates.
C            SGN   = Determines direction of trace.
C            DS    = Step size (arc length in km).
C
C          REFERENCES:
C            Stassinopoulos E. G. , Mead Gilbert D., X-841-72-17 (1971) GSFC,
C            Greenbelt, Maryland
C------------------------------------------------------------------------------
C          HISTORY:
C          Oct 1973: Initial version completed on the 29th by W. Clark, NOAA ERL
C                    Laboratory.
C          Feb 1988: Revised by H. Passi, NCAR.
C          Apr 2004: Replace computed GO TO with IF blocks because some compilers
C                    are threatening to remove this old feature
C
      implicit none
      SAVE
C
C     Statement functions
C
      REAL                RDUS
C
C     Argument variables
C
      INTEGER             IAPX
C
C     Local variables
C
      INTEGER             I,           J
C
      REAL                D,           D12,         D2,          D24
      REAL                D6,          E,           F,           RC
      REAL                RP,          YP(3,4)
C
C     Common variables
C
      REAL                YAPX(3,3)
C
      COMMON / APXIN  /   YAPX
C
C     Common variables
C
      REAL                BB,          BX,          BY,          BZ
C
      COMMON / FLDCOMD/   BX,          BY,          BZ,          BB
C
C     Common variables
C
      INTEGER             NSTP
C
      REAL                DS,          SGN,         Y(3),        YOLD(3)
C
      COMMON / ITRA   /   NSTP,        Y,           YOLD,        SGN
      COMMON / ITRA   /   DS
C          Statement function
      RDUS(D,E,F) = SQRT (D**2 + E**2 + F**2)

      IAPX = 1

C          Cartesian component magnetic field (partial) derivitives steer the trace
      YP(1,4) = SGN*BX/BB
      YP(2,4) = SGN*BY/BB
      YP(3,4) = SGN*BZ/BB

      IF (NSTP .LE. 7) THEN
	DO 10 I=1,3
	IF (NSTP .EQ. 1) THEN
	  D2        = DS/2.
	  D6        = DS/6.
	  D12       = DS/12.
	  D24       = DS/24.
	  YP(I,1)   = YP(I,4)
	  YOLD(I)   = Y(I)
	  YAPX(I,1) = Y(I)
	  Y(I)      = YOLD(I) + DS*YP(I,1)

	ELSE IF (NSTP .EQ. 2) THEN
	  YP(I,2) = YP(I,4)
	  Y(I)    = YOLD(I) + D2*(YP(I,2)+YP(I,1))

	ELSE IF (NSTP .EQ. 3) THEN
	  Y(I) = YOLD(I) + D6*(2.*YP(I,4)+YP(I,2)+3.*YP(I,1))

	ELSE IF (NSTP .EQ. 4) THEN
	  YP(I,2)   = YP(I,4)
	  YAPX(I,2) = Y(I)
	  YOLD(I)   = Y(I)
	  Y(I)      = YOLD(I) + D2*(3.*YP(I,2)-YP(I,1))

	ELSE IF (NSTP .EQ. 5) THEN
	  Y(I) = YOLD(I) + D12*(5.*YP(I,4)+8.*YP(I,2)-YP(I,1))

	ELSE IF (NSTP .EQ. 6) THEN
	  YP(I,3)   = YP(I,4)
	  YOLD(I)   = Y(I)
	  YAPX(I,3) = Y(I)
	  Y(I)      = YOLD(I) + D12*(23.*YP(I,3)-16.*YP(I,2)+5.*YP(I,1))

	ELSE IF (NSTP .EQ. 7) THEN
	  YAPX(I,1) = YAPX(I, 2)
	  YAPX(I,2) = YAPX(I, 3)
	  Y(I)      = YOLD(I) + D24*(9.*YP(I,4) + 19.*YP(I,3) -
     +                               5.*YP(I,2) +     YP(I,1))
	  YAPX(I,3) = Y(I)
	ENDIF
   10   CONTINUE
	IF (NSTP .EQ. 6 .OR. NSTP .EQ. 7) THEN        ! signal if apex passed
	  RC = RDUS (YAPX(1,3), YAPX(2,3), YAPX(3,3))
	  RP = RDUS (YAPX(1,2), YAPX(2,2), YAPX(3,2))
	  IF (RC .LT. RP) IAPX = 2
	ENDIF

      ELSE                 ! NSTP > 7

	DO 30 I=1,3
	YAPX(I,1) = YAPX(I,2)
	YAPX(I,2) = Y(I)
	YOLD(I)   = Y(I)
	Y(I)      = YOLD(I) + D24*(55.*YP(I,4) - 59.*YP(I,3) +
     +                             37.*YP(I,2) -  9.*YP(I,1))
	YAPX(I,3) = Y(I)

	DO 20 J=1,3
   20   YP(I,J) = YP(I,J+1)
   30   CONTINUE
	RC = RDUS (   Y(1),    Y(2),    Y(3))
	RP = RDUS (YOLD(1), YOLD(2), YOLD(3))
	IF (RC .LT. RP) IAPX = 2
      ENDIF

      RETURN
      END

       SUBROUTINE FNDAPX (ALT,ZMAG,A,ALAT,ALON)

C          Find apex coordinates once tracing (in subroutine ITRACE) has
C          signalled that the apex has been passed.
C          INPUTS:
C            ALT  = Altitude of starting point
C            ZMAG = Downward component of geomagnetic field at starting point
C          OUTPUT
C            A    = Apex radius, defined as (Apex height + Req)/Req, where
C                   Req = equatorial Earth radius.
C                   A is analogous to the L value in invariant coordinates.
C            ALAT = Apex Lat. (deg)
C            ALON = Apex Lon. (deg)
C
C          COMMON BLOCKS:
C            COMMON /APXIN/  YAPX(3,3)
C            COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
C
C          APXIN has step locations determined in ITRACE:
C            YAPX  = Matrix of cartesian coordinates (loaded columnwise) of the
C                    three points about the apex.  Set in subroutine ITRACE.
C                                                                               
C          DIPOLE has IGRF variables obtained from routines in magfld.f:
C            COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
C            ELON  = East longitude of geomagnetic dipole north pole (deg)
C            VP    = Magnitude (T-m) of dipole component of magnetic potential at
C                    geomagnetic pole and geocentric radius of 6371.2 km
C            CTP   = cosine of COLAT
C            STP   = sine   of COLAT
C                                                                               
C          EXTERNALS:
C            FINT = Second degree interpolation routine
C------------------------------------------------------------------------------
C          HISTORY:
C          Oct 1973: Initial version completed on the 23rd by Clark, W., NOAA
C                    Boulder.
C          Aug 1994: Revision on the 3rd by A.D. Richmond, NCAR
C          Apr 2004: Repair problem noted by Dan Weimer where the apex location
C                    produced by FINT may still have a non-zero vertical magnetic
C                    field component.

       implicit none
C
C     External functions
C
      EXTERNAL            FINT
C
      REAL                AMAX1,       FINT
C
C     Parameter variables
C
      REAL                RTOD
      PARAMETER           (RTOD = 57.2957795130823)
      REAL                DTOR
      PARAMETER           (DTOR = .01745329251994330)
      REAL                REQ
      PARAMETER           (REQ = 6378.137)
C
C     Argument variables
C
      REAL                A,           ALAT,        ALON,        ALT
      REAL                ZMAG
C
C     Local variables
C
      INTEGER             I,           NITR
C
      REAL                ABDOB,       ANG,         BA,          BD(3)
      REAL                BDA,         BE,          BEA,         BMAG
      REAL                BN,          BNA,         CANG,        CTE
      REAL                GDLN,        GDLT,        HT,          HTA
      REAL                R,           RASQ,        RHO,         SANG
      REAL                STE,         STFCPA,      STFSPA,      XLON
      REAL                Y(3)
C
C     Common variables
C
      REAL                YAPX(3,3)
C
      COMMON / APXIN  /   YAPX
C
C     Common variables
C
      REAL                COLAT,       CTP,         ELON,        STP
      REAL                VP
C
      COMMON / DIPOLE /   COLAT,       ELON,        VP,          CTP
      COMMON / DIPOLE /   STP
C
C          Get geodetic height and vertical (downward) component of the magnetic
C          field at last three points found by ITRACE
      DO 10 I=1,3
      RHO  = SQRT (YAPX(1,I)**2 + YAPX(2,I)**2)
      CALL CONVRT (3,GDLT,HT, RHO,YAPX(3,I))
      GDLN = RTOD*ATAN2 (YAPX(2,I),YAPX(1,I))
   10 CALL FELDG (1,GDLT,GDLN,HT, BN,BE,BD(I),BMAG)

C          Interpolate to where Bdown=0 to find cartesian coordinates at dip equator
      NITR = 0
   20 Y(1) = FINT (BD(1),BD(2),BD(3),YAPX(1,1),YAPX(1,2),YAPX(1,3), 0.)
      Y(2) = FINT (BD(1),BD(2),BD(3),YAPX(2,1),YAPX(2,2),YAPX(2,3), 0.)
      Y(3) = FINT (BD(1),BD(2),BD(3),YAPX(3,1),YAPX(3,2),YAPX(3,3), 0.)

C          Insure negligible Bdown or
C
C            |Bdown/Btot| < 2.E-6
C
C          For instance, Bdown must be less than 0.1 nT at low altitudes where
C          Btot ~ 50000 nT.  This ratio can be exceeded when interpolation is
C          not accurate; i.e., when the middle of the three points interpolated
C          is too far from the dip equator.  The three points were initially
C          defined with equal spacing by ITRACE, so replacing point 2 with the
C          most recently fit location will reduce the interpolation span.
      RHO  = SQRT (Y(1)**2 + Y(2)**2)
      GDLN = RTOD*ATAN2 (Y(2),Y(1))
      CALL CONVRT (3,GDLT,HTA, RHO,Y(3))
      CALL FELDG (1,GDLT,GDLN,HTA, BNA,BEA,BDA,BA)
      ABDOB = ABS(BDA/BA)

      IF (ABDOB .GT. 2.E-6) THEN
	IF (NITR .LT. 4) THEN        ! 4 was chosen because tests rarely required 2 iterations
	  NITR      = NITR + 1
	  YAPX(1,2) = Y(1)
	  YAPX(2,2) = Y(2)
	  YAPX(3,2) = Y(3)
	  BD(2)     = BDA
	  GO TO 20
	ELSE
	  WRITE (0,'(''APEX: Imprecise fit of apex: |Bdown/B| ='',1PE7.1
     +    )') ABDOB
	ENDIF
      ENDIF

C          Ensure altitude of the Apex is at least the initial altitude when
C          defining the Apex radius then use it to define the Apex latitude whose
C          hemisphere (sign) is inferred from the sign of the dip angle at the
C          starting point
      A = (REQ + AMAX1(ALT,HTA)) / REQ
      IF (A .LT. 1.) THEN
	WRITE (0,'(''APEX: A can not be less than 1; A, REQ, HTA: '',1P3
     +E15.7)') A,REQ,HTA
	CALL EXIT (1)
      ENDIF
      RASQ = ACOS (SQRT(1./A))*RTOD
      ALAT = SIGN (RASQ,ZMAG)

C          ALON is the dipole longitude of the apex and is defined using
C          spherical coordinates where
C            GP   = geographic pole.
C            GM   = geomagnetic pole (colatitude COLAT, east longitude ELON).
C            XLON = longitude of apex.
C            TE   = colatitude of apex.
C            ANG  = longitude angle from GM to apex.
C            TP   = colatitude of GM.
C            TF   = arc length between GM and apex.
C            PA   = ALON be geomagnetic longitude, i.e., Pi minus angle measured
C                   counterclockwise from arc GM-apex to arc GM-GP.
C          then, spherical-trigonometry formulas for the functions of the angles
C          are as shown below.  Notation uses C=cos, S=sin and STFCPA = sin(TF) * cos(PA),
C                                                              STFSPA = sin(TF) * sin(PA)
      XLON = ATAN2 (Y(2),Y(1))
      ANG  = XLON-ELON*DTOR
      CANG = COS (ANG)
      SANG = SIN (ANG)
      R    = SQRT (Y(1)**2 + Y(2)**2 + Y(3)**2)
      CTE  = Y(3)/R
      STE  = SQRT (1.-CTE*CTE)
      STFCPA = STE*CTP*CANG - CTE*STP
      STFSPA = SANG*STE
      ALON = ATAN2 (STFSPA,STFCPA)*RTOD
      RETURN
      END

      SUBROUTINE DIPAPX (GDLAT,GDLON,ALT,BNORTH,BEAST,BDOWN, A,ALON)

C          Compute A, ALON from local magnetic field using dipole and spherical
C          approximation.
C
C          INPUTS:
C            GDLAT  = geodetic latitude, degrees
C            GDLON  = geodetic longitude, degrees
C            ALT    = altitude, km
C            BNORTH = geodetic northward magnetic field component (any units)
C            BEAST  = eastward magnetic field component
C            BDOWN  = geodetic downward magnetic field component
C          OUTPUTS:
C            A      = apex radius, 1 + h_A/R_eq
C            ALON   = apex longitude, degrees
C
C          Use spherical coordinates and define:
C            GP    = geographic pole.
C            GM    = geomagnetic pole (colatitude COLAT, east longitude ELON).
C            G     = point at GDLAT,GDLON.
C            E     = point on sphere below apex of dipolar field line passing
C                    through G.
C            TD    = dipole colatitude of point G, found by applying dipole
C                    formula for dip angle to actual dip angle.
C            B     = Pi plus local declination angle.  B is in the direction
C                    from G to E.
C            TG    = colatitude of G.
C            ANG   = longitude angle from GM to G.
C            TE    = colatitude of E.
C            TP    = colatitude of GM.
C            A     = longitude angle from G to E.
C            APANG = A + ANG
C            PA    = geomagnetic longitude, i.e., Pi minus angle measured
C                    counterclockwise from arc GM-E to arc GM-GP.
C            TF    = arc length between GM and E.
C          Then, using notation C=cos, S=sin, COT=cot, spherical-trigonometry
C          formulas for the functions of the angles are as shown below.  Note:
C            STFCPA = sin(TF) * cos(PA)
C            STFSPA = sin(TF) * sin(PA)
C
C          COMMON BLOCKS:
C            COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
C
C          DIPOLE has IGRF variables obtained from routines in magfld.f:
C            COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
C            ELON  = East longitude of geomagnetic dipole north pole (deg)
C            VP    = Magnitude (T-m) of dipole component of magnetic potential at
C                    geomagnetic pole and geocentric radius of 6371.2 km
C            CTP   = cosine of COLAT
C            STP   = sine   of COLAT
C------------------------------------------------------------------------------
C          HISTORY:
C          May 1994:  Completed on the 1st by A. D. Richmond
      implicit none
C
C     Parameter variables
C
      REAL                RTOD
      PARAMETER           (RTOD = 57.2957795130823)
      REAL                RE
      PARAMETER           (RE = 6371.2)
      REAL                DTOR
      PARAMETER           (DTOR = .01745329251994330)
      REAL                REQ
      PARAMETER           (REQ = 6378.137)
C
C     Argument variables
C
      REAL                A,           ALON,        ALT,         BDOWN
      REAL                BEAST,       BNORTH,      GDLAT,       GDLON
C
C     Local variables
C
      REAL                ANG,         BHOR,        CA,          CANG
      REAL                CAPANG,      CB,          COTTD,       CTD
      REAL                CTE,         CTG,         HA,          R
      REAL                SA,          SANG,        SAPANG,      SB
      REAL                STD,         STE,         STFCPA,      STFSPA
      REAL                STG
C
C     Common variables
C
      REAL                COLAT,       CTP,         ELON,        STP
      REAL                VP
C
      COMMON / DIPOLE /   COLAT,       ELON,        VP,          CTP
      COMMON / DIPOLE /   STP
C

      BHOR = SQRT(BNORTH*BNORTH + BEAST*BEAST)
      IF (BHOR .EQ. 0.) THEN
	ALON = 0.
	A    = 1.E34
	RETURN
      ENDIF
      COTTD  = BDOWN*.5/BHOR
      STD    = 1./SQRT(1.+COTTD*COTTD)
      CTD    = COTTD*STD
      SB     = -BEAST /BHOR
      CB     = -BNORTH/BHOR
      CTG    = SIN (GDLAT*DTOR)
      STG    = COS (GDLAT*DTOR)
      ANG    = (GDLON-ELON)*DTOR
      SANG   = SIN(ANG)
      CANG   = COS(ANG)
      CTE    = CTG*STD + STG*CTD*CB
      STE    = SQRT(1. - CTE*CTE)
      SA     = SB*CTD/STE
      CA     = (STD*STG - CTD*CTG*CB)/STE
      CAPANG = CA*CANG - SA*SANG
      SAPANG = CA*SANG + SA*CANG
      STFCPA = STE*CTP*CAPANG - CTE*STP
      STFSPA = SAPANG*STE
      ALON = ATAN2 (STFSPA,STFCPA)*RTOD
      R    = ALT + RE
      HA   = ALT + R*COTTD*COTTD
      A    = 1. + HA/REQ
      RETURN
      END

      FUNCTION FINT (X1,X2,X3,Y1,Y2,Y3, XFIT)
C          Second degree interpolation used by FNDAPX
C          INPUTS:
C            X1   = point 1 ordinate value
C            X2   = point 2 ordinate value
C            X3   = point 3 ordinate value
C            Y1   = point 1 abscissa value
C            Y2   = point 2 abscissa value
C            Y3   = point 3 abscissa value
C            XFIT = ordinate value to fit
C          RETURNS:
C            YFIT = abscissa value corresponding to XFIT
C
C          MODIFICATIONS:
C          Apr 2004: Change from subroutine to function, rename variables and
C                    add intermediates which are otherwise calculated twice
      implicit none
      REAL                FINT
C
C     Argument variables
C
      REAL                X1,          X2,          X3,          XFIT
      REAL                Y1,          Y2,          Y3
C
C     Local variables
C
      REAL                X12,         X13,         X23,         XF1
      REAL                XF2,         XF3
C

      X12 = X1-X2
      X13 = X1-X3
      X23 = X2-X3
      XF1 = XFIT-X1
      XF2 = XFIT-X2
      XF3 = XFIT-X3

      FINT = (Y1*X23*XF2*XF3 - Y2*X13*XF1*XF3 + Y3*X12*XF1*XF2) /
     +                                                     (X12*X13*X23)
      RETURN
      END
