      SUBROUTINE GGRID (NVERT,GLAMN,GLAMX,GLOMN,GLOMX,ALTMN,ALTMX,
     +                 GPLAT,GPLON,GPALT,NLAT,NLON,NALT)

C          Define grid points suitable for preparing interpolation tables for
C          Apex coordinates, Modified Apex Coordinates and Quasi-Dipole
C          coordinates.
C
C          INPUTS:
C            NVERT = spatial resolution parameter must be at least 2 but 30 or
C                    or more is recommended.  Interpolation accuracy increases
C                    with larger NVERT until about 100 when points near the
C                    poles become too close for accurate east-west gradient
C                    determination affecting B(1), G, H, W, Bhat(1), D1(1),
C                    D2(1), D3, E1, E2, E3, F1(2), and F2(2)); this effects all
C                    latitudes when NVERT is 1000 or more.
C            GLAMN = minimum latitude (degrees) must be in the range -90 to 90.
C            GLAMX = maximum latitude (degrees) must be in the range GLAMN to 90.
C            GLOMN = minimum longitude (degrees) must be in the range -180 to 180.
C            GLOMX = maximum longitude (degrees) must be in the range GLOMN to
C                    the smaller of 540 or GLOMN+360.
C            ALTMN = minimum altitude (km) must be at least 0.
C            ALTMX = maximum altitude (km) must be at least ALTMN.
C          RETURNS:
C            GPLAT = Grid point latitudes (degrees north).
C            GPLON = Grid point longitudes (degrees east).
C            GPALT = Grid point altitudes (km above earth).
C            NLAT  = Number of assigned values in GPLAT.
C            NLON  = Number of assigned values in GPLON.
C            NALT  = Number of assigned values in GPALT.
C
C          This is the suggested method for defining grid points used to
C          prepare interpolation tables (by entry APXWRA of apxntrp.f).
C
C          When determining geomagnetic coordinates, the dipole component is
C          separated so it may be computed analytically (to improve accuracy)
C          and the non-dipole component is interpolated from the tables.  As
C          altitude increases, the non-dipole component becomes spatially
C          smoother and relatively less important, so vertical spacing can
C          increase with altitude; hence, the quantity
C
C            Re/(Re+height)
C
C          is divided into NVERT segments, where Re is the mean earth radius
C          (6371.2 km).  The latitudinal grid spacing is chosen to be
C          approximately the same at the Earth's surface as the vertical
C          distance between the two lowest grid points in altitude or
C          180/(3*NVERT).  The longitude spacing is chosen to be roughly
C          equivalent to low and middle latitude spacing or 360/(5*NVERT);
C          however longitude spacing decreases when approaching the poles.
C          Grid values returned are the minimal subset encompassing the
C          specified latitude, longitude and altitude extremes.
C
C          HISTORY:
C          Aug 1995: Adapted from the initial version of apxntrp.f written by
C                    A. D. Richmond.  This grid definition code was separated from
C                    apxntrp.f to allow arbitrary grid definition.  R. Barnes, NCAR.
C          May 2004: Clean up comments; change print diagnostics to unit 0 (stderr);
C                    alias GLOMX to GLONMX s.t. it can be changed w/o affecting the
C                    input argument.
C
      implicit none
C
      DOUBLE PRECISION    DBLE
C
      INTEGER             INT,         MAX0,        MIN0
C
      REAL                AMAX1,       AMIN1
C
C     Parameter variables
C
      DOUBLE PRECISION    DRE
      PARAMETER           (DRE = 6371.2D0)
C
      REAL                RE
      PARAMETER           (RE = 6371.2)
C
C     Argument variables
C
      INTEGER             NALT,        NLAT,        NLON,        NVERT
C
      REAL                ALTMN,       ALTMX,       GLAMN,       GLAMX
      REAL                GLOMN,       GLOMX,       GPALT(*)
      REAL                GPLAT(*),    GPLON(*)
C
C     Local variables
C
      DOUBLE PRECISION    DIHT,        DLAT,        DLON,        DNV
C
      INTEGER             I,           IHT,         IHTMN,       IHTMX
      INTEGER             J,           K,           LAMN,        LAMX
      INTEGER             LOMN,        LOMX
C
      REAL                GLOMXL,      X
C
C          Code from subroutine SETLIM in the previous version of apxntrp
C          which establishes indices defining the geographic grid (plus
C          a few modifications, such as some double precision variables).
      IF (GLAMX     .LT. GLAMN) GO TO 9100
      IF (GLOMX     .LT. GLOMN) GO TO 9200
      IF (ALTMX     .LT. ALTMN) GO TO 9300
      IF (ABS(GLAMN).GT.   90.) GO TO 9400
      IF (ABS(GLOMN).GT.  180.) GO TO 9500
      IF (ALTMN     .LT.    0.) GO TO 9600

      DNV  = DBLE (NVERT)
      DLON = 360.D0 / (5.D0*DNV)
      DLAT = 180.D0 / (3.D0*DNV)
      DIHT = 1.D0 / DNV

      LAMN =  MAX0 (INT((DBLE(GLAMN) +90.D0 )/DLAT)     , 0)
      LAMX =  MIN0 (INT((DBLE(GLAMX) +90.D0 )/DLAT+1.D0), 3*NVERT)

      LOMN  = MAX0 (INT((DBLE(GLOMN) +180.D0)/DLON)  , 0)
      GLOMXL= AMIN1 (GLOMX,GLOMN+360.)
      LOMX  = MIN0 (INT((DBLE(GLOMXL)+180.D0)/DLON+1.D0),10*NVERT)

      X = RE/(RE+ALTMX)/DIHT - 1.E-5
      IHTMN = AMAX1 (X,1.)
      IHTMN = MIN0 (IHTMN,NVERT-1)
      X = RE/(RE+ALTMN)/DIHT + 1.E-5
      I = X + 1.
      IHTMX = MIN0(I, NVERT)

      NLAT = LAMX - LAMN + 1
      NLON = LOMX - LOMN + 1
      NLON = MIN0 (NLON,5*NVERT+1)
      NALT = IHTMX - IHTMN + 1

C          Code from old MAKEXYZV which converts from indices to lat,lon,alt
      DO 110 I=1,NLAT
  110 GPLAT(I) = DLAT*DBLE(LAMN+I-1) - 90.D0
      DO 120 J=1,NLON
  120 GPLON(J) = DLON*DBLE(LOMN+J-1) - 180.D0
      DO 130 K=1,NALT
      IHT = IHTMX - K + 1
  130 GPALT(K) = DRE*(DBLE(NVERT-IHT) - 1.D-05) / (DBLE(IHT)+1.D-5)

C        Adjustments required to match test case results
      IF (GPLON(NLON-1) .GE. GLOMXL) NLON = NLON - 1
      GPALT(1) = AMAX1 (GPALT(1),0.)

      RETURN

 9100 WRITE (0,'(''GGRID:  GLAMX < GLAMN'')')
      CALL EXIT (1)
 9200 WRITE (0,'(''GGRID:  GLOMX < GLOMN'')')
      CALL EXIT (1)
 9300 WRITE (0,'(''GGRID:  ALTMX < ALTMN'')')
      CALL EXIT (1)
 9400 WRITE (0,'(''GGRID:  |GLAMN| > 90.'')')
      CALL EXIT (1)
 9500 WRITE (0,'(''GGRID:  |GLOMN| > 180.'')')
      CALL EXIT (1)
 9600 WRITE (0,'(''GGRID:  ALTMN < 0.'')')
      CALL EXIT (1)
      END
