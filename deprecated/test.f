      PROGRAM TEST
C          Sample program demonstrating calls to various magnetic-field and
C          magnetic-coordinate subroutines.  Questions about these routines
C          should be directed to Roy Barnes (bozo@ncar.ucar.edu 303/497-1230).
C          See file README for compilation instructions.

      implicit none
C
      INTEGER             INT
C
C     Parameter variables
C
C          Declarations needed for APXMKA, APXWRA, or APXRDA
C            MSGUN = Fortran unit number for diagnostics
C            IUN   = Fortran unit number for I/O
C            MLAT,MLON,MALT = Maximum number of grid latitudes, longitudes
C                    and altitudes.
C            NGRF  = Number of epochs in the current DGRF/IGRF; see COFRM in
C                    file magfld.f
      INTEGER             MSGUN
      PARAMETER           (MSGUN = 6)
      INTEGER             IUN
      PARAMETER           (IUN = 12)
      INTEGER             MLAT
      PARAMETER           (MLAT = 91)
      INTEGER             MLON
      PARAMETER           (MLON = 151)
      INTEGER             MALT
      PARAMETER           (MALT = 6)
      INTEGER             NGRF
      PARAMETER           (NGRF = 8)
      INTEGER             LWK
      PARAMETER           (LWK = MLAT*MLON*MALT*5+MLAT+MLON+MALT)
C
C     Local variables
C
C          Dimensions of non-scalar arguments returned by APXMALL:
      REAL B(3),BHAT(3),D1(3),D2(3),D3(3),E1(3),E2(3),E3(3),F1(2),F2(2)

C          FILNAM1 will contain gridded arrays for the epoch given by "DATE."
C          FILNAM2 will contain gridded arrays for all epochs from TGRF(1)
C          to TGRF(NGRF).
      CHARACTER*80 FILNAM1,FILNAM2
C
      INTEGER             I,           IDAY,        IHR,         IMN
      INTEGER             IST,         IYR,         NALT,        NDATE
      INTEGER             NLAT,        NLON,        NVERT
C
      REAL                A,           AI,          ALAT,        ALATI
      REAL                ALON,        ALONI,       ALONJ,       ALT
      REAL                ALTMN,       ALTMX,       BABS
      REAL                BDOWN,       BE3,         BEAST
      REAL                BMAG,        BNRTH,       COLAT,       CSZA
      REAL                D
      REAL                DATE
      REAL                ELON,        F
      REAL                FLDMAG,      GDLAT,       GDLON,       GLAMN
      REAL                GLAMX,       GLAT,        GLOMN,       GLOMX
      REAL                GLON,        GPALT(MALT), GPLAT(MLAT)
      REAL                GPLON(MLON), HR,          QDLAT,       QDLON
      REAL                SBSLLAT,     SBSLLON,     SEC,         SI
      REAL                SIM,         SZA,         TGRF(NGRF),  VA
      REAL                VMP,         VP,          W,           WK(LWK)
      REAL                XLATM,       XLATQD,      XLONM,       XMAG
      REAL                XMLT,        YMAG,        ZMAG
C
      DATA TGRF /1965.,1970.,1975.,1980.,1985.,1990.,1995.,2000./
C          Sample date (corresponding to 1994 July 2)
      DATE  = 1994.50
      NDATE = 1

C          Name both sample write files to be created during tests (
      FILNAM1 = 'Arecibo_and_Millstone_Hill-1994.5'
      WRITE (FILNAM2,'(''Arecibo_and_Millstone_Hill-'',I4,''.'',I4)')
     +                     INT (TGRF(1)) , INT (TGRF(NGRF))

C--------- TEST 1:  HIGH GRID RESOLUTION OVER ARECIBO
      WRITE (6,'(1X,/,''TEST 1:  Compare direct APEX model output with h
     +igh spatial resolution '',/,''         interpolation for one time
     +(DATE='',F8.2'')'')') DATE

C          Set up IGRF coefficients for the date desired
      CALL COFRM (DATE)

C          Geographic (geodetic) coordinates of a point 300 km above Arecibo
      GLAT = 18.3
      GLON = -66.75
      ALT  = 300.
      WRITE (6,'(1X,/,''Choose location above Arecibo: GLAT,GLON,ALT='',
     +3F8.2,'' deg, deg, km'')') GLAT,GLON,ALT

C          Compute magnetic-field components for location:
      CALL FELDG (1,GLAT,GLON,ALT, BNRTH,BEAST,BDOWN,FLDMAG)
      WRITE (6,'(1X,/,''B components (N, E, down, total) from FELDG:'',
     +4F9.5,'' Gauss'')')        BNRTH,BEAST,BDOWN,FLDMAG

C          Compute apex coordinates and magnetic field for this location
C          (note that APEX also calls FELDG to get field components):
      CALL APEX (DATE,GLAT,GLON,ALT, A,ALAT,ALON,BMAG,XMAG,YMAG,ZMAG,VA)
      WRITE (6,'(''B components (N, E, down, total) from APEX:  '',4F9.0
     +,'' nT'')')           XMAG,YMAG,ZMAG,BMAG

C          Set up arrays needed for interpolation routines, first for a box
C          around a single point with NVERT = 600 to get high accuracy
      NVERT = 600
      GLAMN = GLAT
      GLAMX = GLAT
      GLOMN = GLON
      GLOMX = GLON
      ALTMN = ALT
      ALTMX = ALT
      CALL GGRID (NVERT,GLAMN,GLAMX,GLOMN,GLOMX,ALTMN,ALTMX,
     +            GPLAT,GPLON,GPALT,NLAT,NLON,NALT)
      WRITE (6,'(1X,/,''GGRID produced for NVERT='',I3,'' and limits: la
     +titude '',F9.3,'' - '',F9.3,/,40X,''longitude '',F9.3,'' - '',F9.3
     +,/,40X,'' altitude '',F9.3,'' - '',F9.3)') NVERT,GLAMN,GLAMX,
     +                                     GLOMN,GLOMX,ALTMN,ALTMX
      WRITE (6,'(9X,''NLAT ='',I4,''  Latitudes:'',5F9.3,/,(31X,5F9.3))'
     +      )    NLAT,(GPLAT(I),I=1,NLAT)
      WRITE (6,'(9X,''NLON ='',I4,'' Longitudes:'',5F9.3,/,(31X,5F9.3))'
     +      )    NLON,(GPLON(I),I=1,NLON)
      WRITE (6,'(9X,''NALT ='',I4,''  Altitudes:'',5F9.3,/,(31X,5F9.3))'
     +      )    NALT,(GPALT(I),I=1,NALT)

C          Initialize interpolation arrays, but do not write them
      CALL APXMKA (MSGUN, DATE, GPLAT,GPLON,GPALT,NLAT,NLON,NALT,
     +            WK,LWK, IST)
      IF (IST .NE. 0) STOP
      WRITE (6,'(''A, ALAT, ALON, V from APEX: '',4F7.2,'' Re, deg, deg,
     + T-m'')') A,ALAT,ALON, VA

C          Use interpolation routine to get APEX coordinates:
      CALL APXALL (GLAT,GLON,ALT, WK, AI,ALATI,ALONI, IST)
      IF (IST .NE. 0) STOP
      WRITE (6,'(''A, ALAT, ALON from APXALL:  '',4F7.2)')AI,ALATI,ALONI

C          Reference height for modified apex coordinates and integer version
      HR  = 110.
      IHR = INT (HR + .5)

C          Get modified-apex and quasi-dipole coordinates and associated
C          parameters:
      CALL APXMALL (GLAT,GLON,ALT,HR, WK,
     +             B,BHAT,BABS,SI,XLONM,
     +             XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
     +             XLATQD,F,F1,F2 , IST)
      WRITE (6,'(''M('',I3,'') lat, lon, v from APXMALL:   '',3F7.2)')
     +                   IHR,                            XLATM,XLONM,VMP
      WRITE (6,'(''qsi dpl lat, lon  from APXMALL:    '',3F7.2)')
     +                                                      XLATQD,XLONM

C          Rewrite magnetic-field components from APEX to compare with the
C          interpolated values to be obtained from interpolation (APXMALL):
      WRITE (6,'(''B comps (N, E, down, total), sinI from APEX:    '',4F
     +9.0,F6.3)')                          XMAG,YMAG,ZMAG,BMAG,ZMAG/BMAG
      WRITE (6,'(''B comps (N, E, down, total), sinI from APXMALL: '',4F
     +9.0,F6.3)')                                B(2),B(1),-B(3),BABS,SI
      WRITE (6,'(''Be3, sinI_m from APXMALL:
     +                      '',F9.0,F6.3)') BE3,SIM
      WRITE (6,'(''W, D, F from APXMALL:'',E11.4,2F7.4)') W,D,F
      WRITE (6,'(''    East    North    Up'',/,
     +           ''d1'',3F8.4)') D1(1),D1(2),D1(3)
      WRITE (6,'(''e1'',3F8.4)') E1(1),E1(2),E1(3)
      WRITE (6,'(''f1'',3F8.4)') F1(1),F1(2)
      WRITE (6,'(''d2'',3F8.4)') D2(1),D2(2),D2(3)
      WRITE (6,'(''e2'',3F8.4)') E2(1),E2(2),E2(3)
      WRITE (6,'(''f2'',3F8.4)') F2(1),F2(2)
      WRITE (6,'(''d3'',3F8.4)') D3(1),D3(2),D3(3)
      WRITE (6,'(''e3'',3F8.4)') E3(1),E3(2),E3(3)
      WRITE (6,'(''b^'',3F8.4)') BHAT(1),BHAT(2),BHAT(3)

C--------- TEST 2:  CREATE AND SAVE LOW RESOLUTION GRID OVER ARO & MLH
C          Now set up arrays over the volume (lat = 15 to 45; lon = -75 to
C          -60; ht = 0 to 1000), using NVERT=30, and store them in FILNAM1
C          for one time (DATE):
      NVERT = 30
      WRITE (6,'(1X,/,''TEST 2:  Make comparison with low resolution int
     +erpolation for region over'',/,''         Arecibo and Millstone Hi
     +ll and write interpolation arrays'')')

      GLAMN = 15.
      GLAMX = 45.
      GLOMN = -75.
      GLOMX = -60.
      ALTMN = 0.
      ALTMX = 1000.
      CALL GGRID (NVERT,GLAMN,GLAMX,GLOMN,GLOMX,ALTMN,ALTMX,
     +            GPLAT,GPLON,GPALT,NLAT,NLON,NALT)

      WRITE (6,'(1X,/,''GGRID produced for NVERT='',I3,'' and limits: la
     +titude '',F8.3,'' - '',F8.3,/,40X,''longitude '',F8.3,'' - '',F8.3
     +,/,40X,'' altitude '',F8.3,'' - '',F8.3)') NVERT,GLAMN,GLAMX,
     +                                     GLOMN,GLOMX,ALTMN,ALTMX
      WRITE (6,'(9X,''NLAT ='',I4,''  Latitudes:'',5F8.3,/,(31X,5F8.3))'
     +      )    NLAT,(GPLAT(I),I=1,NLAT)
      WRITE (6,'(9X,''NLON ='',I4,'' Longitudes:'',5F8.3,/,(31X,5F8.3))'
     +      )    NLON,(GPLON(I),I=1,NLON)
      WRITE (6,'(9X,''NALT ='',I4,''  Altitudes:'',5F8.3,/,(31X,5F8.3))'
     +      )    NALT,(GPALT(I),I=1,NALT)

      CALL APXWRA (MSGUN, FILNAM1,IUN, DATE,NDATE,
     +            GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,LWK, IST)
      IF (IST .NE. 0) STOP

      WRITE (6,'(1X,/,''Location remains above Arecibo: GLAT,GLON,ALT=''
     +,3F8.2)') GLAT,GLON,ALT
C          Rewrite true apex coordinates for comparison to interpolated values
      WRITE(6,'(1X,/,''A, ALAT, ALON, V from APEX: '',4F7.2)')
     +                                              A,ALAT,ALON,VA

C          Use interpolation routine to get apex coordinates:
      CALL APXALL (GLAT,GLON,ALT, WK, AI,ALATI,ALONI, IST)
      WRITE(6,'(''A, ALAT, ALON from APXALL:  '',4F7.2)')AI,ALATI,ALONI

C          Get modified-apex and quasi-dipole coordinates and associated
C          parameters:
      CALL APXMALL (GLAT,GLON,ALT,HR, WK,
     +             B,BHAT,BABS,SI,XLONM,
     +             XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
     +             XLATQD,F,F1,F2 , IST)
      WRITE (6,'(''M('',I3,'') lat, lon, v from APXMALL:   '',3F7.2)')
     +                   IHR,                            XLATM,XLONM,VMP
      WRITE (6,'(''qsi dpl lat, lon  from APXMALL:    '',3F7.2)')
     +                                                      XLATQD,XLONM

C          Reprint magnetic-field components from APEX to compare with the
C          interpolated values to be obtained from interpolation (APXMALL):
      WRITE (6,'(''B comps (N, E, down, total), sinI from APEX:   '',4F9
     +.0,F6.3)')                           XMAG,YMAG,ZMAG,BMAG,ZMAG/BMAG
      WRITE (6,'(''B comps (N, E, down, total), sinI from APXMALL:'',4F9
     +.0,F6.3)')                           B(2),B(1),-B(3),BABS,SI

C--------- TEST 3:  READ BACK LOW RESOLUTION GRID OVER ARO & MLH
C          Test reacquisition of stored arrays in FILNAM1
      WRITE (6,'(1X,/,''TEST 3:  Read back arrays stored in '',A,/
     +,1X)') FILNAM1(1:33)
      CALL APXRDA (MSGUN, FILNAM1,IUN, DATE, WK,LWK, IST)
      IF (IST .NE. 0) STOP

C          Geographic (geodetic) coordinates of a point 1000 km above
C          Millstone Hill:
      GLAT = 42.6
      GLON = -71.5
      ALT  = 1000.
      WRITE (6,'(1X,/,''Location is now above Millstone: GLAT,GLON,ALT='
     +',3F8.2)') GLAT,GLON,ALT

C          Get "exact" apex coordinates and magnetic field for new location
      CALL APEX (DATE,GLAT,GLON,ALT, A,ALAT,ALON,BMAG,XMAG,YMAG,ZMAG,VA)
      WRITE (6,'(1X,/,''A, ALAT, ALON, V from APEX: '',3F7.2,F8.2)')
     +                  A, ALAT, ALON, VA

C          Use interpolation routine to get apex coordinates:
      CALL APXALL (GLAT,GLON,ALT, WK, AI,ALATI,ALONI, IST)
      WRITE (6,'(''A, ALAT, ALON from APXALL:  '',4F7.2)')
     +             AI,ALATI,ALONI

C          Compute "exact" magnetic-field components for location:
      CALL FELDG (1,GLAT,GLON,ALT, BNRTH,BEAST,BDOWN,FLDMAG)
      WRITE (6,'(1X,/,''B components (N, E, down, total) from FELDG:  ''
     +,4F9.5)')                                 BNRTH,BEAST,BDOWN,FLDMAG

C          Compute and write interpolated magnetic-field components
      CALL APXMALL (GLAT,GLON,ALT,HR,WK,
     +             B,BHAT,BABS,SI, XLONM,
     +             XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
     +             XLATQD,F,F1,F2, IST)
      WRITE (6,'(''B comps (N, E, down, total), sinI from APXMALL:'',4F9
     +.0,F6.3)')                                 B(2),B(1),-B(3),BABS,SI

C--------- TEST 4:  Use APXQ2G to determine geographic coordinates from
C          given quasi-dipole coordinates of a point that lies in the
C          volume spanned by the arrays:
      WRITE (6,'(1X,/,''Test 4:  Use APXQ2G to get geographic coordinate
     +s given quasi-dipole coordinates'',/,
     +''         then reverse the calculation by calling APXMALL'')')
      QDLAT = 45.
      QDLON = 2.
      ALT   = 110.
      WRITE (6,'(1X,/,''Choose quasi-dipole coordinates:           QDLAT
     +,QDLON,ALT='',3F7.2)') QDLAT,QDLON,ALT
      CALL APXQ2G (QDLAT,QDLON,ALT, WK, GDLAT,GDLON, IST)
      WRITE (6,'(''APXQ2G returned geographic coordinates:    GDLAT,GDLO
     +N    ='',2F7.2)') GDLAT,GDLON

C          Use the calculated geographic coordinates to recompute the
C          quasi-dipole coordinates, for comparison with the original input:
      CALL APXMALL (GDLAT,GDLON,ALT,HR, WK,
     +              B,BHAT,BABS,SI, QDLON,
     +              XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
     +              QDLAT,F,F1,F2, IST)
      WRITE (6,'(''Input GDLAT,GDLON,ALT and APXMALL returns: QDLAT,QDLO
     +N,ALT='',3F7.2)') QDLAT,QDLON,ALT

C--------- TEST 5:  CREATE AND WRITE LOW RESOLUTION GRIDS OVER ARO & MLH
C          FOR MULTIPLE TIMES, THEN READ BACK AND INTERPOLATE TO DATE.
C          Set up arrays for the DGRF/IGRF periods and store in FILNAM2
      WRITE (6,'(1X,/,''TEST 5: Interpolate in time too by first setting
     + up arrays for dates: '',/,(7X,5F9.3))') TGRF

      CALL APXWRA (MSGUN, FILNAM2,IUN, TGRF,NGRF,
     +            GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,LWK, IST)
      IF (IST .NE. 0) STOP
      WRITE(6,'(1X,/,''Stored arrays in file "'',A,''"'')')FILNAM2(1:36)

C          Acquire arrays stored in FILNAM2 and interpolate in time to DATE
      WRITE (6,'(''Retrieving arrays, interpolating in time to'',F9.3,
     +'' (call APXRDA)'')') DATE
      CALL APXRDA (MSGUN, FILNAM2,IUN, DATE, WK,LWK, IST)
      IF (IST .NE. 0) STOP

C          Geographic (geodetic) coordinates of a point 1000 km above
C          Millstone Hill:
      GLAT = 42.6
      GLON = -71.5
      ALT  = 1000.
      WRITE (6,'(''Location is again above Millstone: GLAT,GLON,ALT'',
     +           3F8.2)') GLAT,GLON,ALT

C          Use interpolation routine to get apex coordinates:
      CALL APXALL (GLAT,GLON,ALT, WK, AI,ALATI,ALONI, IST)
      WRITE (6,'(''A, ALAT, ALON from APXALL:  '',4F7.2)')AI,ALATI,ALONI

C--------- TEST 6:  TEST RELATED ROUTINES SUBSOL, MAGLOCTM, MLT2ALON
C          Ensure that IGRF coefficients are set up for the correct date
C          (since COFRM may have been called by APXWRA for a different DATE).
      WRITE (6,'(1X,/,''TEST 6: Use related routines SUBSOL, MAGLOCTM an
     +d MLT2ALON to obtain the'',/,''        subsolar point, use it to g
     +et magnetic local time and solar zenith'',/,''        angle, then
     +convert magnetic local time back to APEX longitude'')')

      CALL COFRM (DATE)
      CALL DYPOL (COLAT,ELON,VP)
      WRITE (6,'(1X,/,''First call COFRM and reestablish DATE ='',F10.3,
     +'' (The previous APXWRA call'',/,''set the date of the IGRF coeffi
     +cients at'',F9.3,'').  MAGLOCTM and MLT2ALON inputs'',/,''include
     +the magnetic dipole location, which is obtained calling DYPOL:'',/
     +,''COLAT,ELON,VP='',3F9.4,'' deg, deg, T-m'')') DATE,TGRF(NGRF),
     +   COLAT,ELON,VP

C          Find subsolar latitude, longitude for a specific time:
      IYR  = 1994
      IDAY = 182
      IHR  = 12
      IMN  = 0
      SEC  = 0.
      WRITE (6,'(1X,/,''SUBSOL inputs: IYR,IDAY,IHR,IMN,SEC ='',4I5,F5.0
     +,'' UT'')')                      IYR,IDAY,IHR,IMN,SEC
      CALL SUBSOL (IYR,IDAY,IHR,IMN,SEC,SBSLLAT,SBSLLON)
      WRITE (6,'(''SUBSOL returns: SBSLLAT,SBSLLON ='',2F8.2,'' deg'')')
     +                             SBSLLAT,SBSLLON

      CALL COSSZA (GLAT,GLON,SBSLLAT,SBSLLON,CSZA)
      SZA = ACOS (CSZA)*180./3.1415926535898
      WRITE (6,'(1X,/,''Using this subsolar point and the location above
     + Millstone Hill,'',/,''COSSZA returns solar zenith angle ='',F8.2,
     +'' degrees'')') SZA

      CALL MAGLOCTM (ALONI,SBSLLAT,SBSLLON,COLAT,ELON,XMLT)
      WRITE (6,'(1X,/,''Using this subsolar point, MAGLOCTM returns:  XM
     +LT ='',F7.3,'' hours'')') XMLT

      CALL MLT2ALON (XMLT,SBSLLAT,SBSLLON,COLAT,ELON,ALONJ)
      WRITE (6,'(''Using this subsolar point and XMLT, MLT2ALON returns'
     +',F7.2,'' deg'',/,''In comparision, APXALL originally returned:  A
     +LON = '',F7.2,'' deg'')') ALONJ,ALONI
      
C--------- TEST 7:  TEST OUT OF BOUNDS GRID INPUTS
C          Test bound-checking by deliberately feeding in a point outside
C          the bounds of the precalculated array:
      WRITE (6,'(1X,/,''Test 7:  Demonstrate response to inputting locat
     +ion outside array boundaries'')')

      GLAT = -42.6
      GLON = -71.5
      ALT  = 1000.
      WRITE (6,'(1X,/,''Location is in southern hemisphere: GLAT,GLON,AL
     +T'',3F8.2)') GLAT,GLON,ALT
      CALL APXALL (GLAT,GLON,ALT, WK, A,ALAT,ALON, IST)
      WRITE (6,'(''APXALL returns A,ALAT,ALON,IST='',3F10.2,I4)')
     +                            A,ALAT,ALON,IST

      END

