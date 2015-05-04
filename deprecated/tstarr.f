      PROGRAM TSTARR
C          Sample program demonstrating use of Apex interpolation routines
C          using precalculated gridded arrays from FILNAM. 
C          Requires that you first execute mkglob!
C
C          Dimensions of non-scalar arguments returned by APXMALL:
      DIMENSION B(3),BHAT(3),
     +          D1(3),D2(3),D3(3), E1(3),E2(3),E3(3), F1(2),F2(2)

C          FILNAM contains the lookup table. 
      CHARACTER*80 FILNAM
      PARAMETER (  FILNAM = 'apxntrp_grid2011')

C          Declarations needed for APXMKA, APXWRA, or APXRDA
C            MSGUN = Fortran unit number for diagnostics
C            IUN   = Fortran unit number for I/O
C            MLAT,MLON,MALT = Maximum number of grid latitudes, longitudes
C                    and altitudes.
      PARAMETER (MSGUN=6, MLAT=121,MLON=201,MALT=7, IUN=12,
     +           LWK= MLAT*MLON*MALT*5 + MLAT+MLON+MALT)
      DIMENSION GPLAT(MLAT),GPLON(MLON),GPALT(MALT),WK(LWK)

C--------- TEST 1:  Load into memory global tables for DATE
      WRITE (6,'(1X,/,''Test 1:  Read previously created tables, loading
     + into memory array '',/,''for 2011 January 1, and obtain the grid 
     +coordinates.'')')

      DATE  = 2011.0
      NDATE = 1
      NLAT  = MLAT
      NLON  = MLON
      NALT  = MALT
      WRITE (6,'(1X,/,''Retrieving arrays (with APXRDA) from file '',A)'
     +)                   FILNAM
      CALL APXRDA (MSGUN, FILNAM,IUN, DATE, WK,LWK, IST)
      IF (IST .NE. 0) STOP

      WRITE (6,'(1X,/,''Retrieving grid coordinates (with APXGGC)'')')
      CALL APXGGC (MSGUN, WK,LWK, GPLAT,GPLON,GPALT,NLAT,NLON,NALT,IST)
      WRITE (6,'(9X,''NLAT ='',I4,''  Latitudes:'',5F8.3,/,(31X,5F8.3))'
     +)    NLAT,(GPLAT(I),I=1,NLAT)
      WRITE (6,'(9X,''NLON ='',I4,'' Longitudes:'',5F8.3,/,(31X,5F8.3))'
     +)    NLON,(GPLON(I),I=1,NLON)
      WRITE (6,'(9X,''NALT ='',I4,''  Altitudes:'',5F8.3,/,(31X,5F8.3))'
     +)    NALT,(GPALT(I),I=1,NALT)

C--------- TEST 2:  Obtain Apex, quasi-dipole and modified Apex
C   coordinates for 
      GDLAT = 42.6
      GDLON = -71.5
      ALT = 250.
C HR = reference height for Modified Apex Coordinates, in km.
      HR = 0.
      WRITE (6,'(1X,/,''Obtain Apex, quasi-dipole and modified Apex coor
     +dinates for location '',/,''above Millstone Hill:  GDLAT,GDLON,ALT
     + = '',3F8.2)') GDLAT,GDLON,ALT
      CALL APXMALL (GDLAT,GDLON,ALT,HR, WK,
     +              B,BHAT,BABS,SI, QDLON,
     +              XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
     +              QDLAT,F,F1,F2, IST)
      WRITE (6,'(''Input GDLAT,GDLON,ALT and APXMALL returns:  '',/,
     +''QDLAT,QDLON,ALT='',3F7.2)') QDLAT,QDLON,ALT

C--------- TEST 3:  Use APXQ2G to get geographic coordinates of original
C   point and conjugate point.
      WRITE (6,'(1X,/,''Use APXQ2G to get geographic coordinates of orig
     +inal point.'')')
      CALL APXQ2G (QDLAT,QDLON,ALT, WK, GDLAT,GDLON, IST)
      WRITE (6,'(''APXQ2G returns geographic coordinates: GDLAT,GDLON,AL
     +T = '',3F7.2)') GDLAT,GDLON,ALT
      WRITE (6,'(1X,/,''Use APXQ2G to get geographic coordinates of conj
     +ugate point.'')')
      CALL APXQ2G (-QDLAT,QDLON,ALT, WK, GDLAT,GDLON, IST)
      WRITE (6,'(''APXQ2G returned conjugate geographic coordinates:'',/
     +,''GDLAT,GDLON,ALT ='',3F7.2)') GDLAT,GDLON,ALT

      END
