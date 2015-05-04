      PROGRAM MKGLOB
C          Make a global Apex interpolation tables file; usage:
C
C            mkglob.exe nvert altmx apxntrp_file epoch [epoch2 ...]
C
C          where
C            mkglob.exe   = executable created from this source
C            nvert        = spatial resolution parameter must be at least 2 but 30
C                           or more is recommended.  Resolution increases as nvert
C                           increases until 100 where accuracy degrades close to the
C                           poles in east-west gradients (B(1), G, H, W, Bhat(1),
C                           D1(1), D2(1), D3, E1, E2, E3, F1(2), and F2(2)).  nvert
C                           greater than 1000 degrades accuracy at all latitudes.
C                           Comments in ggrid.f describe the algorithm.
C            altmx        = Maximum altitude (km) required in interpolation tables.
C            apxntrp_file = name of interpolation tables file to be created.
C            epoch        = date (year.fraction); must not be before the earliest
C                           IGRF (currently 1900).  Multiple epochs must be in
C                           chronological order.
C
C          A set of interpolation tables is written to 'apxntrp_file' for each
C          requested epoch.  Tables are set up for a grid created by subroutine
C          GGRID (file ggrid.f) based on nvert, altmx and hard-coded limits of
C          latitude (-90 to 90), longitude (-180 to 180) and minimum altitude
C          (0 km).  The following table lists dimensions and actual maximum
C          altitudes for two suggested nvert values and various altmx:
C
C                 nvert altmx   nlat nlon nalt gpalt(nalt)
C                   30   500.0    91  151    4   707.909
C                   30  1000.0    91  151    6  1274.237
C                   30  2000.0    91  151    9  2316.796
C                   30  3000.0    91  151   11  3185.595
C                   30  4000.0    91  151   13  4247.461
C                   30  5000.0    91  151   15  5574.792
C                   30 10000.0    91  151   20 11004.784
C                   30 20000.0    91  151   24 20933.904
C                   30 30000.0    91  151   26 31855.924
C                   30 40000.0    91  151   27 41412.680
C                   30 50000.0    91  151   28 57340.586
C                   30 60000.0    91  151   29 89196.320
C                   40   500.0   121  201    4   516.582
C                   40  1000.0   121  201    7  1124.327
C                   40  2000.0   121  201   11  2123.730
C                   40  3000.0   121  201   14  3067.611
C                   40  4000.0   121  201   17  4247.462
C                   40  5000.0   121  201   19  5212.795
C                   40 10000.0   121  201   26 10618.655
C                   40 20000.0   121  201   32 21945.213
C                   40 30000.0   121  201   34 30035.605
C                   40 40000.0   121  201   36 44598.297
C                   40 50000.0   121  201   37 57340.641
C                   40 60000.0   121  201   38 78577.852
C
C          The time to execute and resulting file size is proportional to the grid
C          dimensions and number of epochs; it has exceeded an hour:
C
C                   File        --------- Grid ---------                Float.Pt.
C             Time  Size   No.  Dimensions   -- GGRID --               Clock Rate
C             (sec) (MB) Epochs lat lon alt  nvert mxalt  - Hardware -    MHz    ------ OS ------ Hostname   Date
C              293   37    9    121 201  7     40  1000   Intel 686              RedHat Linux V.9  cedar-l  May 2004
C             1356   37    9    121 201  7     40  1000   Sun-blade-100   500    Sun Solaris 5.8   hermes   May 2004
C
C              354   37    9    121 201  7     40  1000   Intel 686              RedHat Linux V.8  cedar-l  Feb 2003
C
C              944   37    9    121 201  7     40  1000   SGI Origin2000  250*   IRIX64 v6.5       dataproc May 2001
C                                                                            *   2 flops/cycle
C             1218   37    9    121 201  7     40  1000   Sun Ultra       400    Sun Solaris 5.7   huron    May 2000
C             1740   33    8    121 201  7     40  1000   Sun 4           248
C             4500   16    8     91 121  7     30  1000   Sun 4            70?
C            16943   74    9    121 201  7     40  1000   Cray J924se            Unicos            ouray    May 2000
C            12700   66    8    121 201  7     40  1000   Cray J924se            Unicos            ouray    May 1999
C             6300   32    8     91 121  7     30  1000   Cray J924se            Unicos            ouray    1998
C
C          Please direct questions to Roy Barnes (bozo@ucar.edu or 303-497-1572).
C
C          INSTALLATION:
C          (1) Make executable mkglob.exe using the Makefile (make mkglob) or
C              do it manually; e.g.,
C
C                f77 -O -o ../bin/mkglob.exe mkglob.f ggrid.f apxntrp.f
C
C          (2) Optional: current implementation uses ../bin/envapex to define
C              environment variables including $APXLIB, the Apex object library
C              and $APXNTRF, the global interpolation tables file.  Also optional
C              is the perl script driver for this program ../bin/mkglob which
C              prepares command line arguments and launches the executable.
C
C          HISTORY
C          Sep 1995: Initially written
C          May 2004: Change hard-coded assignment of file name, nvert, altmax
C                    and epochs to command arguments.  Also change print format
C                    to show three decimal places.

C          FILNAM is the name of the file to be created/overwritten
      CHARACTER*200 FILNAM

C          Declarations needed for APXMKA, APXWRA, or APXRDA
C            MSGUN = Fortran unit number for diagnostics
C            IUN   = Fortran unit number for writing FILNAM
C            MLAT  = Maximum number of grid latitudes
C            MLON  = Maximum number of grid longitudes
C            MALT  = Maximum number of grid altitudes
C            MEPO  = Maximum number of epochs
C            LWK   = Dimension of work array (IWK) must be at least
C                    NLAT*NLON*NALT*5 + NLAT+NLON+NALT
C            GLAMN,GLAMX = Min. and max. required latitudes
C            GLOMN,GLOMX = Min. and max. required longitudes
C            ALTMN,ALTMX = Min. and max. required altitudes
C            NEPO  = Number of epochs

      PARAMETER (MSGUN=0, IUN=12, MLAT=300,MLON=500,MALT=100, MEPO=12,
     +           LWK= MLAT*MLON*MALT*5 + MLAT+MLON+MALT,
     +           GLAMN = -90. , GLOMN = -180. ,
     +           GLAMX =  90. , GLOMX =  180. , ALTMN = 0.)

      DIMENSION GPLAT(MLAT),GPLON(MLON),GPALT(MALT),EPOC(MEPO),WK(LWK)
      CHARACTER STR20*20

C          Get command arguments (apxntrp file name)
      NARG = IARGC()
      IF (NARG .LT. 4) GO TO 9100

      CALL GETARG (1, STR20)
      READ (STR20,'(I20)',ERR=9200) NVERT

      CALL GETARG (2, STR20)
      READ (STR20,'(F20.0)',ERR=9300) ALTMX

      CALL GETARG (3, FILNAM)

      NEPO = NARG-3
      IF (NEPO .GT. MEPO) GO TO 9400
      J = 0
      DO 10 I=4,NARG
      J = J + 1
      CALL GETARG (I, STR20)
   10 READ (STR20,'(F20.0)',ERR=9500) EPOC(J)

      CALL GGRID (NVERT,GLAMN,GLAMX,GLOMN,GLOMX,ALTMN,ALTMX,
     +            GPLAT,GPLON,GPALT,NLAT,NLON,NALT)
      WRITE (6,'(1X,/,''GGRID produced for NVERT='',I3,'' and limits: la
     +titude '',F8.3,'' - '',F8.3,/,40X,''longitude '',F8.3,'' - '',F8.3
     +,/,40X,'' altitude '',F8.3,'' - '',F8.3)') NVERT,GLAMN,GLAMX,
     +                                     GLOMN,GLOMX,ALTMN,ALTMX
      WRITE (6,'(9X,''NLAT ='',I4,''  Latitudes:'',10F9.3,/,(31X,10F9.3)
     +)') NLAT,(GPLAT(I),I=1,NLAT)
      WRITE (6,'(9X,''NLON ='',I4,'' Longitudes:'',10F9.3,/,(31X,10F9.3)
     +)') NLON,(GPLON(I),I=1,NLON)
      WRITE (6,'(9X,''NALT ='',I4,''  Altitudes:'',10F9.1,/,(31X,10F9.1)
     +)') NALT,(GPALT(I),I=1,NALT)

      CALL APXWRA (MSGUN,FILNAM,IUN, EPOC,NEPO,
     +            GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,LWK, IST)
      CALL EXIT (IST)


 9100 WRITE (MSGUN,'(''mkglob.exe usage: mkglob.exe nvert altmx apxntrp_
     +file epoch [epoch2 ...]'')')
      CALL EXIT (1)
 9200 WRITE (MSGUN,'(''mkglob.exe: trouble parsing first command argumen
     +t (nvert)'')')
      CALL EXIT (1)
 9300 WRITE (MSGUN,'(''mkglob.exe: trouble parsing second command argume
     +nt (altmx)'')')
      CALL EXIT (1)
 9400 WRITE (MSGUN,'(''mkglob.exe: increase dimension (MEPO) to at least
     +'',I5)') NEPO
      CALL EXIT (1)
 9500 WRITE (MSGUN,'(''mkglob.exe: expecting date (yyyy.fraction), got '
     +',A)') STR20
      CALL EXIT (1)
      END
