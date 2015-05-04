      PROGRAM XGLOB
C          Example of reading an Apex interpolation tables file; usage:
C
C            xglob.exe apxntrp.file
C
C          where
C            xglob.exe    = executable created from this source
C            apxntrp.file = interpolation tables file created by APXWRA in apxntrp.f
C
C------------------------------------------------------------------------------
C          INSTALLATION:
C          (1) See apxntrp.f and mkglob.f to prepare apxntrp.file
C          (2) Make executable xglob.exe using the Makefile (make xglob) or
C              do it manually; e.g.
C
C                f77 -O -o ../bin/xglob.exe xglob.f apxntrp.f
C
C          (3) Optional: current implementation uses ../bin/envapex to define
C              environment variables, including $APXNTRF, the local full path
C              name of apxntrp.file and it uses ../bin/xglob, a script to drive
C              this executable.
C
C------------------------------------------------------------------------------
C          EXTERNALS:
C
C          apxntrp.f  - Apex interpolation subroutines
C          Makefile   - 'make' file for Unix systems which defines source code
C                       dependencies; i.e. instructions to make executables
C          xapxntrp.f - example program for interpolation routines which creates
C                       a file (FILNAM2) which can be read back here
C
C------------------------------------------------------------------------------
C          MODIFICATIONS:
C          May 2004: Change FILNAM assignment from hard-coded to a command
C                    argument; change name from tstarr.f to xglob.f

C          Dimensions of non-scalar arguments returned by APXMALL:
      DIMENSION B(3),BHAT(3),
     +          D1(3),D2(3),D3(3), E1(3),E2(3),E3(3), F1(2),F2(2)

C          FILNAM = Name of the previously created interpolation tables file
C          CALNM  = Name of current routine being called; used for error diagnostics
      CHARACTER FILNAM*200, CALNM*7

C          Declarations needed for APXMKA, APXWRA, or APXRDA
C            MSGUN = Fortran unit number for diagnostics
C            IUN   = Fortran unit number for I/O
C            MLAT  = Maximum number of grid latitudes (bigger than necessary)
C            MLON  = Maximum number of grid longitudes
C            MALT  = Maximum number of grid altitudes
      PARAMETER (MSGUN=6, MLAT=300,MLON=500,MALT=40, IUN=12,
     +           LWK= MLAT*MLON*MALT*5 + MLAT+MLON+MALT)
      DIMENSION GPLAT(MLAT),GPLON(MLON),GPALT(MALT),WK(LWK)

C          Get command argument (apxntrp file name)
      NARG = IARGC()
      IF (NARG .NE. 1) THEN
	WRITE (MSGUN,'(''xglob requires one command argument: apxntrp fi
     +le name'')')
	CALL EXIT (1)
      ENDIF
      CALL GETARG (1, FILNAM)


      WRITE (6,'(''Test 1:  Read previously created tables, loading into
     + memory tables interpolated'',/,''         in time to 2 Jul 1994 a
     +nd obtain the table grid coordinates'',/)')

      DATE  = 1994.5                                                     ! time UT, yyyy.fraction
      CALL APXRDA (MSGUN, FILNAM,IUN, DATE, WK,LWK, IST)                 ! read back tables, possibly interpolating in time
      IF (IST .LT. 0) WRITE (6,9100) IST
      IF (IST .GT. 0) THEN
	CALNM = 'APXRDA'
	GO TO 9200
      ENDIF
      CALL APXGGC (MSGUN, WK,LWK, GPLAT,GPLON,GPALT,NLAT,NLON,NALT,IST)  ! get interp. table grid coordinates
      IF (IST .GT. 0) THEN
	CALNM = 'APXGGC'
	GO TO 9200
      ENDIF

      WRITE (6,'(''         Retrieved tables by calling APXRDA with:'',/
     +,''         FILNAM = '',A,'' (read file name)'',/,''         DATE
     +  = '',F9.3,'' (time UT, yyyy.fraction)'')') FILNAM(:LENNB(FILNAM
     +)),DATE
      WRITE (6,'(''         APXGGC found: '',/,''         NLAT  = '',I9,
     +'' (no. latitudes)'',/,''         GPLAT = '',10F9.3,/,(17X,10F9.3)
     +))') NLAT,(GPLAT(I),I=1,NLAT)
      WRITE (6,'(''         NLON  = '',I9,  '' (no. longitudes)'',/,''
     +       GPLON = '',10F9.3,/,(17X,10F9.3))')NLON,(GPLON(I),I=1,NLON)
      WRITE (6,'(''         NALT  = '',I9,  '' (no. altitudes)'' ,/,''
     +       GPALT = '',10F9.3,/,(17X,10F9.3))')NALT,(GPALT(I),I=1,NALT)


      WRITE (6,'(/,''Test 2:  Repeat test 4 of xapxntrp.f demonstrating
     +geomagnetic-geographic coordinate'',/,''         conversion; note
     +Apex, Modified Apex and Quasi-Dipole share the same'',/,''
     + longitude coordinate (ALON) and altitude coordinate (ALT):'',/,52
     +X,''----------------------------- Results ------------------------
     +-----'',/,52X,''-------- Geographic -------   --- Apex ----  Mod.
     +Apex Lat.   Quasi-'',/,52X,''Latitude Longitude Altitude   Lat.
     +Long.    Hr=0  Hr=110 Dipole Lat.'',/,44X,''Routine  degrees  degr
     +ees    km.      deg.    deg.     deg.   deg.    deg.'',/,44X,''Cal
     +led    GLAT     GLON      ALT      ALAT    ALON    ALAT0  ALAT1
     +QDLAT'')')

      QDLAT = 45.                                                        ! Quasi-dipole latitude, degrees
      QDLON = 2.                                                         ! Quasi-dipole longitude, degrees
      ALT   = 110.                                                       ! Quasi-dipole altitude, km
      CALL APXQ2G (QDLAT,QDLON,ALT, WK, GLAT,GLON, IST)                  ! Quasi-dipole to geographic
      IF (IST .GT. 0) THEN
	CALNM = 'APXQ2G'
	GO TO 9200
      ENDIF
      CALL APXALL (GLAT ,GLON ,ALT, WK, A,ALAT,ALON, IST)                ! Geographic to Apex
      IF (IST .GT. 0) THEN
	CALNM = 'APXALL'
	GO TO 9200
      ENDIF
      CALL APXA2G (ALAT ,ALON ,ALT, WK, GLATA,GLONA, IST)                ! Apex to geographic
      IF (IST .GT. 0) THEN
	CALNM = 'APXA2G'
	GO TO 9200
      ENDIF
      HR0   = 0.                                                         ! Modified Apex reference altitude, km
      CALNM = 'APXMALL'
      CALL APXMALL (GLATA,GLONA,ALT,HR0, WK,                             ! Geographic to Modified Apex (Hr=0)
     +              B,BHAT,BTAM,SI, ALON0,
     +              ALAT0,VMP,WM,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
     +              QLAT0,F,F1,F2, IST)
      IF (IST .GT. 0) GO TO 9200
      HR1 = 110.                                                         ! Modified Apex reference altitude, km
      CALL APXMALL (GLATA,GLONA,ALT,HR1, WK,                             ! Geographic to Modified Apex (Hr=110 km)
     +              B,BHAT,BTAM,SI, ALON1,
     +              ALAT1,VMP,WM,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
     +              QLAT1,F,F1,F2, IST)
      IF (IST .GT. 0) GO TO 9200

      CALNM = 'APXM2G'
      CALL APXM2G (ALAT0,ALON0,ALT,HR0, WK, GLAT0,GLON0, IST)            ! Modified Apex (Hr=0)   to geographic
      IF (IST .GT. 0) GO TO 9200
      CALL APXM2G (ALAT1,ALON1,ALT,HR1, WK, GLAT1,GLON1, IST)            ! Modified Apex (HR=110) to geographic
      IF (IST .GT. 0) GO TO 9200

      WRITE (6,'(''     choose Quasi-Dipole (Q-D) coordinates:'',27X,F9.
     +3,9X,F8.3,17X,F7.3)') ALT,QDLON,QDLAT
      WRITE (6,'(''                         Q-D to geographic: APXQ2G ''
     +,2F9.3)') GLAT,GLON
      WRITE (6,'(''                        geographic to Apex: APXALL ''
     +,27X,2F9.3)') ALAT,ALON
      WRITE (6,'(''                        Apex to geographic: APXA2G ''
     +,2F9.3)') GLATA,GLONA
      WRITE (6,'(''geographic to Modified Apex (Hr=0)   & Q-D: APXMALL''
     +,37X,2F8.3,8X,F8.3)') ALON0,ALAT0,QLAT0
      WRITE (6,'(''geographic to Modified Apex (Hr=110) & Q-D: APXMALL''
     +,37X,F8.3,8X,2F8.3)') ALON1,ALAT1,QLAT1
      WRITE (6,'(''      Modified Apex (Hr=0)   to geographic: APXM2G ''
     +,2F9.3,/,''      Modified Apex (Hr=110) to geographic: APXM2G '',2
     +F9.3)') GLAT0,GLON0, GLAT1,GLON1
C                                                    ----------------------------- Results -----------------------------
C                                                    -------- Geographic -------   --- Apex ----  Mod. Apex Lat.   Quasi-
C                                                    Latitude Longitude Altitude   Lat.   Long.    Hr=0  Hr=110 Dipole Lat.
C                                            Routine  degrees  degrees    km.      deg.    deg.     deg.   deg.    deg.
C                                            Called    GLAT     GLON      ALT      ALAT    ALON    ALAT0  ALAT1   QDLAT
C     choose Quasi-Dipole (Q-D) coordinates:                           xxxxx.xxx         -xxx.xxx                 -xx.xxx
C                         Q-D to geographic: APXQ2G   -xx.xxx -xxx.xxx
C                        geographic to Apex: APXALL                              -xx.xxx -xxx.xxx
C                        Apex to geographic: APXA2G   -xx.xxx -xxx.xxx
Cgeographic to Modified Apex (Hr=0)   & Q-D: APXMALL                                     -xxx.xxx -xx.xxx         -xx.xxx
Cgeographic to Modified Apex (Hr=110) & Q-D: APXMALL                                     -xxx.xxx         -xx.xxx -xx.xxx
C      Modified Apex (Hr=0)   to geographic: APXM2G   -xx.xxx -xxx.xxx
C      Modified Apex (Hr=110) to geographic: APXM2G   -xx.xxx -xxx.xxx


      WRITE (6,'(/,''Test 3:  Demonstrate consistent read-back interpola
     +tion by repeating that part'',/,''         of xapxntrp.f test 5:''
     +,/)')

      GLAT  =  42.6                                                      ! geographic latitude  of Millstone Hill
      GLON  = -71.5                                                      ! geographic longitude of Millstone Hill
      ALT   = 1000.                                                      ! geographic altitude above ground
      CALNM = 'APXALL'
      CALL APXALL (GLAT,GLON,ALT, WK, AI,ALATI,ALONI, IST)               ! geographic to Apex via interpolation
      IF (IST .GT. 0) GO TO 9200

      WRITE (6,'(''         Time is already set; choose location over Mi
     +llstone Hill:'',/,''         GLAT  = '',F8.3,'' (geographic latitu
     +de, degrees)'',/,''         GLON  = '',F8.3,'' (geographic longitu
     +de, degrees)'',/,''         ALT   = '',F8.3,'' (geographic altitud
     +e, degrees)'')') GLAT,GLON,ALT
      WRITE (6,'(''         Apex coordinates from APXALL: A, ALAT,ALON''
     +,3F11.3,'' (compare this to xapxntrp output)'')') AI,ALATI,ALONI


      WRITE (6,'(/,''Test 4:  Repeat test 3 for dates before and after t
     +hose currently stored'',/,''         in the look-up tables file'')
     +')
      DATE = 1964.5                                                      ! time UT, yyyy.fraction
      WRITE (6,'(''         Choose early time: DATE = '',F9.3)') DATE
      CALL APXRDA (MSGUN, FILNAM,IUN, DATE, WK,LWK, IST)                 ! read back tables and interpolate in time
      IF (IST .LT. 0) WRITE (6,9100) IST
      IF (IST .GT. 0) THEN
	CALNM = 'APXRDA'
	GO TO 9200
      ENDIF
      CALL APXALL (GLAT,GLON,ALT, WK, AI,ALATI,ALONI, IST)               ! geographic to Apex via interpolation
      IF (IST .GT. 0) GO TO 9200
      WRITE (6,'(''         Apex coordinates from APXALL: A, ALAT,ALON''
     +,3F11.3,'' (note effect of time difference)'')') AI,ALATI,ALONI

      DATE = 2010.5                                                      ! time UT, yyyy.fraction
      WRITE (6,'(''         Choose late time:  DATE = '',F9.3)') DATE
      CALL APXRDA (MSGUN, FILNAM,IUN, DATE, WK,LWK, IST)                 ! read back tables and interpolate in time
      IF (IST .LT. 0) WRITE (6,9100) IST
      IF (IST .GT. 0) THEN
	CALNM = 'APXRDA'
	GO TO 9200
      ENDIF
      CALL APXALL (GLAT,GLON,ALT, WK, AI,ALATI,ALONI, IST)               ! geographic to Apex via interpolation
      IF (IST .GT. 0) GO TO 9200
      WRITE (6,'(''         Apex coordinates from APXALL: A, ALAT,ALON''
     +,3F11.3)') AI,ALATI,ALONI
      CALL EXIT (0)

C          Error trap diagnostics
 9100 FORMAT ('         APXRDA negative return status (IST=,',I3,') impl
     +ies extrapolated date')
 9200 WRITE (6,'(9X,A,'' positive return status (IST),'',I3,'', implies
     +fatal error)'')') CALNM(:LENNB(CALNM)),IST
      CALL EXIT (1)

      END

      FUNCTION LENNB (STR)
C          Return the position of the last non-blank in STR or 1 (STR is blank)

      CHARACTER*(*) STR

      LENNB = LEN (STR) + 1
   10 LENNB = LENNB - 1
      IF (STR(LENNB:LENNB) .EQ. ' ' .AND. LENNB .GT. 1) GO TO 10

      RETURN
      END
