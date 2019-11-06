!  Sane wrappers for apxntrp_legacy.f subroutines. 
!
!  See apxntrp_legacy.f comments for information on subroutine arguments.
!
!  apxntrp_legacy.f is legacy code that consistenly follows worst practices
!  in software engineering (ie. GOTO statements, common blocks, entry
!  statements, etc.).
!
!  Author: Peter Schmitt
!  Date: August 23, 2012
!

SUBROUTINE APXMKA(MSGUN, EPOCH, GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,IST)
!          Create tables which are used later by APXALL, APXMALL, APXQ2G,
!          APXA2G or APXM2G.  The tables are computed for a grid of geographic
!          locations (latitude, longitude and altitude) which are inputs; i.e.
!          one must define a suitable grid, a compromise between global extent
!          and resolution.  The smallest possible grid (2,2,2) with small
!          spatial increments bounding the point of interest will insure high
!          interpolation accuracy at the risk of round-off errors in vector
!          quantities which are determined from differences between adjacent
!          grid point values.  Such degraded accuracy has been observed in
!          east-west gradients very close to the poles using a 301x501x15 grid
!          which has increments 0.6x0.72 (deg. lat. by lon.).  This introduces
!          the other extreme, a global array extending to several Earth radii
!          which can produce huge tables.  Examples are given separately in
!          program xapxntrp.f and the global table generator program mkglob.f;
!          for global coverage up to at least 1000 km (ALTMX=1000) we have used
!          two grids, whose resolution is specified by NVERT:
!
!                     Dimension    file size MB   Number of
!             NVERT (latxlonxalt) (32-bit words) Times (epochs)
!               30     91,121,7      15.8            8
!               40    121,201,7      32.7            8
!               40    121,201,7      36.8            9

  implicit none
!            MSGUN  = Fortran unit number to write diagnostic messages.
  INTEGER, intent(in) :: MSGUN
!            EPOCH  = Time formatted as UT year and fraction; e.g., 1990.0
!                     is 0 UT 1 Jan 1990.  For APXMKA EPOCH is single valued;
!                     for APXWRA EPOCH contains NEPOCH values.
  REAL, intent(in) :: EPOCH
!            GPLAT  = Grid point latitudes, an array of NLAT real values in
!                     ascending numerical order with each value in the range
!                     -90. to +90. degrees with positive north.
  REAL, intent(in) :: GPLAT(NLAT)
!            GPLON  = Grid point longitudes, an array of NLON real values in
!                     ascending numerical order with each value in the range
!                     -270. to 270. degrees with positive east.
  REAL, intent(in) :: GPLON(NLON)
!            GPALT  = Grid point altitudes, an array of NALT real values in
!                     ascending numerical order with none less than zero in
!                     units of km above ground.
  REAL, intent(in) :: GPALT(NALT)
!            NLAT   = Number of latitudes  in GPLAT.
  INTEGER, intent(in) :: NLAT
!            NLON   = Number of longitudes in GPLON.
  INTEGER, intent(in) :: NLON
!            NALT   = Number of altitudes  in GPALT.
  INTEGER, intent(in) :: NALT
!            WK     = Work array which should not be altered between
!                     initialization (APXMKA, APXWRA or APXRDA) and use (APXGGC,
!                     APXALL, APXMALL, APXQ2G, APXA2G, APXM2G).
  REAL, intent(inout) :: WK(NLAT*NLON*NALT*6)
!            LWK    = Dimension of WK >=  NLAT*NLON*NALT*5 + NLAT+NLON+NALT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            IST    = Return status: 0 (okay) or > 0 (failed)
  INTEGER, intent(out) :: IST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER :: LWK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LWK = NLAT*NLON*NALT*6
  call APXMKA_legacy(MSGUN, EPOCH, GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,LWK, IST)
END SUBROUTINE APXMKA

SUBROUTINE APXGGC(MSGUN,WK,GPLAT,GPLON,GPALT,NLAT,NLON,NALT,IST)
!          Get grid coordinates for tables currently in memory; i.e., after
!          calling APXMKA, APXWRA or APXRDA.
  implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            MSGUN  = Fortran unit number to write diagnostic messages.
  integer, intent(in) :: MSGUN
!            NLAT   = Number of latitudes  in GPLAT
  INTEGER, intent(in) :: NLAT
!            NLON   = Number of longitudes in GPLON.
  INTEGER, intent(in) :: NLON
!            NALT   = Number of altitudes  in GPALT.
  INTEGER, intent(in) :: NALT
!            WK     = Work array which should not be altered between
!                     initialization (APXMKA, APXWRA or APXRDA) and use (APXGGC,
!                     APXALL, APXMALL, APXQ2G, APXA2G, APXM2G).
  REAL, intent(inout) :: WK(nlat*nlon*nalt*6)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            GPLAT  = Real array of grid point latitudes  in degrees.
  REAL, intent(out) :: GPLAT(NLAT)
!            GPLON  = Real array of grid point longitudes in degrees.
  REAL, intent(out) :: GPLON(NLON)
!            GPALT  = Real array of grid point altitudes  in km.
  REAL, intent(out) :: GPALT(NALT)

!            IST    = Return status: 0 (okay).
  INTEGER, intent(out) :: IST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer nlat_out, nlon_out, nalt_out
!            LWK    = Dimension of WK >=  NLAT*NLON*NALT*5 + NLAT+NLON+NALT
  INTEGER :: LWK  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LWK=NLAT*NLON*NALT*6
  call APXGGC_legacy(MSGUN,WK,LWK, GPLAT,GPLON,GPALT,nlat_out,nlon_out,nalt_out,IST)

  if (nlat_out .ne. nlat) then
     write(*,*) "Error in APXGGC: Expected nlat=",nlat,"got",nlat_out
  endif
  if (nlon_out .ne. nlon) then
     write(*,*) "Error in APXGGC: Expected nlon=",nlon,"got",nlon_out
  endif
  if (nalt_out .ne. nalt) then
     write(*,*) "Error in APXGGC: Expected nalt=",nalt,"got",nalt_out
  endif
END SUBROUTINE APXGGC

SUBROUTINE APXWRA(MSGUN, FILNAM,IUN, EPOCH,NEPOCH,GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,LWK, IST)
!          APXWRA create interpolation tables for one or more times and the tables in a
!          file.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
!            MSGUN  = Fortran unit number to write diagnostic messages.
  INTEGER, intent(in) :: MSGUN
!            FILNAM = file name where arrays are stored.
  CHARACTER, intent(in) :: FILNAM*(*)
!            IUN    = Fortran unit number to be associated with FILNAM.
  INTEGER, intent(in) :: IUN
!            EPOCH  = Time formatted as UT year and fraction; e.g., 1990.0
!                     is 0 UT 1 Jan 1990.  For APXMKA EPOCH is single valued;
!                     for APXWRA EPOCH contains NEPOCH values.
  REAL, intent(in) :: EPOCH(*)
!            NEPOCH = Number of times in EPOCH.
  INTEGER, intent(in) :: NEPOCH
!            GPLAT  = Grid point latitudes, an array of NLAT real values in
!                     ascending numerical order with each value in the range
!                     -90. to +90. degrees with positive north.
  REAL, intent(in) :: GPLAT(*)
!            GPLON  = Grid point longitudes, an array of NLON real values in
!                     ascending numerical order with each value in the range
!                     -270. to 270. degrees with positive east.
  REAL, intent(in) :: GPLON(*)
!            GPALT  = Grid point altitudes, an array of NALT real values in
!                     ascending numerical order with none less than zero in
!                     units of km above ground.
  REAL, intent(in) :: GPALT(*)
!            NLAT   = Number of latitudes  in GPLAT.
  INTEGER, intent(in) :: NLAT
!            NLON   = Number of longitudes in GPLON.
  INTEGER, intent(in) :: NLON
!            NALT   = Number of altitudes  in GPALT.
  INTEGER, intent(in) :: NALT
!            LWK    = Dimension of WK >=  NLAT*NLON*NALT*5 + NLAT+NLON+NALT
  INTEGER, intent(in) :: LWK
!            WK     = Work array which should not be altered between
!                     initialization (APXMKA, APXWRA or APXRDA) and use (APXGGC,
!                     APXALL, APXMALL, APXQ2G, APXA2G, APXM2G).
  REAL, intent(in) :: WK(LWK)

!            IST    = Return status: 0 (okay) or > 0 (failed)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER, intent(out) :: IST
  !call APXWRA_legacy(MSGUN, FILNAM,IUN, EPOCH,NEPOCH,GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,LWK, IST)
  IST=1
  write(*,*) "WARNING: APXWRA is not implemented.  Need to rewrite APXWRA."

END SUBROUTINE APXWRA

SUBROUTINE APXRDA(MSGUN, FILNAM,IUN, inputDATE, WK,LWK, IST)
!          Read back tables previously created by calling APXWRA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
!            MSGUN  = Fortran unit number to write diagnostic messages.
  INTEGER, intent(in) :: MSGUN
!            FILNAM = file name where arrays are stored.
  CHARACTER, intent(in) :: FILNAM*(*)
!            IUN    = Fortran unit number to be associated with FILNAM.
  INTEGER, intent(in) :: IUN
!            DATE   = Time formatted as UT year and fraction; e.g., 1990.0
!                     is 0 UT 1 Jan 1990.  DATE is only used when FILNAM
!                     contains tables for more than one time.
  REAL, intent(in) :: inputDATE
!            LWK    = Dimension of WK >=  NLAT*NLON*NALT*5 + NLAT+NLON+NALT
  INTEGER, intent(in) :: LWK
!            WK     = Work array which should not be altered between
!                     initialization (APXMKA, APXWRA or APXRDA) and use (APXGGC,
!                     APXALL, APXMALL, APXQ2G, APXA2G, APXM2G).
  REAL, intent(in) :: WK(LWK)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            IST    = status: 0 (okay), -1 (minor date problem) or > 0 (failed)
  INTEGER, intent(out) :: IST
  !call APXRDA_legacy(MSGUN, FILNAM,IUN, inputDATE, WK,LWK, IST)
  IST=1
  write(*,*) "WARNING: APXRDA is not implemented.  Need to rewrite APXRDA."

END SUBROUTINE APXRDA

SUBROUTINE  APXMALL(GLAT,GLON,ALT,HR, WK, LWK,                & !Inputs
                    B,BHAT,BMAG,SI,                           & !Mag Fld
                    ALON,                                     & !Apx Lon
                    XLATM,VMP,WM,D,BE3,SIM,D1,D2,D3,E1,E2,E3, & !Mod Apx
                    XLATQD,F,F1,F2 , IST)                       !Qsi-Dpl

!          Determine Modified Apex and Quasi-dipole coordinates using tables
!          currently in memory; i.e. after calling APXMKA_legacy, APXWRA_legacy, or APXRDA_legacy.
!          The above reference (Richmond, 1995) describes HR and all returned
!          variables.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
!            GLAT   = Geographic latitude (degrees); must be within the grid
!                     latitude limits.
  REAL, intent(in) :: GLAT
!            GLON   = Geographic longitude in degrees; must be within one
!                     revolution of the grid longitude limits.
  REAL, intent(in) :: GLON
!            ALT    = Altitude, km above ground.
  REAL, intent(in) :: ALT
!            HR     = Modified Apex coordinates reference altitude, km.
  REAL, intent(in) :: HR
  INTEGER, intent(in) :: LWK
!            WK     = same as entry APXMKA_legacy
  REAL, intent(in) :: WK(LWK)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            B      = magnetic field components (east, north, up), in nT
  REAL, intent(out) :: B(3)
!            BHAT   = components (east, north, up) of unit vector along the
!                     geomagnetic field direction
  REAL, intent(out) :: BHAT(3)
!            BMAG   = magnetic field magnitude, nT
  REAL, intent(out) :: BMAG
!            SI     = sin(I) where I is the angle of inclination of the field
!                     line below the horizontal
  REAL, intent(out) :: SI
!            ALON   = Apex longitude = Modified Apex longitude = Quasi-Dipole
!                     longitude, degrees
  REAL, intent(out) :: ALON
!            XLATM  = Modified Apex latitude, degrees
  REAL, intent(out) :: XLATM
!            VMP    = magnetic potential, in T.m
  REAL, intent(out) :: VMP
!            WM     = Wm of reference above, in km**2 /nT; i.e., 10**15 m**2 /T)
  REAL, intent(out) :: WM
!            D      = D of reference above
  REAL, intent(out) :: D
!            BE3    = B_e3 of reference above (= Bmag/D), in nT
  REAL, intent(out) :: BE3
!            SIM    = sin(I_m) of reference above
  REAL, intent(out) :: SIM
!            D1,D2,D3,E1,E2,E3 = Modified Apex coordinates base vectors, each with
!                     three components (east, north, up) as described in reference above
  REAL, intent(out) :: D1(3)
  REAL, intent(out) :: D2(3)
  REAL, intent(out) :: D3(3)
  REAL, intent(out) :: E1(3)
  REAL, intent(out) :: E2(3)
  REAL, intent(out) :: E3(3)
!            XLATQD = Quasi-Dipole latitude, degrees
  REAL, intent(out) :: XLATQD
!            F      = F for Quasi-Dipole coordinates described in reference above
  REAL, intent(out) :: F
!            F1,F2  = Quasi-Dipole coordinates base vectors with components
!                     (east, north) described in reference above
  REAL, intent(out) :: F1(2)
  REAL, intent(out) :: F2(2)
!            IST    = Return status: 0 (okay) or 1 (failed).
  INTEGER, intent(out) :: IST

  call APXMALL_legacy(GLAT,GLON,ALT,HR, WK,                     & !Inputs
                      B,BHAT,BMAG,SI,                           & !Mag Fld
                      ALON,                                     & !Apx Lon
                      XLATM,VMP,WM,D,BE3,SIM,D1,D2,D3,E1,E2,E3, & !Mod Apx
                      XLATQD,F,F1,F2 , IST)                       !Qsi-Dpl
END SUBROUTINE APXMALL


SUBROUTINE APXALL(GLAT,GLON,ALT, WK, A,ALAT,ALON, NLAT,NLON,NALT, IST)
!          Determine Apex coordinates using tables currently in memory; i.e.,
!          after calling APXMKA_legacy, APXWRA_legacy, or APXRDA_legacy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
!            GLAT = Geographic latitude (degrees) must be within grid latitude limits.
  REAL, intent(in) :: GLAT
!            GLON = Geographic longitude in degrees; must be within one
!                   revolution of the grid longitude limits.
  REAL, intent(in) :: GLON
!            ALT  = Altitude, km.
  REAL, intent(in) :: ALT
!            WK   = same as entry APXMKA_legacy
  REAL, intent(inout) :: WK(NLAT*NLON*NALT)

  INTEGER, intent(in) :: NLAT
  INTEGER, intent(in) :: NLON
  INTEGER, intent(in) :: NALT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            A    = Apex radius, normalized by Req
  REAL, intent(out) :: A
!            ALAT = Apex latitude, degrees
  REAL, intent(out) :: ALAT
!            ALON = Apex longitude, degrees
  REAL, intent(out) :: ALON
!            IST  = Return status:  okay (0); or failure (1).
  INTEGER, intent(out) :: IST

  call APXALL_legacy(GLAT,GLON,ALT, WK, A,ALAT,ALON, IST)
END SUBROUTINE APXALL


SUBROUTINE APXQ2G(QDLAT,QDLON,ALT, WK, LWK, GDLAT,GDLON, IST)
!          Convert from Quasi-Dipole to geographic coordinates using tables
!          currently in memory; i.e. after calling APXMKA_legacy, APXWRA_legacy, or APXRDA_legacy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
!            QDLAT = Quasi-Dipole latitude in degrees
  REAL, intent(in) :: QDLAT
!            QDLON = Quasi-Dipole longitude in degrees
  REAL, intent(in) :: QDLON
!            ALT   = altitude in km
  REAL, intent(in) :: ALT
!            WK    = same as entry APXMKA_legacy
  INTEGER, intent(in) :: LWK
  REAL, intent(in) :: WK(lwk)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            GDLAT = geodetic latitude in degrees
  REAL, intent(out) :: GDLAT
!            GDLON = geodetic longitude in degrees
  REAL, intent(out) :: GDLON
!            IST   = status: 0 (okay), -1 (results are not as close as PRECISE),
!                    or > 0 (failed)
  INTEGER, intent(out) :: IST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call  APXQ2G_legacy(QDLAT,QDLON,ALT, WK, GDLAT,GDLON, IST)
END SUBROUTINE APXQ2G


SUBROUTINE APXA2G(ALAT,ALON,ALT, WK, LWK, GDLAT,GDLON, IST)
!          Convert from Quasi-Dipole to geographic coordinates using tables
!          currently in memory; i.e. after calling APXMKA_legacy, APXWRA_legacy, or APXRDA_legacy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
!            ALAT  = Apex latitude in degrees
  REAL, intent(in) :: ALAT
!            ALON  = Apex longitude in degrees
  REAL, intent(in) :: ALON
!            ALT   = altitude in km
  REAL, intent(in) :: ALT
!            WK    = same as entry APXMKA_legacy
  INTEGER, intent(in) :: LWK
  REAL, intent(in) :: WK(LWK)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            GDLAT = geodetic latitude in degrees
  REAL, intent(out) :: GDLAT
!            GDLON = geodetic longitude in degrees
  REAL, intent(out) :: GDLON
!            IST   = status: 0 (okay), -1 (results are not as close as PRECISE),
!                    or > 0 (failed)
  INTEGER, intent(out) :: IST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call APXA2G_legacy(ALAT,ALON,ALT, WK, GDLAT,GDLON, IST)
END SUBROUTINE APXA2G


SUBROUTINE APXM2G(XLATM,ALON,ALT,HR, WK,LWK, GDLAT,GDLON, IST)
!          Convert from Modified Apex to geographic coordinates using tables
!          currently in memory; i.e. after calling APXMKA_legacy, APXWRA_legacy, or APXRDA_legacy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none
!            XLATM = Modified Apex latitude in degrees
  REAL, intent(in) :: XLATM
!            ALON  = Modified Apex longitude (same as Apex longitude) in degrees
  REAL, intent(in) :: ALON
!            ALT   = altitude in km
  REAL, intent(in) :: ALT
!            HR   = Reference altitude, km (used only for Modified Apex coords)
  REAL, intent(in) :: HR
!            WK    = same as entry APXMKA_legacy
  INTEGER, intent(in) :: LWK
  REAL, intent(in) :: WK(LWK)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            GDLAT = geodetic latitude in degrees
  REAL, intent(out) :: GDLAT
!            GDLON = geodetic longitude in degrees
  REAL, intent(out) :: GDLON
!            IST   = status: 0 (okay), -1 (results are not as close as PRECISE),
!                    or > 0 (failed)
  INTEGER, intent(out) :: IST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call APXM2G_legacy(XLATM,ALON,ALT,HR, WK, GDLAT,GDLON, IST)
END SUBROUTINE APXM2G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  The following are helper subroutines used by the Apex subroutines   !
!  (APXMKA, APXGGC, APXWRA, APXRDA, APXALL, APXQ2G, APXA2G, APXM2G),   !
!  defined above                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INTRP(GLAT,GLON,ALT, X,Y,Z,V, NLAT,NLON,NALT, &
                 GPLAT,GPLON,GPALT, FX,FY,FZ,FV, &
                 DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN, &
                 DFVDLN,DFXDH,DFYDH,DFZDH,DFVDH, CALNM,IST)
  implicit none
  REAL, intent(in) :: GLAT
  REAL, intent(in) :: GLON
  REAL, intent(in) :: ALT
  REAL, intent(in) :: X(NLAT,NLON,NALT)
  REAL, intent(in) :: Y(NLAT,NLON,NALT)
  REAL, intent(in) :: Z(NLAT,NLON,NALT)
  REAL, intent(in) :: V(NLAT,NLON,NALT)
  INTEGER, intent(in) :: NLAT
  INTEGER, intent(in) :: NLON
  INTEGER, intent(in) :: NALT
  REAL, intent(in) :: GPLAT(NLAT)
  REAL, intent(in) :: GPLON(NLON)
  REAL, intent(in) :: GPALT(NALT)
  CHARACTER, intent(in) :: CALNM*(*)

  REAL, intent(out) :: FX
  REAL, intent(out) :: FY
  REAL, intent(out) :: FZ
  REAL, intent(out) :: FV
  REAL, intent(out) :: DFXDTH,DFYDTH,DFZDTH,DFVDTH
  REAL, intent(out) :: DFXDLN,DFYDLN,DFZDLN,DFVDLN
  REAL, INTENT(out) :: DFXDH,DFYDH,DFZDH,DFVDH
  INTEGER, intent(out) :: IST
  
  call  INTRP_legacy (GLAT,GLON,ALT, X,Y,Z,V, NLAT,NLON,NALT, &
                      GPLAT,GPLON,GPALT, FX,FY,FZ,FV, &
                      DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN, &
                      DFVDLN,DFXDH,DFYDH,DFZDH,DFVDH, CALNM,IST)
END SUBROUTINE INTRP

SUBROUTINE INTRPSC(GLAT,GLON,ALT, X,Y,Z, NLAT,NLON,NALT, &
                   GPLAT,GPLON,GPALT, FX,FY,FZ , CALNM,IST)
  implicit none
  REAL, intent(in) :: GLAT
  REAL, intent(in) :: GLON
  REAL, intent(in) :: ALT
  REAL, intent(in) :: X(NLAT,NLON,NALT)
  REAL, intent(in) :: Y(NLAT,NLON,NALT)
  REAL, intent(in) :: Z(NLAT,NLON,NALT)
  INTEGER, intent(in) :: NLAT
  INTEGER, intent(in) :: NLON
  INTEGER, intent(in) :: NALT
  REAL, intent(in) :: GPLAT(NLAT)
  REAL, intent(in) :: GPLON(NLON)
  REAL, intent(in) :: GPALT(NALT)
  CHARACTER, intent(in) :: CALNM*(*)

  REAL, intent(out) :: FX
  REAL, intent(out) :: FY
  REAL, intent(out) :: FZ
  INTEGER, intent(out) :: IST

  call  INTRPSC_legacy(GLAT,GLON,ALT, X,Y,Z, NLAT,NLON,NALT, &
                       GPLAT,GPLON,GPALT, FX,FY,FZ , CALNM,IST)
END SUBROUTINE INTRPSC

SUBROUTINE TRILIN(U,NLAT,NLON,XI,YJ,ZK,FU,DFUDX,DFUDY,DFUDZ)
  implicit none
  REAL, intent(in) :: U(NLAT,NLON,2)
  INTEGER, intent(in) :: NLAT
  INTEGER, intent(in) :: NLON
  REAL, intent(in) :: XI
  REAL, intent(in) :: YJ
  REAL, intent(in) :: ZK

  REAL, intent(out) :: FU
  REAL, intent(out) :: DFUDX
  REAL, intent(out) :: DFUDY
  REAL, intent(out) :: DFUDZ
  call TRILIN_legacy(U,NLAT,NLON,XI,YJ,ZK,FU,DFUDX,DFUDY,DFUDZ)
END SUBROUTINE TRILIN

SUBROUTINE TRILINS(U,NLAT,NLON,XI,YJ,ZK,FU)
  implicit none
  REAL, intent(in) :: U(NLAT,NLON,2)
  INTEGER, intent(in) :: NLAT
  INTEGER, intent(in) :: NLON
  REAL, intent(in) :: XI
  REAL, intent(in) :: YJ
  REAL, intent(in) :: ZK

  REAL, intent(out) :: FU
  call  TRILINS_legacy(U,NLAT,NLON,XI,YJ,ZK,FU)
END SUBROUTINE TRILINS

SUBROUTINE ADPL(GLAT,GLON,CTH,STH,FX,FY,FZ,FV, DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
  implicit none
  REAL, intent(in) :: GLAT
  REAL, intent(in) :: GLON

  REAL, intent(inout) :: FX
  REAL, intent(inout) :: FY
  REAL, intent(inout) :: FZ
  REAL, intent(inout) :: FV
  REAL, intent(inout) :: DFXDTH,DFYDTH,DFZDTH,DFVDTH
  REAL, intent(inout) :: DFXDLN,DFYDLN,DFZDLN,DFVDLN

  REAL, intent(out) :: cth
  REAL, intent(out) :: sth
  call ADPL_legacy(GLAT,GLON,CTH,STH,FX,FY,FZ,FV, DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
END SUBROUTINE ADPL

SUBROUTINE  ADPLSC(GLAT,GLON,FX,FY,FZ)
  implicit none
  REAL, intent(in) :: GLAT
  REAL, intent(in) :: GLON

  REAL, intent(inout) :: FX
  REAL, intent(inout) :: FY
  REAL, intent(inout) :: FZ
  
  call ADPLSC_legacy(GLAT,GLON,FX,FY,FZ)
END SUBROUTINE ADPLSC

SUBROUTINE GRADXYZV(ALT,CTH,STH, &
                    DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN, &
                    DFXDH,DFYDH,DFZDH,DFVDH,GRADX,GRADY,GRADZ,GRADV)

  implicit none
  REAL, intent(in) :: ALT
  REAL, intent(in) :: CTH
  REAL, intent(in) :: STH
  REAL, intent(in) :: DFXDTH,DFYDTH,DFZDTH,DFVDTH
  REAL, intent(in) :: DFXDLN,DFYDLN,DFZDLN,DFVDLN
  REAL, intent(in) :: DFXDH,DFYDH,DFZDH,DFVDH

  REAL, intent(out) :: GRADX(3)
  REAL, intent(out) :: GRADY(3)
  REAL, intent(out) :: GRADZ(3)
  REAL, intent(out) :: GRADV(3)

  call GRADXYZV_legacy(ALT,CTH,STH, &
                       DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN, &
                       DFXDH,DFYDH,DFZDH,DFVDH,GRADX,GRADY,GRADZ,GRADV)
END SUBROUTINE GRADXYZV

SUBROUTINE  GRAPXYZV(ALT,CTH,STH,DFXDLN,DFYDLN,DFZDLN,DFVDLN,GRADX,GRADY,GRADZ,GRADV)
  implicit none
  REAL, intent(in) :: ALT
  REAL, intent(in) :: CTH
  REAL, intent(in) :: STH
  REAL, intent(in) :: DFXDLN,DFYDLN,DFZDLN,DFVDLN

  REAL, intent(out) :: GRADX(3)
  REAL, intent(out) :: GRADY(3)
  REAL, intent(out) :: GRADZ(3)
  REAL, intent(out) :: GRADV(3)

  call GRAPXYZV_legacy(ALT,CTH,STH,DFXDLN,DFYDLN,DFZDLN,DFVDLN,GRADX,GRADY,GRADZ,GRADV)
END SUBROUTINE GRAPXYZV

SUBROUTINE GRADLPV(HR,ALT,FX,FY,FZ,FV,GRADX,GRADY,GRADZ,GRADV, &
                   XLATM,XLONM,VMP,GRCLM,CLMGRP,QDLAT,RGRLP,B,CLM,R3_2)
  implicit none
  REAL, intent(in) :: HR
  REAL, intent(in) :: ALT
  REAL, intent(in) :: FX
  REAL, intent(in) :: FY
  REAL, intent(in) :: FZ
  REAL, intent(in) :: FV
  REAL, intent(in) :: GRADX(3)
  REAL, intent(in) :: GRADY(3)
  REAL, intent(in) :: GRADZ(3)
  REAL, intent(in) :: GRADV(3)

  REAL, intent(out) :: XLATM
  REAL, intent(out) :: XLONM
  REAL, intent(out) :: VMP
  REAL, intent(out) :: GRCLM(3)
  REAL, intent(out) :: CLMGRP(3)
  REAL, intent(out) :: QDLAT
  REAL, intent(out) :: RGRLP(3)
  REAL, intent(out) :: B(3)
  REAL, intent(out) :: CLM
  REAL, intent(out) :: R3_2

  call GRADLPV_legacy(HR,ALT,FX,FY,FZ,FV,GRADX,GRADY,GRADZ,GRADV, &
                      XLATM,XLONM,VMP,GRCLM,CLMGRP,QDLAT,RGRLP,B,CLM,R3_2)
END SUBROUTINE GRADLPV

SUBROUTINE XYZ2APX(ALT,FX,FY,FZ,A,ALAT,ALON,IERR)
  implicit none
  REAL, intent(in) :: ALT
  REAL, intent(in) :: FX
  REAL, intent(in) :: FY
  REAL, intent(in) :: FZ

  REAL, intent(out) :: A
  REAL, intent(out) :: ALAT
  REAL, intent(out) :: ALON
  INTEGER, intent(out) :: IERR

  call XYZ2APX_legacy(ALT,FX,FY,FZ,A,ALAT,ALON,IERR)
END SUBROUTINE XYZ2APX

SUBROUTINE BASVEC(HR,XLATM,GRCLM,CLMGRP,RGRLP,B,CLM,R3_2,BMAG,SIM,SI,F,D,W,BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
  implicit none

  REAL, intent(in) :: HR
  REAL, intent(in) :: XLATM
  REAL, intent(in) :: GRCLM(3)
  REAL, intent(in) :: CLMGRP(3)
  REAL, intent(in) :: RGRLP(3)
  REAL, intent(in) :: B(3)
  REAL, intent(in) :: CLM
  REAL, intent(in) :: R3_2

  REAL, intent(out) :: BMAG
  REAL, intent(out) :: SIM
  REAL, intent(out) :: SI
  REAL, intent(out) :: F
  REAL, intent(out) :: D
  REAL, intent(out) :: W
  REAL, intent(out) :: BHAT(3)
  REAL, intent(out) :: D1
  REAL, intent(out) :: D2
  REAL, intent(out) :: D3
  REAL, intent(out) :: E1
  REAL, intent(out) :: E2
  REAL, intent(out) :: E3
  REAL, intent(out) :: F1
  REAL, intent(out) :: F2
  
  call BASVEC_legacy(HR,XLATM,GRCLM,CLMGRP,RGRLP,B,CLM,R3_2,BMAG,SIM,SI,F,D,W,BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
END SUBROUTINE BASVEC

SUBROUTINE CKGP(CALNM,MSGUN,NLAT,NLON,NALT,GPLAT,GPLON,GPALT,IST)
  implicit none
  CHARACTER, intent(in) :: CALNM*(*)
  INTEGER, intent(in) :: MSGUN
  INTEGER, intent(in) :: NLAT
  INTEGER, intent(in) :: NLON
  INTEGER, intent(in) :: NALT
  REAL, intent(IN) :: GPLAT(NLAT)
  REAL, intent(IN) :: GPLON(NLON)
  REAL, intent(IN) :: GPALT(NALT)
  INTEGER, intent(out) :: IST

  call CKGP_legacy(CALNM,MSGUN,NLAT,NLON,NALT,GPLAT,GPLON,GPALT,IST)
END SUBROUTINE CKGP

SUBROUTINE MAKEXYZV(EPOCH,NLAT,NLON,NALT,GPLAT,GPLON,GPALT, X,Y,Z,V)
  implicit none
  REAL, intent(in) :: EPOCH
  INTEGER, intent(in) :: NLAT
  INTEGER, intent(in) :: NLON
  INTEGER, intent(in) :: NALT
  REAL, intent(in) :: GPLAT(NLAT)
  REAL, intent(in) :: GPLON(NLON)
  REAL, intent(in) :: GPALT(NALT)

  REAL, intent(out) :: X
  REAL, intent(out) :: Y
  REAL, intent(out) :: Z
  REAL, intent(out) :: V

  call MAKEXYZV_legacy(EPOCH,NLAT,NLON,NALT,GPLAT,GPLON,GPALT, X,Y,Z,V)
END SUBROUTINE MAKEXYZV

SUBROUTINE SETMISS(XMISS, &
                   XLATM,XLONM,VMP,B,BMAG,BE3,SIM,SI,F,D,W, &
                   BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
  implicit none

  REAL, intent(in) :: XMISS

  REAL, intent(out) :: XLATM
  REAL, intent(out) :: XLONM
  REAL, intent(out) :: VMP
  REAL, intent(out) :: B(3)
  REAL, intent(out) :: BMAG
  REAL, intent(out) :: BE3
  REAL, intent(out) :: SIM
  REAL, intent(out) :: SI
  REAL, intent(out) :: F
  REAL, intent(out) :: D
  REAL, intent(out) :: W
  REAL, intent(out) :: BHAT(3)
  REAL, intent(out) :: D1(3)
  REAL, intent(out) :: D2(3)
  REAL, intent(out) :: D3(3)
  REAL, intent(out) :: E1(3)
  REAL, intent(out) :: E2(3)
  REAL, intent(out) :: E3(3)
  REAL, intent(out) :: F1(3)
  REAL, intent(out) :: F2(3)

  call SETMISS_legacy(XMISS, &
                      XLATM,XLONM,VMP,B,BMAG,BE3,SIM,SI,F,D,W, &
                      BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
END SUBROUTINE SETMISS

SUBROUTINE GM2GC(GMLAT,GMLON,GCLAT,GCLON)
  implicit none
    REAL, intent(in) :: GMLAT
    REAL, intent(in) :: GMLON

    REAL, intent(out) :: GCLAT
    REAL, intent(out) :: GCLON
  call GM2GC_legacy(GMLAT,GMLON,GCLAT,GCLON)
END SUBROUTINE GM2GC
