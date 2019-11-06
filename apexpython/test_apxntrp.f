      ! This program was written to demonstrate that the Python bindings
      ! to the Apex codebase work properly.  Compare the results of 
      ! this code to the Python unittest at test_apxntrp.py
      program test

      integer nlat
      integer nnlat
      parameter (nnlat=13)
      integer nlon
      integer nnlon
      parameter (nnlon=21)
      integer nalt
      integer nnalt
      parameter (nnalt=2)
      real workArray(nnlat*nnlon*nnalt*6)
      real gplat(nnlat)
      real gplon(nnlon)
      real gpalt(nnalt)
      real gdlat, gdlon
      integer ist
      
      real a,alat,alon

      write(*,*) "Hello world"

      nlat=13
      nlon=21
      nalt=2

      call GGRID(4, 
     $     -90.0, 90.0,
     $     -180.0, 180.0,
     $     300.0, 300.0,
     $     gplat, gplon, gpalt, 
     $     nlat, nlon, nalt)
!      call GGRID(30, 
!     $     -19.0,-18.0,
!     $     -67.0,-66.0,
!     $     299.0,301.0, 
!     $     gplat, gplon, gpalt, 
!     $     nlat, nlon, nalt)
!
      write(*,*) "nlat:",nlat,"gplat:",gplat
      write(*,*) "nlon:",nlon,"gplon:",gplon
      write(*,*) "nalt:",nalt,"gpalt:",gpalt

      call APXMKA(6, 1994.5, 
     $     gplat,gplon,gpalt,
     $     nlat,nlon,nalt, 
     $     workArray,
     $     ist)
      write(*,*) "APXMKA ist", ist

      call APXGGC(6, workArray, gplat,gplon,gpalt, nlat,nlon,nalt, ist)
      write(*,*) "APXGGC gplat:",gplat
      write(*,*) "APXGGC gplon:",gplon
      write(*,*) "APXGGC gpalt:",gpalt
      write(*,*) "APXGGC ist", ist

      call APXALL(-18.3, -66.75, 300.0, 
     $     workArray, 
     $     a,alat,alon, 
     $     nlat,nlon,nalt, ist)
      write(*,*) "APXALL a:",a
      write(*,*) "APXALL alat:",alat
      write(*,*) "APXALL alon:",alon
      write(*,*) "APXALL ist",ist

      call APXA2G(alat, alon, 300.0, workArray, nlat*nlon*nalt*6,
     $     gdlat,gdlon,ist)
      write(*,*) "APXA2G gdlat:",gdlat
      write(*,*) "APXA2G gdlon:",gdlon
      write(*,*) "APXA2G ist:",ist

      end program test
