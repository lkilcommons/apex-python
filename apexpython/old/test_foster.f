      ! This program was written to demonstrate that the Python bindings
      ! to the Apex codebase work properly.  Compare the results of 
      ! this code to the Python unittest at test_apxntrp.py
      program test

      integer nlat, nlatp1
      parameter (nlat=36)
      parameter (nlatp1=nlat+1)
      integer nlon, nlonp2
      parameter (nlon=72)
      parameter (nlonp2=nlon+2)
      integer nalt
      parameter (nalt=2)
      real workArray(nlat*nlon*nalt*6)
      real gplat(nlatp1)
      real gplon(nlonp2)
      real gpalt(nalt)
      real gdlat, gdlon
      integer ist
      
      integer i,j
      real dellat, dellon

      write(*,*) "Hello world"

      gpalt(1) = 90.0
      gpalt(2) = 170.0

      dellat = 180.0 / float(nlat)
      do j= 1,nlatp1
         gplat(j) = (j-1)*dellat - 90.0
      enddo

      dellon = 360.0 / float(nlon)
      do i=1,nlonp2
         gplon(i) = (float(i) - 1.5)*dellon - 180.0
      enddo

      call APXMKA(6, 2002.0, 
     $     gplat,gplon,gpalt,
     $     nlat,nlon,nalt, 
     $     workArray,
     $     ist)
      write(*,*) "APXMKA ist", ist
      
      !call APXGGC(6, workArray, gplat,gplon,gpalt, nlat,nlon,nalt, ist)
      !write(*,*) "APXGGC gplat:",gplat
      !write(*,*) "APXGGC gplon:",gplon
      !write(*,*) "APXGGC gpalt:",gpalt
      !write(*,*) "APXGGC ist", ist

      call APXQ2G(-29.9158214580855, 18.0, 130.0, workArray, 
     $     nlat*nlon*alt*6,gdlat, gdlon, ist)
      write(*,'(''gdlat='',F16.12,'' gdlon='',F16.12,'' ist='',I3)')
     $     gdlat,gdlon,ist

      end program test
