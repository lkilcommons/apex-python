
#builtin
import sys
import os.path
import datetime

# custom
from apexpy import apex

# 3rd-party
import numpy
#import matplotlib.pyplot as pp
#import matplotlib as mpl

class apex_converter:

	def __init__(self,glatmin=-90.,glatmax=90.,glonmin=-180.,glonmax=180.,altmin=300., altmax=900.,
		nvert=50.,epoch=2010.0):
		"""	
		Bulk conversion of geographic coordinates to Modified Apex

		Uses numpy-style docstrings for clarity. 
		Parameters described are for class constructor

		Parameters
		----------	
		glatmin : float, optional
			Grid boundary: minimum geographic latitude (default=-90.)
		glatmax : float, optional
			Grid boundary: maximum geographic latitude (default=-90.)
		glonmin : float, optional
			Grid boundary: minimum geographic longitude (default=-180.)
		glonmax : float, optional
			Grid boundary: maximum geographic longitude (default=180.)
		altmin : float, optional
			Grid boundary: minimum altitude in km (default=300.)
		altmax : float, optional
			Grid boundary: minimum altitude in km (default=900.)
		nvert : float, optional
			Grid spacing factor. Determines vertical and lat,lon resolution.
			Recommended values are between 30 and 100. (default=50.)
		epoch : float, optional
			Epoch for input into IGRF; the year 
			
		Altitude is in km
		Latitude and longitude are in degrees. Negative longitude instead of lon>180.
		nvert is the grid spacing factor. Recommended values are between 30 and 100. 
		Larger values will cause the program
		to use more memory to build larger interpolation tables, 
		but interpolation accuracy will increase up to about nvert of 100
		See ggrid.f for more information.
		"""
		#Much of the following code is ripped out of test_apxntrp, written by Peter Schmidt at NCAR HAO

		# due to a shortcoming in f2py (documented here http://cens.ioc.ee/pipermail/f2py-users/2008-December/001764.html), 
		# we must know dimensions apriori.   
		
		#Code to figure out what nlat, nlon and nalt are from the nvert spacing factor. Translated from fortran code below:
		#Code from ggrid.f:
		#DNV  = DBLE (NVERT)
		#DLON = 360.D0 / (5.D0*DNV)
		#DLAT = 180.D0 / (3.D0*DNV)
		#DIHT = 1.D0 / DNV
		#
		#LAMN =  MAX0 (INT((DBLE(GLAMN) +90.D0 )/DLAT)     , 0)
		#LAMX =  MIN0 (INT((DBLE(GLAMX) +90.D0 )/DLAT+1.D0), 3*NVERT)
		#
		#LOMN  = MAX0 (INT((DBLE(GLOMN) +180.D0)/DLON)  , 0)
		#GLOMXL= AMIN1 (GLOMX,GLOMN+360.)
		#LOMX  = MIN0 (INT((DBLE(GLOMXL)+180.D0)/DLON+1.D0),10*NVERT)
		#
		#X = RE/(RE+ALTMX)/DIHT - 1.E-5
		#IHTMN = AMAX1 (X,1.)
		#IHTMN = MIN0 (IHTMN,NVERT-1)
		#X = RE/(RE+ALTMN)/DIHT + 1.E-5
		#I = X + 1.
		#IHTMX = MIN0(I, NVERT)
		#
		#NLAT = LAMX - LAMN + 1
		#NLON = LOMX - LOMN + 1
		#NLON = MIN0 (NLON,5*NVERT+1)
		#NALT = IHTMX - IHTMN + 1
		
		#Set local variables to their attribute counterparts
		self.glatmin = glatmin
		self.glatmax = glatmax
		self.glonmin = glonmin
		self.altmin = altmin
		self.altmax = altmax
		self.nvert = nvert
		self.epoch = epoch

		re = 6371.2
		dlon = 360./(5.*nvert)
		dlat = 180./(3.*nvert)
		diht = 1./nvert

		lamn = max([(glatmin+90)/dlat,0.])
		lamx = min([(glatmax+90)/dlat+1.,3.*nvert])

		lomn = max([(glonmin+180.)/dlon,0])
		glomxl = min([glonmax,glonmin+360.])
		lomx = min([(glonmax+180.)/dlon+1.,10*nvert])
		
		x1 = re/(re+altmax)/diht - 1.0e-5
		ihtmn = max([x1,1.])
		ihtmn = min([ihtmn,nvert-1.])
		x2 = re/(re+altmin)/diht + 1.0e-5
		I = x2 + 1.	 
		ihtmx = min([I,nvert])

		#print "debug:\n"
		#print "lamn,lamx = %f,%f" % (lamn,lamx)
		#print "lomn,lomx = %f,%f" % (lomn,lomx)
		#print "ihtmn,ihtmx = %f,%f" % (ihtmn,ihtmx)

		#Use floor to make sure we have an integer number of points in the array
		self.nlat = numpy.ceil(lamx - lamn + 1.)
		self.nlon = numpy.ceil(lomx - lomn + 1.)
		self.nlon = min([self.nlon,5*nvert+1.])
		self.nalt = numpy.ceil(ihtmx - ihtmn + 1.)

		self.workArray = numpy.empty(int(self.nlat)*int(self.nlon)*int(self.nalt)*6, numpy.float32)
		# Define grid points suitable for preparing interpolation
		# tables for Apex coordinates, Modified Apex Coordinates and
		# Quasi-Dipole coordinates
		# LMK - nvert = 40 was setting used in Tomoko Matsuo's get_apexB_mod.f fortran code, which was used for conjunction analysis
		print "Initializing lat,lon,alt grid with following parameters:\n"
		print "Epoch, Nvert = %f, %f\n" % (epoch,nvert)
		(epoch,nvert)
		print "Lat min,max,npts = %f,%f,%f\n" % (glatmin,glatmax,self.nlat)
		print "Lon min,max,npts = %f,%f,%f\n" % (glonmin,glonmax,self.nlon)
		print "Alt min,max,npts = %f,%f,%f\n" % (altmin,altmax,self.nalt)
		
		(gplat,gplon,gpalt) = apex.ggrid(nvert=nvert,
										 glamn=glatmin, glamx=glatmax,
										 glomn=glonmin,   glomx=glonmax,
										 altmn=altmin,   altmx=altmax,
										 nlat=self.nlat, nlon=self.nlon,  nalt=self.nalt)
		self.grid = dict()
		self.grid['lat'] = gplat
		self.grid['lon'] = gplon
		self.grid['alt'] = gpalt

		print "Preparing interpolation tables..."
		# Create interpolation tables for a single time
		ist = apex.apxmka(msgun=6, #msgun is message unit?
						  epoch=epoch,
						  gplat=gplat, gplon=gplon, gpalt=gpalt,
						  wk=self.workArray,
						  nlat=self.nlat,   nlon=self.nlon,   nalt=self.nalt) 
		if ist != 0:
			raise RuntimeError('Unable to make the interpolation tables: Return code was %d' % (ist))

	def identicalInputs(self,lat,lon,alt,hr):
		"""
		Test if the self.lastrun represents the results obtained by transforming inputs time,lat,lon,alt,hr
		"""
		retval = False
		if hasattr(self,'lastrun'): #Check for lastrun attribute existance 
			if all(inputval in self.lastrun for inputval in ['lat','lon','alt','hr']): #Check for valid lastrun dictionary keys
				if len(lat) == len(self.lastrun['lat']): #Check for equal size
					equality_checks = [numpy.array_equal(lat,self.lastrun['lat']),
									   numpy.array_equal(lon,self.lastrun['lon']), 
									   numpy.array_equal(alt,self.lastrun['alt']),
									   hr  == self.lastrun['hr']]
					if all(equality_checks): 
						#Check for equal values
						retval = True
		return retval

	def solarParams(self,year,dayofyear,utseconds):
		"""
		Gets the geographic latitude and longitude of the 
		northern geomagnetic pole and the subsolar point(s)
		for a particular year and day of year.

		Note that the geomagnetic north pole coordinate returned
		are determined by the epoch setting of this apex_converter
		instance and not explicitly by the year/day-of-year argument.
		If you are not using a year/day-of-year that agrees 
		with the epoch setting (self.epoch), you'll get incorrect
		results. 

		Parameters
		----------
			year : int or numpy.array, length 1
				The year for which to get the subsolar point(s)
			dayofyear : int or numpy.array, length 1
				The day of the year for which to get the subsolar point(s)
			utseconds : numpy.array
				The Universal Time second of the day for which
				to obtain the subsolar points.
		
		Returns
		-------
			subsllat : numpy.array len(utseconds)
				geographic latitude(s) of subsolar point
			subsllon : numpy.array len(utseconds)
				geographic lons(s) of subsolar point
			nplat : numpy.array len(1)
				north geographic latitude of magnetic north pole
			nplon : numpy.array len(1)
				east geographic longitude of magnetic north pole
		"""

		if not isinstance(year,numpy.ndarray):
			year = numpy.array(year)
		if not isinstance(dayofyear,numpy.ndarray):
			dayofyear = numpy.array(dayofyear)
		if not isinstance(utseconds,numpy.ndarray):
			utseconds = numpy.array(utseconds)

		hour = numpy.floor(utseconds/3600)
		minute = numpy.floor((utseconds - hour*3600)/60)
		second = utseconds - hour*3600 - minute*60

		epoch = self.epoch
		#Set IGRF coeffcients to current epoch
		apex.cofrm(epoch)
		clatp,polon,vp = apex.dypol()

		nplat = 90.-clatp
		nplon = polon

		sbsllat = numpy.zeros_like(utseconds)
		sbsllon = numpy.zeros_like(utseconds)
		if utseconds.size > 1:
			for k in range(utseconds.size):
				sbsllat[k],sbsllon[k] = apex.subsol(year,dayofyear,hour[k],minute[k],second[k])
		else: 
			sbsllat,sbsllon = apex.subsol(year,dayofyear,hour,minute,second)

		return sbsllat,sbsllon,nplat,nplon

	def mlt2alon(self,mlt,year,dayofyear,utseconds):
		"""
		Converts Magnetic Local Time to Apex Longitudes

		Calls apex subroutine apex.mlt2alon (mag local time) 
		after setting the IGRF coefficient global variables using
		the same epoch as was used to initially 
		generate the interpolation tables. 
		
		Parameters
		----------
		mlt : numpy.array
			a numpy array of length n, apex longitude
		year : int or numpy.array
			The year for which to compute the MLT
		dayofyear : int or numpy.array
			The day-of-year(s) for which MLT is to be found
		utseconds : int or numpy.array 
			the second-of-the-day(s) in UT
			of the time for which the MLT is to be found
		
		Options for time argument (year,dayofyear,utseconds) length:
		1. All have length 1 for MLT of all alons for a fixed time.
		2. Year,and Day of year can have length 1, 
			and then UT seconds can have length n 
		 	for different times on same day
		3. All can have length n to find the MLT 
			of each alon with a unique date and time 
			year is year(s) for which MLT is to be found.

		Returns
		-------
		mlt : numpy.array
			The magnetic local time corresponding to alon for
			the times specified
		"""
		#Make sure everything that is supposed to be an array is one
		#Had to do this to make sure everything has a length, since we check that later
		if not isinstance(mlt,numpy.ndarray):
			mlt = numpy.array(mlt)
		if not isinstance(year,numpy.ndarray):
			year = numpy.array(year)
		if not isinstance(dayofyear,numpy.ndarray):
			dayofyear = numpy.array(dayofyear)
		if not isinstance(utseconds,numpy.ndarray):
			utseconds = numpy.array(utseconds)
		
		print "Computing %d Magnetic Local Time values to Apex Longitude...\n" % ( len(mlt) )
		hour = numpy.floor(utseconds/3600)
		minute = numpy.floor((utseconds - hour*3600)/60)
		second = utseconds - hour*3600. - minute*60.

		#Set IGRF coeffcients to current epoch
		apex.cofrm(self.epoch)
		clatp,polon,vp = apex.dypol()
		if year.size==1 and dayofyear.size==1 and utseconds.size==1: 
			sbsllat1,sbsllon1 = apex.subsol(year,dayofyear,hour,minute,second)
			sbsllat = numpy.ones_like(mlt)*sbsllat1
			sbsllon = numpy.ones_like(mlt)*sbsllon1
		elif year.size==mlt.size and dayofyear.size==mlt.size and utseconds.size==mlt.size:
			sbsllat = numpy.zeros_like(mlt)
			sbsllon = numpy.zeros_like(mlt)
			for k,lon in enumerate(mlt):
				sbsllat[k],sbsllon[k] = apex.subsol(year[k],dayofyear[k],hour[k],minute[k],second[k])
		elif year.size==1 and dayofyear.size==1 and utseconds.size==mlt.size:
			sbsllat = numpy.zeros_like(mlt)
			sbsllon = numpy.zeros_like(mlt)
			for k,lon in enumerate(mlt):
				sbsllat[k],sbsllon[k] = apex.subsol(year,dayofyear,hour[k],minute[k],second[k])
		else:
			raise ValueError("Length of time arguments must either be 1 or equal to length of longitude argument!\nyear.size=%d\ndayofyear.size=%d\nutseconds.size=%d" % (year.size,dayofyear.size,utseconds.size))
		alon = numpy.zeros_like(mlt)
		for k in xrange(len(mlt)):
			#Call signitude subroutine mlt2alon(xmlt,sbsllat,sbsllon,clatp,polon,alonx)
			alon[k] = apex.mlt2alon(mlt[k],sbsllat[k],sbsllon[k],clatp,polon)
		#Sanity check magnetic local time
		alon[alon<0] = alon[alon<0]+360.
		alon[alon>360.] = numpy.mod(alon[alon>360.],360.)
		return alon

	def alon2mlt(self,alon,year,dayofyear,utseconds):
		"""
		Converts Apex Longitudes to Magnetic Local Time

		Calls apex subroutine apex.magloctm (mag local time) 
		after setting the IGRF coefficient global variables using
		the same epoch as was used to initially 
		generate the interpolation tables. 
		
		Parameters
		----------
		alon : numpy.array
			a numpy array of length n, apex longitude
		year : int or numpy.array
			The year for which to compute the MLT
		dayofyear : int or numpy.array
			The day-of-year(s) for which MLT is to be found
		utseconds : int or numpy.array 
			the second-of-the-day(s) in UT
			of the time for which the MLT is to be found
		
		Options for time argument (year,dayofyear,utseconds) length:
		1. All have length 1 for MLT of all alons for a fixed time.
		2. Year,and Day of year can have length 1, 
			and then UT seconds can have length n 
		 	for different times on same day
		3. All can have length n to find the MLT 
			of each alon with a unique date and time 
			year is year(s) for which MLT is to be found.

		Returns
		-------
		mlt : numpy.array
			The magnetic local time corresponding to alon for
			the times specified
		"""
		#Make sure everything that is supposed to be an array is one
		#Had to do this to make sure everything has a length, since we check that later
		if not isinstance(alon,numpy.ndarray):
			alon = numpy.array(alon)
		if not isinstance(year,numpy.ndarray):
			year = numpy.array(year)
		if not isinstance(dayofyear,numpy.ndarray):
			dayofyear = numpy.array(dayofyear)
		if not isinstance(utseconds,numpy.ndarray):
			utseconds = numpy.array(utseconds)
		
		print "Computing %d Apex Longitude values to Magnetic Local Time...\n" % ( len(alon) )
		hour = numpy.floor(utseconds/3600)
		minute = numpy.floor((utseconds - hour*3600)/60)
		second = utseconds - hour*3600. - minute*60.

		epoch = self.epoch
		#Set IGRF coeffcients to current epoch
		apex.cofrm(epoch)
		#Obtain the dipole position
		#Original comments:
		# 		      SUBROUTINE DYPOL (COLAT,ELON,VP)
		# C          Computes parameters for dipole component of geomagnetic field.
		# C          COFRM must be called before calling DYPOL!
		# C          940504 A. D. Richmond
		# C
		# C          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(196),GV(196),ICHG
		# C            NMAX = Maximum order of spherical harmonic coefficients used
		# C            GB   = Coefficients for magnetic field calculation
		# C            GV   = Coefficients for magnetic potential calculation
		# C            ICHG = Flag indicating when GB,GV have been changed
		# C
		# C          RETURNS:
		# C            COLAT = Geocentric colatitude of geomagnetic dipole north pole
		# C                    (deg)
		# C            ELON  = East longitude of geomagnetic dipole north pole (deg)
		# C            VP    = Magnitude, in T.m, of dipole component of magnetic
		# C                    potential at geomagnetic pole and geocentric radius
		# C                    of 6371.2 km
		clatp,polon,vp = apex.dypol()

		# 		      SUBROUTINE SUBSOL (IYR,IDAY,IHR,IMN,SEC,SBSLLAT,SBSLLON)
		# C          Find subsolar geographic latitude and longitude given the
		# C          date and time (Universal Time).
		# C
		# C          This is based on formulas in Astronomical Almanac for the
		# C          year 1996, p.  C24. (U.S.  Government Printing Office,
		# C          1994).  According to the Almanac, results are good to at
		# C          least 0.01 degree latitude and 0.025 degree longitude
		# C          between years 1950 and 2050.  Accuracy for other years has
		# C          not been tested although the algorithm has been designed to
		# C          accept input dates from 1601 to 2100.  Every day is assumed
		# C          to have exactly 86400 seconds; thus leap seconds that
		# C          sometimes occur on June 30 and December 31 are ignored:
		# C          their effect is below the accuracy threshold of the algorithm.
		# C
		# C          961026 A. D. Richmond, NCAR
		# C
		# C          INPUTS:
		# C            IYR  = Year (e.g., 1994). IYR must be in the range: 1601 to 2100.
		# C            IDAY = Day number of year (e.g., IDAY = 32 for Feb 1)
		# C            IHR  = Hour of day    (e.g., 13 for 13:49)
		# C            IMN  = Minute of hour (e.g., 49 for 13:49)
		# C            SEC  = Second and fraction after the hour/minute.
		# C          Note:  While IYR is bounds tested, there is no constraint
		# C                 placed on values: IDAY,IHR,IMN,SEC; e.g., IHR=25 is
		# C                 properly interpreted.
		# C          RETURNS:
		# C            SBSLLAT = geographic latitude of subsolar point (degrees)
		# C            SBSLLON = geographic longitude of subsolar point (degrees,
		# C                      between -180 and +180)
		if year.size==1 and dayofyear.size==1 and utseconds.size==1: 
			sbsllat1,sbsllon1 = apex.subsol(year,dayofyear,hour,minute,second)
			sbsllat = numpy.ones_like(alon)*sbsllat1
			sbsllon = numpy.ones_like(alon)*sbsllon1
		elif year.size==alon.size and dayofyear.size==alon.size and utseconds.size==alon.size:
			sbsllat = numpy.zeros_like(alon)
			sbsllon = numpy.zeros_like(alon)
			for k,lon in enumerate(alon):
				sbsllat[k],sbsllon[k] = apex.subsol(year[k],dayofyear[k],hour[k],minute[k],second[k])
		elif year.size==1 and dayofyear.size==1 and utseconds.size==alon.size:
			sbsllat = numpy.zeros_like(alon)
			sbsllon = numpy.zeros_like(alon)
			for k,lon in enumerate(alon):
				sbsllat[k],sbsllon[k] = apex.subsol(year,dayofyear,hour[k],minute[k],second[k])
		else:
			raise ValueError("Length of time arguments must either be 1 or equal to length of longitude argument!\nyear.size=%d\ndayofyear.size=%d\nutseconds.size=%d" % (year.size,dayofyear.size,utseconds.size))
		#Original Comments from magloctm.f:
		#SUBROUTINE MAGLOCTM (ALON,SBSLLAT,SBSLLON,CLATP,POLON,MLT)
		# C  Computes magnetic local time from magnetic longitude, subsolar coordinates,
		# C   and geomagnetic pole coordinates.
		# C  950302 A. D. Richmond, NCAR
		# C  Algorithm:  MLT is calculated from the difference of the apex longitude,
		# C   alon, and the geomagnetic dipole longitude of the subsolar point.
		# C
		# C   Inputs:
		# C    ALON    = apex magnetic longitude of the point (deg)
		# C    SBSLLAT = geographic latitude of subsolar point (degrees)
		# C    SBSLLON = geographic longitude of subsolar point (degrees)
		# C    CLATP   = Geocentric colatitude of geomagnetic dipole north pole (deg)
		# C    POLON   = East longitude of geomagnetic dipole north pole (deg)
		# C
		# C   Output:
		# C    mlt (real) = magnetic local time for the apex longitude alon (hours)
		# C
		# C To go from mlt to alon (see comments following Entry mlt2alon for definition 
		# C  of variables), use:
		# C
		# C     CALL MLT2ALON (MLT,SBSLLAT,SBSLLON,CLATP,POLON,ALON)
		# C
		# C  NOTE: If the calling routine uses subroutine magloctm in conjunction with 
		# C   file magfld.f (which is used by subroutine APEX), then clatp and polon can 
		# C   be found by invoking
		# C
		# C     CALL DYPOL (CLATP,POLON,VP)
		# C
		# C   where vp is an unneeded variable.  (Note that subroutine COFRM must have
		# C   been called before DYPOL, in order to set up the coefficients for the
		# C   desired epoch.)  Alternatively, if subroutine apxntrp is
		# C   used to get alon from previously computed arrays, then
		# C   clatp and polon can be obtained for use in magloctm by adding
		# C
		# C     COMMON /APXDIPL/ CLATP,POLON,DUM1,DUM2,DUM3
		# C
		# C   to the calling routine (where DUM1,DUM2,DUM3 are unneeded dummy variables).
		# C
		mlt = numpy.zeros_like(alon)
		for k in xrange(len(alon)):
			mlt[k] = apex.magloctm(alon[k],sbsllat[k],sbsllon[k],clatp,polon)
		#Sanity check magnetic local time
		mlt[mlt<0] = mlt[mlt<0]+24
		mlt[mlt>24] = numpy.mod(mlt[mlt>24],24)
		return mlt

	def getTransformationResults(self,lat,lon,alt,hr=110.):
		"""
		Converts a set of points to Apex, also saves the results from APXMALL
		into a class parameter so they can be reused.

		Parameters
		----------
		lat : numpy.array
			Geographic latitude in degrees with N positive (length n)
		lon : numpy.array 
			Geographic longitude in degress with E positive (length n)
		alt : numpy.array 
			Altitude of measurement in km (length n)
		hr  : float 
			Modified Apex reference altitude in km (see Richmond,1995 for description of coordinates)

		Saves the following variables to the dictionary self.lastrun
		-Before '=' indicates dictionary key
		-All values are numpy arrays of dimensions [n,1], [n,2] or [n,3]
		b      = magnetic field components (east, north, up), in nT
		bhat   = components (east, north, up) of unit vector along the geomagnetic field direction
		bmag   = magnetic field magnitude, nT
		si     = sin(I) where I is the angle of inclination of the field
		             line below the horizontal
		alon   = Apex longitude = Modified Apex longitude = Quasi-Dipole
		             longitude, degrees
		xlatm  = Modified Apex latitude, degrees
		vmp    = magnetic potential, in T.m
		wm     = Wm of reference above, in km**2 /nT; i.e., 10**15 m**2 /T)
		d     = D of reference above
		be3    = B_e3 of reference above (= Bmag/D), in nT
		sim    = sin(I_m) of reference above
		d1,d2,d3,e1,e2,e3 = Modified Apex coordinates base vectors, each with
		        three components (east, north, up) as described in reference above
		xlatqd = Quasi-Dipole latitude, degrees
		f      = F for Quasi-Dipole coordinates described in reference above
		f1,f2  = Quasi-Dipole coordinates base vectors with components
		             (east, north) described in reference above
		ist    = Return status: 0 (okay) or 1 (failed).
		"""

		lon[lon>180.] = lon[lon>180.]-360.
		lon[lon<-180.] = lon[lon<-180.]+360.
		#if hr < self.altmin or hr > self.altmax:
		#	raise ValueError('Reference height %f is not within altitude limits (%f-%f)!' % (hr,self.altmin,self.altmax)) 
		if not self.identicalInputs(lat,lon,alt,hr):
			npts = len(lat)
			#Init outputs
			# real intent(in) :: glat
			#       real intent(in) :: glon
			#       real intent(in) :: alt
			#       real intent(in) :: hr
			  #       real dimension(lwk),intent(in) :: wk
			  #       integer intent(in) :: lwk
			  #       real dimension(3),intent(out) :: b
			  #       real dimension(3),intent(out) :: bhat
			  #       real intent(out) :: bmag
			  #       real intent(out) :: si
			  #       real intent(out) :: alon
			  #       real intent(out) :: xlatm
			  #       real intent(out) :: vmp
			  #       real intent(out) :: wm
			  #       real intent(out) :: d
			  #       real intent(out) :: be3
			  #       real intent(out) :: sim
			  #       real dimension(3),intent(out) :: d1
			  #       real dimension(3),intent(out) :: d2
			  #       real dimension(3),intent(out) :: d3
			  #       real dimension(3),intent(out) :: e1
			  #       real dimension(3),intent(out) :: e2
			  #       real dimension(3),intent(out) :: e3
			  #       real intent(out) :: xlatqd
			  #       real intent(out) :: f
			  #       real dimension(2),intent(out) :: f1
			  #       real dimension(2),intent(out) :: f2

			#Init scalar variables
			bmag,si,alon,xlatm,vmp,wm,d,be3,sim,xlatqd,f,ist = numpy.zeros([npts,1]),numpy.zeros([npts,1]),numpy.zeros([npts,1]),numpy.zeros([npts,1]),\
														numpy.zeros([npts,1]),numpy.zeros([npts,1]),numpy.zeros([npts,1]),numpy.zeros([npts,1]),\
														numpy.zeros([npts,1]),numpy.zeros([npts,1]),numpy.zeros([npts,1]),numpy.zeros([npts,1])
			b,bhat,d1,d2,d3,e1,e2,e3 = numpy.zeros([npts,3]),numpy.zeros([npts,3]),numpy.zeros([npts,3]),\
									numpy.zeros([npts,3]),numpy.zeros([npts,3]),numpy.zeros([npts,3]),\
									numpy.zeros([npts,3]),numpy.zeros([npts,3])
			f1 = numpy.zeros([npts,2])
			f2 = numpy.zeros([npts,2])
			
			print "Transforming %d points from lat,lon,alt to apex..." % (npts)

			for i in xrange(npts):
				if ist[i] == 0: #Return status == success
					(b[i,:],bhat[i,:],bmag[i],
					 si[i],alon[i],xlatm[i],vmp[i],wm[i],d[i],be3[i],sim[i],
					 d1[i,:],d2[i,:],d3[i,:],e1[i,:],e2[i,:],e3[i,:],
					 xlatqd[i],f[i],f1[i,:],f2[i,:],ist[i]) = apex.apxmall(glat=lat[i],
														glon=lon[i],
														alt=alt[i],
														hr=hr,
														wk=self.workArray)
				else: 
					raise RuntimeError('Call to APXMALL failed at point #%d for unknown reasons (Fortran code does not introspect)' % i)
			
			varnames = ['lat','lon','alt','hr','b','bhat','bmag',
			 'si','alon','xlatm','vmp','wm','d','be3','sim',
			 'd1','d2','d3','e1','e2','e3',
			 'xlatqd','f','f1','f2','ist']

			#Make a dictionary of the output from this run and store it in lastrun, 
			#that way, we can avoid re-running when inputs are identical  
			self.lastrun = dict()
			for name in varnames: 
				self.lastrun[name] = eval(name)
		else:
			print "Inputs identical to last run, results already in attribute lastrun\n"
			

	def geo2apex(self,lat,lon,alt,hr=110.):
		"""
		Does a simple transformation of observation positions from geographic to apex
		
		Parameters
		----------
		lat : numpy.array 
			Observation Geographic latitudes
		lon : numpy.array 
			Observation Geographic longitudes
		alt : numpy.array 
			Observation Altitudes in km
		hr=110. : float
			Modified Apex reference height in km

		Returns
		-------
		alat : numpy.array
			Modified Apex Latitude
		alon : numpy.array
			Apex Longitude
		qdlat : numpy.array
			Quasi-Dipole Latitude
		"""
		self.getTransformationResults(lat,lon,alt,hr=hr)

		alat = self.lastrun['xlatm']
		alon = self.lastrun['alon']
		qdlat = self.lastrun['xlatqd']
		
		return alat,alon,qdlat

	def apex2geo(self,alat,alon,alt,hr=110.):
		"""
		Does a simple transformation of observation positions from modified apex to geodetic
		
		Parameters
		----------
		alat : numpy.array 
			Modified Apex latitudes
		alon : numpy.array 
			Modified Apex longitudes
		alt : numpy.array 
			Altitude to get coordinates for in km
		hr=110. : float
			Modified Apex reference height in km

		Returns
		-------
		glat : numpy.array
			Modified Apex Latitude
		glon : numpy.array
			Apex Longitude
		"""
		gdlats,gdlons = numpy.zeros_like(alat),numpy.zeros_like(alon)
		try:
			n_alts = len(alt)
		except:
			n_alts = 1

		if n_alts == 1:
			alt = numpy.ones_like(alat)*alt

		for i in range(len(alat)):
			(this_gdlat, this_gdlon, this_ist) = apex.apxm2g(xlatm=alat[i],
	                                          alon=alon[i],
	                                          alt=alt[i],
	                                          hr=hr,
	                                          wk=self.workArray,
	                                          lwk=len(self.workArray))
			gdlats[i],gdlons[i] = this_gdlat,this_gdlon

		return gdlats,gdlons

	def measurement2apex(self,lat,lon,alt,v,hr=110.):
		#import pdb
		"""
		Converts a vector measurement and it's associated location to modified apex coordinates
		
		Transforms measurement vector v [m x 3] in geo east, north, up
		to apex. Lat, lon, and alt are [m x 1] vectors of locations association with 
		Measurements in vector v, hr is reference height
		
		Parameters
		----------
		lat : numpy.array 
			geographic latitude
		lon : numpy.array
			geographic longitude
		alt : numpy.array
			geographic altitude
		v : numpy.array 
			[nx3] vector to transform
		hr : float, optional
			Apex reference height in km (default=110.)

		Returns
		-------
		alat : numpy.array
			Modified Apex Latitude of measurement
		alon : numpy.array
			Apex Longitude of measurement
		v_d : numpy.array
			Vector in Apex 'd' basis 

		"""
		self.getTransformationResults(lat,lon,alt,hr=hr)

		alat = self.lastrun['xlatm']
		alon = self.lastrun['alon']
		e1 = self.lastrun['e1']
		e2 = self.lastrun['e2']
		e3 = self.lastrun['e3']

		#Make array with same dimension as lat
		v_d1 = numpy.zeros_like(lat)
		v_d2 = numpy.zeros_like(lat)
		v_d3 = numpy.zeros_like(lat)
		
		for c in numpy.arange(len(v[:,0])):
			v_d1[c] = numpy.dot(v[c,:],e1[c,:])
			v_d2[c] = numpy.dot(v[c,:],e2[c,:])
			v_d3[c] = numpy.dot(v[c,:],e3[c,:])

		v_d = numpy.column_stack((v_d1,v_d2,v_d3))
		#pdb.set_trace()
		return alat,alon,v_d

		