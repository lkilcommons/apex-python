"""
apexVisualizer.py
Tool to visualize modified apex coordinates.
"""     

from mpl_toolkits.basemap import Basemap # Maps
#from mayavi import mlab #3-d Visualization

import sys
sys.path.append('/home/liamk/mirror/superdarn/davitpy')
import models.igrf as igrf #use the davitpy wrapper on the IGRF
import numpy as np
from numpy import cos,sin

sys.path.append('/home/liamk/mirror/Projects/conjunction/src') #Add the conjunction search source to the system path
													#This is needed so the code can find satplottools
sys.path.append('/home/liamk/mirror/Projects/apex-python/src') #This is so the apex converter and sathat code can be found

import apex_converter
import datetime
import pdb



def drawApexContours(year,dayofyear,secofday,hemi='N'):
	"""
	Display apex latitude and longitude lines on an orthographic
	projection of the globe. 

	Parameters
	----------
		year : int
		dayofyear : int
		secofday : int
			Date and time for apex transform.
		satloc : [lat,lon] - int
			Geographic Latitude and Longitude for POV of projection
			(perspective of looking down from a satellite at satloc)

	Returns
	-------
		map - Basemap instance for the orthographic projection	
	"""

	#Initialize Apex Conversion
	#--------------------------

	#Figure out the epoch for the apex converter
	if np.mod(year,4)==0:
		epoch = year+dayofyear/366.
	else:
		epoch = year+dayofyear/365.

	#Apex converter defaults:
	#glatmin=-90.,glatmax=90.,glonmin=-180.,glonmax=180.,altmin=300., altmax=900.,
	#	nvert=50.,epoch=2010.0
	converter = apex_converter.apex_converter(epoch=epoch)

	#Get subsolar point and lat/lon of the geomagnetic dipole north pole
	subsol_lat,subsol_lon,dnp_lat,dnp_lon = converter.solarParams(year,dayofyear,secofday)

	#Compute the datetime for this year, dayofyear, secofday
	#Used by nightshade
	# set up orthographic map projection with
	# perspective of satellite looking down at 50N, 100W.
	# use low resolution coastlines.
	doy2datetime = lambda doy,year,secofday: datetime.datetime(year,1,1,0,0,0)+datetime.timedelta(days=doy)+datetime.timedelta(seconds=secofday)
	time = doy2datetime(year,dayofyear,secofday)


	#Set up the map
	#--------------

	#map = Basemap(projection='geos',lat_0=45,lon_0=-100,resolution='l')
	if hemi=='N':
		satloc = [60,subsol_lon+180.]
	elif hemi=='S':
		satloc = [-60,subsol_lon+180.]

	map = Basemap(projection='ortho',lat_0=satloc[0],lon_0=satloc[1],resolution='l')
	
	# draw coastlines, country boundaries, fill continents.
	map.drawcoastlines(linewidth=0.25)
	map.drawcountries(linewidth=0.25)
	map.fillcontinents(color='coral',lake_color='aqua')
	# draw the edge of the map projection region (the projection limb)
	map.drawmapboundary(fill_color='aqua')
	# draw lat/lon grid lines every 30 degrees.
	map.drawmeridians(np.arange(0,360,30))
	map.drawparallels(np.arange(-90,90,30))
	# make a regular lat/lon grid.
	nlats = 73; nlons = 145; delta = 2.*np.pi/(nlons-1)
	lats = (0.5*np.pi-delta*np.indices((nlats,nlons))[0,:,:])*180./np.pi
	lons = (delta*np.indices((nlats,nlons))[1,:,:])*180./np.pi
	alts = np.ones_like(lats)*800. #altitudes of 800 km

	#Get the apex latitude and longitude for all grid points at an altitude of 800 km for a reference height of 110 km
	alat,alon,qdlat = converter.geo2apex(lats.flatten(),lons.flatten(),alts.flatten(),hr=110.)

	#Get the magnetic local time for the grid
	mlt = converter.alon2mlt(alon,year,dayofyear,secofday)

	#Create grid-shaped versions
	alats = np.reshape(alat,(nlats,nlons))
	alons = np.reshape(alon,(nlats,nlons))
	mlts = np.reshape(mlt,(nlats,nlons))

	# compute native map projection coordinates of lat/lon grid.
	x, y = map(lons, lats)

	# draw the terminator
	map.nightshade(time)

	# contour apex latitude data over the map in blue
	if hemi=='N':
		lat_lines = np.arange(0,90,10.)
	elif hemi=='S':
		lat_lines = np.arange(-80,0,10.)

	lt_lines = np.arange(0,24,6.)
	
	cs_lat = map.contour(x,y,alats,lat_lines,linewidths=1.5,colors='b')
	pp.clabel(cs_lat, fontsize=9, inline=1)
	
	# contour apex longitude data over the map in red
	#cs_lon = map.contour(x,y,mlts,lt_lines,linewidths=1.5,colors='r')
	#pp.clabel(cs_lon, fontsize=9, inline=1)

	# Plot the geomagnetic north pole
	#map.plot(dnp_lat,dnp_lon,'go',latlon=True,label='Geomangetic North Pole')

	pp.show()
	return map

def gridIGRF(year,dayofyear,altmin,altmax,altstep=20.,geodetic=False):
	#import pdb
	#-----------------------------------------------
	#Make a spherical grid of IGRF main field values
	#-----------------------------------------------

	#Some stuff I stole from a DavitPy tutorial
	if geodetic: 
		itype = 1 # Geodedic coordinates
	else: 
		itype = 2 # Geocentric coordinates

	stp = 5. 
	xlti, xltf, xltd = -90.,90.,stp # latitude start, stop, step
	xlni, xlnf, xlnd = -180.,180.,stp # longitude start, stop, step
	
	#All altitudes
	alts = np.arange(altmin,altmax,altstep)

	#Dimensions so we can define array size
	nlats = int((xltf-xlti)/stp)+1
	nlons = int((xlnf-xlni)/stp)+1
	nalts = int((altmax-altmin)/stp)

	ifl = 0 # Main field

	#Figure out the epoch for the IGRF
	if np.mod(year,4)==0:
		epoch = year+dayofyear/366.
	else:
		epoch = year+dayofyear/365.

	#Initialize 
	gridshape = (nlats,nlons,nalts)
	
	#Not using spherical values currently
	#lats = np.zeros(gridshape)
	#lons = np.zeros(gridshape)
	#alts = np.zeros(gridshape)
	#B_east = np.zeros(gridshape)
	#B_north = np.zeros(gridshape)
	#B_up = np.zeros(gridshape)

	#Initialize cartesian coordinate arrays
	X,Y,Z = np.zeros(gridshape),np.zeros(gridshape), np.zeros(gridshape)
	BX,BY,BZ,B = np.zeros(gridshape),np.zeros(gridshape), np.zeros(gridshape), np.zeros(gridshape)

	for k,alt in enumerate(alts):
		#Explaination of IGRF inputs/output from original fortran source

		# C Edits to switch IGRF to subroutine
		# C   INPUTS:
		# C       - ITYPE:
		# C           - 1 - geodetic (shape of Earth is approximated by a spheroid)
		# C           - 2 - geocentric (shape of Earth is approximated by a sphere)
		# C       - DATE: date in years A.D. (ignored if IOPT=2)
		# C       - ALT: altitude or radial distance in km (depending on ITYPE)
		# C       - XLTI,XLTF,XLTD: latitude (initial, final, increment) in decimal degrees
		# C       - XLNI,XLNF,XLND: longitude (initial, final, increment) in decimal degrees
		# C       - IFL: value for MF/SV flag:
		# C           - 0 for main field (MF)
		# C           - 1 for secular variation (SV)
		# C           - 2 for both
		# C   OUTPUTS:
		# C       - aLat is the latitude of each point
		# C       - aLon is the longitude of each point
		# C       - D is declination in degrees (+ve east)
		# C       - I is inclination in degrees (+ve down)
		# C       - H is horizontal intensity in nT
		# C       - X is north component in nT
		# C       - Y is east component in nT
		# C       - Z is vertical component in nT (+ve down)
		# C       - F is total intensity in nT
			
		# Call fortran subroutine
		lat,lon,d,s,h,be,bn,bu,f = igrf.igrf11(itype,epoch,alt,ifl,xlti,xltf,xltd,xlni,xlnf,xlnd)
		
		print('Finished iteration %d.' % k)
		#Assign values to master arrays
		#lats[:,:,k]=lat
		#lons[:,:,k]=lon
		#B_east[:,:,k]=be
		#B_north[:,:,k]=bn
		#B_up[:,:,k]=bu

		#Everything comes out 1 dimensional from igrf11
		#so reshape to be the size of the grid
		be = np.reshape(be,(nlats,nlons))
		bn = np.reshape(bn,(nlats,nlons))
		bu = np.reshape(bu,(nlats,nlons))
		lon = np.reshape(lon,(nlats,nlons))
		lat = np.reshape(lat,(nlats,nlons))
		
		R = alt
		theta = lon/180.*np.pi
		phi = (90.-lat)/180*np.pi

		#Grab the main field intesity(nT)
		B[:,:,k] = np.reshape(f,(nlats,nlons))

		#Transform lat,lon,alt to XYZ, (spherical to cartesian)
		X[:,:,k] = R*cos(theta)*sin(phi)
		Y[:,:,k] = R*sin(theta)*sin(phi)
		Z[:,:,k] = R*cos(phi) 

		#Main Field to Cartesian Using standard spherical unit vector transformation
		#mathworld.wolfram.com/SphericalCoordinates.html
		BX[:,:,k] = cos(theta)*sin(phi)*bu + -1.*sin(theta)*be + cos(theta)*cos(phi)*bn
		BY[:,:,k] = sin(theta)*sin(phi)*bu + cos(theta)*be + sin(theta)*cos(phi)*bn
		BZ[:,:,k] = cos(phi)*bu + 0.*be + -1.*sin(phi)*bn

 	
	return X,Y,Z,BX,BY,BZ,B
	
	
	
def plotGlobe3D():
	"""
	Generates a Mayavi 3d model of the globe with the continents
	drawn on the surface.

	Parameters
	----------

	Returns
	-------

	"""