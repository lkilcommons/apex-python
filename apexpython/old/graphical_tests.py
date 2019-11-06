
def graphical_test(self,npts=50,hemisphere='N',testalt=400.,UT=12.,hr=110.,saveit=False,plotdir='/home/liamk/mirror/Projects/apex-python'):
	"""
	Generates a series of 'rings' in geographic coordinates 
	(fixed geographic latitude and altitude, 
	evenly spaced geographic longitude around the globe). 
	By setting the time of 'aquisition' in UT hours, 
	and the Modified Apex reference height hr using named 
	parameter-type arguments, the user can test how data 
	from different geographic locations and times would appear
	in Apex Coordinates. It visualizes the data in a 'dial plot' 
	format (polar coordinates centered on the Modified 
	Apex coordinates pole (plus-minus 90 degrees)),
	with the magnetic local time as the 
	azimuthal coordinate, and the absolute value 
	of latitude as radial. 
	
	Parameters
	----------
	npts : int, optional
		Number of points per ring (default:50)
	hemisphere : str, optional
		Sets sign of latitudes used for rings. N for >0, S for < 0, (default='N')
	testalt : float, optional
		Altitude in km of simulated geographical data (default=400.)
	UT : float, optional
		Universal time of aquisition of simulated 
		geographical data in hours (default=12.)
	hr : float, optional
		Modified Apex reference height in km (default=110.) 
	saveit : bool, optional 
		Save a copy of the plot as a png file (defualt=False)
	plotdir : str, optional
		Where you would like the png to save.

	Returns
	-------
	ax : matplotlib.axes
		Axes object on which plot was displayed

	"""
	import satplottools
	UT=numpy.mod(UT,24.)

	lon_range = numpy.linspace(-180,180,num=npts)

	#ring = numpy.transpose(numpy.vstack((numpy.ones(npts),lon_range)))
	#Above is longhand version of:
	ring = numpy.column_stack((numpy.ones(npts),lon_range))
	
	ringlats = [80.,70.,60.,50.]

	lats = []
	lons = []

	if hemisphere == 'S':
		ringlats = [-1*lat for lat in ringlats]

	for r in ringlats:
		thisring = ring.copy()
		thisring[:,0] = thisring[:,0]*r
		lats.append(thisring[:,0])
		lons.append(thisring[:,1])

	lats = numpy.hstack(lats)
	lons = numpy.hstack(lons)

	year = numpy.floor(self.epoch)
	dayofyear = 1

	alts = numpy.ones_like(lats)*testalt
	alats,alons,qdlats = self.geo2apex(lats,lons,alts,hr=hr)
	mlts = self.alon2mlt(alons,year,dayofyear,UT*3600.)

	r,theta = satplottools.latlt2polar(alats,mlts,hemisphere=hemisphere)
	X = r*numpy.cos(theta)
	Y = r*numpy.sin(theta)

	f = pp.figure()
	ax = pp.axes()
	satplottools.draw_dialplot(ax,padding=0.)
	ax.plot(X,Y,'b.')
	pp.show()
	f.suptitle('%s Hemisphere\nApex Latitude And MLT of rings of constant GEO lat\n Observation Altitude %dkm\n Reference Altitude %dkm\n Time: Jan 1, %d, %d:00:00 UT' % (hemisphere,testalt,hr,year,UT),fontsize=9)
	
	aspectratio=1.0
	ratio_default=(ax.get_xlim()[1]-ax.get_xlim()[0])/(ax.get_ylim()[1]-ax.get_ylim()[0])
	ax.set_aspect(ratio_default*aspectratio)

	if saveit:
		f.savefig(os.path.join(plotdir,'graphical_test_%s_%dkm_%dUT.png' % (hemisphere,testalt,UT)))

	return ax

#TEST CODE TO VALIDATE TRANSFORMATION WORKS, REQUIRES OUTSIDE DATA AND LIBRARIES
def testAgainstOldDMSP(self):
	import satplottools
	import readDMSPAmpr as readr
	import pdb
	#we will use minutes 61-181 of satellite F15 data from may 29 2010 as our test set of data
	#!! The apex_coverter needs to be initialized with an epoch date of 2010.0 and an nvert value of 40 in order for this test,
	# to be valid !!
	year = 2010.
	dayofyear = 149.

	geo_data = readr.readGeoDMSP(60,180,satnum=15,dayofmay=29)
	#Columns = UT second, Glat, Glon, Alt, dB_east, dB_north
	
	ref_apx_data = readr.readDMSPRange(60,180)
	ref_apx_data = ref_apx_data['DMSP-F15']
	ref_cols = [1,2,3,4,5,6,7] #time,lat,lon,mlt,db_d1,db_d2,db_d3
	ref_apx_data = ref_apx_data[:,ref_cols]


	dB_up = numpy.zeros_like(geo_data[:,1])
	dB_geo = numpy.column_stack((geo_data[:,4],geo_data[:,5],dB_up))

	alat,alon,dB_apx = self.measurement2apex(geo_data[:,1],geo_data[:,2],geo_data[:,3],dB_geo)
	#amlt = geo_data[:,0]/3600.+(alon-70.)/15 #Approximation
	amlt = self.alon2mlt(alon,year,dayofyear,geo_data[:,0])

	apx_data = numpy.column_stack((geo_data[:,0],alat,alon,amlt,dB_apx))

	f1 = pp.figure()
	ax1 = pp.axes()
	ax1.plot(apx_data[:,0],apx_data[:,4],'ro',label='PythonApex dB_d1')
	ax1.plot(ref_apx_data[:,0],ref_apx_data[:,4],'b.',label='MatsuoFortranApex dB_d1')
	ax1.set_title('Python Apex Validation, dB_d1')
	ax1.legend()
	
	f2 = pp.figure()
	ax2 = pp.axes()
	ax2.plot(apx_data[:,0],apx_data[:,5],'ro',label='PythonApex dB_d2')
	ax2.plot(ref_apx_data[:,0],ref_apx_data[:,5],'b.',label='MatsuoFortranApex dB_d2')
	ax2.set_title('Python Apex Validation, dB_d2')
	ax2.legend()

 	f3 = pp.figure()
	ax3 = pp.axes()
	ax3.plot(apx_data[:,0],apx_data[:,6],'ro',label='PythonApex dB_d3')
	ax3.plot(ref_apx_data[:,0],ref_apx_data[:,6],'b.',label='MatsuoFortranApex dB_d3')
	ax3.set_title('Python Apex Validation, dB_d3')
	ax3.legend()

	f4 = pp.figure()
	ax4 = pp.axes()
	ax4.plot(apx_data[:,0],apx_data[:,1],'ro',label='PythonApex apex_latitude')
	ax4.plot(ref_apx_data[:,0],ref_apx_data[:,1],'b.',label='MatsuoFortranApex apex latitude')
	ax4.set_title('Python Apex Validation, apex latitude')
	ax4.legend()

	f5 = pp.figure()
	ax5 = pp.axes()
	ax5.plot(apx_data[:,0],apx_data[:,2],'ro',label='PythonApex apex_longitude')
	ax5.plot(ref_apx_data[:,0],ref_apx_data[:,2],'b.',label='MatsuoFortranApex apex longitude')
	ax5.set_title('Python Apex Validation, apex longitude')
	ax5.legend()

	f6 = pp.figure()
	ax6 = pp.axes()
	ax6.plot(apx_data[:,0],apx_data[:,3],'ro',label='PythonApex apex_MLT')
	ax6.plot(ref_apx_data[:,0],ref_apx_data[:,3],'b.',label='MatsuoFortranApex apex MLT')
	ax6.set_title('Python Apex Validation, apex longitude')
	ax6.legend()


	pp.show()
	pdb.set_trace()
	return ax1, ax2