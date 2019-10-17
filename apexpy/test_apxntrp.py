# custom
import apex

# 3rd-party
import numpy

# standard
import unittest

class ApxntrpTest(unittest.TestCase):
    """
    Unit tests for primary Apex routines (implemented in apxntrp.f90).
    Validated Python bindings by comparing results to test_apxntrp.f
    """
    def setUp(self):
        # due to a shortcoming in f2py (documented here http://cens.ioc.ee/pipermail/f2py-users/2008-December/001764.html), 
        # we must know dimensions apriori.   Hard-coding them here:
        self.nlat=13
        self.nlon=21
        self.nalt=2

        self.workArray = numpy.empty(self.nlat*self.nlon*self.nalt*6, numpy.float32)
        # Define grid points suitable for preparing interpolation
        # tables for Apex coordinates, Modified Apex Coordinates and
        # Quasi-Dipole coordinates
        (gplat,gplon,gpalt) = apex.ggrid(nvert=4,
                                         glamn=-90.0, glamx=90.0,
                                         glomn=-180.0,   glomx=180.0,
                                         altmn=300.0,   altmx=300.0,
                                         nlat=self.nlat, nlon=self.nlon,  nalt=self.nalt)

        # Create interpolation tables for a single time
        ist = apex.apxmka(msgun=6,
                          epoch=1994.5,
                          gplat=gplat, gplon=gplon, gpalt=gpalt,
                          wk=self.workArray,
                          nlat=self.nlat,   nlon=self.nlon,   nalt=self.nalt)
        self.assertEqual(ist, 0)

    def test_APXGGC(self):
        """ Get grid coordinates for tables currently in memory """
        (gplat, gplon, gpalt,ist) = apex.apxggc(6, self.workArray, self.nlat,self.nlon,self.nalt,)
        
        self.assertEquals(0, ist)
        
        for i,val in enumerate(gplat):
            self.assertAlmostEqual( -90.0 + 15.0*i, val )

        for i,val in enumerate(gplon):
            self.assertAlmostEqual( -180.0 + 18.0*i, val, places=4 )

        self.assertAlmostEqual(0.0, gpalt[0], places=2)
        self.assertAlmostEqual(2123.70507812, gpalt[1], places=2)


    def test_APXALL(self):
        """ Geographic to Apex """
        (a,alat,alon, ist) = apex.apxall(glat=-18.3,
                                         glon=-66.75,
                                         alt=300.0,
                                         wk=self.workArray,
                                         nlat=self.nlat,
                                         nlon=self.nlon,
                                         nalt=self.nalt)

        self.assertAlmostEqual(1.0594660043716431, a, places=5)
        self.assertAlmostEqual(-13.704476356506348, alat, places=4)
        self.assertAlmostEqual(3.9752681255340576, alon, places=5)
        self.assertEqual(0, ist)

    def test_APXA2G(self):
        """ Apex to Geographic (works with global grid only) """
        #FIXME:  A2G results don't look quite right.  Perhaps the grid is too low resolution?
        
        (gdlat, gdlon, ist) = apex.apxa2g(alat=-13.840052604675293,
                                          alon=3.7327430248260498,
                                          alt=300.0,
                                          wk=self.workArray,
                                          lwk=len(self.workArray))
        
        self.assertAlmostEqual(-5.1449408531188965, gdlat, places=4)
        self.assertAlmostEqual(-68.195671081542969, gdlon, places=4)
        self.assertEqual(0, ist)


    def test_APXM2G(self):
        """ Modified Apex to Geographic (works with global grid only) """
        # Modified Apex coordinates obtained via APXMALL on (glat=-19.0, glon=-67.0)
        (gdlat, gdlon, ist) = apex.apxm2g(xlatm=-13.840052604675293,
                                          alon=3.7327430248260498,
                                          alt=300,
                                          hr=8.5,
                                          wk=self.workArray,
                                          lwk=len(self.workArray))

        self.assertAlmostEqual(-4.825462818145752, gdlat, places=4)
        self.assertAlmostEqual(-68.226242065429688, gdlon, places=4)
        self.assertEqual(0, ist)
                                          

    def test_APXQ2G(self):
        """ Quasi-Dipole Apex to Geographic (works with global grid only) """
        # Quasi-Dipole Apex coordinates obtained via APXMALL on (glat=-19.0, glon=-67.0)
        (gdlat,gdlon, ist) = apex.apxq2g(qdlat=-6.8306035995483398,
                                         qdlon=3.7327430248260498,
                                         alt=300.0,
                                         wk=self.workArray,
                                         lwk=len(self.workArray))

        self.assertAlmostEqual(-18.999721527099609, gdlat, places=3)
        self.assertAlmostEqual(-67.000106811523438, gdlon, places=4)

        self.assertEqual(0, ist)


    def test_APXMALL(self):
        """ Determine Modified Apex and Quasi-dipole coordinates using tables currently in memory """
        (b,bhat,bmag,
         si,alon,xlatm,vmp,wm,d,be3,sim,
         d1,d2,d3,e1,e2,e3,
         xlatqd,f,f1,f2,ist) = apex.apxmall(glat=-19.0,
                                            glon=-67.0,
                                            alt=300,
                                            hr=8.5,
                                            wk=self.workArray)

        self.assertAlmostEqual( -2117.1929, b[0],places=1)
        self.assertAlmostEqual( 20533.721, b[1], places=1)
        self.assertAlmostEqual( 4550.9702, b[2], places=1)

        self.assertAlmostEqual( -0.10015912, bhat[0], places=6)
        self.assertAlmostEqual( 0.97139913, bhat[1], places=6)
        self.assertAlmostEqual( 0.21529505, bhat[2], places=5)

        #FIXME: calcualte this?
        self.assertAlmostEqual( 21138.294921875, bmag, places=1)

        self.assertAlmostEqual(-0.21529504656791687, si, places=5)
        self.assertAlmostEqual(3.7327430248260498, alon, places=5)
        self.assertAlmostEqual(-13.840052604675293, xlatm, places=4)
        self.assertAlmostEqual(18.23223876953125, vmp, places=4)
        self.assertAlmostEqual(1101.16796875, wm, places=2)
        self.assertAlmostEqual(0.7504001259803772, d, places=3)
        self.assertAlmostEqual(28169.365234375, be3, places=1)
        self.assertAlmostEqual(-0.44198882579803467, sim, places=5)

        self.assertAlmostEqual(0.84624624, d1[0], places=5)
        self.assertAlmostEqual(0.08271649, d1[1], places=5)
        self.assertAlmostEqual(0.02047677, d1[2], places=5)

        self.assertAlmostEqual(-0.02326523, d2[0], places=5)
        self.assertAlmostEqual(0.1886366, d2[1], places=5)
        self.assertAlmostEqual(-0.8619411, d2[2], places=5)

        self.assertAlmostEqual(-0.13347428, d3[0], places=5)
        self.assertAlmostEqual(1.29450822, d3[1], places=5)
        self.assertAlmostEqual(0.28690699, d3[2], places=5)

        self.assertAlmostEqual(1.16991103, e1[0], places=5)
        self.assertAlmostEqual(0.12172192, e1[1], places=5)
        self.assertAlmostEqual(-0.0049389, e1[2], places=5)

        self.assertAlmostEqual(0.00277541, e2[0], places=5)
        self.assertAlmostEqual(0.24552709, e2[1], places=5)
        self.assertAlmostEqual(-1.10651326, e2[2], places=5)

        self.assertAlmostEqual(-0.07515942, e3[0], places=5)
        self.assertAlmostEqual(0.72893804, e3[1], places=5)
        self.assertAlmostEqual(0.16155744, e3[2], places=5)

        self.assertAlmostEqual(-6.8306035995483398, xlatqd, places=3)

        self.assertAlmostEqual(0.87059986591339111, f, places=5)

        self.assertAlmostEqual(0.94297397, f1[0], places=5)
        self.assertAlmostEqual(0.19558513, f1[1], places=4)

        self.assertAlmostEqual(-0.08844995, f2[0], places=5)
        self.assertAlmostEqual(0.90490347, f2[1], places=5)
        
        self.assertEqual(0,ist)

#class ApxntrpHelpersTest(unittest.TestCase):
#    def test_intrpsc(self):
#        (fx,fy,fz,ist) = apex.intrpsc(glat=-19.0,
#                                      glon=-67.0,
#                                      alt=300.0,
#                                      x=[[[0.0]]],
#                                      y=[[[0.0]]],
#                                      z=[[[0.0]]],
#                                      nlat=1,
#                                      nlon=1,
#                                      nalt=1,
#                                      gplat=[0.0],
#                                      gplon=[0.0],
#                                      gpalt=[0.0],
#                                      calnm='intrpsc')        
#        self.assertAlmostEqual(0.0, fx)
#        self.assertAlmostEqual(0.0, fy)
#        self.assertAlmostEqual(0.0, fz)
#
#        self.assertAlmostEqual(1.0, ist)

class FosterTest(unittest.TestCase):
    def test_apxq2g_legacy(self):        
        nalt=2
        gpalt=numpy.empty(nalt, numpy.float32)
        gpalt[0] = 90.0
        gpalt[1] = 170.0

        nlat=36
        nlatp1=nlat+1
        gplat = numpy.empty(nlatp1, numpy.float32)
        dellat = 180.0 / float(nlat)
        for j in range(nlatp1):
            gplat[j] = j*dellat - 90.0

        nlon = 72
        nlonp2 = nlon+2
        dellon = 360.0 / float(nlon)
        gplon = numpy.empty(nlonp2, numpy.float32)
        for i in range(nlonp2):
            gplon[i] = (float(i)-0.5) * dellon - 180.0

        workArray = numpy.empty(nlatp1*nlonp2*nalt*6, numpy.float32)

        ist = apex.apxmka(msgun=6, epoch=2002,
                          gplat=gplat, gplon=gplon, gpalt=gpalt,
                          wk=workArray,
                          nlat=nlatp1, nlon=nlonp2, nalt=nalt)
        self.assertEqual(ist, 0)

        # gdlat, gdlon, alt values obtained from Ben Foster example
        # which q2g would print error:
        #
        #     Warning 10 iterations only reduced the angular difference to
        #              0.00055 degrees ( 0.00050 degrees is the test criterion)
        (gdlat, gdlon, ist) = apex.apxq2g(qdlat = -29.9158214580855,
                                          qdlon=18.0000000000000,
                                          alt=130.0,
                                          wk=workArray,
                                          lwk=len(workArray))

        print(gdlat, gdlon)
        print('lat:',gplat.min(), gplat.max())
        print('lon:',gplon.min(), gplon.max())
        print('alt:',gpalt.min(), gpalt.max())

def main():
  #Added by LMK
  #unittest.main()
  suite = unittest.TestLoader().loadTestsFromTestCase(ApxntrpTest)
  unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == '__main__':
    main()

