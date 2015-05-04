#sathat.py 
#estimates the direction of the along and cross track unit vectors
#of a spacecraft in geocentric (earth-centered spherical) coordinates 
#given only the timestamps in UTC seconds and two latitude and longitude pairs
#and assuming no altitude change 
#The algorithm is based on spherical trig
#The returned vectors are estimated for time ut2
from numpy import * 
def sathat(ut1,ut2,lat1,lat2,lon1,lon2):

  #Correct for earth rotation
  corrlon1 = lon1 - (ut2-ut1)/240.

  arc_a = (90.-lat1)*pi/180. #colat of point 1 in rad
  arc_b = (90 -lat2)*pi/180. #colat of point 2 in rad
  ang_C = (lon2-corrlon1)*pi/180. #angle subtended by arc connecting points 1 and 2 
  
  #Do some spherical trig
  cos_arc_c = cos(arc_a)*cos(arc_b) + sin(arc_a)*sin(arc_b)*cos(ang_C)
  sin_arc_c = sqrt(1.-cos_arc_c**2)
  cos_ang_A = (cos(arc_a)*sin(arc_b) - sin(arc_a)*cos(arc_b)*cos(ang_C))/sin_arc_c
  sin_ang_A = sin(arc_a)*(sin(ang_C)/sin_arc_c)
  cos_ang_B = (cos(arc_b)*sin(arc_a)-sin(arc_b)*cos(arc_a)*cos(ang_C))/sin_arc_c
  sin_ang_B = sin(arc_b)/sin_arc_c*sin(ang_C)

  s_along = array([sin_ang_B,cos_ang_B,0]) #Along track
  s_cross = array([-1*cos_ang_B,sin_ang_B,0]) #Cross track (+ to the left)

  return s_along,s_cross

def sathats(ut_times,lats,lons):
  s_alongs = zeros([len(lats),3])
  s_crosses = zeros([len(lats),3])
  for k in arange(len(lats)-1):
    s_alongs[k+1,:],s_crosses[k+1,:] = sathat(ut_times[k],ut_times[k+1],lats[k],lats[k+1],lons[k],lons[k+1])
  return s_alongs, s_crosses

# Original fortran code
#       subroutine sathat(ut1,ut2,glat1,glon1,glat2,glon2, 
#      +                   s1,s2)
# C     Input/Output Argument Initializations
#       CHARACTER LOGFILE*40
#       REAL ut1, ut2
#       REAL pi,glat1,glat2,glon1,glon2,s1(3),s2(3)

# C     Internal variable initializations	
#       REAL sin_ang_A, cos_ang_A, sin_ang_B, cos_ang_B, 
#      + ang_C, sin_ang_C, cos_ang_C
#      + arc_a, sin_arc_a, cos_arc_a, 
#      + arc_b, sin_arc_b, cos_arc_b
#      + sin_arc_c, cos_arc_c,
#      + corr_glon1

# C     calculates unit vectors along-track (s1) and cross-track (s2) at point 2
# C     (where point 2 = glat2,glon2 at time ut2)
#       corr_glon1 = glon1 - (ut2-ut1)/240.
# C      WRITE(6,'(''glat1,   glat2,  glon1, corr_glon1,  glon2'')')
# C      WRITE(6,*)glat1,glat2,glon1,corr_glon1,glon2   
      
#       PI = 3.141592653
        
#       arc_a = (90.-glat1)*PI/180.
#       arc_b = (90.-glat2)*PI/180.
#       ang_C = (glon2-corr_glon1)*PI/180.
# C      WRITE(6,'(''arc_a   ,arc_b   ang_C'')')
# C      WRITE(6,*) arc_a, arc_b, ang_C
      
# C     Sines and Cosines of inital quantities
#       sin_arc_a = sin(arc_a)
#       cos_arc_a = cos(arc_a)
#       sin_arc_b = sin(arc_b)
#       cos_arc_b = cos(arc_b)
#       sin_ang_C = sin(ang_C)
#       cos_ang_C = cos(ang_C)

# C      WRITE(6,'(''sin_arc_a, cos_arc_a, sin_arc_b, cos_arc_b'')')
# C      WRITE(6,*) sin_arc_a,cos_arc_a,sin_arc_b,cos_arc_b
# C      WRITE(6,'(''sin_ang_C,  cos_ang_C'')')
# C      WRITE(6,*) sin_ang_C, cos_ang_C

# C     Derived Quantities
#       cos_arc_c = cos_arc_a*cos_arc_b + sin_arc_a*sin_arc_b*cos_ang_C
#       sin_arc_c = sqrt(1.-cos_arc_c*cos_arc_c)
#       cos_ang_A = (cos_arc_a*sin_arc_b - sin_arc_a*cos_arc_b*cos_ang_C)/
#      +sin_arc_c
#       sin_ang_A = sin_arc_a*(sin_ang_C/sin_ang_c)
#       cos_ang_B = (cos_arc_b*sin_arc_a-sin_arc_b*cos_arc_a*cos_ang_C)/
#      +sin_arc_c
#       sin_ang_B = sin_arc_b/sin_arc_c*sin_ang_C
# C      WRITE(6,'(''cos_arc_c, sin_arc_c,cos_B,sin_B,cos_A,sin_A'')')
# C      WRITE(6,*) cos_arc_c, sin_arc_c, cos_ang_B, sin_ang_B,
# C     +cos_ang_A,sin_ang_A

# C     Unit Vectors
#       s1(1) =  sin_ang_B
#       s1(2) =  cos_ang_B
#       s1(3) = 0 !assuming no altitude change
#       s2(1) = -1*cos_ang_B
#       s2(2) =  sin_ang_B
#       s2(3) = 0 !assuming no altitude change

# C      WRITE(6,'(''UT Times: '',I6,'','',I6)') ut1,ut2
# C      WRITE(6,'(''s1 = <'',F10.4,'','',F10.4,''>'')') s1(1), s1(2)
# C      WRITE(6,'(''s2 = <'',F10.4,'','',F10.4,''>'')') s2(1), s2(2)
#       return
 
#       end