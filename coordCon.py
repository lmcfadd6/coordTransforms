# Imports
import numpy as np

# Constants
JUL_DAYS_BEFORE_JAN_1_2000 = 2451545.0
DAYS_PER_YEAR = 365.25
YEARS_PER_CENTURY = 100
T2ANG_COEFFS = [280.46061837, 360.98564736629, 0.0003879332, 38710000]
DEG_PER_HOUR = 15

# Earth constants
EARTH_a = 6378137.0
EARTH_b = 6356752.314245


class Angle:
    """
    Angle object to easilt convert between radians and degrees

    Input:
    ang [float] - angle in radians
    deg [boolean] - if True, angle is given in degrees. If False (default), angle is given in radians
    Example:

        a = Angle(np.pi/2)
        or
        a = Angle(90, deg=True)
    """

    def __init__(self, ang, deg=False):

        if deg:
            self.deg = ang%360
            self.rad = self.deg/180*np.pi
        else:
            self.rad = ang%(2*np.pi)
            self.deg = self.rad*180/np.pi

    def __str__(self):

        return "Angle: {:.2f}".format(self.deg)

    def __add__(self, other):

        return Angle(self.rad + other.rad)

    def __sub__(self, other):

        return Angle(self.rad - other.rad)

    def unmod(self, deg=False):
        """ Retruns angle in range -180 < x < 180 instead of 0 < x < 360
        """

        if deg:
            if self.deg < 180:
                return self.deg
            return self.deg - 360
        else:
            if self.rad < np.pi:
                return self.rad
            return self.rad - 2*np.pi


class RightAsc:

    """ Quick object to convert an angle in degrees to a right ascension

    input: 
    angle [float] - angle in degrees
    Example:

        a = RightAsc(90)
        print(a)
        >> Right Ascension: 6.0h 0.0m 0.00s

    """
    def __init__(self, angle):

        self.angle = angle

        total_hours = angle/DEG_PER_HOUR

        self.hour, r = divmod(total_hours, 1)
        self.min, r = divmod(r*60, 1)
        self.sec = r*60


    def __str__(self):

        return "Right Ascension: {:}h {:}m {:.2f}s".format(self.hour, self.min, self.sec)

    def asFloat(self):
        return self.angle

class Cart:
    """
    A Cart(esian) object takes a geocentric vector defined as [x, y, z] [in meters], 
    and calculates various parameters from it
    """

    def __init__(self, x, y, z):

        try:
            self.x = float(x)
            self.y = float(y)
            self.z = float(z)
        except ValueError:
            print("[WARNING] Cartesian values must be a float")
            return None

        self.xyz = [self.x, self.y, self.z]
        self.h = (x*x + y*y)**0.5
        self.r = (x*x + y*y + z*z)**0.5

        self.phi =   Angle(np.arctan2(self.y, self.x))
        self.theta = Angle(np.arctan2(self.h, self.z))

        self.lat = Angle(np.arctan2(EARTH_a**2*self.z, EARTH_b**2*self.h))
        self.lon = Angle(np.arctan2(self.y, self.x))

        rho = self.r/EARTH_a
        phiprime = Angle(np.arctan2(self.z, self.h))
        u = Angle(np.arctan2(EARTH_a*self.z, EARTH_b*self.h))

        self.alt = EARTH_a*(rho*np.cos(phiprime.rad) - np.cos(u.rad))/np.cos(self.lat.rad)

    def __str__(self):

        return "Cart Obj: x={:.2f} y={:.2f} z={:.2f}".format(self.x, self.y, self.z)

    def geoPos(self, pnt=True):
        """ returns the position in lat/lon/alt instead of x/y/z
        """

        if pnt:
            return "Cart Obj: lat={:.4f}N lon={:.4f}E alt={:.2f}km".format(self.lat.deg, self.lon.deg, self.alt/1000)
        else:
            return [self.lat.unmod(deg=True), self.lon.unmod(deg=True), self.alt]

    def rotate(self, ang, axis):
        """ Rotates vector <ang> degrees around an axis
            inputs:
            ang [Angle Obj] - angle to rotate coordinate system by
            axis ["x", "y", or "z"] - axis to rotate vector around 
        """
        
        if axis == "x":
            M = np.array([[1,          0,               0     ], \
                          [0, np.cos(ang.rad), np.sin(ang.rad)], \
                          [0, -np.sin(ang.rad), np.cos(ang.rad)]])
        elif axis == "y":
            M = np.array([[np.cos(ang.rad), 0, -np.sin(ang.rad)], \
                          [0, 1, 0], \
                          [np.sin(ang.rad), 0, np.cos(ang.rad)]])
        elif axis == "z":
            M = np.array([[np.cos(ang.rad), np.sin(ang.rad), 0], \
                          [-np.sin(ang.rad), np.cos(ang.rad), 0], \
                          [0, 0, 1]])
        else:
            print("Unrecognized Axis")
            return None

        vect = np.inner(M, self.xyz)

        return Cart(vect[0], vect[1], vect[2])

def jd2Angle(jd):

    """ Converts julian date to the angle between greenwich and the First Point of Aries
    input:
    jd [float] - Julian Date

    returns:
    theta [Angle Object] - angle between greenwich and the First Point of Aries
    """
    try:
        jd = float(jd)
    except ValueError:
        print("[WARNING] Julian Date is not a float")
        return None

    t = jd - JUL_DAYS_BEFORE_JAN_1_2000
    t_cent = t/DAYS_PER_YEAR/YEARS_PER_CENTURY
    theta = T2ANG_COEFFS[0] + T2ANG_COEFFS[1]*t + T2ANG_COEFFS[2]*t_cent**2 - t_cent**3/T2ANG_COEFFS[3]
    theta = Angle(theta, deg=True)

    return theta

def cart2Radec(cart):

    """
    Extension to Cart object to correctly convert theta and phi to Right Ascension and Declination
    inputs:
    cart [Cart Object] - cart object to convert

    returns:
    ra [Right Ascension Obj] - Right ascension of cart
    dec [float] - Declination of cart
    """
    if not hasattr(cart, "theta") or not hasattr(cart, "phi"):
        print("[WARNING] Cartesian object does not have angles!")
        return None, None

    dec = 90 - cart.theta.deg
    ra  = cart.phi.deg

    ra = RightAsc(ra)

    return ra, dec

def aries2greenwich(cart, JD):

    """
    Rotates the coordinate system of cart from the First Point of Aries to the location of 
    Greenwich at a given julian date

    inputs:
    cart [Cart Obj] - vector to rotate
    JD [float] - Julian date (to find the position of Greenwich)
    """

    # angle is negative in the notes
    shift = Angle(360 - jd2Angle(JD).deg, deg=True)

    return cart.rotate(shift, axis="z")

if __name__ == "__main__":

    # 1A
    print("############ 1A Testing")
    print("Angle from JD 2456293.520833...")
    print(jd2Angle(2456293.520833))
    print("Angle from JD 2451545.0... (should be 280.46... by equations)")
    print(jd2Angle(2451545.0))
    print("Angle from half a JD later...")
    print(jd2Angle(2451545.5))
    print("Angle from a full JD later...")
    print(jd2Angle(2451546.0))
    print("JD where Greenwich = Aries...")
    print(jd2Angle(2451546.2176))
    print("Angle from a non-float JD")
    print(jd2Angle("not a float"))
    print("")

    # 1B
    print("############ 1B Testing")
    print("Using a non-float cartesian")
    r = Cart(1, "A", 1)
    print("RA and DEC from (1, -0.0001, 1) - expected: 23h 59m - 24h 0m and ~45 deg")
    r, d = cart2Radec(Cart(1, -0.0001, 1))
    print(r)
    print(d)
    print("RA and DEC from (0, 1, sqrt(3)) - expected: 8h and 60 deg")
    r, d = cart2Radec(Cart(0, 1, np.sqrt(3)))
    print(r)
    print(d)
    print("")

    print("############ 1C Testing")
    print("(1, 0, 1) at JD = 2451545.0, should be rotated around Z axis by 280.46... deg")
    print(aries2greenwich(Cart(1, 0, 1), 2451545.0))
    print("(1, 0, 1) at JD = 2451546.2176, should be rotated around Z axis by 0 deg")
    print(aries2greenwich(Cart(1, 0, 1), 2451546.2176))
    print("")

    print("############ 1D Testing")
    print("Position of north pole at the center of earth")
    print(Cart(0, 0, 1).geoPos())
    print("Position of y=0 (0 Longitude)")
    print(Cart(600000, 0, 600000).geoPos())
    print("Equator position, 90 degrees away from before")
    print(Cart(0, 600000, 0).geoPos())
    print("")

    ### CODE BELOW IS FOR Q2
    print("Q2")
    print("######################")
    ### Q2
    ISS_rad= 6680000
    tilt = 51.6 #deg
    T = 92*60 #sec

    
    iss = Cart(ISS_rad, 0, 0)
    ra, dec = cart2Radec(iss)
    
    #ISS rotates horizontally over period
    # x = cos(360/T)
    # y = sin(360/T)

    points_x = []
    points_y = []
    points_z = []
    points_ra = []
    points_dec = []
    points_lat = []
    points_lon = []
    points_alt = []

    JD_const = 2456674.5
    for i in range(T*3):
        if i%60 == 0:

            JD = JD_const + i/24/60/60

            # Circle which sweeps out orbit
            iss = Cart(ISS_rad*np.cos(i/T*2*np.pi), ISS_rad*np.sin(i/T*2*np.pi), 0)

            # Rotate to inclination
            iss = iss.rotate(Angle(51.6, deg=True), "x")

            if i < T:
                ra, dec = cart2Radec(iss)
                points_ra.append(ra.asFloat())
                points_dec.append(dec)

            # JD shift
            iss = aries2greenwich(iss, JD)
            
            points_x.append(iss.x)
            points_y.append(iss.y)
            points_z.append(iss.z)

            lat, lon, alt = iss.geoPos(pnt=False)
            points_lat.append(lat)
            points_lon.append(lon)
            points_alt.append(alt)


    points_x = np.array(points_x)
    points_y = np.array(points_y)
    points_z = np.array(points_z)

    import matplotlib.pyplot as plt

    from mpl_toolkits.mplot3d import Axes3D


    fig = plt.figure()
    ax = Axes3D(fig)
    l = len(points_x)//3
    ax.scatter(points_x[0:l]/1000, points_y[0:l]/1000, points_z[0:l]/1000, 'r')
    ax.scatter(points_x[l:2*l]/1000, points_y[l:2*l]/1000, points_z[l:2*l]/1000, 'b')
    ax.scatter(points_x[2*l:3*l]/1000, points_y[2*l:3*l]/1000, points_z[2*l:3*l]/1000, 'g')
    ax.scatter(ISS_rad/1000, 0, 0, 'k', marker="*")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.show()

    plt.scatter(np.array(points_ra[0:l])/15, points_dec[0:l])
    plt.scatter(np.array(points_ra[l:2*l])/15, points_dec[l:2*l])
    plt.scatter(np.array(points_ra[2*l:3*l])/15, points_dec[2*l:3*l])
    plt.xlabel("Right Ascension [h]")
    plt.ylabel("Declination [deg]")
    plt.title("RA and DEC of the ISS")
    plt.show()

    # plt.scatter(points_lon, points_lat)
    # plt.show()
    print(points_lon[0], points_lat[0], points_alt[0])
    plt.scatter(points_lon, points_lat)
    plt.xlabel("Longitude (deg E)")
    plt.ylabel("Latitude (deg N)")
    plt.title("Orbits of the ISS")
    plt.show()

    # The following basemap parameters are copied from my seismic project - which is mostly taken from StackOverflow
    from mpl_toolkits.basemap import Basemap
    m = Basemap(projection='robin',lon_0=0,resolution='c')
    # m.drawmapscale(13, 47.75, 13, 47.5, 10, barstyle='fancy')
    # m.bluemarble(scale=1)   # full scale will be overkill
    # m.drawcoastlines()  # add coastlines
    m.fillcontinents(color='grey',lake_color='aqua')


    x, y = m(points_lon, points_lat)  # transform coordinates
    plt.scatter(x, y, 32, color='red', zorder=3)
    plt.title("Orbits of the ISS")
    plt.show()