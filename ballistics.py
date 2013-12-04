from numpy import *
from scipy.optimize import fmin

class Ballistics( object ) :

    """
    x = x0 + v*delt + .5 * delt**2 * a

    v = v0 + a * delt

    """


    Cd = 0.42  # http://en.wikipedia.org/wiki/Drag_coefficient

    g = array([0,-9.81])   # gravity (m/s^2)

    dt = .0005  # sec)

    ft2m = 0.3048
    in2m = .0254
    m2in = 1/.0254
    grain2kg = 6.479891e-5
    kg2lb = 2.20462
    J2ftlb = 1/1.3558179483314

    MV = 2500*ft2m

    v0 = array([MV,0])*ft2m

    m = 200*grain2kg    # mass (kg)

    rho = 1.1839   # kg/m^3  
    # Variable should make adjustable with atmospheric
    # http://en.wikipedia.org/wiki/Density_of_air

    A = pi*(.3*in2m)**2/4

    sightHeight = 1.5*in2m  


    # def __init__(self) :
    #     self.setLaunchAngle(0)

    def zeroSight(self,range,theta0_guess=.1*pi/180) :
        """
        Zero the sights at specified range (in meters).
        """

        # Iteratively adjust launch angle until on target (within some
        # tolerance at desired range.

        def minFunc(angle) :
            self.setLaunchAngle(angle[0])
            self.solve()
            idx = find(self.pos[:,0] > range)
            return self.pos[idx[0],1]**2

        X = fmin(minFunc,theta0_guess)
        self.setLaunchAngle(X[0])
        return X[0]

    def solve(self) :
        
        pos = [array([0,-self.sightHeight])]
        vel = [self.v0]
        t = [0]

        while pos[-1][1] > -2.0 :  # go until 2m drop

            a = self.g + self.dragAccel(vel[-1])

            pos.append(pos[-1] + vel[-1]*self.dt + .5*a*self.dt**2)
            # update vel
            vel.append(vel[-1] + a*self.dt)
            t.append(t[-1] + self.dt)

        self.pos = array(pos)
        self.vel = array(vel)
        self.energy = .5*self.m*( self.vel[:,0]**2 + self.vel[:,1]**2 )
        self.t = array(t)

    def plot(self) :
        figure(1); 

        plot(self.pos[:,0],self.pos[:,1])
        xlabel('Pos (m)')
        ylabel('Pos (m)')

        figure(2); 
        subplot(211)
        plot(self.t,self.vel[:,0])
        ylabel('x vel (m/s)')

        subplot(212)
        plot(self.t,self.vel[:,1])
        ylabel('y vel (m/s)')

        xlabel('Time (sec)');

        figure(3);
        plot(self.pos[:,0],self.energy*self.J2ftlb)
        xlabel('x pos (m)')
        ylabel('Energy (ft-lb)')
        


    def dragAccel(self,v) :
        """
        Return drag acceleration 
        Fd = .5*rho*v**2*Cd*A = m*a
        """
        return -.5*self.rho*v**2*self.Cd*self.A/self.m

    def setLaunchAngle(self,angle) :
        self.theta0 = angle
        self.v0 = array([ self.MV*cos(self.theta0),
                          self.MV*sin(self.theta0) ])

    def BC2Cd(self,BC) :
        """
        Convert ballistics coefficient to drag coefficient

        http://en.wikipedia.org/wiki/Ballistic_coefficient
        
        BC_bullets = SD/i = SD*Cg/Cb
        
        SD = sectional density lb/in^2; mass in lbs / caliber squared
             in inches

        Cg = 0.5191

        Cb = SD*Cg/BC
        """

        SD = self.m*Ballistics.kg2lb / (self.cal**2)
        Cg = 0.5191

        self.Cd = SD*Cg/BC
        self.BC = BC

        return self.Cd

    def calcAirDensity(self,Tfarenheight,h_ft) :
        """
        See: http://en.wikipedia.org/wiki/Density_of_air

        """
        
        # convert T to Kelvin
        T = (Tfarenheight-32)*5./9 + 273.15
        h = h_ft*self.ft2m

        p0 = 101.325 
        T0 = 288.15
        g = 9.81
        M = .0289644
        R = 8.31447
        L = .0065
        Rs = 287.058

        p = p0*(1-L*h/T0)**(g*M/R/L)


        return (p*1000)*M/R/T
        

class Bullet30_06( Ballistics ) :
    cal = .308   # in inches
    A = pi*(cal*Ballistics.in2m)**2/4
    MV = 2500*Ballistics.ft2m
    m = 200*Ballistics.grain2kg

class Bullet22LR( Ballistics ) :
    cal = .223   # in inches
    A = pi*(cal*Ballistics.in2m)**2/4
    MV = 1750*Ballistics.ft2m
    m = 30*Ballistics.grain2kg
    
class Bullet300WinMag( Ballistics ) :
    cal = .308   # in inches
    A = pi*(cal*Ballistics.in2m)**2/4
    MV = 2950*Ballistics.ft2m
    m = 180*Ballistics.grain2kg
    

if __name__ == '__main__' :
    
    close('all')
    b22 = Bullet22LR()
    b22.solve()

    b30_06 = Bullet30_06()
    b30_06.solve()

    b300winmag = Bullet300WinMag()
    b300winmag.BC2Cd(0.509)
    b300winmag.solve()

    b22.plot()
    b30_06.plot()
    b300winmag.plot()
    show()
