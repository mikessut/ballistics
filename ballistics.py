from numpy import *
from scipy.optimize import fmin, newton

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
    m2yrd = 1.09361
    yrd2m = 1/1.09361
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

        X = fmin(minFunc,theta0_guess,ftol=.001**2)
        self.setLaunchAngle(X[0])
        return X[0]

    def endCond(self,pos,vel) :
        return pos[-1][1] < -2.0 # go until 2m drop

    def solve(self) :
        
        pos = [array([0,-self.sightHeight])]
        vel = [self.v0]
        t = [0]

        while not self.endCond(pos,vel) :
        
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

        plot(self.pos[:,0]*self.m2yrd,self.pos[:,1]*self.m2in,label=self.name)
        xlabel('Pos (yrd)')
        ylabel('Pos (in)')
        legend(loc=0)
        grid('on')

        savefig("Trajectory.png")

        figure(2); 
        subplot(211)
        plot(self.t,self.vel[:,0]/Ballistics.ft2m,label=self.name)
        ylabel('x vel (fps)')
        legend(loc=0)
        grid('on')

        subplot(212)
        plot(self.t,self.vel[:,1]/Ballistics.ft2m,label=self.name)
        ylabel('y vel (fps)')

        xlabel('Time (sec)');
        legend(loc=0)
        grid('on')

        savefig("Velocities.png")

        figure(3);
        plot(self.pos[:,0]*self.m2yrd,self.energy*self.J2ftlb,label=self.name)
        xlabel('x pos (yrd)')
        ylabel('Energy (ft-lb)')
        legend(loc=0)
        grid('on')

        savefig("Energy.png")

        figure(4)
        subplot(211)
        plot(self.t,self.pos[:,0])
        ylabel('x position (m)')
        grid('on')

        subplot(212)
        plot(self.t,self.pos[:,1])
        ylabel('y position (m)')
        grid('on')
        xlabel('Time (sec)')
        

        figure(5)
        plot(self.pos[:,0]*self.m2yrd,
             hstack([0, diff(self.energy*self.J2ftlb)/diff(self.pos[:,0]*self.m2yrd)]),
             label=self.name)
        xlabel('x pos (yrd)')
        ylabel('dEnergy/dx (ft-lb/yrd)')
        grid('on')
        legend(loc=0)
        
        savefig("dEnergy.png")

        # figure(4);
        # # 1./60 deg = 1 MOA
        # idx = self.pos[:,0] > 100
        # plot(self.pos[idx,0]*self.m2yrd,arctan(self.pos[idx,1]/self.pos[idx,0])*180/pi*60,label=self.name)
        # xlabel('x pos (yrd)')
        # ylabel('MOA')
        # grid('on')
        # legend(loc=0)
        


    def dragAccel(self,v) :
        """
        Return drag acceleration 
        Fd = .5*rho*v**2*Cd*A = m*a
        """

        vmag = norm(v)
        theta = arctan2(v[1],v[0])
        dragMag = .5*self.rho*vmag**2*self.Cd*self.A/self.m
        return -dragMag*array([cos(theta),sin(theta)])

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

    @staticmethod
    def calcAirDensity(Tfarenheight,h_ft) :
        """
        See: http://en.wikipedia.org/wiki/Density_of_air

        """
        
        # convert T to Kelvin
        T = (Tfarenheight-32)*5./9 + 273.15
        h = h_ft*Ballistics.ft2m

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
    name = '.30-06'
    cal = .308   # in inches
    A = pi*(cal*Ballistics.in2m)**2/4
    MV = 2700*Ballistics.ft2m
    m = 180*Ballistics.grain2kg  # 150, 180 fairly common

class Bullet300WinMag( Ballistics ) :
    name = '300 Winchester Magnum'
    cal = .308   # in inches
    A = pi*(cal*Ballistics.in2m)**2/4
    MV = 2950*Ballistics.ft2m
    m = 180*Ballistics.grain2kg

class Bullet7mmRemMag( Ballistics ) :
    """ Billy Odonnel's choise.
    """
    name = '7mm Remington Magnum'
    cal = .284   # in inches
    A = pi*(cal*Ballistics.in2m)**2/4
    MV = 3110*Ballistics.ft2m
    m = 150*Ballistics.grain2kg

class Bullet308Winchester( Ballistics ) :

    name = ".308 Winchester"
    cal = .308
    A = pi*(cal*Ballistics.in2m)**2/4
    MV = 2700*Ballistics.ft2m
    m = 165*Ballistics.grain2kg  # could go lighter

class Bullet270Winchester( Ballistics ) :
    name = ".270 Winchester"
    cal = .277
    A = pi*(cal*Ballistics.in2m)**2/4
    MV = 2850*Ballistics.ft2m
    m = 150*Ballistics.grain2kg


class Bullet65mmLapua( Ballistics ) :
    """ Tyler's
    """
    name = '6.5x47mm Lapua'
    cal = .264   # in inches
    A = pi*(cal*Ballistics.in2m)**2/4
    MV = 2900*Ballistics.ft2m
    m = 123*Ballistics.grain2kg

class Bullet22LR( Ballistics ) :
    name = '.22LR'
    cal = .223   # in inches
    A = pi*(cal*Ballistics.in2m)**2/4
    MV = 1750*Ballistics.ft2m
    m = 30*Ballistics.grain2kg
    
    
"""
Desired characteristics:

1) Straight trajectory
  a) High initial velocity
  b) Low drag force
     i) High mass
    ii) Small diameter

2) Low recoil
  a) Small mass
  b) Low initial velocity

3) Take down power

Clearly there are contradictions between these.  One thing that is
clear is we want small diameter.  That seems to be a good thing no
matter what. 

"""

#if __name__ == '__main__' :
def bulletComparisons() :
    
    bullets = [#Bullet22LR(),
               Bullet300WinMag(),
               Bullet30_06(),
               Bullet308Winchester(),
        Bullet270Winchester(),
               Bullet7mmRemMag()]
               #Bullet65mmLapua()]


    close('all')

    legtxt = []
    for b in bullets :
        b.zeroSight(150*Ballistics.yrd2m)
        b.solve()
        b.plot()
        legtxt.append(b.name)

    legend()

    show()

def tempComparison() :
    tempF = linspace(15,60,5)

    b = Bullet270Winchester()
    b.rho = b.calcAirDensity(50,0)
    b.zeroSight(100*Ballistics.yrd2m)
    basename = b.name
    for T in tempF :
        b.rho = b.calcAirDensity(T,0)
        b.name = basename + ' T = ' + str(T) + ' F'
        b.solve()
        b.plot()

def altComparison() :
    alt = linspace(5000,14000,5)

    b = Bullet270Winchester()
    b.rho = b.calcAirDensity(50,0)
    b.zeroSight(100*Ballistics.yrd2m)
    basename = b.name
    for a in alt :
        b.rho = b.calcAirDensity(50,a)
        b.name = basename + ' MSL = ' + str(a) + ' ft'
        b.solve()
        b.plot()

def straightDown() :
    class BulletDown( Bullet270Winchester ) :
        def __init__(self) :
            super(BulletDown,self).__init__()
            self.setLaunchAngle(-pi/2)
        
        def endCond(self, pos, vel) :
            return pos[-1][1] < -5000

    b = BulletDown()

    b.solve()
    b.plot()

    return b
