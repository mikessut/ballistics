from numpy import *

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
    grain2kg = 6.479891e-5

    MV = 2500*ft2m

    v0 = array([MV,0])*ft2m

    m = 200*grain2kg    # mass (kg)

    rho = 1.1839   # kg/m^3  
    # Variable should make adjustable with atmospheric
    # http://en.wikipedia.org/wiki/Density_of_air

    A = pi*(.3*in2m)**2/4

    theta0 = 1*pi/180


    def __init__(self) :
        self.setLaunchAngle(1)

    def solve(self) :
        
        pos = [array([0,0])]
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
        plot(self.pos[:,0],self.energy)
        xlabel('x pos (m)')
        ylabel('Energy (J)')
        


    def dragAccel(self,v) :
        """
        Return drag acceleration 
        Fd = .5*rho*v**2*Cd*A = m*a
        """
        return -.5*self.rho*v**2*self.Cd*self.A/self.m

    def setLaunchAngle(self,deg) :
        self.v0 = array([ self.MV*cos(deg*pi/180),
                          self.MV*sin(deg*pi/180) ])
        

class Bullet30_06( Ballistics ) :

    A = pi*(.308*Ballistics.in2m)**2/4
    MV = 2500*Ballistics.ft2m
    m = 200*Ballistics.grain2kg

class Bullet22LR( Ballistics ) :
    A = pi*(.223*Ballistics.in2m)**2/4
    MV = 1750*Ballistics.ft2m
    m = 30*Ballistics.grain2kg
    
class Bullet300WinMag( Ballistics ) :
    A = pi*(.308*Ballistics.in2m)**2/4
    MV = 3029*Ballistics.ft2m
    m = 200*Ballistics.grain2kg
    

if __name__ == '__main__' :
    
    close('all')
    b22 = Bullet22LR()
    b22.solve()

    b30_06 = Bullet30_06()
    b30_06.solve()

    b300winmag = Bullet300WinMag()
    b300winmag.solve()

    b22.plot()
    b30_06.plot()
    b300winmag.plot()
    show()
