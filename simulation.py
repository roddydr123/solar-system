"""Final code for Compsim project with no satellite
David Roddy
"""




""" changed the velocity and position methods to Beeman
stores the past two acceleration values in ahist[0,1,2,3]
added potential energy calculator and total energy calculator
added a function to graph the total energy over time
added a function to order self.bodies by mass
(added a satellite in a seperate function
added a function to print the time when the satellite comes within a given radius of mars
modified mars' (and satellite's) position to make the satellite come close to it
gives initial velocity of satellite relative to Earth)
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from utils import Body

G = 6.673E-11
    
    
class Simulation(object):

    
    def __init__(self, satellite):
        
        # Add new bodies here, remember to include them in self.bodies.  The
        # third block calculates the initial velocity of each body 
        # if calcinit == True.  Change the velocity components and set calcinit
        # to False to start with your own conditions.  Second boolean arg is
        # True if you want to print the period of the body at intervals.  List
        # 'nrghist' holds total energy values as they're calculated to be
        # displayed graphically.  List 'time' holds the time each energy value
        # was calculated for the graphing function 'nrggraph'.  The list of
        # bodies is sorted by mass so that the first body is the most massive
        # and all periods are calculated around it.  Boolean 'self.reached' is
        # used to print the time the satellite reaches Mars only once.
        
        sol = Body([0,0],[0,0],1E10,'y',1.9884E30,'Sol',False,False)
        mercury = Body([0,5.79E10],[0,0],5E9,'grey',3.3011E23,'Mercury',True,True)
        venus = Body([0,1.082E11],[-3E4,0],7E9,'g',4.8675E24,'Venus',True,True)
        earth = Body([0,1.496E11],[3E5,0],7E9,'b',5.972E24,'Earth',True,True)
        mars = Body([-1.5E11,171576251270.38998],[0,0],5E9,'r',6.417E23,'Mars',True,True)
        
        self.bodies = [sol,mercury,venus,earth,mars]
        self.bodies.sort(key=lambda x: x.mass, reverse=True)
        self.nrghist = []
        self.time = []
        self.reached = False
        self.satellite = satellite
        
        for a,body in enumerate(self.bodies):
            body = self.bodies[a]
            if body.calcinit == True:
                body.v_init(self.bodies)
                if body.name == 'Mars':         # makes sure Mars is orbiting in the correct direction
                    body.velocity[0] = -1 *body.velocity[0]
                    body.velocity[1] = -1 *body.velocity[1]

        if self.satellite is True:
            self.addsat()
        
    def animation_init(self):
        
        # This function is called once before animate to put the patches in place.
        
        return self.patches
    
    def addsat(self):
        
        # A new body called satellite is appended to self.bodies and given similar
        # starting coordinates and velocity.  The numbers on the right can be
        # changed to vary the initial position and velocity relative to Earth.
        # Last section prints the initial velocity relative to Earth.
        
        satellite = Body([0,0],[0,0],4E9,'grey',3E2,'Satellite',False,True)
        
        for x,body in enumerate(self.bodies):
            body = self.bodies[x]
            if body.name == 'Earth':
                sx = -875
                sy = 0
                satellite.position[0] = body.position[0] - 6E7
                satellite.position[1] = body.position[1] - 0
                satellite.velocity[0] = body.velocity[0] + sx
                satellite.velocity[1] = body.velocity[1] + sy
        
        vinit = np.sqrt(sx**2 + sy**2)
        print("inital satellite speed is " + str(vinit) + " m/s")
        self.bodies.append(satellite)
    
    def timesat(self,i):
        
        # The first block finds the distance between Mars and the satellite.
        # If this is less than a given number, the time taken to get to this
        # point is printed.
        
        for a,body in enumerate(self.bodies):
            body = self.bodies[a]
            if body.name == 'Satellite':
                q = body.position
            if body.name == 'Mars':
                n = body.position
        s = np.sqrt((q[0] - n[0])**2 + (q[1] - n[1])**2)
        
        if s <= 7E9 and self.reached == False:
            w = round((self.dt * i / 31536000*12), 2)
            self.reached = True
            print("reached Mars in " + str(w) + " months")

    def eperiod(self,i):
        
        # Finds the orbital period using Kepler's 3rd law and prints this
        # to the terminal every so often
        
        for a,body in enumerate(self.bodies):
            body = self.bodies[a]
            if body.period == True and i % 200 == 0:
                r = np.sqrt((body.position[0] - self.bodies[0].position[0])**2 + (body.position[1] - self.bodies[0].position[1])**2)
                f = np.sqrt(r**3 / (G * self.bodies[0].mass))
                t = 2 * np.pi * f
                h = t / 31536000        # conversion to years
                print(str(body.name) + " : " + str(round(h,2)) + " years")
    
    def animate(self, i, satellite):
        
        # self.dt sets the speed of the animation.  Then p.step() updates the
        # positions of all bodies.  The third block calls the function to
        # find the total kinetic energy once every 100 animation frames.
        # The fourth block updates the positions of the patches based on the
        # new locations of the bodies.
        
        self.dt = 22000
        
        p = self.bodies[0]
        p.step(self.bodies,self.dt)

        if satellite is True:
            self.timesat(i)
        else:
            self.eperiod(i)
        
        if i % 10 == 0:
            K = p.kinetic(self.bodies)
            U = p.potential(self.bodies)
            T = K + U
            self.nrghist.append(T)
            self.time.append(i * self.dt / 31536000)
            
        if i % 500 == 0 and i != 0:
            # self.nrggraph()
            pass

        for a in range(len(self.bodies)):
            body = self.bodies[a]
            self.patches[a].center = (body.position[0], body.position[1])
        return self.patches
        
    def sim(self, satellite):
        
        # This func adds patches in the locations of the bodies and
        # sets the size of the axes.  Then the ani variable repeatedly calls
        # the animate function using FuncAnimation.
        
        fig1= plt.figure()
        ax1 = fig1.add_subplot(111)
        self.patches = []
        
        for body_index, body in enumerate(self.bodies):
            body = self.bodies[body_index]
            self.patches.append(plt.Circle((body.position[0], body.position[1]), body.radius, color = body.colour, animated = False))
            ax1.add_patch(self.patches[body_index])
        
        ax1.axis('scaled')
        ax1.set_xlim(-2.5e11,2.5e11)
        ax1.set_ylim(-2.5e11,2.5e11)
        
        anim = FuncAnimation(fig1, self.animate, fargs=[satellite], init_func = self.animation_init, frames = 10000, repeat = False, interval = 20, blit = True)
        plt.show()
        
    def nrggraph(self):
        
        # Creates a plot of time vs total energy which is displayed every 500
        # frames (I believe).
        
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        x = np.array(self.time)
        y = np.array(self.nrghist)
        plt.plot(x,y)
        plt.xlabel('time / years')
        plt.ylabel('total energy / Joules')
        plt.show()

def main():

    satellite = False
    
    try:
        if sys.argv[1] == "sat":
            satellite = True
    except IndexError:
        pass
    
    test = Simulation(satellite)
    test.sim(satellite)


main()