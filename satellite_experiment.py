"""
Final code for satellite experiment Compsim project
David Roddy
"""

# Satellite only launches properly if dt = 22000 and its x position is
# 60 000km (6E7) from Earth, otherwise it flies off to the right.




""" writing Beeman algorithm to use for solar system
i think the acceleration calculation stays the same as it isnt part of the algo
need to change the velocity and position methods
need to store the past two acceleration values
make it ahist[0,1,2,3]
added potential energy calculator and total energy calculator
added a function to graph the total energy over time (nulled in this code)
added a function to order self.bodies by mass
added a satellite in a seperate function
added a function to print the time when the satellite comes within a given radius of mars
modified mars' (and satellite's) position to make the satellite come close to it
gives initial velocity of satellite relative to Earth
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

G = 6.673E-11


class Body(object):
    
    
    def __init__(self,position,velocity,radius,colour,mass,name,calcinit,period):
        
        # Constructor method to create instance variables for each body.
        
        self.position = position
        self.velocity = velocity
        self.radius = radius
        self.colour = colour
        self.mass = mass
        self.name = name
        self.calcinit = calcinit
        self.period = period
        self.ahist = [0,0,0,0]
        
    def v_init(self,bodies):
        
        # Using v = sqrt(GM/r) finds the orbital velocity.  This
        # is then split into components in the second code block so that the
        # direction of motion is always perpendicular to the orbital radius.
        
        r = (self.position[0])**2 + (self.position[1])**2
        r = math.sqrt(r)
        v = math.sqrt(G*bodies[0].mass/r)
        
        if self.position[0] == 0:
            self.position[0] = 1    
        d = (self.position[1]/self.position[0])
        self.velocity[1] = math.sqrt((v**2) / (1 + d**2))
        self.velocity[0] = - d * self.velocity[1]
        
    def calc_a(self,bodies,body_index):
        
        # This function is called on each body in turn.  The acceleration of
        # this body due to each other body is then calculated and added to
        # a total. This is returned and passed to calc_v.
        
        acceleration = [0,0]
        target = bodies[body_index]
        for index, otherbod in enumerate(bodies):
            if index != body_index:
                r = (target.position[0]-otherbod.position[0])**2 + (target.position[1]-otherbod.position[1])**2
                r = math.sqrt(r)
                force = G * otherbod.mass / r**3
                acceleration[0] += force * (otherbod.position[0] - target.position[0])
                acceleration[1] += force * (otherbod.position[1] - target.position[1])
        return acceleration

    def calc_v(self,bodies,dt):
        
        # Calls calc_a on each body in turn and finds the new velocity of each
        # body based on its output.  Stores the previous acceleration values
        # in body.ahist so they can be used by Beeman algo.
        
        for a,body in enumerate(bodies):
            acceleration = self.calc_a(bodies,a) 
            body.ahist.append(acceleration[0])
            body.ahist.append(acceleration[1])
            body.velocity[0] += 1/6 * ((2 * body.ahist[4]) + (5 * body.ahist[2]) - (body.ahist[0])) * dt   
            body.velocity[1] += 1/6 * ((2 * body.ahist[5]) + (5 * body.ahist[3]) - (body.ahist[1])) * dt
    def calc_p(self,bodies,dt):
        
        # For every body in the simulation, its position is updated using the
        # Beeman algo and previous acceleration values stored in ahist (which
        # is a six-element list).  After the calculations are done, the first
        # two elements of the list are deleted corresponding with the oldest
        # values of acceleration.
        
        for num,body in enumerate(bodies):
            body.position[0] += (body.velocity[0] * dt) + 1/6 * (4 * body.ahist[2] - body.ahist[0]) * (dt)**2
            body.position[1] += (body.velocity[1] * dt) + 1/6 * (4 * body.ahist[3] - body.ahist[1]) * (dt)**2
            del body.ahist[0:2]
                
    def step(self,bodies,dt):
        
        # Called by animate to perform a single step of size dt.
        
        self.calc_v(bodies,dt)
        self.calc_p(bodies,dt)
        
    def kinetic(self,bodies):
        
        # Finds the total kinetic energy of the system by
        # adding up the contribution from each body using 1/2 mv**2.
        
        K = 0
        for a,body in enumerate(bodies):
            body = bodies[a]
            vsq = body.velocity[0]**2 + body.velocity[1]**2
            K += 0.5 * body.mass * vsq
        return K
    
    def potential(self,bodies):
        
        # Goes through each pair of bodies and finds the potential between
        # them.  Adds them up and returns to animate function to be added
        # to kinetic.
        
        P = 0
        for a,body in enumerate(bodies):
            body = bodies[a]
            for num,target in enumerate(bodies):
                target = bodies[num]
                if a != num:
                    r = math.sqrt((body.position[0] - target.position[0])**2 + (body.position[1] - target.position[1])**2)
                    U = G * body.mass * target.mass / r
                    P += -1/4 * U
        return P
    
    
class Simulation(object):

    
    def __init__(self):
        
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
        mercury = Body([0,5.79E10],[0,0],5E9,'grey',3.3011E23,'Mercury',True,False)
        venus = Body([0,1.082E11],[-3E4,0],7E9,'g',4.8675E24,'Venus',True,False)
        earth = Body([0,1.496E11],[3E5,0],7E9,'b',5.972E24,'Earth',True,False)
        mars = Body([-1.5E11,171576251270.38998],[0,0],5E9,'r',6.417E23,'Mars',True,False)
        
        self.bodies = [sol,mercury,venus,earth,mars]
        self.bodies.sort(key=lambda x: x.mass, reverse=True)
        self.nrghist = []
        self.time = []
        self.reached = False
        
        for a,body in enumerate(self.bodies):
            body = self.bodies[a]
            if body.calcinit == True:
                body.v_init(self.bodies)
                if body.name == 'Mars':         # makes sure Mars is orbiting in the correct direction
                    body.velocity[0] = -1 *body.velocity[0]
                    body.velocity[1] = -1 *body.velocity[1]
        self.addsat()
        
    def init(self):
        
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
        
        vinit = math.sqrt(sx**2 + sy**2)
        print("inital speed is " + str(vinit) + " ms^-1")
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
        s = math.sqrt((q[0] - n[0])**2 + (q[1] - n[1])**2)
        
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
                r = math.sqrt((body.position[0] - self.bodies[0].position[0])**2 + (body.position[1] - self.bodies[0].position[1])**2)
                f = math.sqrt(r**3 / (G * self.bodies[0].mass))
                t = 2 * math.pi * f
                h = t / 31536000        # conversion to years
                print(str(body.name) + " : " + str(round(h,2)) + " years")
    
    def animate(self,i):
        
        # self.dt sets the speed of the animation.  Then p.step() updates the
        # positions of all bodies.  The third block calls the function to
        # find the total kinetic energy once every 100 animation frames.
        # The fourth block updates the positions of the patches based on the
        # new locations of the bodies.
        
        self.dt = 22000
        
        p = self.bodies[0]
        p.step(self.bodies,self.dt)
        self.eperiod(i)
        self.timesat(i)
        """
        if i % 10 == 0:
            K = p.kinetic(self.bodies)
            U = p.potential(self.bodies)
            T = K + U
            self.nrghist.append(T)
            self.time.append(i * self.dt / 31536000)
            
        if i % 500 == 0 and i != 0:
            self.nrggraph()
        """

        for a in range(len(self.bodies)):
            body = self.bodies[a]
            self.patches[a].center = (body.position[0], body.position[1])
        return self.patches
        
    def sim(self):
        
        # This func adds patches in the locations of the bodies and
        # sets the size of the axes.  Then the ani variable repeatedly calls
        # the animate function using FuncAnimation.
        
        fig1= plt.figure()
        ax1 = fig1.add_subplot(111)
        self.patches = []
        
        for body_index, body in enumerate(self.bodies):
            body = self.bodies[body_index]
            self.patches.append(plt.Circle((body.position[0], body.position[1]), body.radius, color = body.colour, animated = True))
            ax1.add_patch(self.patches[body_index])
        
        ax1.axis('scaled')
        ax1.set_xlim(-5e11,5e11)
        ax1.set_ylim(-5e11,5e11)
        
        ani = FuncAnimation(fig1, self.animate, init_func = self.init, frames = 10000, repeat = False, interval = 20, blit = True)
            
        
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
    
    
    test = Simulation()
    test.sim()


main()