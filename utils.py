import math


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