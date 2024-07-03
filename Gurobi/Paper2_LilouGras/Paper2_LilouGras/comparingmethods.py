#######################################################
# Author: Lilou Gras
# Inspired from: Temporal vs. spatial formulation of autonomous overtaking algorithms, by Johan Karlsson, Nikolce Murgovski and Jonas Sjoberg
# Mail : lilou.gras@gmail.com
# Date: 06/2023
#######################################################

from gurobipy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from celluloid import Camera



# Function doing the overtaking in the time domain (the one we are used to)
def time_domain():
    mod = Model() # Creating the model

    deltaT = 27/110 # It takes 27s (from the paper)
    N = 110   # defining number of steps, here 90
    print(deltaT*N) 

    # Setting initial values for our vehicles
    x_0 = 0
    y_0 = 1.05
    x_f = 600
    vx_0 = 70/3.6

    # Setting initial values for the leading vehicle
    xL_0 = 70
    yL = 0.9
    vL = 50/3.6

    # Setting the displacement of the leading vehicle
    xL = np.zeros(N)
    xL[0] = xL_0
    for i in range (1,N):
        xL[i] = xL[i-1] + vL*deltaT
    Lo = 30
    Wo= 0.3
    

    # Creating variables, position, velocities, acceleration and jerk
    x = mod.addVars(N, lb=0, ub=2000, vtype=GRB.CONTINUOUS, name="x")
    y = mod.addVars(N, lb =0, ub=5, vtype=GRB.CONTINUOUS, name="y")
    vx = mod.addVars(N, lb=0, ub=40, vtype=GRB.CONTINUOUS, name= 'vx' )
    vy = mod.addVars(N, lb=-2, ub =2, vtype=GRB.CONTINUOUS, name="vy")
    ax = mod.addVars(N, lb=-4, ub=3, vtype=GRB.CONTINUOUS, name="ax")
    ay = mod.addVars(N, lb=-1, ub=1, vtype=GRB.CONTINUOUS, name="ay")
    jx = mod.addVars(N, lb=-5, ub=5, vtype=GRB.CONTINUOUS, name="jx")
    jy = mod.addVars(N, lb=-5, ub=5, vtype=GRB.CONTINUOUS, name="jy")

    # Creating the binary variables for the big-M method 
    d1 = mod.addVars(N, name = 'd1', vtype = GRB.BINARY)
    d2 = mod.addVars(N, name = 'd2', vtype = GRB.BINARY) 
    d3 = mod.addVars(N, name = 'd3', vtype = GRB.BINARY) 
    d4 = mod.addVars(N, name = 'd4', vtype = GRB.BINARY) 

    # Defining the cost
    cost = LinExpr()
    for i in range (N):
        if i == 0:
            cost+= 5*(x[N-1]-x_f)**2
        cost +=  3*(y[i]-y_0)**2 + (vx[i]-vx_0)**2 + 10*ay[i]**2 + 10*jy[i]**2 + ax[i]**2 + 4*jx[i]**2
        # cost += ax[i]**2 + 4*jx[i]**2 + 3*(y[i]-y_0)**2 + 4*(vy[i]**2) + 5*ay[i]**2  + 20*jy[i]**2
    mod.setObjective( cost , GRB.MINIMIZE)

    # Set initial constraints
    mod.addConstr(x[0] == x_0,  name="initx")
    mod.addConstr(y[0] == y_0,  name="inity")
    mod.addConstr(vx[0] == vx_0, name="initvx")

    # Set constraints for position, velocities and acceleration
    mod.addConstrs((x[i] + deltaT*(vx[i]) + 0.5*(deltaT**2)*ax[i] + (1/6)*(deltaT**3)*jx[i] == x[i+1] for i in range(N-1)),"c1")
    mod.addConstrs((y[i] + deltaT*vy[i] + 0.5*(deltaT**2)*ay[i] + (1/6)*(deltaT**3)*jy[i] == y[i+1] for i in range(N-1)),"c2")
    mod.addConstrs((vx[i] + deltaT*ax[i] + 0.5*(deltaT**2)*jx[i]  == vx[i+1] for i in range(N-1)),"c3")
    mod.addConstrs((vy[i] + deltaT*ay[i] + 0.5*(deltaT**2)*jy[i] == vy[i+1] for i in range(N-1)),"c4")
    mod.addConstrs((ax[i] + deltaT*jx[i] == ax[i+1] for i in range(N-1)),"j3")
    mod.addConstrs((ay[i] + deltaT*jy[i] == ay[i+1] for i in range(N-1)),"j4")


    M=10000 # take a large value to make sure

    # Obstacle avoidance 
    mod.addConstrs((x[i] <= (xL[i] - Lo) + M*(1- d1[i]) for i in range(N)), "c5")
    mod.addConstrs((x[i] >= (xL[i] - Lo) - M*d1[i] for i in range(N)), "c6")

    mod.addConstrs((x[i] >= (xL[i] + Lo) - M*(1-d2[i]) for i in range(N)), "c7")
    mod.addConstrs((x[i] <= (xL[i] + Lo) + M*d2[i] for i in range(N)), "c8")

    mod.addConstrs((y[i] >= (yL + 1.3*Wo) - M*(1-d3[i]) for i in range(N)), "c9")
    mod.addConstrs((y[i] <= (yL + 1.3*Wo) + M*d3[i] for i in range(N)), "c10")

    mod.addConstrs((y[i] <= (yL - 1.3*Wo) + M*(1-d4[i]) for i in range(N)), "c11")
    mod.addConstrs((y[i] >= (yL - 1.3*Wo) - M*d4[i] for i in range(N)), "c12")

    mod.addConstrs((d1[i] + d2[i] + d3[i] + d4[i] == 1 for i in range(N)), "c13")

    # Launching the optimizer
    mod.optimize()
    
    if mod.Status == GRB.OPTIMAL: 
        fig = plt.figure(1)
        camera = Camera(fig)
        xr, yr, vxr, vyr, axr, ayr = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
        for i in range(N-1):
            # !! to keep if you remove plotting !! #########
            xr[i] = x[i].X
            yr[i] = y[i].X
            vxr[i] = vx[i].X
            vyr[i] = vy[i].X
            axr[i] = ax[i].X
            ayr[i] = ay[i].X
            ##############################################
            plt.title('Obstacle avoidance video in the time domain')
            plt.plot([x[i].X, x[i].X],[y[i].X, y[i].X], 'ko-')

            plt.plot([xL[i]-(Lo/2), xL[i]-(Lo/2), xL[i]+(Lo/2), xL[i]+(Lo/2), xL[i]-(Lo/2)], [yL-(Wo/1.5), yL+(Wo/1.5), yL+(Wo/1.5), yL-(Wo/1.5), yL-(Wo/1.5)], 'r-', linewidth=1.5)
            camera.snap() #plotting an animated graph
            
        anim = camera.animate()
        plt.show()

    return(xr, yr, vxr, vyr, axr, ayr)

# Function doing the same overtaking but this time in the spatial domain
def spatial_domain():
    N=110 # nb of steps
    xe = np.arange(0, 601, 600/(N-1)) # Setting up for the sampling
    xref = np.zeros(N)
    deltaT = 27/N # We still need this to calculate velocities in m/s 

    # Setting up values for leading vehicle
    xL = 75
    vL = 50/3.6
    yL = 0.9

    Lo = 20
    Wo = 0.3

    # Creating the samples for the spatial domain
    for i in range(N):
        xref[i] = xe[i] - (50/3.6)*i*deltaT

    print(xref)


    deltaX = xref[1]-xref[0] # the increment is useful to get the velocites


    mod = Model() # Creating the model

    # Initial values for the ego vehicle
    vx_0 = 70/3.6 #from km/h in m/s
    vr = 70/3.6
    y_0 = 1.05


    # Setting up the variables, here we don't have x since it is in the spatial domain, we use x~ instead of t (refer to paper)
    y = mod.addVars(N, lb =0, ub=5, vtype=GRB.CONTINUOUS, name="y")
    vx = mod.addVars(N, lb=-20, ub=40,  vtype=GRB.CONTINUOUS, name="vx")
    vy = mod.addVars(N, lb=-4, ub =4, vtype=GRB.CONTINUOUS, name="vy")
    ax = mod.addVars(N, lb=-4, ub=3,  vtype=GRB.CONTINUOUS, name="ax")
    ay = mod.addVars(N, lb=-1, ub=1, vtype=GRB.CONTINUOUS, name="ay")
    jx = mod.addVars(N, lb=-1, ub=1, vtype=GRB.CONTINUOUS, name="jx")
    jy = mod.addVars(N, lb=-1, ub=1, vtype=GRB.CONTINUOUS, name="jy")


    # Same as before, for big-M method
    d1 = mod.addVars(N, name = 'd1', vtype = GRB.BINARY)
    d2 = mod.addVars(N, name = 'd2', vtype = GRB.BINARY) 
    d3 = mod.addVars(N, name = 'd3', vtype = GRB.BINARY) 
    d4 = mod.addVars(N, name = 'd4', vtype = GRB.BINARY) 

    # Cost function
    cost = LinExpr()
    for i in range (N):

        cost +=  10*(y[i]-y_0)**2 + (vx[i]-vr)**2 + 2*ax[i]**2 + 2*jx[i]**2 + 0.5*ay[i]**2 + 0.5*jy[i]**2
    mod.setObjective( cost , GRB.MINIMIZE)

    # Set initial constraints
    mod.addConstr(y[0] == y_0,  name="inity")
    mod.addConstr(vx[0] == vx_0, name="initvx")

    # set constraints for position, velocity, and acceleration, here note that in the y domain we have deltaT/deltaX, that's because we have dy/dx~= dy/dt*dt/dx~
    mod.addConstrs((vx[i] + (deltaT/deltaX)*ax[i] + 1/2*deltaX**2*jx[i] == vx[i+1] for i in range(N-1)),"c1")
    mod.addConstrs((ax[i] + (deltaT/deltaX)*jx[i] == ax[i+1] for i in range(N-1)),"cx")
    mod.addConstrs((y[i] + (deltaT/deltaX)*vy[i] + 1/2*(deltaT/deltaX)**2*ay[i] + 1/6*(deltaT/deltaX)**3*jy[i] == y[i+1] for i in range(N-1)),"c2")
    mod.addConstrs((vy[i] + (deltaT/deltaX)*ay[i] + 1/2*(deltaT/deltaX)**2*jy[i] == vy[i+1] for i in range(N-1)),"c3")
    mod.addConstrs((ay[i] + (deltaT/deltaX)*jy[i]  == ay[i+1] for i in range(N-1)),"c4")

    M=10000

        # Obstacle avoidance 

    mod.addConstrs((float(xref[i]) <= (xL - Lo) + M*(1- d1[i]) for i in range(N)), "c5")
    mod.addConstrs((float(xref[i]) >= (xL - Lo) - M*d1[i] for i in range(N)), "c6")

    mod.addConstrs((float(xref[i]) >= (xL + Lo) - M*(1-d2[i]) for i in range(N)), "c7")
    mod.addConstrs((float(xref[i]) <= (xL + Lo) + M*d2[i] for i in range(N)), "c8")

    mod.addConstrs((y[i] >= (yL + 1.3*Wo) - M*(1-d3[i]) for i in range(N)), "c9")
    mod.addConstrs((y[i] <= (yL + 1.3*Wo) + M*d3[i] for i in range(N)), "c10")

    mod.addConstrs((y[i] <= (yL - 1.3*Wo) + M*(1-d4[i]) for i in range(N)), "c11")
    mod.addConstrs((y[i] >= (yL - 1.3*Wo) - M*d4[i] for i in range(N)), "c12")

    mod.addConstrs((d1[i] + d2[i] + d3[i] + d4[i] == 1 for i in range(N)), "c13")


    mod.optimize()

    yr, vxr, vyr, ayr, axr = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
    for i in range(N):
    # !! to keep if you remove plotting !! #########
        yr[i] = y[i].X
        vxr[i] = vx[i].X
        axr[i] = ax[i].X
        vyr[i] = vy[i].X
        ayr[i] = ay[i].X



    return(yr, vyr, ayr, xref, vxr, axr, deltaX)



# Launching the functions
ys, vys, ays, xs, vxs, axs, dx = spatial_domain()

xv, yv, vxv, vyx, axv, ayv = time_domain()

# Plotting
xv[-1] = xv[-2]
yv[-1] = yv[-2]
fig, ax = plt.subplots(3)

fig.suptitle('Overtaking of the leading vehicle in time domain')

ax[0].set(xlabel='Relative Longitudinal distance (m)', ylabel='Lateral position (m)')
ax[0].add_patch(Rectangle((70,0.9), 20, 0.3))
ax[0].plot(xs, yv)

ax[1].set(xlabel='Relative Longitudinal distance (m)', ylabel='Lateral velocity (m/s)')
ax[1].plot(xs, vyx*dx)

ax[2].plot(xs, ayv*(dx**2))
ax[2].set(xlabel='Relative Longitudinal distance (m)', ylabel='Lateral acceleration (m/s^2)')

plt.show()



fig2, ax2 = plt.subplots(3)
fig2.suptitle('Overtaking of the lead vehicle in the spatial domain')
    
ax2[0].plot(xs, ys)
ax2[0].add_patch(Rectangle((70,0.9), 20, 0.3))
ax2[0].set(xlabel='Relative longitudinal distance (m)', ylabel = 'Lateral position (m)')


ax2[1].plot(xs,vys)
ax2[1].set(xlabel='Relative longitudinal distance (m)', ylabel = 'Lateral speed (m/s)')

ax2[2].plot(xs,ays)
ax2[2].set(xlabel='Relative longitudinal distance (m)', ylabel = 'Lateral acceleration (m/s^2)')
plt.show()

fig3, ax3 = plt.subplots(3)
fig3.suptitle('Comparison of spatial and time domain')

    
ax3[0].plot(xs, ys, label='spatial domain')
ax3[0].plot(xs, yv, label='time domain')
ax3[0].add_patch(Rectangle((70,0.9), 20, 0.3))
ax3[0].set(xlabel='Relative longitudinal distance (m)', ylabel = 'Lateral position (m)')


ax3[1].plot(xs,vys, label='spatial domain')
ax3[1].plot(xs, vyx*dx, label='time domain')
ax3[1].set(xlabel='Relative longitudinal distance (m)', ylabel = 'Lateral speed (m/s)')

ax3[2].plot(xs,ays, label='spatial domain')
ax3[2].plot(xs, ayv*(dx**2), label='time domain')
ax3[2].set(xlabel='Relative longitudinal distance (m)', ylabel = 'Lateral acceleration (m/s^2)')

ax3[0].legend(loc="upper right")
ax3[1].legend(loc="upper right")
ax3[2].legend(loc="upper right")
plt.show()