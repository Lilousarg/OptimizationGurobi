#######################################################
# Author: Lilou Gras
# Inspired from: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7795682
# Mail : lilou.gras@gmail.com
# Date: 05/2023
#
# Aim: Taking values from obstacles and finding an optimized path for our vehicle
# Inputs: -xo: array of x initial position of the obstacles
#         -yo: array of y initial position of the obstacles
#         -vo: array of constant speed of the obstacles
#              (array length is the number of obstacles)
#
# Outputs: -xv : array of size N for x position of our vehicle
#          -yv : array of size N for y position of our vehicle
#          -vxv : array of size N for x velocities of our vehicle
#          -vyv : array of size N for y velocities of our vehicle
#          -axv : array of size N for x acceleration of our vehicle
#          -ayv : array of size N for y acceleration of our vehicle
#######################################################


from gurobipy import *
import numpy as np


def avoiding_2obstacles(xo, yo, vo):
    mod = Model() # Creating the model

    deltaT = 0.2  # defining time step, here 1s
    N = 35    # defining number of steps, here 50
    # print(deltaT*N) 

    # Setting initial values for our vehicles
    x_0 = 0
    y_0 = 2.5
    x_f = 100
    vx_0 = 15

    # Creating obstacles dictionnary
    obs = {}
    print(len(xo))
    for n in range(len(xo)):
        obs['xo0' + str(n)] = xo[n]
        obs['yo0' + str(n)] = yo[n]
        obs['vo' + str(n)] = vo[n]
        xoi = np.ones(N)
        for i in range(1,N):
            xoi[0] = obs['xo0' + str(n)]
            xoi[i] = xoi[i-1] + obs['vo' + str(n)]*deltaT
        obs['xo' +str(n)] =  xoi
    
    print(obs)
    Lo = 10
    Wo= 1
    

    # Creating variables, position, velocities, acceleration and jerk
    x = mod.addVars(N, lb=0, ub=2000, vtype=GRB.CONTINUOUS, name="x")
    y = mod.addVars(N, lb =0, ub=5, vtype=GRB.CONTINUOUS, name="y")
    vx = mod.addVars(N, lb=0, ub=20,  vtype=GRB.CONTINUOUS, name="vx")
    vy = mod.addVars(N, lb=-2, ub =2, vtype=GRB.CONTINUOUS, name="vy")
    ax = mod.addVars(N, lb=-4, ub=3, vtype=GRB.CONTINUOUS, name="ax")
    ay = mod.addVars(N, lb=-1, ub=1, vtype=GRB.CONTINUOUS, name="ay")
    jx = mod.addVars(N, lb=-5, ub=5, vtype=GRB.CONTINUOUS, name="jx")
    jy = mod.addVars(N, lb=-5, ub=5, vtype=GRB.CONTINUOUS, name="jy")

    # Creating the binary arrays for big-M method depending on the number of obstacles
    d = {}
    for i in range(len(xo)):
        d['d1_' + str(i)] = mod.addVars(N, name = 'd1_'+str(i), vtype = GRB.BINARY)
        d['d2_' + str(i)] = mod.addVars(N, name = 'd2_'+str(i), vtype = GRB.BINARY) 
        d['d3_' + str(i)] = mod.addVars(N, name = 'd3_'+str(i), vtype = GRB.BINARY) 
        d['d4_' + str(i)] = mod.addVars(N, name = 'd4_'+str(i), vtype = GRB.BINARY)  

    # Defining the cost
    cost = LinExpr()
    for i in range (N):
        if i == 0:
            cost+= 100*(x[N-1]-x_f)**2
        cost += ax[i]**2 + 4*jx[i]**2 + 3*(y[i]-2.5)**2 + 0.05*(vy[i]**2) + 0.5*ay[i]**2 + (vx[i]-12)**2
    mod.setObjective( cost , GRB.MINIMIZE)

    # Set initial constraints
    mod.addConstr(x[0] == x_0,  name="initx")
    mod.addConstr(y[0] == y_0,  name="inity")
    mod.addConstr(vx[0] == vx_0, name="initvx")

    # Set constraints for position, velocities and acceleration
    mod.addConstrs((x[i] + deltaT*vx[i] + 0.5*(deltaT**2)*ax[i] + (1/6)*(deltaT**3)*jx[i] == x[i+1] for i in range(N-1)),"c1")
    mod.addConstrs((y[i] + deltaT*vy[i] + 0.5*(deltaT**2)*ay[i] + (1/6)*(deltaT**3)*jy[i] == y[i+1] for i in range(N-1)),"c2")
    mod.addConstrs((vx[i] + deltaT*ax[i] + 0.5*(deltaT**2)*jx[i]  == vx[i+1] for i in range(N-1)),"c3")
    mod.addConstrs((vy[i] + deltaT*ay[i] + 0.5*(deltaT**2)*jy[i] == vy[i+1] for i in range(N-1)),"c4")
    mod.addConstrs((ax[i] + deltaT*jx[i] == ax[i+1] for i in range(N-1)),"j3")
    mod.addConstrs((ay[i] + deltaT*jy[i] == ay[i+1] for i in range(N-1)),"j4")

    M=10000

    # Obstacle avoidance for every obstacle
    for n in range(len(yo)):
        mod.addConstrs((x[i] <= (obs['xo' + str(n)][i] - Lo) + M*(1- d['d1_' + str(n)][i]) for i in range(N)), "c5")
        mod.addConstrs((x[i] >= (obs['xo' + str(n)][i] - Lo) - M*d['d1_' + str(n)][i] for i in range(N)), "c6")

        mod.addConstrs((x[i] >= (obs['xo' + str(n)][i] + Lo) - M*(1-d['d2_' + str(n)][i]) for i in range(N)), "c5")
        mod.addConstrs((x[i] <= (obs['xo' + str(n)][i] + Lo) + M*d['d2_' + str(n)][i] for i in range(N)), "c6")

        mod.addConstrs((y[i] >= (obs['yo0' +str(n)] + Wo) - M*(1-d['d3_' + str(n)][i]) for i in range(N)), "c5")
        mod.addConstrs((y[i] <= (obs['yo0' +str(n)] + Wo) + M*d['d3_' + str(n)][i] for i in range(N)), "c6")

        mod.addConstrs((y[i] <= (obs['yo0' +str(n)] - Wo) + M*(1-d['d4_' + str(n)][i]) for i in range(N)), "c5")
        mod.addConstrs((y[i] >= (obs['yo0' +str(n)] - Wo) - M*d['d4_' + str(n)][i] for i in range(N)), "c6")

        mod.addConstrs((d['d1_' + str(n)][i] + d['d2_' + str(n)][i] + d['d3_' + str(n)][i] + d['d4_' + str(n)][i] == 1 for i in range(N)), "c13")

    # Launching the optimizer
    mod.optimize()


    # Plotting to see if it works, not necessary while in use  
    # !!! I also put the values in an array, it is required to give back the array of position, velocity and acceleration !!!

    import matplotlib.pyplot as plt
    from celluloid import Camera

    
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
            c = ['b-', 'g-', 'y-', 'k-', 'r-', 'b-', 'g-', 'y-', 'k-' ]
            
            plt.plot([x[i].X, x[i].X],[y[i].X, y[i].X], 'ko-')
            for n in range(len(yo)):
                plt.plot([obs['xo' + str(n)][i]-(Lo/2), obs['xo' + str(n)][i]-(Lo/2), obs['xo' + str(n)][i]+(Lo/2), obs['xo' + str(n)][i]+(Lo/2), obs['xo' + str(n)][i]-(Lo/2)], [obs['yo0' +str(n)]-(Wo/1.5), obs['yo0' +str(n)]+(Wo/1.5), obs['yo0' +str(n)]+(Wo/1.5), obs['yo0' +str(n)]-(Wo/1.5), obs['yo0' +str(n)]-(Wo/1.5)], c[n], linewidth=1.5)
            camera.snap()
            
        anim = camera.animate()
        plt.show()
            # plt.axis('equal')
        for i in range(N-1):    
            plt.figure(2)
            plt.plot([x[i].X, x[i+1].X],[vx[i].X, vx[i+1].X], 'ko-')
            plt.xlabel("x (m)")
            plt.ylabel("vx (m/s)")

            plt.figure(3)
            plt.plot([x[i].X, x[i+1].X],[vy[i].X, vy[i+1].X], 'ko-')
            plt.xlabel("x (m)")
            plt.ylabel("vy (m/s)")

        #     plt.figure(4)
        #     plt.plot([x[i].X, x[i+1].X],[ax[i].X, ax[i+1].X], 'ko-')
        plt.show()

    return(xr, yr, vxr, vyr, axr, ayr)


# To test the function, should be removed when using with the simulation 
x_obstacles = [100, 140]
y_obstacles = [2, 3]
v_obstacles = [-10, -10]
xv, yv, vxv, vyx, axv, ayv = avoiding_2obstacles(x_obstacles, y_obstacles, v_obstacles)

print(xv[:5], yv[:5], vxv[:5], vyx[:5], axv[:5], ayv[:5])