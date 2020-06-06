"""---------------------------------------------------------------
1D HEAT DIFFUSION SIMULATION OF A ROD
------------------------------------------------------------------

------------------------------------------------------------------
Consider a 50m rod made of silver.            
Thermal diffusivity of the Silver (kappa) = 1.65 (cm)^2/sec. 
------------------------------------------------------------------"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

#total time taken for the simulation (sec)
total_time = 1000

#total length of the rod (meters)
length_X = 50

#unit step of length considered
delta_X = 1

#unit step of time considered
delta_time = 1

#thermal diffusivitty of the material (here Silver)
kappa = 1.65

#defining the temperature at the right and left ends of the rod
T_right = 50
T_left = 100

#total number of time steps
len_time = int((total_time)/delta_time)

#total number of length steps
lenX = int((length_X/delta_X) + 1)

#defining alpha
alpha = (kappa * delta_X)/(delta_time**2)

T = np.empty((len_time, lenX))
T.fill(0)

T[:,(lenX-1):] = T_right
T[:,:1] = T_left

A = np.empty((lenX-2,lenX-2))
for i in range(lenX-2):
    for j in range(lenX-2):
        if i == j :
            A[i][j] = (1+2*alpha)
        elif i == j+1 or i == j-1:
            A[i][j] = -alpha
        else :
            A[i][j] = 0
            
for l in range(0,len_time-1):
    
    #Defining B matrix
    B = np.empty(lenX-2)               
    for i in range(lenX-2):
        B[i] = T[l][i+1]       
        if i == 0:
            B[0] += alpha * T[l][0]
        if i == lenX-3:
            B[lenX-3] += alpha * T[l][lenX-1]
    
    #Solution after increment
    x = np.empty(lenX-2)
    x = np.linalg.solve(A,B)          
    
    #Appending it in the temperature array
    for i in range(1,lenX-1):
        T[l+1][i] = x[i-1]
        
#plt.rcParams['animation.ffmpeg_path']='Your file path to ffmpeg.exe'
plt.rcParams['animation.ffmpeg_path']='C:\\Users\\karth\\.spyder-py3\\ffmpeg\\bin\\ffmpeg.exe'
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

fig_1 = plt.figure()
ax = plt.axes(xlim = (0,lenX-1), ylim= (0,100))
line, = ax.plot([], [], lw=2)
plt.xlabel('length (m)')
plt.ylabel('temperature (C)')

def init():
    line.set_data([], [])
    return line,

def animate(i):
    line.set_xdata(np.arange(0,lenX,delta_X))
    line.set_ydata(T[i])
    plt.title("Temperature at time = %i sec"  % (i*delta_time))
    return line,

anim = animation.FuncAnimation(fig_1, animate, frames=len_time,
                               init_func=init, interval=10,
                               blit=True)

anim.save('1D_heat_diffusion.mp4', fps = 20)


#code for 3d plot begins
fig_2 = plt.figure(figsize = (12,8))
ax =fig_2.add_subplot(111, projection = '3d')
x = np.arange(0,lenX)
t = np.arange(0,len_time)
Time,X = np.meshgrid(t,x)


colourMap = plt.cm.jet

colourMap = plt.cm.jet

ax.plot_wireframe(X,Time,T[:,:len_time].transpose(), cmap=colourMap, rstride=10, cstride=50)

ax.set_xlabel('Length (m)')
ax.set_ylabel('Time (sec)')
ax.set_zlabel('Temperature (C)')
ax.set_title(r"$\bf{Silver}$" + "\n" +
             " Thermal diffusivity $\kappa$ = {} $(cm)^2/s$".format(kappa),
             fontsize=16, fontname="Times New Roman Bold")
ax.view_init(azim=-50, elev = 20)
fig_2.savefig('1D_heat_diffusion.png', dpi=300, bbox_inches='tight')