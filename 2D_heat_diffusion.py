"""---------------------------------------------------------------
2D HEAT DIFFUSION SIMULATION OF A SQUARE PLATE
------------------------------------------------------------------

------------------------------------------------------------------
Consider a 250m * 250m square plate made of silver.            
Thermal diffusivity of the Silver (kappa) = 1.65 (cm)^2/sec. 
------------------------------------------------------------------"""



import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


#total time taken for the simulation (sec)
total_time = 2000

# total length of the plate in X-direction
length_X = 250

#total length of the plate in Y-direction
length_Y = 250

#unit step of length in X-direction 
delta_X = 10

#unit step of length in Y-direction
delta_Y = 10

#unit step of time
delta_time = 10

#thermal diffusivity of the material (here Silver)
kappa = 1.65

#defining the boundary temperatures
T_top = 100
T_bottom = 0
T_right = 50
T_left = 75


#total number of time divisions
len_time = int((total_time)/delta_time)        

#total number of divisions in x direction 
lenX = int((length_X/delta_X) + 1)   

#total number of divisions in y direction
lenY = int((length_Y/delta_Y) + 1)   

#defining alpha
alpha = (kappa * delta_X)/(delta_time**2)

T = np.empty((len_time, lenY, lenX))
T.fill(0)

T[:,(lenY-1):,:] = T_top
T[:,:1,:] = T_bottom
T[:,:,(lenX-1):] = T_right
T[:,:,:1] = T_left


A = np.empty((lenY-2,lenX-2))
for i in range(lenX-2):
    for j in range(lenY-2):
        if i == j :
            A[i][j] = 2*(1+alpha)
        elif i == j+1 or i == j-1:
            A[i][j] = -alpha
        else :
            A[i][j] = 0
            

for l in range(0,len_time-2,2):

    D = np.empty((lenY-2,lenX-2))   
    
    #Defining D matrix (eqn 29)
    for i in range(lenY-2):
        for j in range(lenX-2):
            D[i][j] = alpha * T[l][i][j+1] + 2*(1-alpha)* T[l][i+1][j+1] + alpha * T[l][i+2][j+1]        
            if j == 0:
                D[i][0] += alpha * T[l][i+1][0]
            if j == lenX-3:
                D[i][lenX-3] += alpha * T[l][i+1][lenX-1]
    
    #Solution after half increment
    x = np.empty((lenY-2,lenX-2))
    for i in range(lenX-2):
        x[i] = np.linalg.solve(A,D[i])          
     
    #Appending it in the odd temperature array  
    for i in range(lenX-2):
        for j in range(lenY-2):                  
            T[l+1][i+1][j+1] = x[i][j]          
        
    #Defining B matrix (eqn 36)   
    B = np.empty((lenY-2,lenX-2))               
    for i in range(lenY-2):
        for j in range(lenX-2):
            B[i][j] = alpha * T[l+1][j+1][i] + 2*(1-alpha)* T[l+1][j+1][i+1] + alpha * T[l+1][j+1][i+2]        
            if j == 0:
                B[i][0] += alpha * T[l+1][0][i+1]
            if j == lenX-3:
                B[i][lenX-3] += alpha * T[l+1][lenX-1][i+1]
    
    #Solution after full increment
    y = np.empty((lenY-2,lenX-2))              
    for i in range(lenX-2):
        y[i] = np.linalg.solve(A,B[i])
        
    #Appending it in the even temperature array
    for i in range(lenX-2):
        for j in range(lenY-2):                
            T[l+2][i+1][j+1] = y[j][i]      


if len_time%2 == 0:
    T = T[:len_time-1]

Temp = [T[0]]
for i in range(len(T)-2):
    Temp.append(T[i+2])


#plt.rcParams['animation.ffmpeg_path']='Your file path to ffmpeg.exe'
plt.rcParams['animation.ffmpeg_path'] ='C:\\Users\\karth\\.spyder-py3\\ffmpeg\\bin\\ffmpeg.exe'
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)


X, Y =np.meshgrid(np.arange(0,lenX), np.arange(0,lenY))
colorinterpolation = 50
colourMap = plt.cm.jet
fig = plt.figure()
plt.contourf(X,Y,Temp[0],colorinterpolation, cmap = colourMap)
plt.colorbar()
plt.xlabel('X (10m)')
plt.ylabel('Y (10m)')

def animate(i):
    cont = plt.contourf(X,Y,Temp[i],colorinterpolation, cmap = colourMap)
    plt.title(r"$\bf{Silver}$" + "\n" + "$\kappa$ = 1.65 $(cm)^2/s$   time = %i sec"  % (i*delta_time))
    return cont

anim = animation.FuncAnimation(fig,animate,frames=len(Temp),repeat=True)

anim.save('2D_heat_diffusion.mp4', fps = 10)
