from fenics import *
#from mshr import *
from poisson_rectangle import solvePoisson_rectangle
import numpy as np
import scipy.io as sio
from fenics2nparray import *
import matplotlib.pyplot as plt

resolution = 900

corners = [[0, 0], [2, 2]]
centre = 1
hole_coor = [[1, 1]]    	

R_star = 141.0	
M_star = 1.0e-3

r_ves = 6/R_star
p_ves = [80/(M_star*R_star**2)]
M = Expression("1-0.5*(x[1]<1)", degree=4)
         
p_solution, mesh = solvePoisson_rectangle(corners, hole_coor, r_ves, p_ves, M, resolution)
mesh_coor =  mesh.coordinates()

meshfig = plot(mesh)
plt.show()
#meshfig.write_png("myMesh_twovessels")
fig = plot(p_solution, title="Ground truth pO2 values")
#fig = plot(p_solution, interactive=True, title="Ground truth pO2 values")
plt.show()
#fig.write_png("mymesh_p")

d = 0.007
filename = 'groundTruth_varyingM'
x = np.arange(0, 2.0001, d)
y = np.arange(0, 2.0001, d)

p_grid, r1= fenics2nparray(p_solution, p_ves[0], x, y, hole_coor)

plt.imshow(p_grid)
plt.show()

sio.savemat(filename, {'P':p_grid, 'r':r1, 'd':d, 'M_star':M_star, 'Hx':x, 'Hy':y, 'res':resolution})
