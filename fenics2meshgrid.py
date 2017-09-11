from fenics import *
import numpy as np

def fenics2meshgrid(data, boundary_value, x, y):
	"""

	"""
	Nx = len(x)
	Ny = len(y)	
	X,Y = np.meshgrid(x,y)
	data_grid = np.zeros([Nx, Ny])

	for i in range(Nx):
	    for j in range(Ny):
		x_val = X[i,j]
		y_val = Y[i,j]
		point = Point(x_val, y_val)
		try:
	            data_val = data(point)
		    data_grid[i,j] = data_val
		except:
		    data_val = boundary_value
		    data_grid[i,j] = boundary_value

	r = np.sqrt(X**2 + Y**2)

	return data_grid, r
