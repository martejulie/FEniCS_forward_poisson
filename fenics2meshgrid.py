from fenics import *

def fenics2meshgrid(data, boundary_value, x, y):
	"""

	"""
	Nx = len(x)
	Ny = len(y)	
	X,Y = meshgrid(x,y)
	p_grid = zeros([Nx, Ny])

	for i in range(Nx):
	    for j in range(Ny):
		x_val = X[i,j]
		y_val = Y[i,j]
		point = Point(x_val, y_val)
		try:
		    p_val = p_solution(point)
		    p_grid[i,j] = p_val
		except:
		    p_val = boundary_value
		    p_grid[i,j] = boundary_value

	r = np.sqrt(X**2 + Y**2)

	return p_grid, r
