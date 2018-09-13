from fenics import *
import numpy as np

def fenics2nparray(data, boundary_value, r_ves, x, y, center):
	"""
	Writes fenics solution to numpy array.

	Arguments: 
		data (fenics solution): input data 
		boundary_value (int): if point from numpy array not found on fenics mesh, value is set to boundary_value
		x (array): x values of numpy array
		y (array): y values of numpy array  

	Returns: 
		data_grid (array): output data
		r (array): sqrt(X**2 + Y**2)
	"""
        Nx = len(y)
        Ny = len(x)	
        X,Y = np.meshgrid(x,y)
        data_grid = np.zeros([Nx, Ny])

        if len(center) == 1:
            r0 = np.sqrt((X-center[0][0])**2 + (Y-center[0][1])**2)
            r1 = np.ones([Nx,Ny])*1e3
            r2 = np.ones([Nx,Ny])*1e3
        elif len(center) == 2:
            r0 = np.sqrt((X-center[0][0])**2 + (Y-center[0][1])**2)
            r1 = np.sqrt((X-center[1][0])**2 + (Y-center[1][1])**2)
            r2 = np.ones([Nx,Ny])*1e3
        elif len(center) == 3:
            r0 = np.sqrt((X-center[0][0])**2 + (Y-center[0][1])**2)
            r1 = np.sqrt((X-center[1][0])**2 + (Y-center[1][1])**2)
            r2 = np.sqrt((X-center[2][0])**2 + (Y-center[2][1])**2)
        else:
            print "Mesh must only have three holes or less."
            return 

        for i in range(Nx):
            for j in range(Ny):
                x_val = X[i,j]
                y_val = Y[i,j]
                point = Point(x_val, y_val)
                try:
                    data_val = data(point)
                    data_grid[i,j] = data_val
                except:
                    if r0[i,j] < r_ves:
                        data_val = boundary_value[0]
                        data_grid[i,j] = boundary_value[0]
                        print "Point found within hole. Value set to", boundary_value[0]
                    elif r1[i,j] < r_ves:
                        data_val = boundary_value[1]
                        data_grid[i,j] = boundary_value[1]
                        print "Point found within hole. Value set to", boundary_value[1]
                    elif r2[i,j] < r_ves:
                        data_val = boundary_value[2]
                        data_grid[i,j] = boundary_value[2]
                        print "Point found within hole. Value set to", boundary_value[2]
        
        if len(center) == 1:
            return data_grid, r0
        elif len(center) == 2:
            return data_grid, r0, r1
        elif len(center) == 3:
            return data_grid, r0, r1, r2
