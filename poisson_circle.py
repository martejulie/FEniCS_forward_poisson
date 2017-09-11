from pylab import *
from fenics import *
from mshr import *
import scipy.io as sio

def solvePoisson_circle(r_domain, r_hole, hole_boundary_value, M, resolution):
	"""
	Docstring
	"""

	big_circle = Circle(Point(0,0), r_domain)
	small_circle = Circle(Point(0,0), r_hole)

	domain = big_circle - small_circle
	mesh = generate_mesh(domain, resolution)

	V = FunctionSpace(mesh, 'CG', 1)

	p = TrialFunction(V)
	M = Constant(M)
	v = TestFunction(V)
	form = (inner(nabla_grad(p), nabla_grad(v)) + M*v )*dx
	(a,L) = system(form)

	def boundary(x, on_boundary):
	    r = np.sqrt((x[0])**2 + (x[1])**2)
	    b = ((r < r_hole+5) and on_boundary)
	    return b

	bc = DirichletBC(V, hole_boundary_value, boundary)

	p_solution = Function(V)	
	solve(a==L, p_solution, bc)

	return p_solution


if __name__ == "__main__":

	plot(p_solution, interactive=True)

	Nx = 251
	Ny = 251
	x = linspace(-125, 125, Nx)
	y = linspace(-125, 125, Ny)

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
	import matplotlib.pyplot as plt
	plt.imshow(p_grid)
	plt.show()

	filename = "oneVesselReady2Go.mat"
	sio.savemat(filename, {'P':p_grid, 'R_ves':r_vessel, 'M_true':C, 'P_ves':boundary_value, 'Hx':x, 'Hy':y, 'r':r})
