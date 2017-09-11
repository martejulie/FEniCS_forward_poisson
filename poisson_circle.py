from pylab import *
from fenics import *
from mshr import *
import scipy.io as sio
from fenics2meshgrid import *
import matplotlib.pyplot as plt

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

	p_solution = solvePoisson_circle(200, 6, 80, 1.14e-3, 200) 
	plot(p_solution, interactive=True)
	Nx = 251
	Ny = 251
	x = linspace(-125, 125, Nx)
	y = linspace(-125, 125, Ny)
	p_grid, r = fenics2meshgrid(p_soltuion, 80, x, y)

	plt.imshow(p_grid)
	plt.show()

#	sio.savemat(filename, {'P':p_grid, 'R_ves':r_vessel, 'M_true':C, 'P_ves':boundary_value, 'Hx':x, 'Hy':y, 'r':r})

