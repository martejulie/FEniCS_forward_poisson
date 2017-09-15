from pylab import *
from fenics import *
from mshr import *
import scipy.io as sio
from fenics2meshgrid import *
import matplotlib.pyplot as plt

def solvePoisson_circle(r_domain, r_hole, hole_boundary_value, M, resolution):
	"""
	Solves the Poisson equation 
	nabla**2 p = M
	on a circular mesh with a centered hole
	
	Boundary conditions: 
		dp/dr = 0 when r = r_domain
		p = hole_boundary_value when r = r_hole

	Arguments:
		r_domain (int): mesh radius
		r_hole (int): hole radius
		hole_boundary_value (int): value of p at r = r_hole
		M (float): constant
		resolution (int): resolution of fenics mesh

	Returns:
		p_solution (fenics solution)
		mesh (fenics mesh)
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

	return p_solution, mesh


if __name__ == "__main__":

	filename = "circleMesh_res200_d1"	
	M = 1.14e-3
	resolution = 200
	p_solution, mesh = solvePoisson_circle(200, 6, 80, M, resolution) 
 
	meshfig = plot(mesh, interactive=True)
	meshfig.write_png("firstMesh")	
	
	fig = plot(p_solution, interactive=True, title="Ground truth pO2 values")
	fig.write_png("po2fenics_firstMesh")


	d = 1	
	N = 251

#	d = 5
#	N = 51

#	d = 10
#	N = 26

	x = linspace(-125, 125, N)
	y = linspace(-125, 125, N)
	p_grid, r = fenics2numpyArray(p_solution, 80, x, y)

	plt.imshow(p_grid)
	plt.show()

	sio.savemat(filename, {'P':p_grid, 'r':r, 'd':d, 'M_true':M, 'Hx':x, 'Hy':y, 'res':resolution})

