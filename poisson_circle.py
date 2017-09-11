from pylab import *
from fenics import *
from mshr import *
import scipy.io as sio

r_domain = 200
r_vessel = 6

big_circle = Circle(Point(0,0), r_domain)
small_circle = Circle(Point(0,0), r_vessel)

domain = big_circle - small_circle

resolution = 800
mesh = generate_mesh(domain, resolution)

V = FunctionSpace(mesh, 'CG', 1)

p = TrialFunction(V)
C = 1.14e-3
M = Constant(C)
v = TestFunction(V)
form = (inner(nabla_grad(p), nabla_grad(v)) + M*v )*dx
(a,L) = system(form)

def boundary1(x, on_boundary):
    r = np.sqrt((x[0])**2 + (x[1])**2)
    #b = ((r < hole_radius+DOLFIN_EPS) and on_boundary)
    b = ((r < r_vessel+5) and on_boundary)
    return b

boundary_value = 80
bc1 = DirichletBC(V, boundary_value, boundary1)
bcs = [bc1]

p_solution = Function(V)
solve(a==L, p_solution, bcs)
plot(p_solution, interactive=True)

if __name__ == "__main__":

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
