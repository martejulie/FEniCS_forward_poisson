from fenics import *
from mshr import *
import numpy as np
import scipy.io as sio

filename = 'po2grid.mat'

x0 = 0.
x1 = 200.
y0 = 0.
y1 = 200.
corners = [x0, y0, x1, y1]
r = Rectangle(Point(corners[0],corners[1]), Point(corners[2],corners[3]))

R_ves = 6.
p_ves = [70., 60.]
C = 5.33e-3

center1 = [75., 75.]
center2 = [125., 160.]

resolution = 150

# Stop making changes
# *******************************
circle1 = Circle(Point(center1[0], center1[1]), R_ves)
circle2 = Circle(Point(center2[0], center2[1]), R_ves)

domain = r - circle1 - circle2

mesh = generate_mesh(domain, resolution)

meshfig = plot(mesh, interactive=True)
#meshfig.write_png("po2fenics_firstMesh")

V = FunctionSpace(mesh, 'CG', 1)

p = TrialFunction(V)
M = Constant(C)
v = TestFunction(V)
form = (inner(nabla_grad(p), nabla_grad(v)) + M*v )*dx
(a,L) = system(form)

def boundary1(x, on_boundary):
    r = np.sqrt((x[0]-center1[0])**2 + (x[1]-center1[1])**2)
    b = ((r < R_ves+DOLFIN_EPS) and on_boundary)
    return b

def boundary2(x, on_boundary):
    r = np.sqrt((x[0]-center2[0])**2 + (x[1]-center2[1])**2)
    b = ((r < R_ves+DOLFIN_EPS) and on_boundary)
    return b


bc1 = DirichletBC(V, p_ves[0], boundary1)
bc2 = DirichletBC(V, p_ves[1], boundary2)
bcs = [bc1, bc2]

p_solution = Function(V)
solve(a==L, p_solution, bcs)

fig = plot(p_solution, interactive=True, title="Ground truth pO2 values")
#fig.write_png("po2fenics_firstExample")

mesh_coor =  mesh.coordinates()
p_array = p_solution.vector().array()

p_array = p_array[vertex_to_dof_map(V)]

sio.savemat(filename, {'xvec':mesh_coor[:,0], 'yvec':mesh_coor[:,1], 'xyvec':mesh_coor, 'cvec':p_array, 'center1':center1, 'center2':center2, 'r_ves':R_ves, 'M_true':C, 'p_ves':p_ves, 'corners':corners})
