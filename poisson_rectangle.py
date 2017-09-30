from fenics import *
from mshr import *
import numpy as np
import scipy.io as sio
from fenics2nparray import *
import matplotlib.pyplot as plt

def solvePoisson_rectangle(corners, hole_coor, r_hole, hole_boundary_value, M, resolution):
    """
    Solves the Poissons equtaion 
    nabla**2 p = M
    on a rectangle mesh with one or two holes.

    Boundary conditions: 
	dp/dr = 0 at outer boundaries
	p = hole_boundary_value at hole boundaries
    
    Arguments:
	corners (array): mesh corners coordinates
	r_hole (int): hole radius
	hole_boundary_value (array): value of p at hole boundaries 
	M (float): constant
	resolution (int): resolution of fenics mesh

    Returns:
	p_solution (fenics solution)
	mesh (fenics mesh)
    """
    
    r = Rectangle(Point(corners[0][0],corners[0][1]), Point(corners[1][0],corners[1][1]))  
    domain = r
    for i in range(len(hole_coor)):
        hole = Circle(Point(hole_coor[i][0], hole_coor[i][1]), r_hole)
	domain = domain - hole

    mesh = generate_mesh(domain, resolution)

    V = FunctionSpace(mesh, 'CG', 1)

    p = TrialFunction(V)
    M = Constant(M)
    v = TestFunction(V)
    form = (inner(nabla_grad(p), nabla_grad(v)) + M*v )*dx
    (a,L) = system(form)

    def boundary1(x, on_boundary):
        r = np.sqrt((x[0]-hole_coor[0][0])**2 + (x[1]-hole_coor[0][1])**2)
        b = ((r < r_hole+5) and on_boundary)
        return b
    bc1 = DirichletBC(V, hole_boundary_value[0], boundary1)
    bcs = [bc1]	

    if len(hole_coor)>1:
        def boundary2(x, on_boundary):
	    r = np.sqrt((x[0]-hole_coor[1][0])**2 + (x[1]-hole_coor[1][1])**2)
	    b = ((r < r_hole+5) and on_boundary)
	    return b
    	bc2 = DirichletBC(V, hole_boundary_value[1], boundary2)
        bcs.append(bc2)

    p_solution = Function(V)
    solve(a==L, p_solution, bcs)

    return p_solution, mesh


if __name__ == "__main__":

    filename = 'rectangleMesh_res900_d1.mat'
    resolution = 900

    corners = [[0, 0], [1000, 1000]]
    hole_coor = [[500., 500.]]

    r_ves = 6.
    p_ves = [80.]
    M = 1.33e-4
             
    p_solution, mesh = solvePoisson_rectangle(corners, hole_coor, r_ves, p_ves, M, resolution)
    mesh_coor =  mesh.coordinates()
    
#    meshfig = plot(mesh, interactive=True)
#    meshfig.write_png("po2fenics_firstMesh_rectangular")
#    fig = plot(p_solution, interactive=True, title="Ground truth pO2 values")
#    fig.write_png("po2fenics_firstExample_rectangular")

    d = 1
    N = 251
    	
#    d = 10
#    N = 26	

    x = np.linspace(375, 625, N)
    y = np.linspace(375, 625, N)
    p_grid, r = fenics2nparray(p_solution, 80, x, y, hole_coor)

    plt.imshow(p_grid)
    plt.show()
   
    sio.savemat(filename, {'P':p_grid, 'r':r, 'd':d, 'M_true':M, 'Hx':x, 'Hy':y, 'res':resolution})
