from fenics import *
from mshr import *
import numpy as np
import scipy.io as sio

def solvePoisson_rectangle(corners, hole_coor, hole_radius, hole_boundary_value, M):
    """
    Solves the Poissons equtaion 
    nabla**2 p = M
    on a rectangle mesh with two holes.
    """
    
    r = Rectangle(Point(corners[0][0],corners[0][1]), Point(corners[1][0],corners[1][1]))
    circle1 = Circle(Point(hole_coor[0][0], hole_coor[0][1]), hole_radius)
    circle2 = Circle(Point(hole_coor[1][0], hole_coor[1][1]), hole_radius)

    domain = r - circle1 - circle2

    resolution = 200
    mesh = generate_mesh(domain, resolution)

    V = FunctionSpace(mesh, 'CG', 1)

    p = TrialFunction(V)
    M = Constant(M)
    v = TestFunction(V)
    form = (inner(nabla_grad(p), nabla_grad(v)) + M*v )*dx
    (a,L) = system(form)

    def boundary1(x, on_boundary):
        r = np.sqrt((x[0]-hole_coor[0][0])**2 + (x[1]-hole_coor[0][1])**2)
        b = ((r < hole_radius+5) and on_boundary)
        return b

    def boundary2(x, on_boundary):
        r = np.sqrt((x[0]-hole_coor[1][0])**2 + (x[1]-hole_coor[1][1])**2)
        b = ((r < hole_radius+5) and on_boundary)
        return b

    bc1 = DirichletBC(V, hole_boundary_value[0], boundary1)
    bc2 = DirichletBC(V, hole_boundary_value[1], boundary2)
    bcs = [bc1, bc2]

    p_solution = Function(V)
    solve(a==L, p_solution, bcs)

    p_array = p_solution.vector().array()
    p_array = p_array[vertex_to_dof_map(V)]

    return mesh, p_solution, p_array


if __name__ == "__main__":

    filename = 'po2FenicsSolution.mat'

    #corners = [[0, 0], [184, 216]]
    #hole_coor = [[72., 112.], [128., 56.]]

    #r_ves = 2.
    #p_ves = [81., 85.]
    #M = 4.86e-3

    corners = [[0, 0], [800, 800]]
    hole_coor = [[300., 300.], [500., 500.]]
    #corners = [[-400, -400], [400, 400]]
    #hole_coor = [[-100., -100.], [100., 100.]]

    r_ves = 6.
    p_ves = [80., 80.]
    M = 3.54e-4
    #M = 3.54e-6
             
    mesh, p_solution, p_array = solvePoisson(corners, hole_coor, r_ves, p_ves, M)
    mesh_coor =  mesh.coordinates()
    
    meshfig = plot(mesh, interactive=True)
    #meshfig.write_png("po2fenics_firstMesh")
    fig = plot(p_solution, interactive=True, title="Ground truth pO2 values")
    #fig.write_png("po2fenics_firstExample")

    sio.savemat(filename, {'mesh_coor':mesh_coor, 'p_array':p_array, 'corners':corners, 'hole_coor':hole_coor, 'r_ves':r_ves, 'M_true':M, 'p_ves':p_ves, 'corners':corners})
