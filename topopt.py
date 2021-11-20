from dolfin import *
import numpy as np, sklearn.metrics.pairwise as sp
# A 55 LINE TOPOLOGY OPTIMIZATION CODE ---------------------------------
def main(nelx, nely, volfrac, penal, rmin):
    sigma = lambda _u: 2.0 * mu * sym(grad(_u)) + lmbda * tr(sym(grad(_u))) * Identity(len(_u))
    psi = lambda _u: lmbda / 2 * (tr(sym(grad(_u))) ** 2) + mu * tr(sym(grad(_u)) * sym(grad(_u)))
    xdmf = XDMFFile("output/density.xdmf")
    mu, lmbda = Constant(0.3846), Constant(0.5769)
    # PREPARE FINITE ELEMENT ANALYSIS ----------------------------------
    mesh = RectangleMesh(Point(0, 0), Point(nelx, nely), nelx, nely, "right/left")
    U = VectorFunctionSpace(mesh, "P", 1)
    D = FunctionSpace(mesh, "DG", 0)
    u, v = TrialFunction(U), TestFunction(U)
    u_sol, density, density_old, density_new = Function(U), Function(D, name="density"), Function(D), Function(D)
    density.vector()[:] = volfrac
    V0 = assemble(1 * TestFunction(D) * dx)
    # DEFINE SUPPORT ---------------------------------------------------
    support = CompiledSubDomain("near(x[0], 0.0, tol) && on_boundary", tol=1e-14)
    bcs = [DirichletBC(U, Constant((0.0, 0.0)), support)]
    # DEFINE LOAD ------------------------------------------------------
    load_marker = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    CompiledSubDomain("x[0]==l && x[1]<=1", l=nelx).mark(load_marker, 1)
    ds = Measure("ds")(subdomain_data=load_marker)
    F = dot(v, Constant((0.0, -1.0))) * ds(1)
    # SET UP THE VARIATIONAL PROBLEM AND SOLVER ------------------------
    K = inner(density ** penal * sigma(u), grad(v)) * dx
    solver = LinearVariationalSolver(LinearVariationalProblem(K, F, u_sol, bcs))
    # PREPARE DISTANCE MATRICES FOR FILTER -----------------------------
    midpoint = [cell.midpoint().array()[:] for cell in cells(mesh)]
    distance_mat = np.maximum(rmin - sp.euclidean_distances(midpoint, midpoint), 0)
    distance_sum = distance_mat.sum(1)
    # START ITERATION --------------------------------------------------
    loop, change = 0, 1
    while change > 0.01 and loop < 100:
        loop = loop + 1
        density_old.assign(density)
        # FE-ANALYSIS --------------------------------------------------
        solver.solve()
        # OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS ------------------
        objective = density ** penal * psi(u_sol)
        sensitivity = project(-diff(objective, density), D).vector()[:]
        # FILTERING/MODIFICATION OF SENSITIVITIES ----------------------
        sensitivity = np.divide(distance_mat @ np.multiply(density.vector()[:], sensitivity), np.multiply(density.vector()[:], distance_sum))
        # DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD --------------
        l1, l2, move = 0, 100000, 0.2
        while l2 - l1 > 1e-4:
            l_mid = 0.5 * (l2 + l1)
            density_new.vector()[:] = np.maximum(0.001, np.maximum(density.vector()[:] - move, np.minimum(1.0, np.minimum(density.vector()[:] + move, density.vector()[:] * np.sqrt(-sensitivity / V0 / l_mid)))))
            current_vol = assemble(density_new * dx)
            l1, l2 = (l_mid, l2) if current_vol > volfrac * V0.sum() else (l1, l_mid)
        # PRINT RESULTS ------------------------------------------------
        change = max(density_new.vector()[:] - density_old.vector()[:])
        print("it.: {0} , obj.: {1:.3f} Vol.: {2:.3f}, ch.: {3:.3f}".format(loop, project(objective, D).vector().sum(), current_vol / V0.sum(), change))
        density.assign(density_new)
        xdmf.write(density, loop)
