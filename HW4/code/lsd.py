import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
import time

from tri import solve_tridiagonal

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

def main(condition=0, show_plots=True, show_prints=True):
    # -------------------------------------------------------------------------
    # Parameters (translated from lsd.m)
    # -------------------------------------------------------------------------
    thickness = 0.10                         # Airfoil thickness

    # SLOR parameters
    omega = 1.97
    itmax = 1000
    itplot = 200
    
    gamma = 1.4

    # Mesh input parameters
    
    # case from hw3
    if condition == 1:
        j_le = 33     # leading edge index
        j_te = 63     # trailing edge index
        jmax = 95     # total x points
        kmax = 33     # total y points
        xsf = 1.18    # x stretching factor
        ysf = 1.18    # y stretching factor
        airfoil_type = 0   # airfoil type - 0 for NACA0010 and 1 for biconvex airfoil
        minf = 0.0
    elif condition == 7:
        j_le = 63
        j_te = 123
        jmax = 185
        kmax = 63
        xsf = 1.1
        ysf = 1.0
        airfoil_type = 0   # airfoil type - 0 for NACA0010 and 1 for biconvex airfoil
        minf = 0.8
    else:
        j_le = 33
        j_te = 73
        jmax = 105
        kmax = 43
        xsf = 1.2
        ysf = 1.2
        if condition == 2:
            minf = .75
            airfoil_type = 0
        elif condition == 3:
            minf = 0.75
            airfoil_type = 1
        elif condition == 4:
            minf = 0.8
            airfoil_type = 0
        elif condition == 5:
            minf = 0.8
            airfoil_type = 1
        elif condition == 6:
            minf = 0.75
            airfoil_type = 0
            ysf = 1.0


    
    kconst = 3      # number of constant-spaced mesh points
    dxdy = 1.0      # ratio of dx to dy at airfoil surface

    # Upper wall boundary?
    if condition == 6:
        iwall = 1
    else:
        iwall = 0
    if iwall == 1:
        ysf = 1.0
    ibc = 0         # Use full bc if 0, or simplified if 1

    # -------------------------------------------------------------------------
    # Mesh Generation
    # -------------------------------------------------------------------------
    # Generate x-array
    dx1 = 1.0 / (j_te - j_le + 1)      # equal spacing along airfoil
    dy1 = dx1 / dxdy
    # Create x coordinates (using MATLAB-style definition but converting indices)
    x = np.zeros(jmax)
    for j in range(jmax):
        # MATLAB: x(j) = (j - j_le)*dx1 + 0.5*dx1, j=1...jmax => Python index: j = 0...jmax-1, j+1 used
        x[j] = ((j + 1) - j_le) * dx1 + 0.5 * dx1

    # Upstream stretching (MATLAB: for ji=j_le-kconst-1:-1:1)
    # Convert: for i from (j_le - kconst - 1) down to 1 (MATLAB) => Python indices: from (j_le - kconst - 1 - 1) down to 0.
    for i in range(j_le - kconst - 2, -1, -1):
        x[i] = x[i+1] + (x[i+1] - x[i+2]) * xsf


    # Downstream stretching (MATLAB: for ji=j_te+kconst+1:1:jmax)
    for i in range(j_te + kconst, jmax):
        x[i] = x[i-1] + (x[i-1] - x[i-2]) * xsf

    # Generate y-array
    y = np.zeros(kmax)
    for k in range(kmax):
        # MATLAB: y(k) = (k-1)*dy1 - 0.5*dy1, k=1...kmax
        y[k] = (k) * dy1 - 0.5 * dy1

    # Stretching in y-direction for interior points (MATLAB: for ki=kconst+1:kmax-kconst)
    for k in range(kconst, kmax - kconst):
        y[k] = y[k-1] + (y[k-1] - y[k-2]) * ysf
    # Outermost points: uniform spacing
    for k in range(kmax - kconst, kmax):
        y[k] = y[k-1] + (y[k-1] - y[k-2])

    # Create 2D mesh arrays
    xmesh = np.zeros((jmax, kmax))
    ymesh = np.zeros((jmax, kmax))
    for k in range(kmax):
        for j in range(jmax):
            xmesh[j, k] = x[j]
            ymesh[j, k] = y[k]

    # Calculate spacing in x-direction
    dx = np.full(jmax, dx1)
    for j in range(1, jmax-1):
        dx[j] = 0.5 * (x[j+1] - x[j-1])
    dx[0] = 2 * dx[1] - dx[2]
    dx[-1] = 2 * dx[-2] - dx[-3]

    # Calculate arrays for mesh scaling in x-direction
    dxp2 = 1.0 / (dx ** 2)
    dxm2 = 1.0 / (dx ** 2)
    for j in range(1, jmax-1):
        dxp2[j] = 1 / (dx[j]* (x[j+1] - x[j]))
        dxm2[j] = 1 / (dx[j]* (x[j] - x[j-1]))
        
    # set initial things
    dxp2[0] = dxp2[1]
    dxm2[0] = dxm2[1]

    # Calculate spacing in y-direction
    dy = np.full(kmax, dy1)
    dy[0] = y[1] - y[0]
    for k in range(1, kmax-1):
        dy[k] = 0.5 * (y[k+1] - y[k-1])
    dy[-1] = 2 * dy[-2] - dy[-3]

    # Calculate arrays for mesh scaling in y-direction
    dyp2 = 1.0 / (dy ** 2)
    dym2 = 1.0 / (dy ** 2)
    for k in range(1, kmax-1):
        dyp2[k] = 1 / (dy[k] * (y[k+1] - y[k]))
        dym2[k] = 1 / (dy[k] * (y[k] - y[k-1]))

    # -------------------------------------------------------------------------
    # Initialization
    # -------------------------------------------------------------------------
    phi = np.zeros((jmax, kmax))
    # Incompressible potential data for NACA0010 (Abbott and von Doenhoff)
    potx = np.array([0.0, 0.005, 0.0125, 0.025, 0.050, 0.075, 0.1, 0.15, 0.2, 0.25,
                     0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0])
    potcp = np.array([1.0, 0.282, -0.061, -0.237, -0.325, -0.341, -0.341, -0.341,
                      -0.329, -0.309, -0.284, -0.237, -0.190, -0.138, -0.094, -0.040,
                      0.04, 0.075, 1.0]) 
    if show_prints:
        print("DX1 is", dx1, "DY1 is", dy1)
        print("XMAX is", x[-1], "YMAX is", y[-1])
    
    # Plot the initial mesh
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.plot_wireframe(xmesh, ymesh, phi, rstride=2, cstride=2)
    ax1.set_title('Mesh')
    ax1.set_xlabel('x/c')
    ax1.set_ylabel('y/c')
    

    # Pre-allocate arrays for SLOR iterations
    a = np.zeros(kmax)
    b = np.zeros(kmax)
    c = np.zeros(kmax)
    f_arr = np.zeros(kmax)
    cpg = np.zeros((jmax, kmax))
    cp = np.zeros(jmax)
    cpu = np.zeros(jmax)
    res = np.zeros((jmax, kmax))
    l2reshist = np.zeros(itmax)
    
    istop = False
    l2res1 = None
    supersonic_points = np.zeros(itmax) 
    # -------------------------------------------------------------------------
    # SLOR Iterations
    # -------------------------------------------------------------------------
    for it in range(1, itmax+1):
        if istop:
            break

        # --- Apply boundary condition along airfoil (k=0 corresponds to k=1 in MATLAB) ---
        bc = np.zeros(jmax)
        # Loop from j_le to j_te (MATLAB indices) => Python indices: j_le-1 to j_te-1
        for j in range(j_le - 1, j_te):
            # NACA00xx airfoil
            if airfoil_type == 0:
                x_int = 1.008930411365  # NACA00xx airfoil slope constant
                # Avoid division by zero in sqrt; x[j] is positive in the airfoil region
                dphidy = 5 * thickness * (0.2969 * 0.5 * sqrt(x_int / x[j])
                                        - 0.126 * x_int
                                        - 0.3516 * 2 * x_int**2 * x[j]
                                        + 0.2843 * 3 * x_int**3 * x[j]**2
                                        - 0.1015 * 4 * x_int**4 * x[j]**3)
            # Biconvex airfoil
            elif airfoil_type == 1:
                dphidy = 2 * thickness * (1 - 2 * x[j])

            # Compute velocity in x-direction using central differences
            velx = 0.5 * (
                ((phi[j+1, 0] - phi[j, 0]) / dx[j]) * (dxp2[j] / dxm2[j]) +
                ((phi[j, 0] - phi[j-1, 0]) / dx[j]) * (dxm2[j] / dxp2[j])
            )
            if ibc == 1:
                velx = 0.0  # simplified boundary condition
            bc[j] = -dphidy * (1.0 + velx) * dy[0]

        # Apply airfoil boundary condition (k=0)
        for j in range(jmax):
            phi[j, 0] = phi[j, 1] + bc[j]

        # Upper wall boundary condition if applicable (kmax-1 corresponds to kmax in MATLAB)
        if iwall == 1:
            for j in range(jmax):
                phi[j, kmax-1] = phi[j, kmax-2]

        # Compute A matrix
        A = np.zeros((jmax, kmax))
        for j in range(1, jmax-1):
            for k in range(1, kmax-1):
                A[j, k] = 1 - minf**2 - minf**2*(gamma + 1) * 1 * 1/(x[j+1] - x[j-1]) * ( (phi[j+1, k] - phi[j, k]) * (x[j] - x[j-1]) / (x[j+1] - x[j]) \
                    + (phi[j, k] - phi[j-1, k]) * (x[j+1] - x[j]) / (x[j] - x[j-1]))

                
        # --- Calculate residual and its L2 norm ---
        l2res = 0.0
        for k in range(1, kmax-1):
            for j in range(1, jmax-1):
                res[j, k] = (A[j,k] *
                             ((phi[j+1, k] - phi[j, k]) * dxp2[j] - (phi[j, k] - phi[j-1, k]) * dxm2[j])
                             + ((phi[j, k+1] - phi[j, k]) * dyp2[k] - (phi[j, k] - phi[j, k-1]) * dym2[k]))
                
                l2res += res[j, k] ** 2
        l2res = sqrt(l2res / (jmax * kmax))
        if it == 1:
            l2res1 = l2res
        
        # Stop condition if residual changes too much
        if l2res1 / l2res >= 10000 or l2res1 / l2res <= 1.0 / 1000:
            istop = True
        
        if show_prints:
            print(f"Iteration {it} : L2(RES) = {l2res}")
        
        mu = np.where(A >= 0, 0, 1)
        omega_jk = np.ones((jmax, kmax)) * omega
        # omega_jk[A < 0] = 1
        
        
        # --- SLOR update using tridiagonal solver along each interior x-row ---
        dphi = np.zeros((jmax, kmax))
        for j in range(1, jmax-1):
            # Set up the tridiagonal system for row j
            # k = 0 (airfoil boundary)
            a[0] = 0.0
            b[0] = 1.0
            c[0] = -1.0
            f_arr[0] = 0.0
            # Interior points k = 1 to kmax-2
            for k in range(1, kmax-1):
                a[k] = dym2[k]
                b[k] = -(1  - mu[j,k]) *(dym2[k] + dyp2[k] + A[j,k] * (dxm2[j] + dxp2[j])) + mu[j-1,k] * A[j-1,k] * dxp2[j-1]
                c[k] = dyp2[k]
                f_arr[k] = -omega_jk[j,k] * (1  - mu[j,k]) *  A[j,k] * (res[j, k] + dphi[j-1, k] * dxm2[j]) - mu[j-1,k] * omega_jk[j,k] * A[j-1,k] * ((dxp2[j-1] + dxm2[j-1]) * dphi[j-1, k] - dxm2[j-1] * dphi[j-2, k])
            # k = kmax-1 (upper boundary)
            a[kmax-1] = 0.0
            b[kmax-1] = 1.0
            c[kmax-1] = 0.0
            f_arr[kmax-1] = 0.0
            if iwall == 1:
                a[kmax-1] = -1.0
            # Solve tridiagonal system for this row
            dphi[j, :] = solve_tridiagonal(a.copy(), b.copy(), c.copy(), f_arr.copy())
        
        # Update solution
        phi = phi + dphi
        l2reshist[it-1] = l2res

        # --- Calculate pressure coefficients along the airfoil and upper boundary ---
        cp[0] = -2 * ((phi[1, 0] - phi[0, 0]) / (x[1] - x[0]))
        cpu[0] = -2 * ((phi[1, kmax-1] - phi[0, kmax-1]) / (x[1] - x[0]))
        for j in range(1, jmax-1):
            cp[j] = -2 * 0.5 * (
                ((phi[j+1, 0] - phi[j, 0]) / dx[j]) * (dxp2[j] / dxm2[j]) +
                ((phi[j, 0] - phi[j-1, 0]) / dx[j]) * (dxm2[j] / dxp2[j])
            )
            cpu[j] = -2 * 0.5 * (
                ((phi[j+1, kmax-1] - phi[j, kmax-1]) / dx[j]) * (dxp2[j] / dxm2[j]) +
                ((phi[j, kmax-1] - phi[j-1, kmax-1]) / dx[j]) * (dxm2[j] / dxp2[j])
            )
        cp[-1] = -2 * ((phi[-1, 0] - phi[-2, 0]) / (x[-1] - x[-2]))
        cpu[-1] = -2 * ((phi[-1, kmax-1] - phi[-2, kmax-1]) / (x[-1] - x[-2]))
        
        # --- Compute number of supersonic points ---
        supersonic_point = np.sum(A < 0)
        # supersonic_points.append(supersonic_point)
        supersonic_points[it-1] = supersonic_point
        
        if show_prints:
            print(f"Number of supersonic points: {supersonic_point}")

        
        # Plot Cp contours if required
        if it < 10 or it % itplot == 0:
            for k in range(kmax):
                # k=0 boundary
                cpg[0, k] = -2 * ((phi[1, k] - phi[0, k]) / (x[1] - x[0]))
                for j in range(1, jmax-1):
                    cpg[j, k] = -2 * 0.5 * (
                        ((phi[j+1, k] - phi[j, k]) / dx[j]) * (dxp2[j] / dxm2[j]) +
                        ((phi[j, k] - phi[j-1, k]) / dx[j]) * (dxm2[j] / dxp2[j])
                    )
                cpg[-1, k] = -2 * ((phi[-1, k] - phi[-2, k]) / (x[-1] - x[-2]))
            if show_plots:
                plt.figure(2)
                plt.clf()
                cp_levels = 50
                cp_contours = plt.contour(xmesh, ymesh, cpg, cp_levels)
                plt.axis([-0.25, 1.25, 0.0, 1.5])
                plt.title(f'Cp contours at iteration {it}')
                plt.xlabel('x/c')
                plt.ylabel('y/c')
                # plt.pause(0.001)
    
    # Plot pressure coefficient along airfoil and upper wall
    if show_plots:
        plt.figure(3)
        plt.clf()
        plt.plot(x, -cp, 'm', marker='o', linewidth=2, markersize=10, label='Airfoil')
        plt.plot(x, -cpu, 'r', marker='+', linewidth=2, markersize=10, label='Upper wall')
        plt.plot(potx, -potcp/np.sqrt(1 - minf**2), 'gx', linewidth=2, markersize=10, label='External Data')
        plt.axis([-0.25, 1.25, -1.4, 1.0])
        plt.title(f'-Cp as a function of x/c at iteration {it}')
        plt.xlabel('x/c')
        plt.ylabel('-Cp')
        plt.grid(True)
        plt.legend()
        plt.savefig(f"images/pressure_coefficient-{condition}-no-omega.png", dpi=300)
        # plt.pause(0.001)
    
    # Plot the residual history
    if show_plots:
        plt.figure(4)
        # plt.clf()
        plt.semilogy(range(1, it+1), l2reshist[:it], 'm-', color=colors[condition-1], linewidth=2, label=f'case {condition-1}')
        plt.title('Log of L2-norm as a function of iteration')
        plt.xlabel('Iteration')
        plt.ylabel('Log of L2-norm')
        plt.legend()
        

        plt.grid(True)
        plt.figure(4)
        plt.savefig(f"images/residual_history-{condition}-no-omega.png", dpi=300)
                    
    # # plot supersonic points per iteration
    # if show_plots:
    #     plt.figure(6)
    #     plt.clf()
    #     plt.plot(range(1, it+1), supersonic_points[:it], 'm-', linewidth=2)
    #     plt.title('Supersonic points as a function of iteration')
    #     plt.xlabel('Iteration')
    #     plt.ylabel('Supersonic points')
    #     plt.grid(True)
    #     plt.savefig(f"images/supersonic_points-{condition}.png", dpi=300)
    
    # # plot the final Cp contours
    # if show_plots:
    #             plt.figure(2)
    #             plt.clf()
    #             cp_levels = 50
    #             cp_contours = plt.contour(xmesh, ymesh, cpg, cp_levels)
    #             plt.axis([-0.5, 1.5, 0.0, 2])
    #             plt.title(f'Cp contours at iteration {it}')
    #             plt.xlabel('x/c')
    #             plt.ylabel('y/c')
    #             plt.savefig(f"images/cp_contours-{condition}.png", dpi=300)
  
    # # Plot streamlines
    # if show_plots:
    #         # Compute velocity components from potential phi
    #     u = np.zeros((jmax, kmax))
    #     v = np.zeros((jmax, kmax))

    #     # Compute derivatives using central differences
    #     for j in range(1, jmax - 1):
    #         for k in range(1, kmax - 1):
    #             u[j, k] = (phi[j+1, k] - phi[j-1, k]) / (x[j+1] - x[j-1])
    #             v[j, k] = (phi[j, k+1] - phi[j, k-1]) / (y[k+1] - y[k-1])
                
    #     plt.figure()
    #     plt.quiver(xmesh, ymesh, u, v, color='b', scale=10)
    #     plt.title('Velocity Vector Field')
    #     plt.xlabel('x/c')
    #     plt.ylabel('y/c')
    #     plt.axis([-0.25, 1.25, 0.0, 1.5])
    #     plt.grid(True)
    #     plt.savefig(f"images/velocity_vectors-{condition}.png", dpi=300)
    
    return supersonic_point, minf, airfoil_type, it


def collect_times():
    for i in range(2, 8):
        print("====================================")
        print(f"Running condition {i}")
        start_time = time.time()
        supersonic_points, minf, airfoil_type, final_iteration = main(i, False, False)
        end_time = time.time()
        airfoil_name = "NACA0010" if airfoil_type == 0 else "Biconvex"
        print(f"Airfoil: {airfoil_name} with M_inf = {minf}, and {supersonic_points} supersonic points found in {final_iteration} iterations")
        print("Elapsed time:", end_time - start_time, "seconds")
        print("====================================")

def collect_figures():
    for i in range(5, 6):
        print(f"Running condition {i}")
        supersonic_points, minf, airfoil_type, final_iteration = main(i, True, False)

if __name__ == '__main__':
    # collect_times()
    # collect_figures()    
    main(4, True, True)
